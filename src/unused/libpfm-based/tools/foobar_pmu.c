/*
 * foobar_mpu.c - modified version of task_cpu.c from the examples shipped with
 *                the libpfm library.
 *
 * Original copyright notice for task_cpu.c
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright (c) 2010 Google, Inc
 * Contributed by Stephane Eranian <eranian@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
#include <sys/types.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <signal.h>
#include <sys/wait.h>
#include <locale.h>
#include <sys/ioctl.h>
#include <err.h>

#include "perf_util.h"

#define MAX_GROUPS	16

typedef struct {
	const char *events[MAX_GROUPS];
	const char *filename;
	int to_file;
	int append;
	int num_groups;
	int inherit;
	int print;
	int pin;
	pid_t pid;
	FILE *file;
} options_t;

static options_t options;
static volatile int quit;

int
child(char **arg)
{
	/*
	 * execute the requested command
	 */
	execvp(arg[0], arg);
	errx(1, "cannot exec: %s\n", arg[0]);
	/* not reached */
}

static void
read_groups(perf_event_desc_t *fds, int num)
{
	uint64_t *values = NULL;
	size_t new_sz, sz = 0;
	int i, evt, ret;

	/*
	 * 	{ u64		nr;
	 * 	  { u64		time_enabled; } && PERF_FORMAT_ENABLED
	 * 	  { u64		time_running; } && PERF_FORMAT_RUNNING
	 * 	  { u64		value;
	 * 	    { u64	id;           } && PERF_FORMAT_ID
	 * 	  }		cntr[nr];
	 * 	} && PERF_FORMAT_GROUP
	 *
	 * we do not use FORMAT_ID in this program
	 */

	for (evt = 0; evt < num; ) {
		int num_evts_to_read;

		
		num_evts_to_read = 1;
		new_sz = sizeof(uint64_t) * 3;

		if (new_sz > sz) {
			sz = new_sz;
			values = realloc(values, sz);
		}

		if (!values)
			err(1, "cannot allocate memory for values\n");

		ret = read(fds[evt].fd, values, new_sz);
		if (ret != new_sz) { /* unsigned */
			if (ret == -1)
				err(1, "cannot read values event %s", fds[evt].name);

			/* likely pinned and could not be loaded */
			warnx("could not read event %d, tried to read %zu bytes, but got %d",
				evt, new_sz, ret);
		}

		/*
		 * propagate to save area
		 */
		for (i = evt; i < (evt + num_evts_to_read); i++) {
			/*
			 * scaling because we may be sharing the PMU and
			 * thus may be multiplexed
			 */
			fds[i].prev_value = fds[i].value;
			fds[i].value = perf_scale(values);
			fds[i].enabled = values[1];
			fds[i].running = values[2];
		}
		evt += num_evts_to_read;
	}
	if (values)
		free(values);
}

static void
print_counts(perf_event_desc_t *fds, int num)
{
	int i;
	char* buffer = (char*) malloc(512);
	if (buffer==NULL)
		err(1, "cannot allocate memory for sprintf buffer\n");

	read_groups(fds, num);

	for(i=0; i < num; i++) {
		double ratio;
		uint64_t val;

		val = fds[i].value - fds[i].prev_value;

		ratio = 0.0;
		if (fds[i].enabled)
			ratio = 1.0 * fds[i].running / fds[i].enabled;

		/* separate groups */
		//if (perf_is_group_leader(fds, i))
		//	putchar('\n');

		if (fds[i].value < fds[i].prev_value) {
			sprintf(buffer, "%llu %llu %llu %f %s (Error: Inconsistent scaling! -1/value/prev_value/0/name)\n",
			(long long unsigned int)(-1),
			(long long unsigned int)(fds[i].value),
			(long long unsigned int)(fds[i].prev_value),
			(double)(0),
			fds[i].name);
			continue;
		}
		sprintf(buffer, "%llu %llu %llu %f %s (value/enabled/running/scaling/name)\n",
			(long long unsigned int)val,
			(long long unsigned int)(fds[i].enabled),
			(long long unsigned int)(fds[i].running),
			(double)((1.0-ratio)*100.0),
			fds[i].name);
		if(options.to_file){
			fprintf(options.file,"%s",buffer);
		}else{
			printf("%s",buffer);
		}
	}
	free(buffer);
}

static void sig_handler(int n)
{
	quit = 1;
}

static int error(perf_event_desc_t *fds, pid_t pid){
	free(fds);
	if (!options.pid)
		kill(SIGKILL, pid);

	/* free libpfm resources cleanly */
	pfm_terminate();

	return -1;
}

int
parent(char **arg)
{
	perf_event_desc_t *fds, *fds_cpus[1];
	int status, ret, i, num_fds = 0, grp, group_fd;
	int ready[2], go[2];
	char buf;
	pid_t pid;

	if (pfm_initialize() != PFM_SUCCESS)
		errx(1, "libpfm initialization failed");

	memset(fds_cpus, 0, sizeof(fds_cpus));


	for (grp = 0; grp < options.num_groups; grp++) {
		num_fds = 0;
		ret = perf_setup_list_events(options.events[grp], &fds_cpus[0], &num_fds);
		if (ret || !num_fds)
			exit(1);
	}
	fds = fds_cpus[0];

	pid = options.pid;
	if (!pid) {
		ret = pipe(ready);
		if (ret)
			err(1, "cannot create pipe ready");

		ret = pipe(go);
		if (ret)
			err(1, "cannot create pipe go");


		/*
		 * Create the child task
		 */
		if ((pid=fork()) == -1)
			err(1, "Cannot fork process");

		/*
		 * and launch the child code
		 *
		 * The pipe is used to avoid a race condition
		 * between for() and exec(). We need the pid
		 * of the new tak but we want to start measuring
		 * at the first user level instruction. Thus we
		 * need to prevent exec until we have attached
		 * the events.
		 */
		if (pid == 0) {
			close(ready[0]);
			close(go[1]);

			/*
			 * let the parent know we exist
			 */
			close(ready[1]);
			if (read(go[0], &buf, 1) == -1)
				err(1, "unable to read go_pipe");


			exit(child(arg));
		}

		close(ready[1]);
		close(go[0]);

		if (read(ready[0], &buf, 1) == -1)
			err(1, "unable to read child_ready_pipe");

		close(ready[0]);
	}

	for(i=0; i < num_fds; i++) {
		int is_group_leader; /* boolean */

		is_group_leader = perf_is_group_leader(fds, i);
		if (is_group_leader) {
			/* this is the group leader */
			group_fd = -1;
		} else {
			group_fd = fds[fds[i].group_leader].fd;
		}

		/*
		 * create leader disabled with enable_on-exec
		 */
		if (!options.pid) {
			fds[i].hw.disabled = is_group_leader;
			fds[i].hw.enable_on_exec = is_group_leader;
		}

		fds[i].hw.read_format = PERF_FORMAT_SCALE;
		/* request timing information necessary for scaling counts */

		if (options.inherit)
			fds[i].hw.inherit = 1;

		if (options.pin && is_group_leader)
			fds[i].hw.pinned = 1;
		fds[i].fd = perf_event_open(&fds[i].hw, pid, -1, group_fd, 0);
		if (fds[i].fd == -1) {
			warn("cannot attach event%d %s", i, fds[i].name);
			return error(fds, pid);
		}
	}

	if (!options.pid)
		close(go[1]);

	if (options.print) {
		if (!options.pid) {
			while(waitpid(pid, &status, WNOHANG) == 0) {
				sleep(1);
				print_counts(fds, num_fds);
			}
		} else {
			while(quit == 0) {
				sleep(1);
				print_counts(fds, num_fds);
			}
		}
	} else {
		if (!options.pid)
			waitpid(pid, &status, 0);
		else {
			pause();
			for(i=0; i < num_fds; i++) {
				ioctl(fds[i].fd, PERF_EVENT_IOC_DISABLE, 0);
			}
		}
			print_counts(fds, num_fds);
	}

	for(i=0; i < num_fds; i++){
		close(fds[i].fd);
	}
	free(fds);

	/* free libpfm resources cleanly */
	pfm_terminate();

	return 0;
}

static void
usage(void)
{
	printf("usage: foobar_pmu [-h] [-i] [-g] [-p] [-P] [-t pid] [-e event1,event2,...] cmd\n"
		"-h\t\tget help\n"
		"-i\t\tinherit across fork\n"
		"-p\t\tprint counts every second\n"
		"-P\t\tpin events\n"
		"-t pid\tmeasure existing pid\n"
		"-f filename\twrite pmu values to file\n"
		"-a append to file instead of overwrite (only valid with -f flag)\n"
		"-e ev,ev\tgroup of events to measure (multiple -e switches are allowed)\n"
		"(this program is based on the libpfm example program task_cpu)"
		);
}

int
main(int argc, char **argv)
{
	int c;

	setlocale(LC_ALL, "");

	while ((c=getopt(argc, argv,"+he:ipPt:f:b:a")) != -1) {
		switch(c) {
		case 'a':
			options.append = 1;
			break;
		case 'e':
			if (options.num_groups < MAX_GROUPS) {
				options.events[options.num_groups++] = optarg;
			} else {
				errx(1, "you cannot specify more than %d groups.\n",
					MAX_GROUPS);
			}
			break;
		case 'p':
			options.print = 1;
			break;
		case 'P':
			options.pin = 1;
			break;
		case 'i':
			options.inherit = 1;
			break;
		case 't':
			options.pid = atoi(optarg);
			break;
		case 'f':
			options.to_file = 1;
			options.filename = optarg;
			break;
		case 'h':
			usage();
			exit(0);
		case '?':
			if(optopt == 'f')
				errx(1, "you must specify a file name");
		default:
			errx(1, "unknown error");
		}
	}

	if (options.num_groups == 0) {
		options.events[0] = "PERF_COUNT_HW_INSTRUCTIONS";
		options.num_groups = 1;
	}
	if (!argv[optind] && !options.pid)
		errx(1, "you must specify a command to execute or a thread to attach to\n");
	if(options.to_file){
		if (options.append) {
			options.file = fopen(options.filename,"a"); //append
		} else {
			options.file = fopen(options.filename,"w"); //overwrite
		}
		if(!options.file)
			errx(1, "could not open file for writing");
	}

	signal(SIGINT, sig_handler);

	return parent(argv+optind);
}
