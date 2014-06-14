#!/usr/bin/python

import subprocess
import sys
import os

my_config = {
    'versions'             : [
                               "01ref_matlab",
                               "consolidated.naive_options",
                               "consolidated.naive",
                               "consolidated.store_grey",
                               "consolidated.onestep",
                               "consolidated.blocking",
                               "consolidated.inline2",
                               "consolidated.inline2x2",
                               "consolidated.inline2x4",
                               "consolidated.avx"
                             ],
    'do_benchmark'         : True,
    'do_develop'           : False,
    'do_gprof'             : False,
    'logtofile'            : False,
    'perfunc'              : True,
    'no_pyramids'          : False,
    'openmode'             : 'a',
    'optimization_flags'   : "-O3 -m64 -march=corei7-avx",
#    'optimization_flags'   : "-O3 -Ofast -m64 -march=corei7-avx -ffast-math",
#    'optimization_flags'  : "-O3 -m64 -march=native -mno-abm -fno-tree-vectorize",
#    'optimization_flags'  : "-O0 -m64 -march=native",
    'warmup_count'         : 1,
    'warmup_dev'           : 0,
    'warmup_benchmark'     : 2,
     'dev_driver_args'     : "--q --s out --v ../testdata/house_out/dbg_benchmark-3-1.0-1.0-1.0.tif --w 1024 --h 1024 --t 0.1 "
                              "1024:1:1024 1024:1:1024 "
                              "1.0 1.0 1.0 ../testdata/srcImages/dbg_benchmark.0.tif ../testdata/srcImages/dbg_benchmark.1.tif ../testdata/srcImages/dbg_benchmark.2.tif ../testdata/srcImages/dbg_benchmark.3.tif",
#     'driver_args'         : "--q --v ../testdata/house_out/benchmark-3-1.0-1.0-1.0.tif --w 4096 --h 4096 --t 0.1 "
#                              "16:1:4096 16:1:4096 "
#                              "1.0 1.0 1.0 ../testdata/srcImages/benchmark.0.tif ../testdata/srcImages/benchmark.1.tif ../testdata/srcImages/benchmark.2.tif ../testdata/srcImages/benchmark.3.tif"
     'driver_args'         : "--q "
                              "128:1:4096 128:1:4096 "
                              "1.0 1.0 1.0 ../testdata/srcImages/benchmark.0.tif ../testdata/srcImages/benchmark.1.tif ../testdata/srcImages/benchmark.2.tif ../testdata/srcImages/benchmark.3.tif"

}


#    'driver_args'          : "--q --v ../testdata/house_out/A-3-1-1-1.tif --w 752 --h 500 --t 0.1 "
#                             "752:1:752 500:1:500 "
#                             "1.0 1.0 1.0 ../testdata/srcImages/A.0.tif ../testdata/srcImages/A.1.tif ../testdata/srcImages/A.2.tif ../testdata/srcImages/A.3.tif",
#    'dev_driver_args'      : "--q --s out --v ../testdata/house_out/A-3-1-1-1.tif --w 752 --h 500 --t 0.1 "
#                             "752:1:752 500:1:500 "
#                             "1.0 1.0 1.0 ../testdata/srcImages/A.0.tif ../testdata/srcImages/A.1.tif ../testdata/srcImages/A.2.tif ../testdata/srcImages/A.3.tif",
#     'driver_args'          : "--q "
#                              "10:200:3220 10:200:3220 "
#                              "1.0 1.0 1.0 ../testdata/srcImages/huge_House.0.tif ../testdata/srcImages/huge_House.1.tif ../testdata/srcImages/huge_House.2.tif ../testdata/srcImages/huge_House.3.tif",
#     'dev_driver_args'      : "--b "
#                              "1:1:1 1:1:1 "
#                              "1.0 1.0 1.0 ../testdata/srcImages/huge_House.0.tif ../testdata/srcImages/huge_House.1.tif ../testdata/srcImages/huge_House.2.tif ../testdata/srcImages/huge_House.3.tif",
#     'dev_driver_args'      : "--q --s out --v ../testdata/house_out/benchmark-3-1.0-1.0-1.0.tif --w 1024 --h 1024 --t 0.1 "
#                              "16:1:4096 16:1:4096 "
#                              "1.0 1.0 1.0 ../testdata/srcImages/benchmark.0.tif ../testdata/srcImages/benchmark.1.tif ../testdata/srcImages/benchmark.2.tif ../testdata/srcImages/benchmark.3.tif"
#    'dev_driver_args'      : "--q --s out --v ../testdata/house_out/A-3-1-1-1.tif --w 752 --h 500 --t 0.1 "
#                             "752:1:752 500:1:500 "
#                             "1.0 1.0 1.0 ../testdata/srcImages/House.0.tif ../testdata/srcImages/House.1.tif ../testdata/srcImages/House.2.tif ../testdata/srcImages/House.3.tif"

def add_config(key, val, f):
	print >> f, "%-10s = %s" % (key, val)

def update_for_gprof(config):
	config['cost_measure'] = False
	config['debug'] = False
	config['gprof'] = True
	config['warmup_count'] = 0
	config['read_flops'] = False

def update_for_development(config):
	config['cost_measure'] = True
	config['logtofile'] = False
	config['debug'] = True
	config['gprof'] = False
	config['warmup_count'] = config['warmup_dev']
	config['read_flops'] = False
	config['driver_args'] = config['dev_driver_args']

def update_for_benchmark_cost(config):
	config['cost_measure'] = True
	config['debug'] = False
	config['gprof'] = False
	config['warmup_count'] = 0
	config['read_flops'] = False

def update_for_benchmark_performance(config):
	config['cost_measure'] = False
	config['debug'] = False
	config['gprof'] = False
	config['warmup_count'] = config['warmup_benchmark']
	config['read_flops'] = True

def write_config(version, config):
	with open("Make.sysconfig", "w") as f:
		versionparts = version.split('.')
		if len(versionparts) > 1:
			add_config("VERSION", versionparts[0], f)
			add_config("VERSION_SUFFIX", "."+versionparts[1], f)
		else:
			add_config("VERSION", version, f)
			add_config("VERSION_SUFFIX", "", f)


		if 'cost_measure' in config and config['cost_measure']:
			cost = "-DCOST_MEASURE"
		else:
			cost = ""
		add_config("CF_COST", cost, f)

		add_config("CF_WARMUP", "-DWARMUP_COUNT="+str(config['warmup_count']), f)
		add_config("CF_OPT", config['optimization_flags'], f)

		debug = ""
		if 'debug' in config and config['debug']:
			# enable debugging output
			debug += " -DDEBUG -g"
		else:
			#disable all debugging output
			debug += " -DNDEBUG"
		if 'gprof' in config and config['gprof']:
			debug += " -pg"
		if 'read_flops' in config and config['read_flops']:
			debug += " -DREADFLOPS"
		if 'perfunc' in config and config['perfunc']:
			debug += " -DCOST_MODEL_PERFUNC"
		if 'no_pyramids' in config and config['no_pyramids']:
			debug += " -DNO_PYRAMIDS"
		add_config("CF_DEBUG", debug, f)

		add_config("CF_CONFIG", "$(CF_OPT) $(CF_COST) $(CF_WARMUP) $(CF_DEBUG)", f)

def title(t):
	print 70*'-'
	print t
	print 70*'-'
	sys.stdout.flush()

def build_and_run(config):
	build_error = False
	title("Doing Builds")
	for version in config['versions']:
		print "%s: writing build configuration" % version
		write_config(version, config)
		logfile = version + ".build.log"
		print "%s: cleaning " % version
		with open(logfile, config['openmode']) as log:
			if config['logtofile']:
				ret = subprocess.call(["make", "-fMake.system", "clean"], stdout=log)
			else:
				ret = subprocess.call(["make", "-fMake.system", "clean"])
			sys.stdout.flush()
			log.flush()
			if ret == 0:
				print "%s: clean successful. output in %s" % (version, logfile)
			else:
				print "%s: ERROR cleaning. see %s for details" % (version, logfile)
				build_error = True
			sys.stdout.flush()

	for version in config['versions']:
		print "%s: writing build configuration" % version
		write_config(version, config)
		logfile = version + ".build.log"
		print "%s: building" % version
		with open(logfile, config['openmode']) as log:
			if config['logtofile']:
				ret = subprocess.call(["make", "-fMake.system"], stdout=log)
			else:
				ret = subprocess.call(["make", "-fMake.system"])
			sys.stdout.flush()
			log.flush()
			if ret == 0:
				print "%s: build successful. output in %s" % (version, logfile)
			else:
				print "%s: ERROR building. see %s for details" % (version, logfile)
				build_error = True
			sys.stdout.flush()
	if not build_error:
		title("Running Builds")
		for version in config['versions']:
			print "%s: running" % version
			binary = "./bin/driver_" + version
			runfile = version + ".run.log"
			with open(runfile, config['openmode']) as rf:
				arg = binary + " " + config['driver_args']
				#print arg
				#arg = arg.split()
				if config['logtofile']:
					rf.write(arg+'\n')
					ret = subprocess.call(arg, stdout=rf, shell=True)
				else:
					print arg
					ret = subprocess.call(arg, shell=True)
				sys.stdout.flush()
				rf.flush()
				if ret == 0:
					print "%s: run successful. output in %s" % (version, runfile)
					if 'gprof' in config and config['gprof']:
						title("Profiling Builds")
						for version in config['versions']:
							print "%s: profiling" % version
							binary = "./bin/driver_" + version
							profilefile = version + ".profile.log"
							with open(profilefile, config['openmode']) as pf:
								arg = "gprof " + binary + " gmon.out"
								if config['logtofile']:
									pf.write(arg+'\n')
									ret = subprocess.call(arg, stdout=pf, shell=True)
								else:
									print arg
									ret = subprocess.call(arg, shell=True)
								sys.stdout.flush()
								pf.flush()
								if ret == 0:
									print "%s: profiling successful. output in %s" % (version, profilefile)
								else:
									print "%s: ERROR profiling. see %s for details" % (version, profilefile)
									sys.stdout.flush()
				else:
					print "%s: ERROR running. see %s for details" % (version, runfile)
	else:
		sys.exit("\n ===> BUILD ERROR <=== \n")
	title("Done")

if __name__ == "__main__":
	config = my_config
	title("Configuration")
	print "Building and Running: %s" % config['versions']
	print "Build Config:"
	for key in config:
		print "  %-20s : %s" % (key, config[key])
	sys.stdout.flush()

	if config['do_develop']:
		cfg = config.copy()
		title("DEV: DEBUGGING")
		update_for_development(cfg)
		build_and_run(cfg);
		title("DEV: DONE")
	else:
		if config['do_benchmark']:
			try:
				os.remove("flops.out")
			except:
				None
			cfg = config.copy()
			update_for_benchmark_cost(cfg)
			build_and_run(cfg)
			title("BENCHMARK: MEASURING PERFORMANCE")
			update_for_benchmark_performance(cfg)
			build_and_run(cfg);
			title("BENCHMARK: DONE")
		else:
			if config['do_gprof']:
				cfg = config.copy()
				update_for_gprof(cfg)
				build_and_run(cfg)
			else:
				cfg = config.copy()
				build_and_run(cfg)

