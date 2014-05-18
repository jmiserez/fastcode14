#!/usr/bin/python

import subprocess

example_config = {
    'versions'             : [ "01ref_matlab" ],
    'logtofile'            : True,
	'performance_counters' : True,
	'optimization_flags'   : "-O3 -m64 -march=native -mno-abm -fno-tree-vectorize -g",
	'debug'                : True,
	'ndebug'               : False,
	'driver_args'          : "--store zzz --val ../testdata/house_out/A-3-1-1-1.tif --threshold 1.0 752:1:752 500:1:500 1.0 1.0 1.0 ../testdata/srcImages/A.0.tif ../testdata/srcImages/A.1.tif"
}

def add_config(key, val, f):
	print >> f, "%-10s = %s" % (key, val)

def write_config(version, config):
	with open("Make.sysconfig", "w") as f:
		add_config("VERSION", version, f)

		if config['performance_counters']:
			cost = "-DCOST_MEASURE"
		else:
			cost = ""
		add_config("CF_COST", cost, f)

		add_config("CF_OPT", config['optimization_flags'], f)

		debug = ""
		if config['debug']:
			debug += "-DDEBUG"
		if config['ndebug']:
			debug += "-DNDEBUG"
		add_config("CF_DEBUG", debug, f)

		add_config("CF_CONFIG", "$(CF_OPT) $(CF_COST) $(CF_DEBUG)", f)

def title(t):
	print 70*'-'
	print t
	print 70*'-'

if __name__ == "__main__":
	config = example_config
	title("Configuration")
	print "Building and Running: %s" % config['versions']
	print "Build Config:"
	for key in config:
		print "  %-20s : %s" % (key, config[key])
	
	title("Doing Builds")
	for version in config['versions']:
		print "%s: writing build configuration" % version
		write_config(version, config)
		logfile = version + ".build.log"
		print "%s: cleaning " % version
		with open(logfile, "w") as log:
			if config['logtofile']:
				ret = subprocess.call(["make", "-fMake.system", "clean"], stdout=log)
			else:
				ret = subprocess.call(["make", "-fMake.system", "clean"])
			if ret == 0:
				print "%s: clean successful. output in %s" % (version, logfile)
			else:
				print "%s: ERROR cleaning. see %s for details" % (version, logfile)

	for version in config['versions']:
		print "%s: writing build configuration" % version
		write_config(version, config)
		logfile = version + ".build.log"
		print "%s: building" % version
		with open(logfile, "a") as log:
			if config['logtofile']:
				ret = subprocess.call(["make", "-fMake.system"], stdout=log)
			else:
				ret = subprocess.call(["make", "-fMake.system"])
			if ret == 0:
				print "%s: build successful. output in %s" % (version, logfile)
			else:
				print "%s: ERROR building. see %s for details" % (version, logfile)

	title("Running Builds")
	for version in config['versions']:
		print "%s: running" % version
		binary = "./bin/driver_" + version
		runfile = version + ".run.log"
		with open(runfile, "w") as rf:
			arg = binary + " " + config['driver_args']
			#print arg
			#arg = arg.split()
			#print arg
			if config['logtofile']:
				ret = subprocess.call(arg, stdout=rf, shell=True)
			else:
				ret = subprocess.call(arg, shell=True)
			if ret == 0:
				print "%s: run successful. output in %s" % (version, runfile)
			else:
				print "%s: ERROR running. see %s for details" % (version, runfile)

	title("Done")
