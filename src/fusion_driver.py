#!/usr/bin/python

import subprocess

example_config = {
	'versions'             : [ "01ref_matlab" ],
	'performance_counters' : True,
	'optimization_flags'   : "-O3 -fno-tree-vectorize",
	'debug'                : True,
	'ndebug'               : False,
	'driver_args'          : "100:25:150 100:25:300 1.0 1.0 1.0 ../testdata/srcImages/A.0.tif"
	#'driver_args'          : ["100:25:150", "100:25:300", "1.0", "1.0", "1.0", "../testdata/srcImages/A.0.tif"]
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
		print "%s: building" % version
		logfile = version + ".build.log"
		with open(logfile, "w") as log:
			ret = subprocess.call(["make", "-fMake.system"], stdout=log)
		if ret == 0:
			print "%s: build successful. output in %s" % (version, logfile)
		else:
			print "%s: error building. see %s for details" % (version, logfile)
	
	title("Running Builds")
	for version in config['versions']:
		print "%s: running" % version
		binary = "./bin/driver_" + version
		runfile = version + ".run.log"
		with open(runfile, "w") as rf:
			#subprocess.call([[binary] + config['driver_args']], stdout=run)
			arg = binary + " " + config['driver_args']
			#print arg
			#arg = arg.split()
			#print arg
			subprocess.call(arg, stdout=rf, shell=True)
		print "%s: finished. output in %s" % (version, runfile)
	
	title("Done")
