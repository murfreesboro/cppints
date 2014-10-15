"""
main module
"""
__author__ = "Fenglai Liu"
import sys
import os
import headGen

# we have some parameters defined
# module name is the working folder you want to handled
# like twobodyoverlap
# all meaning we will generate entry file as well as
# head file for all working modules
# auxL and maxL is only for some modules, see the headGen.py
# function of useThisCPPFile for more information
module = "all"
maxL = 100
maxAuxL = 100

# now get the dir
has_Dir = False
if len(sys.argv) == 1:
	print "We need at least one argument, which is the input cpp files dir!\n"
	sys.exit()
else:
	for arg in sys.argv:
		if arg == sys.argv[0]:
			continue
		tmp = arg.split("=")
		if len(tmp) != 2:
			print arg
			print "the input argument can not be parsed"
			sys.exit()
		key = tmp[0].strip().upper()
		val = tmp[1].strip()
		if key == "DIR":
			topDir = val
			has_Dir = True
		elif key == "MOD":
			module = val
		elif key == "MAXL":
			maxL = int(val)
		elif key == "AUXL":
			maxAuxL = int(val)
		else:
			print key
			print "the key got from argument can not match anything"
			sys.exit()

# check that whether dir has been defined or not
if not has_Dir:
    print "You did not define the working dir, so we have nothing to work with"
    sys.exit()

# on the top of dir, we should have
# three sub-dir whose name like this
dirList = ["energy", "first_deriv", "second_deriv"]
for iDir in dirList:
    d = topDir + "/" + iDir
    if os.path.exists(d):
        projects = os.listdir(d)
        get_module = False
        for iProj in projects:

            # check with the module
            if module != "all":
                if iProj != module:
                    continue
            else:
                get_module = True

            # now get into real work
            get_module = True
            workDir = d + "/" + iProj
            if os.path.isdir(workDir):
                headGen.filesCreation(workDir, maxL, maxAuxL)

# whether we have found it?
        if not get_module:
            print module
            print "we did not find the matching module"
            sys.exit()
