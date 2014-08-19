"""
main module 
"""
__author__  = "Fenglai Liu"
import sys
import os
import re
import evaluation

# now get the dir
if len(sys.argv) == 1:
    print "We need at least one argument, which is the input cpp files dir!\n"
    sys.exit()
else:
    topDir = sys.argv[1]

# whether user want to specify a dir?
whichDir = "ALL"
if len(sys.argv) > 2:
    whichDir = sys.argv[2]
    print whichDir

# on the top of dir, we should have 
# three sub-dir whose name like this
dirList = ["energy","first_deriv","second_deriv"]
for iDir in dirList:
    d = topDir + "/" + iDir
    if os.path.exists(d):
        if whichDir == "ALL":
            for iProj in projects:
                workDir = d + "/" + iProj

                # check whether it exists
                if not os.path.exists(workDir):
                    print workDir
                    print "the given project dir is wrong"
                    sys.exit()

                if os.path.isdir(workDir):
                    files = os.listdir(workDir)
                    for iFile in files:
                        f = workDir + "/" + iFile
                        if os.path.isfile(f):

                            # we only care about cpp file
                            extension = os.path.splitext(f)[1]
                            if extension != ".cpp":
                                continue

                            # for the file split mode, we will report an error 
                            # here
                            if re.search(r"_vrr", f) is not None or re.search(r"_hrr", f) is not None: 
                                print "we are unable to handle the file split situation"
                                print "please re-generate your cpp files"
                                print "with the definition of single_cpp_max_l = 100"
                                print "so that to prevent the file is splitting into vrr/hrr files"
                                sys.exit()

                            # initilize the evaluation
                            evaluation.initilize(f)

                            # do VRR parse
                            evaluation.vrr_analyze(f)

                            # if we have HRR, do HRR
                            if evaluation.with_hrr():
                                evaluation.hrr_analyze(f)

                            # now print out to screen
                            evaluation.infor_print()

                            # print out the shell quartet information
                            evaluation.shell_quartets_infor(f)

        else:
            workDir = d + "/" + whichDir

            # check whether it exists
            if not os.path.exists(workDir):
               print workDir
               print "the given project dir is wrong"
               sys.exit()

            if os.path.isdir(workDir):
               files = os.listdir(workDir)
               for iFile in files:
                   f = workDir + "/" + iFile
                   if os.path.isfile(f):

                       # we only care about cpp file
                       extension = os.path.splitext(f)[1]
                       if extension != ".cpp":
                           continue

                       # for the file split mode, we will report an error 
                       # here
                       if re.search(r"_vrr", f) is not None or re.search(r"_hrr", f) is not None: 
                           print "we are unable to handle the file split situation"
                           print "please re-generate your cpp files"
                           print "with the definition of single_cpp_max_l = 100"
                           print "so that to prevent the file is splitting into vrr/hrr files"
                           sys.exit()

                       # initilize the evaluation
                       evaluation.initilize(f)

                       # do VRR parse
                       evaluation.vrr_analyze(f)

                       # if we have HRR, do HRR
                       if evaluation.with_hrr():
                           evaluation.hrr_analyze(f)

                       # now print out to screen
                       evaluation.infor_print()

                       # print out the shell quartet information
                       evaluation.shell_quartets_infor(f)


