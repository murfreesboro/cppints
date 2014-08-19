"""
main module 
"""
__author__  = "Fenglai Liu"
import sys
import os
import makefilesGen
import infor

# now get the dir
if len(sys.argv) == 3:
    topDir = sys.argv[1]
    settingFile = sys.argv[2]
else:
    print "We need two arguments, one is the input cpp files dir\n"
    print "This one is necessary, and does not have default value\n"
    print "The other is the setting file for generating makefiles\n"
    sys.exit()

# parse the information from setting file
infor.set_project(settingFile)

# on the top of dir, we should have 
# three sub-dir whose name like this
dirList = ["energy","first_deriv","second_deriv"]
for iDir in dirList:
    d = topDir + "/" + iDir
    if os.path.exists(d):
        projects = os.listdir(d)
        for iProj in projects:
            workDir = d + "/" + iProj
            if os.path.isdir(workDir):

                # does the work dir is in the list?
                if infor.hasJob(workDir):
                    makefilesGen.doMakefile(workDir)

