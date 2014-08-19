"""
main module 
"""
__author__  = "Fenglai Liu"
import sys
import os
import headGen

# now get the dir
if len(sys.argv) == 1:
    print "We need at least one argument, which is the input cpp files dir!\n"
    sys.exit()
else:
    topDir = sys.argv[1]

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
                headGen.filesCreation(workDir,d)

