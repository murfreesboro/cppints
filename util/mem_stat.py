#!/usr/bin/python
"""
this file is used to generate memory statistics for the cpp files
in terms of given module. For the cpp files, it tries to find out
two things:

1 whether it use memory on heap?
2 how much memory on heap it really uses?
"""
__author__  = "Fenglai Liu"
import sys
import os
import re

def getL(s):
    """
    whether the given s is shell symbol
    if it is, we return it's L Code
    else return -1
    """
    if s == "S" or s == "s":
        return 0
    elif s == "P" or s == "p":
        return 1
    elif s == "SP" or s == "sp":
        return 100
    elif s == "D" or s == "d":
        return 2
    elif s == "F" or s == "f":
        return 3
    elif s == "G" or s == "g":
        return 4
    elif s == "H" or s == "h":
        return 5
    elif s == "I" or s == "i":
        return 6
    elif s == "K" or s == "k":
        return 7
    elif s == "L" or s == "l":
        return 8
    elif s == "M" or s == "m":
        return 9
    elif s == "N" or s == "n":
        return 10
    elif s == "O" or s == "o":
        return 11
    elif s == "Q" or s == "q":
        return 12
    elif s == "R" or s == "r":
        return 13
    elif s == "T" or s == "t":
        return 14
    elif s == "U" or s == "u":
        return 15
    elif s == "V" or s == "v":
        return 16
    elif s == "W" or s == "w":
        return 17
    elif s == "X" or s == "x":
        return 18
    elif s == "Y" or s == "y":
        return 19
    elif s == "Z" or s == "z":
        return 20
    else:

        # if L>20, then L symbol
        # is defined as L+number
        # we will see whether this
        # is the case
        if s[0] == "L" and len(s) > 1:
            n = s[1:]
            if n.isdigit():
                return int(n)

        # now this is not the
        # angular momentum
        return -1

def getLCode(filename):
    """
    get the angular momentum code from the given file name
    """

    # drop the extension first
    fname = os.path.basename(filename)
    name  = os.path.splitext(fname)[0]
    nameComponents = name.split("_")

    # now let's get LCode
    L_list = [ ]
    for i in nameComponents:
        lcode = getL(i)
        if lcode >= 0:
            L_list.append(lcode)

    # now return the L list
    return L_list

def checkMem(filename):
    """
    for the given file name, check the memory usage
    """
    nMem = 0

    # whether the cpp file has _vrr and _hrr files?
    name = os.path.splitext(filename)[0]
    vrrFile = name + "_vrr.cpp"
    hrrFile = name + "_hrr.cpp"
    if os.path.isfile(vrrFile):
        cpp = open(vrrFile,"r")
        while True:
            line = cpp.readline()
            if not line: break
            if re.search(r"(?i)getNewMemPos", line) is not None:
		begin = line.index("(")
		end   = line.index(")")
		num   = line[begin+1:end]
		n     = int(num)
                nMem  = nMem + n
        cpp.close()
    if os.path.isfile(hrrFile):
        cpp = open(hrrFile,"r")
        while True:
            line = cpp.readline()
            if not line: break
            if re.search(r"(?i)getNewMemPos", line) is not None:
		begin = line.index("(")
		end   = line.index(")")
		num   = line[begin+1:end]
		n     = int(num)
                nMem  = nMem + n
        cpp.close()

    # now let's open the main cpp file
    cpp = open(filename,"r")
    while True:
        line = cpp.readline()
        if not line: break
        if re.search(r"(?i)getNewMemPos", line) is not None:
	    begin = line.index("(")
	    end   = line.index(")")
	    num   = line[begin+1:end]
	    n     = int(num)
            nMem  = nMem + n
    cpp.close()

    # finally, let's return the nMem
    return nMem

#################################################
# now this is main script for memory statistics #
#################################################
if len(sys.argv) == 3:
    topDir  = sys.argv[1]
    wDir    = sys.argv[2]
else:
    print "We need two arguments, one is the input cpp files dir\n"
    print "This one is necessary, and does not have default value\n"
    print "The other is the work dir name, like twobodyoverlap\n"
    sys.exit()

# on the top of dir, we should have 
# three sub-dir whose name like this
dirList = ["energy","first_deriv","second_deriv"]
for iDir in dirList:
    d = topDir + "/" + iDir
    if os.path.exists(d):
        projects = os.listdir(d)
        gotIt = False
        for iProj in projects:
            workDir = d + "/" + iProj
            if os.path.isdir(workDir) and wDir == iProj:
                gotIt = True

                # now go over the cpp files
                files = os.listdir(workDir)
                for iFile in files:

                    # we only concentrate on cpp file
                    if iFile.find("cpp") < 0:
                        continue

                    # the *_vrr.cpp and *_hrr.cpp will be consider
                    # together with the main cpp file
                    if iFile.find("_vrr") > 0:
                        continue
                    if iFile.find("_hrr") > 0:
                        continue

                    # now do the real job
                    fname = workDir + "/" + iFile
                    nMem = checkMem(fname)
                    if nMem == 0:
                        continue

                    # report what we get
                    L_list = getLCode(iFile)
                    if len(L_list) == 1:
                        print "%d  %d" % (L_list[0], nMem)
                    elif len(L_list) == 2:
                        print "%d  %d  %d" % (L_list[0], L_list[1], nMem)
                    elif len(L_list) == 3:
                        print "%d  %d  %d  %d" % (L_list[0], L_list[1], L_list[2], nMem)
                    elif len(L_list) == 4:
                        print "%d  %d  %d  %d  %d" % (L_list[0], L_list[1], L_list[2], L_list[3], nMem)
                    else:
                        print "something error after we getLCode"
                        sys.exit()


        # if we did not get the work dir, report error
        if not gotIt:
            print workDir
            print "We did not find anything matching the work dir"
            sys.exit()


