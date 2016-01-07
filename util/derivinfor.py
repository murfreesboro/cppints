#!/usr/bin/python
"""
this file is used to extract the derivatives information
which contained in the comments part of cpp files
for the given folder
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

def getFileName(filename):
    """
    get the main file name from the whole dir+filename
    """

    # drop the extension first
    fname = os.path.basename(filename)
    name  = os.path.splitext(fname)[0]
    return name

def isMainFile(filename,jobOrder):
    """
    check that whether the given file is main cpp file
    """
    # drop the extension first
    fname = os.path.basename(filename)
    name  = os.path.splitext(fname)[0]
    nameComponents = name.split("_")

    # whether this is cpp file
    if filename.find("cpp") < 0:
        return False

    # get the last components
    component = nameComponents[-1]
    if jobOrder > 0:
        if component == "d1" or component == "d2":
            return True
        else:
            return False
    else:
        lcode = getL(component)
        if lcode >= 0:
            return True
        else:
            return False


def formLCode(filename):
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

    # now let's form the L code
    code = 0
    length = len(L_list)
    for i in range(length):
        n = 1
        if i == 1:
            n = 1000
        elif i == 2:
            n = 1000000
        elif i == 3:
            n = 1000000000
        code = code + n*L_list[i]

    # now return the code
    return code

#################################################
#         now this is main script               #
#################################################
if len(sys.argv) == 4:
    topDir  = sys.argv[1]
    module  = sys.argv[2]
    jobOrder= sys.argv[3]
else:
    print "We need threearguments:"
    print "1  the top cpp code directory name, like hgp_os"
    print "2  the module name, like esp, twobodyoverlap or eri etc."
    print "3  the derivative order, right now it's only 1 or 2"
    sys.exit()

# generate the output text file
if int(jobOrder) == 1:
    suffix = "_1"
elif int(jobOrder) == 2:
    suffix = "_2"
else:
    print "wrong derivative order given, only 1 or 2 is supported here"
    sys.exit()
out = "deriv_infor_" + module + suffix + ".txt"
print "the output of deriv information will be direct to file: %s" % out

# now let's create the work folder
if int(jobOrder) == 1:
    workDir = topDir + "/" + "first_deriv" + "/" + module
else:
    workDir = topDir + "/" + "second_deriv" + "/" + module
if not os.path.exists(workDir):
    print "working directory does not exist"
    print workDir
    sys.exit()

# now go over the cpp files
output = open(out,"w")
files = os.listdir(workDir)
for iFile in files:

    # whether this is main cpp file
    if not isMainFile(iFile,jobOrder):
        continue

    # get the file name
    name = getFileName(iFile)
    LCode = formLCode(name)
    output.write("\n")
    output.write(str(LCode))
    output.write("\n")

    # now prepare to get the diratives information
    f = workDir + "/" + iFile
    cpp = open(f,"r")
    beginToRead = False
    while True:

        # target the place we need to read in
        line = cpp.readline()
        if not line: break
        if re.search(r"(?i)@@@@", line) is not None:
            beginToRead = True
            continue

        # check whether we should stop reading
        if re.search(r"(?i)//", line) is None and beginToRead:
            break

        # now if begin to read, let's prepare
        if beginToRead:

            # split the line information
            # we get rid of the // sign
            # and print out the rest of
            # the stuff
            line2 = line.replace('//',' ')
            line3 = line2.strip()
            output.write(line3)
            output.write("\n")
    cpp.close()

output.close()
