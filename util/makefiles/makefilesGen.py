"""
This is used to generate the makefiles etc. for integral module
"""
__author__  = "Fenglai Liu"
import infor
import sys
import os
import re

def printCompilation(f,cppfile,inParentDir=False,withLowFlag=False):
    """
    for the given cpp file, print out the line for compilation
    """

    # generate the full path name of cpp file etc.
    name = os.path.splitext(cppfile)[0]

    # form the binary part
    binary = "${BIN}/" + name + ".o"
    if inParentDir:
        cpp    = "${PARSRC}/" + name + ".cpp"
    else:
        cpp    = "${SRC}/" + name + ".cpp"

    # first line
    line = binary + ":  " + cpp
    f.write(line)
    f.write("\n")

    # second line
    f.write("\t\t")
    if withLowFlag:
        line = "$(CC) $(LOWCFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@"
    else:
        line = "$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@"
    f.write(line)
    f.write("\n")

def getL(s):
    """
    whether the given s is shell symbol
    if it is, we return it's Lmax (not LCode)
    else print error message
    """
    if s == "S" or s == "s":
        return 0
    elif s == "P" or s == "p":
        return 1
    elif s == "SP" or s == "sp":
        return 1
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
        print "wrong angular momentum in the cpp integral file name"
        print "can not identify it"
        print s
        sys.exit()

def checkFile(f):
    """
    check the input file to see whether we can list it 
    """
    # get the requirement
    auxL = infor.getAuxMaxL()
    maxL = infor.getMaxL()

    # we need to drop the extension
    name = os.path.splitext(f)[0]

    # split the file names
    infors = name.split("_")

    # get the L information
    operator = infors[2]
    LBra1    = getL(infors[3])
    LBra2    = -1
    if len(infors) >= 5:
        LBra2  = getL(infors[4])
    LKet1 = -1
    if len(infors) >= 6:
        LKet1  = getL(infors[5])
    LKet2 = -1
    if len(infors) >= 7:
        LKet2  = getL(infors[6])

    # check with maxL
    if LBra1<= maxL and LBra2<= maxL and LKet1<= maxL and LKet2<= maxL:
        return True
    else:

        # does this file has auxL?
        oper = operator.lower()
        if oper == "eri":

            # (as|bs)
            if LBra2 == 0 and LKet2 == 0:
                if LBra1 <= auxL and LKet1 <= auxL:
                    return True
                else:
                    return False

            # (ab|cs)
            if LBra2 > 0 and LKet2 == 0:
                if LBra1 <= maxL and Lbra2<= maxL and LKet1 <= auxL:
                    return True
                else:
                    return False

    # print out error message
    print f
    print "something wrong inside the function of checkCPPFiles"
    sys.exit()

def grabCPPFiles(workDir,onLargeFile):
    """
    grab all of cpp files that satisfy the requirement in the 
    infor 
    """
    cppList = [ ]
    cppfiles = os.listdir(workDir)
    for i in cppfiles:

        # is this file what we need?
        if checkFile(i):
           
            # if this is on large file check
            # we only concentrate on large files
            # else we only do small files
             f = workDir + "/" + i
             if onLargeFile and  infor.largeFile(f):
                 cppList.append(i)

             if not onLargeFile and not infor.largeFile(f):
                 cppList.append(i)


    return cppList

def formOBJ(cppList):
    """
    change the cpp file into obj
    """
    objList = [ ]
    for i in cppList:
        name = os.path.splitext(i)[0]
        binary = "${BIN}/" + name + ".o"
        objList.append(binary)
    return objList

def doMakefile(folder):
    """
    generate makefile for the given folder
    """

    # makefile name
    oper = os.path.basename(folder)
    name = "Makefile_" + oper
    name = open(name, "w")

    # first part, comments
    line = "# \n"
    name.write(line)
    line = "# makefile for compiling the module of " + oper + "\n"
    name.write(line)
    line = "# \n"
    name.write(line)
    name.write("\n")
    name.write("\n")

    # now include file
    line = "# include compiler flags\n"
    name.write(line)
    line = "include compile.cfg\n"
    name.write(line)
    name.write("\n")

    # define the src
    line = "# define source code path \n"
    name.write(line)
    p = os.path.abspath(folder)
    parentDir = os.path.abspath(os.path.join(p, os.pardir))
    line    = "PARSRC  = " + parentDir + "\n"
    name.write(line)
    line    = "SRC  = " + "${PARSRC}/" + oper + "\n"
    name.write(line)
    name.write("\n")

    # get the entry file
    line = "# define objs \n"
    name.write(line)
    cppfiles = os.listdir(parentDir)
    getEntryFile = False
    for i in cppfiles:
        if os.path.splitext(i)[1] == ".cpp":
            cpp = i.upper()
            if cpp.find(oper.upper()) >= 0:
                entryFile = os.path.splitext(os.path.basename(i))[0]
                binary = "${BIN}/" + entryFile + ".o"
                line = "OBJ  =  " + binary + "\n"
                name.write(line)
                getEntryFile = True
    if not getEntryFile:
        print "We did not get entry cpp file for module " + folder
        print "Something wrong with it"
        sys.exit()

    # now list all of OBJS, firstly for normal objs
    cppFiles = grabCPPFiles(folder,False)
    objFiles = formOBJ(cppFiles)
    for obj in objFiles:
        line = "OBJ  +=  " + obj + "\n"
        name.write(line)
    name.write("\n")

    # now list all of lower OBJS, if it has
    line = "# define large objs - they need lower compiler opt option\n"
    name.write(line)
    hasLowOBJ = False
    cppFiles = grabCPPFiles(folder,True)
    if len(cppFiles) > 0:
        hasLowOBJ = True
        objFiles = formOBJ(cppFiles)
        for obj in objFiles:
            if obj == objFiles[0]:
                line = "LOWOBJ   =  " + obj + "\n"
                name.write(line)
            else:
                line = "LOWOBJ  +=  " + obj + "\n"
                name.write(line)
        name.write("\n")

    # now step into the real code
    line = "# create lib file \n"
    name.write(line)
    oper = os.path.basename(folder).upper()
    tarName  = "LIB" + oper
    line = "ALL: " + "${OBJ}"
    if hasLowOBJ:
        line = line + " ${LOWOBJ}"
    line = line + "\n"
    name.write(line)
    line = "\t\t"
    name.write(line)
    line = "$(AR) -r ${" + tarName + "}"  + " ${OBJ}"
    if hasLowOBJ:
        line = line + " ${LOWOBJ}"
    line = line + "\n"
    name.write(line)
    name.write("\n")

    # firstly, print out the entry file
    line = "# compile all of objs - entry cpp file first\n"
    name.write(line)
    InParentDir = True
    entryFile = entryFile + ".cpp"
    printCompilation(name,entryFile,InParentDir)
    name.write("\n")
    name.write("\n")

    # now for normal obj
    line = "# these are normal objs whose size are small\n"
    name.write(line)
    cppFiles = grabCPPFiles(folder,False)
    for cpp in cppFiles:
        printCompilation(name,cpp)
    name.write("\n")
    name.write("\n")

    # now for large objs
    InParentDir = False
    withLowCompilerFlag = True
    line = "# these are normal objs whose size are large\n"
    name.write(line)
    if hasLowOBJ:
        cppFiles = grabCPPFiles(folder,True)
        for cpp in cppFiles:
            printCompilation(name,cpp,InParentDir,withLowCompilerFlag)
        name.write("\n")
        name.write("\n")
        
    # final ends
    name.write("\n")
    name.write("\n")
    name.write("\n")

    name.close()
