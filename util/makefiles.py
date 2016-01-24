"""
This is used to set up the jobs for generating makefiles
or CMakeList.txt
"""
import re
import sys
import os
__author__  = "Fenglai Liu"

# functional related infor
makefile_choice = "makefile"  # makefile is to generate makefile for debug, cmake is to generate CMakeList.txt
job             = " "         # what kind of integral module we would like to include?
maxL            = 5           # maximum L for integral
auxMaxL         = 5           # maximum L for the aux basis sett dimension of integrals
file_size       = 5242880     # this is the default file size that we may keep a caution (5M)
                              # file larger than this size should go with lower compiler flags
topDir          = " "         # this is the source code dir


#########################################################
#  initialize the information from setting file
#########################################################
def set_project():

    global makefile_choice
    global job
    global maxL
    global auxMaxL
    global file_size
    global topDir

    # initialize the value to be default
    job = " "
    makefile_choice = "makefile"
    maxL = 5
    auxMaxL = 5
    file_size = 5242880
    topDir = " "

    # for sys.argv length is 0, print out help
    if len(sys.argv) == 1:
        print "help on the command line parameters, format is key=value:\n"
        print "topDir(key): this is the source code location you need to give\n"
        print "maxL(key): set maximum L for cpp files that taking into makefile etc."
        print "maxL the default value is 5"
        print "maxL is used to characterizes the normal angular momentum for cpp file,"
        print "if the cpp file has angular momentum larger than maxL, the cpp file"
        print "are not included into the compilation systems\n"
        print "auxMaxL(key): set aux maximum L for cpp files that taking into makefile etc."
        print "aux maximum L should be same or larger than the maximum L value"
        print "this is used for three body ERI or two body ERI etc."
        print "it's default value is 5\n"
        print "file_size(key): this means any cpp file larger than 5M will be handled with"
        print "low optimization flag as you assigned"
        print "it's default value is 5242880, denotes 5M(1024*1024*5)\n"
        print "job: we generate makefiles etc. for the given job like two body overlap integrals\n"
        print "choice(key): we generate CMakeLists.txt for cmake compilation rather than makefile"
        print "if you want makefile, set makefile instead of cmake, makefile is default choice\n"
        sys.exit()

    # now let's parsing the parameters
    for arg in sys.argv:

        # the first argument is the program name
        if arg.find("py") > 0:
            continue

        # see whether we have the "=", this is necessary
        if arg.find("=") < 0:
            print arg
            print "this is a wrong augument, we did not find = inside"
            print "= connects the left side keyword and right side value"
            print "you should have it in your command line arguments"
            sys.exit()
        
        # now get the key and value
        arg_list = arg.split("=")
        if len(arg_list) == 2:
            oriKey = arg_list[0]
            oriVal = arg_list[1]
        else:
            print arg
            print "this is a wrong augument, we expect key=value form of arguments"
            sys.exit()

        # now transform key
        key = oriKey.lower()
        val = oriVal.lower()

        # read in things
        if key == "maxl":
            maxL = int(val)
            if maxL <= 0:
                print line
                print "Illegal value for setting maxL"
                sys.exit()
        elif key == "auxmaxl":
            auxMaxL = int(val)
            if auxMaxL <= 0:
                print line
                print "Illegal value for setting auxMaxL"
                sys.exit()
        elif key == "job":
            job = val
        elif key == "topdir":
            topDir = val
        elif key == "choice":
            if val != "makefile" and val != "cmake":
                print "Illegal choice for generating the makefiles/cmakelist.txt"
                sys.exit()
            makefile_choice = val
        elif key == "file_size":
            size = int(val)
            if size <= 0:
                print line
                print "Illegal value for setting file size"
                sys.exit()
            file_size = size
        else:
            print key
            print "Illegal key in argument, can not be recognized!"
            sys.exit()

    # final check
    if auxMaxL< maxL:
        print "maxL should be same or larger than the aux max L"
        sys.exit()

    # job need to have something
    print "job is:", job
    if job.isspace():
        print "job is not defined"
        sys.exit()

    # whether we have the top dir defined
    if topDir.isspace():
        print "topDir is not defined"
        sys.exit()
    if topDir[-1] == "/":
        tmp = topDir[0:-1]
        topDir = tmp

def doMakefile():
    global makefile_choice
    if makefile_choice == "makefile":
        return True
    else:
        return False

def getMaxL():
    global maxL
    return maxL

def getAuxMaxL():
    global auxMaxL
    return auxMaxL

def largeFile(name):
    global file_size

    # get the size of this file
    size = os.path.getsize(name)
    if size > file_size:
        return True
    else:
        return False

#########################################################
#  generate the Makefile
#########################################################
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
        # we just return -1
        return -1

def checkFile(f):
    """
    check the input file to see whether we can list it
    """
    # get the requirement
    auxL = getAuxMaxL()
    maxL = getMaxL()

    # we need to drop the extension
    name = os.path.splitext(f)[0]
    ext  = os.path.splitext(f)[1]
    if ext != ".cpp":
        return False

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
                if LBra1 <= maxL and LBra2<= maxL and LKet1 <= auxL:
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
            if onLargeFile and largeFile(f):
                cppList.append(i)

            if not onLargeFile and not largeFile(f):
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

    # append the deriv information
    if folder.find("first_deriv") >=0:
        name = name + "_d1"
    elif folder.find("second_deriv") >=0:
        name = name + "_d2"
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
        print "This is just a warning message, we will continue"

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

    # firstly, print out the entry file if we have one
    if getEntryFile:
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

##############################
# now let's generate the stuff
##############################
# parse the information from the input parameters
set_project()

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
                if iProj == job:
                    doMakefile(workDir)
