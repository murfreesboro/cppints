"""
This is used to generate the head files as well as the
global cpp function for the cppints
"""
__author__ = "Fenglai Liu"
__date__ = "Dec, 2013"
import sys
import os
import re

def printCode(f, line, indentLength=0):
    """
    indenting the code and print the given line
    """
    if indentLength > 0:
        for i in range(indentLength):
            f.write(" ")
    f.write(line)
    f.write("\n")

def printEmptyLine(f):
    """
    print out empty lines
    """
    f.write("\n")

def getL(s, withLMax=False):
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
		  if withLMax:
			  return 1
		  else:
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

def sigCheck(workDir):
    """
    for the given type of work (specified by workdir),
    do we really do significance check
    this should be in accordance with cppints code
    """
    if workDir.find("eri") >= 0:
        return True
    else:
        return False

def getLCode(filename):
    """
    get the angular momentum code from the given file name
    """

    # drop the extension first
    flist = os.path.splitext(filename)
    name = flist[0]
    nameComponents = name.split("_")

    # now let's get LCode
    L_list = [ ]
    for i in nameComponents:
        lcode = getL(i)
        if lcode >= 0:
            L_list.append(lcode)

    # now form the LCode
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
        code = code + n * L_list[i]

    # now return the code
    return code

def getProjectName(workDir):
    """
    from the input dir name, we get the project name
    """

    # we only need the last three sub dir name
    dirList = workDir.split("/")
    proj = dirList[-1]
    order = dirList[-2]
    method = dirList[-3]

    # now get the deriv order
    jobOrder = 0
    if order == "first_deriv":
        jobOrder = 1
    if order == "second_deriv":
        jobOrder = 2

    # finally get the project name
    if jobOrder == 0:
        project = method + "_" + proj
    elif jobOrder == 1:
        project = method + "_" + proj + "_d1"
    elif jobOrder == 2:
        project = method + "_" + proj + "_d2"
    return project

def getArgumentList(fname):
    """
    get the function argument list by reading it from the file
    """
    f = open(fname, "r")
    while True:
        line = f.readline()
        if not line: break
        if re.search(r"(?i)void", line) is not None or re.search(r"(?i)bool", line) is not None:

            # double check
            # for file split situation, the vrr and hrr functions
            # will be declared in front of the file
            if re.search(r"(?i)_vrr", line) is not None:
                continue
            if re.search(r"(?i)_hrr", line) is not None:
                continue

            begin = line.index("(")
            end = line.index(")")
            arg = line[begin + 1:end]
            args = arg.strip()
            return args
    f.close()

    # finally, we need to watch out the case that
    # we did not find anything inside
    print "failed to search the void in function name, wrong in gettArgumentList\n"
    sys.exit()

def getFuncArgumentList(arglist):
    """
    get the function argument list without type information
    """
    args = arglist.split()
    arg = ""
    for i in args:
        if re.search(r"(?i)Double", i) is not None:
            continue
        if re.search(r"(?i)const", i) is not None:
            continue
        if re.search(r"(?i)UInt", i) is not None:
            continue
        if re.search(r"(?i)Int", i) is not None:
            continue
        if re.search(r"(?i)\*", i) is not None:
            continue
        if re.search(r"(?i)\&", i) is not None:
            continue
        arg = arg + i
    return arg

def useLocalMemScrObject(workDir, maxL, maxAuxL):
    """
    whether the cpp files in this project uses the localmemscr object
    we can do it by searching the argument list from the function declaration
    """
    cppfiles = os.listdir(workDir)
    for i in cppfiles:

        # we will omit the cpp files with _vrr and _hrr
        if i.find("_vrr") >= 0:
            continue
        if i.find("_hrr") >= 0:
            continue

        # we only match the cpp file
	extension = os.path.splitext(i)[1]
        if extension != ".cpp":
            continue

        # do we use this cpp file
        if not useThisCPPFile(i, maxL, maxAuxL):
            continue

        # read in the argument list from the cpp file
        workFile = workDir + "/" + i
        arglist = getArgumentList(workFile)
        if re.search(r"(?i)LocalMemScr", arglist) is not None:
            return True

    # finally, we get nothing
    return False

def useThisCPPFile(name, maxL, auxL):
    """
    by parsing the name of the file, we will see that whether
    it needs to satisfy the requirement for maxL and auxL
    """
    # currently we only concentrate on ERI
    if re.search(r"(?i)eri", name) is not None:

        # drop the extension first
        f = os.path.splitext(name)[0]
        nameComponents = f.split("_")

        # now let's get L
        L_list = [ ]
        withLMax = True
        for i in nameComponents:
            lcode = getL(i, withLMax)
            if lcode >= 0:
                L_list.append(lcode)
        if len(L_list) != 4:
            print name
            print "in useThisCPPFile if it's eri file name, L list should be 4 data field"
            sys.exit()

	# firstly, whether it's three body ERI?
	if L_list[3] == 0 and L_list[1] != 0:
	    L = L_list[2]
	    L1 = L_list[0]
	    L2 = L_list[1]
	    if L <= auxL and L1 <= maxL and L2 <= maxL:
		return True

	# now let's see whether it's two body ERI?
	elif L_list[3] == 0 and L_list[1] == 0:
	    L1 = L_list[0]
            L2 = L_list[2]
	    if L1 <= auxL and L2 <= auxL:
		return True

	# now all of normal four body ERI
	else:
	    L1 = L_list[0]
	    L2 = L_list[1]
	    L3 = L_list[2]
	    L4 = L_list[3]
	if L1 <= maxL and L2 <= maxL and L3 <= maxL and L4 <= maxL:
	    return True

        # now we know that we need to ignore this file
	return False

    # all of other situations
    else:
	return True

def filesCreation(workDir, maxL, maxAuxL):
    """
    print out the head files as well as cpp entry files to call the function
    """

    # get the parent working dir
    fileDir = os.path.split(os.path.abspath(workDir))[0]

    # firstly, we need to search the whole folder to see that
    # whether we have localmemscr class used in the cpp files
    hasLocalMemScr = useLocalMemScrObject(workDir, maxL, maxAuxL)

    # get the project name from the dir
    name = getProjectName(workDir)
    incfile = fileDir + "/" + name + ".h"
    head = open(incfile, "w")

    # top macros
    macro = name.upper()
    line = "#ifndef  " + macro + "_H"
    printCode(head, line)
    line = "#define " + macro + "_H"
    printCode(head, line)
    printEmptyLine(head)

    # head file, we have to make up the name
    line = "#include <cstddef>"
    printCode(head, line)
    line = "#include <cassert>"
    printCode(head, line)
    line = "#include <cstdio>"
    printCode(head, line)
    if hasLocalMemScr:
        line = "#include \"localmemscr.h\""
        printCode(head, line)
        line = "using namespace localmemscr;"
        printCode(head, line)
    printEmptyLine(head)

    # typedef information
    # only use it when we do not include localmemscr
    if not hasLocalMemScr:
        line = "typedef int             Int; "
        printCode(head, line)
        line = "typedef long long       LInt; "
        printCode(head, line)
        line = "typedef size_t          UInt; "
        printCode(head, line)
        line = "#ifdef WITH_SINGLE_PRECISION "
        printCode(head, line)
        line = "typedef float           Double;"
        printCode(head, line)
        line = "#else"
        printCode(head, line)
        line = "typedef double          Double;"
        printCode(head, line)
        line = "#endif"
        printCode(head, line)
    printEmptyLine(head)

    # before we create the entry file, firstly let's get
    # the entry function argument list
    # we note, that all of cpp functions share the same
    # function list, except that some of them may have
    # additional localmemscr objeect pass in
    # now let's do it
    entryFunArgList = " "
    cppfiles = os.listdir(workDir)
    for i in cppfiles:

        # we will omit the cpp files with _vrr and _hrr
        if i.find("_vrr") >= 0:
            continue
        if i.find("_hrr") >= 0:
            continue

        # we only match the cpp file
        extension = os.path.splitext(i)[1]
        if extension != ".cpp":
            continue

		# do we use this cpp file
        if not useThisCPPFile(i, maxL, maxAuxL):
            continue

        # read in the argument list from the cpp file
        workFile = workDir + "/" + i
        entryFunArgList = getArgumentList(workFile)
        break

    # set up the entry file - cpp file
    cppfile = fileDir + "/" + name + ".cpp"
    cpp = open(cppfile, "w")

    # include header file into cpp file
    inc = name + ".h"
    line = "#include " + "\"" + inc + "\""
    printCode(cpp, line)
    printEmptyLine(cpp)

    # generating the function names
    # then we are done with this
    # we note, that for each of the function it may have the namespace
    # around it
    # you can choose what you want, modify here
    entryFunArgList = "const LInt& LCode, " + entryFunArgList
    if hasLocalMemScr:
        if re.search(r"(?i)LocalMemScr", entryFunArgList) is None:
            entryFunArgList = entryFunArgList + ", LocalMemScr& scr"
    returnType = "void "
    funcname = returnType + name + "(" + entryFunArgList + ")"
    printCode(cpp, funcname)

    # begin to print the body of entry cpp file
    line = "{"
    printCode(cpp, line)
    printEmptyLine(cpp)

    # step into the main body,waiting for the further printing
    line = "switch (LCode) {"
    printCode(cpp, line, 2)

    # now loop over the dir to get each file name
    cppfiles = os.listdir(workDir)
    for i in cppfiles:

        # we will omit the cpp files with _vrr and _hrr
        if i.find("_vrr") >= 0:
            continue
        if i.find("_hrr") >= 0:
            continue

        # we only match the cpp file
        extension = os.path.splitext(i)[1]
        if extension != ".cpp":
            continue

	# do we use this cpp file
        if not useThisCPPFile(i, maxL, maxAuxL):
            continue

        # read in the argument list from the cpp file
        workFile = workDir + "/" + i
        arglist = getArgumentList(workFile)

        # drop the extension
        # print it in the head file
        flist = os.path.splitext(i)
        fname = flist[0]
        line = returnType + fname + "(" + arglist + ");"
        printCode(head, line)
        printEmptyLine(head)

        # now print the entry cpp file
        # we get the LCode
        LCode = getLCode(fname)
        line = "case " + str(LCode) + ":"
        printCode(cpp, line, 4)

        # now it's function call
        funcArgList = getFuncArgumentList(arglist)
        if sigCheck(workDir):
            line = fname + "(" + funcArgList + ");"
        else:
            line = fname + "(" + funcArgList + ");"
        printCode(cpp, line, 6)
        line = "break;"
        printCode(cpp, line, 6)

    # finish the head file printingg
    line = "#endif"
    printCode(head, line)
    head.close()

    # finish the switch code
    line = "default:"
    printCode(cpp, line, 4)
    line = "#ifdef DEBUG"
    printCode(cpp, line, 0)
    line = "printf(\"%s %lld\\n\","
    line = line + "\"Un-recognized LCode in the integrals calculation \", "
    line = line + "LCode);"
    printCode(cpp, line, 6)
    line = "assert(0);"
    printCode(cpp, line, 6)
    line = "#endif"
    printCode(cpp, line, 0)
    line = "break;"
    printCode(cpp, line, 6)
    line = "}"
    printCode(cpp, line, 2)
    printEmptyLine(cpp)

    # finish the whole function code
    line = "}"
    printCode(cpp, line)
    cpp.close()




