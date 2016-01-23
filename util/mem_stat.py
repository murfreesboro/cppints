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
    get the angular momentum list from the given file name
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

def formLCode(L_list):
    """
    get the angular momentum code from the L list
    """
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

def sameCPPCodes(filename,code):
    """
    search for the working files (not the top driver files),
    whether the given file is same with the L Code
    """
    L_list  = getLCode(filename)
    newcode = formLCode(L_list)
    if newcode == code:
        return True
    return False

def isMainCPPFile(fname):
    """
    check whether the given file is main cpp file
    """
    flist = os.path.splitext(fname)
    name = flist[0]
    nameComponents = name.split("_")

    # we only match the cpp file
    extension = flist[1]
    if extension != ".cpp":
        return False

    # whether this is the derivatives case
    if nameComponents[-1] == "d1" or nameComponents[-1] == "d2":
        return True

    # if not the derivatives, the last one should
    # be the angular momentum
    ang = nameComponents[-1]
    lcode = getL(ang)
    if lcode >= 0:
        return True

    # finally return false
    return False

def checkMem(filename,workDir):
    """
    for the given file name, check the memory usage
    """
    nMem = 0

    # check whether this is the main cpp file?
    if not isMainCPPFile(filename):
        return -1

    # get the code of this cpp work 
    L_list  = getLCode(filename)
    code = formLCode(L_list)
    #print "the main cpp file", filename

    # let's go through the main file to have a look
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
    #print "nmem is", nMem

    # now let's go to it's sub working files to have a look
    files = os.listdir(workDir)
    for iFile in files:

        # we only concentrate on cpp file
        if iFile.find("cpp") < 0:
            continue

        # omit all of main cpp files
        if isMainCPPFile(iFile):
            continue

        # check whether they give the same L codes
        if not sameCPPCodes(iFile,code):
            continue

        # now let's count the memory usage
        fname = workDir + "/" + iFile
        #print fname
        cpp = open(fname,"r")
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

def getWorkPath(topDir,order):
    """
    return the work path through top dir
    """
    projectDir = topDir
    if order == "0":
        projectDir = topDir + "/" + "energy"
    elif order == "1":
        projectDir = topDir + "/" + "first_deriv"
    elif order == "2":
        projectDir = topDir + "/" + "second_deriv"
    else:
        print "illegal order here, we only support order from 0 to 2"
        sys.exit()
    return projectDir

def checkJob(topDir,job,order):
    """
    check the input job, to see that whether
    they are valid
    also return the information about maxL
    and auxMaxL
    """
    # firstly let's check the job
    # also define the work dir
    workDir = job
    if job == "eri" or job == "2eri" or job == "3eri":
        workDir = "eri"
    elif job == "twobodyoverlap":
        workDir = "twobodyoverlap"
    elif job == "threebodyoverlap":
        workDir = "threebodyoverlap"
    elif job == "kinetic":
        workDir = "kinetic"
    elif job == "esp":
        workDir = "esp"
    elif job == "nai":
        workDir = "nai"
    elif job == "mom":
        workDir = "mom"
    else:
        print job
        print "illegal job here, please check the function of checkJob to see th supported operator list"
        sys.exit()

    # set the maxL and auxL
    # in default they are 0, since they are not used
    # however, for some given jobs, they need to be 
    # examined
    maxL = 0
    auxL = 0
    if job == "eri" or job == "2eri" or job == "3eri":

        # we need to loop over the dir
        projectDir = getWorkPath(topDir,order)
        wDir = projectDir + "/eri" 
        if not os.path.isdir(wDir):
            print "you set the eri job, however there are no source code?"
            sys.exit()

        # now let's loop over the files to see L information
        files = os.listdir(wDir)
        for iFile in files:

            # we only concentrate on main cpp file
            if not isMainCPPFile(iFile):
                continue

            # get the L list
            L_list = getLCode(iFile)

            # now see the job
            # for all of eri job
            # we just check the ket2 position
            # this is for setting maxL
            # because ket2 is never referred by 2eri and 3eri
            #
            # for 2eri/3eri, we need to check KET1
            # this is the direction of aux basis set
            L = L_list[3]

            # whether this is SP shell?
            # then we need to reset it
            ang = L
            if L == 100:
                ang = 1

            # now get value
            if ang>maxL:
                maxL = ang

            # for aux basis set
            if job == "3eri" or job == "2eri":
                L = L_list[2]
                ang = L
                if L == 100:
                    ang = 1
                if ang>auxL:
                    auxL = ang

    # now let's return
    print "maxL is", maxL
    print "auxMaxL is", auxL
    return (maxL,auxL,workDir)

def checkERIJob(job,maxL,auxL,filename):
    """
    for the given ERI job, through the information
    of maxL and auxL, we check that whether we 
    should include the given file
    """

    # we just test ERI job
    if job == "eri" or job == "2eri" or job == "3eri":

        # get the L information
        # we also need to correct it
        L = getLCode(filename)
        for i in range(len(L)):
            ang = L[i]
            if ang == 100:
                L[i] = 1

        # now let's see
        if job == "eri":
            if L[0] > maxL or L[1] > maxL or L[2] > maxL or L[3] > maxL:
                return False
        elif job == "2eri":
            if L[0] > auxL or L[1] > 0    or L[2] > auxL or L[3] > 0:
                return False
        elif job == "3eri":
            if L[0] > maxL or L[1] > maxL or L[2] > auxL or L[3] > 0:
                return False

    # now everything is good
    # let's go on
    return True

#################################################
# now this is main script for memory statistics #
#################################################
order = 0
if len(sys.argv) == 4:
    topDir  = sys.argv[1]
    job     = sys.argv[2]
    order   = sys.argv[3]
else:
    print "We need three arguments, one is the input cpp files dir"
    print "This one is necessary, and does not have default value"
    print "The second one is the work dir name, like twobodyoverlap"
    print "The third one is the deriv order, order = 0 is for energy"
    print "order = 1 for derivatives, order = 2 for second derivatives"
    sys.exit()

# check the last symbol of topDir
# with "/" will cause trouble
if topDir[-1] == "/":
    tmp = topDir[0:-1]
    topDir = tmp

# further check the job
# for some job, for example; the 
# eri, 2eri and 3eri job; we need to 
# search the maxL and auxMaxL
maxL    = 0
auxMaxL = 0
wDir    = job
(maxL, auxMaxL, wDir) = checkJob(topDir,job,order)

# on the top of dir, we should have 
# three sub-dir whose name like this
projectDir = getWorkPath(topDir,order)
suffix = "_0"
if order == "1":
    suffix = "_1"
elif order == "2":
    suffix = "_2"

# we should create file to store the result
out = "mem_infor_" + job + suffix + ".txt"

# now the real work
projects = os.listdir(projectDir)
gotIt = False
for iProj in projects:
    workDir = projectDir + "/" + iProj
    if os.path.isdir(workDir) and wDir == iProj:
        gotIt = True
        break

# now go over the cpp files
if gotIt:
    output = open(out,"w")
    files = os.listdir(workDir)
    for iFile in files:

        # we only concentrate on cpp file
        if iFile.find("cpp") < 0:
            continue

        # for ERI job, we need a double check
        if not checkERIJob(job,maxL,auxMaxL,iFile):
            continue

        # now do the real job
        fname = workDir + "/" + iFile
        nMem = checkMem(fname,workDir)
        if nMem < 0:
            continue

        # report what we get
        L_list = getLCode(iFile)
        code   = formLCode(L_list)
        output.write( "%-12d  %-12d " % (code, nMem))
        if len(L_list) == 1:
            output.write( " %-3d\n" % (L_list[0]))
        elif len(L_list) == 2:
            output.write( " %-3d  %-3d\n" % (L_list[0], L_list[1]))
        elif len(L_list) == 3:
            output.write( " %-3d  %-3d  %-3d\n" % (L_list[0], L_list[1], L_list[2]))
        elif len(L_list) == 4:
            output.write( " %-3d  %-3d  %-3d  %-3d\n" % (L_list[0], L_list[1], L_list[2], L_list[3]))
        else:
            print "something error after we getLCode"
            sys.exit()
    output.close()
else:
    print workDir
    print "We did not find anything matching the work dir"
    sys.exit()


