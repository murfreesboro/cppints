"""
this file is used to evaluate the efficiency of the
cpp code generated by the cppints. Including the
memory usage and FLOPS operations etc.
"""
__author__  = "Fenglai Liu"
import sys
import os
import re

#
#  data section
#  all of the data will be re-evaluated for each cpp file parsing
#  we note, that derivatives calculation part, non-RR part are also
#  count into the hrr part evaluation
#
hrr_var_num  = 0        # the counting of local HRR/non-RR/deriv section memory use (only in variable form)
hrr_vec_nvar = 0        # the counting of local HRR/non-RR/deriv section memory use (only in vector)
vrr_var_num  = 0        # the counting of local VRR memory use (only in variable form)
vrr_vec_nvar = 0        # the counting of local VRR memory use (only in vector)
n_ignore_vrr = 0        # number of ignore integrals in VRR part
n_ignore_hrr = 0        # number of ignore integrals in HRR/non-RR/deriv section part
n_add_k2     = 0        # number of adding/substracting operations for K2 loop
n_mul_k2     = 0        # number of multiplication/dividing operations for K2 loop
n_add_k4     = 0        # number of adding/substracting operations for K4 loop
n_mul_k4     = 0        # number of multiplication/dividing operations for K4 loop
n_add_hrr    = 0        # number of adding/substracting operations for HRR/non-RR/deriv sections
n_mul_hrr    = 0        # number of multiplication/dividing operations for HRR/non-RR/deriv sections

def initilize():
    """
    initilize the variables for the running
    """
    global hrr_var_num
    global vrr_var_num
    global hrr_vec_nvar
    global vrr_vec_nvar
    global n_add_k2
    global n_mul_k2
    global n_add_k4
    global n_mul_k4
    global n_add_hrr
    global n_mul_hrr
    global n_ignore_vrr
    global n_ignore_hrr

    hrr_var_num  = 0
    vrr_var_num  = 0
    hrr_vec_nvar = 0
    vrr_vec_nvar = 0
    n_add_k2     = 0
    n_mul_k2     = 0
    n_add_k4     = 0
    n_mul_k4     = 0
    n_add_hrr    = 0
    n_mul_hrr    = 0
    n_ignore_hrr = 0
    n_ignore_vrr = 0

def goto_beginning(f):
    """
    for the given file handler, we go to
    the begining of function
    """
    # rewind the file to the beginning
    f.seek(0)

    # now go to the beginning of function
    while True:
        line = f.readline()
        newline = line.strip()
        if re.search(r"^{", newline) is not None:
            break

def is_comment(line):
    """
    whether this line is comment
    """
    newline = line.strip()
    if re.search(r"^//", newline) is not None:
        return True
    if re.search(r"^\*", newline) is not None:
        return True
    if re.search(r"^/\*", newline) is not None:
        return True
    return False

def is_vrr_file(fname):
    """
    whether this is sub file which do VRR/VRR contraction work
    """
    flist = os.path.splitext(fname)
    name0 = flist[0]
    name1 = os.path.basename(name0)
    if re.search(r"vrr", name1) is not None:
        return True
    return False

def only_vrr(workDir,fname):
    """
    whether the file only contains the VRR part of code?
    we assume that the pass in file is the main file,
    we will check it; and if not we will report an error
    """

    # check whether it's main cpp file?
    if not isMainCPPFile(fname):
        print fname
        print "you should pass in main cpp file in only_vrr function"
        sys.exit()

    # if this is derivative file, we must have at least derivative code
    flist = os.path.splitext(fname)
    name = flist[0]
    nameComponents = name.split("_")
    if nameComponents[-1] == "d1" or nameComponents[-1] == "d2":
        return False

    # now check whether we have HRR
    has_HRR = False
    filename = workDir + "/" + fname
    f = open(filename, "r")
    while True:
        line = f.readline()
        if not line:
            break
        if re.search(r"initilize the HRR steps", line) is not None:
            has_HRR = True
            break
    f.close()
    if has_HRR:
        return False

    # if we reach here, that means it really just
    # has the VRR part
    return True

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

def formFileList(workDir,fname):
    """
    for the given file name(this should be the main cpp file),
    we collect the sub working file list, and finally make
    the main file as well as it's sub files together for 
    a list
    """
    fList = [ ]
    fList.append(fname)

    # get the LCode
    LCode = getLCode(fname)

    # now let's loop over all of file
    # and file those whose LCode is same
    # with the main file one
    cppfiles = os.listdir(workDir)
    for iFile in cppfiles:

        # we only concentrate on cpp file
        if iFile.find("cpp") < 0:
            continue

        # omit the main file, because we 
        # already count it
        if isMainCPPFile(iFile):
            continue

        # get code
        code = getLCode(iFile)
        if code == LCode:
            fList.append(iFile)
    
    # return the file list
    return fList

def hrr_analyze(workDir,filelist):
    """
    analyze the HRR/non-RR/derivatives part of codes
    """
    global hrr_var_num
    global hrr_vec_nvar
    global n_add_hrr
    global n_mul_hrr
    global n_ignore_hrr 

    # check the main file, to see whether only VRR
    # is contained
    mainCPPFile = filelist[0]
    if only_vrr(workDir,mainCPPFile):
        return 

    # now loop over all of non-VRR files
    # for the main file, we need to be careful
    # since we need to direct the code section to
    # the non-VRR places
    for iFile in filelist:

        # open the file
        workFile = workDir + "/" + iFile
        name = open(workFile, "r")
        goto_beginning(name)

        # whether this is the main cpp file?
        isMain = False
        if mainCPPFile == iFile:
            isMain = True

        # if this is VRR working file, continue
        if is_vrr_file(iFile):
            continue

        # for main cpp file, in terms of esp; the top loop
        # is over grid point, we need to pass that
        if isMain:
            if re.search(r"esp", iFile) is not None:
                while True:
                    line = name.readline()
                    newline = line.strip()
                    if re.search(r"^for", newline) is not None and re.search(r"nGrids", newline) is not None:
                        break
                    if not line:
                        print "in hrr_analyze, we did not get the the beginning of VRR for esp??\n"
                        sys.exit()

        # for the main cpp file, we need to bypass the VRR section
        if isMain:

            # read in possible VRR result delcare
            while True:
                line = name.readline()
                newline = line.strip()
                if re.search(r"^for", newline) is not None and re.search(r"{", newline) is not None:
                    break
                if is_comment(line) or line.isspace():
                    continue

                # if we reach the end of line, report an error
                if not line:
                    print "in hrr_analyze, we did not get the vrr beginning??\n"
                    sys.exit()

                # let's count how many variables defined here
                # the variables are defined in two ways:
                # one is with vector, else must be double var
                if re.search(r"DoubleVec", line) is not None:
                    begin = line.index("(")
                    end   = line.index(",")
                    num   = line[begin+1:end]
                    n     = int(num)
                    hrr_vec_nvar = hrr_vec_nvar + n
                elif re.search(r"getNewMemPos", line) is not None:
                    begin = line.index("(")
                    end   = line.index(")")
                    num   = line[begin+1:end]
                    n     = int(num)
                    hrr_vec_nvar = hrr_vec_nvar + n
                elif re.search(r"Double", line) is not None and re.search(r"=", line) is not None:
                    hrr_var_num = hrr_var_num + 1

            # here we need to stop by where VRR really begins
            while True:
                line = name.readline()

                # this is for VRR is in the main file
                if re.search(r"shell quartet name", line) is not None:
                    break

                # this is for VRR part is in a independent file
                if re.search(r"hgp_os", line) is not None and re.search(r"vrr", line) is not None:
                    break

                if not line:
                    print iFile
                    print "in hrr_analyze, failed to get to the part that VRR begins it's sq calculation\n"
                    sys.exit()

            # now let's search where the VRR ends
            while True:
                line = name.readline()
                if re.search(r"}", line) is not None:
                    break
                if not line:
                    print "in hrr_analyze, failed to get to the end of VRR\n"
                    sys.exit()

        # now let's begin to parse
        while True:
            line = name.readline()
            newline = line.strip()

            # omit the non-code lines
            if line.isspace():
                continue

            # jump out the end of file
            if not line: break

            # count the number of ignore integrals
            if is_comment(line): 
	        if re.search(r"integrals are omitted", newline) is not None: 
		    tmp = newline.split(" ")
		    pos = tmp.index("totally")
		    n = tmp[pos+1]
		    ignore = int(n)
                    n_ignore_hrr += ignore
            else:

                # whether we declare a vector?
                if re.search(r"DoubleVec", line) is not None:
                    begin = newline.index("(")
                    end   = newline.index(",")
                    num   = newline[begin+1:end]
                    n     = int(num)
                    hrr_vec_nvar = hrr_vec_nvar + n
                elif re.search(r"getNewMemPos", line) is not None:
                    begin = newline.index("(")
                    end   = newline.index(")")
                    num   = newline[begin+1:end]
                    n     = int(num)
                    hrr_vec_nvar = hrr_vec_nvar + n

                # whether this is a new double var?
                if re.search(r"^Double", newline) is not None:
                    hrr_var_num = hrr_var_num + 1

                # counting the number of adding/substracting
                n = newline.count('+')
                n_add_hrr = n_add_hrr + n
                n = newline.count('-')
                n_add_hrr = n_add_hrr + n

                # counting the number of multiplication/dividing
                n = newline.count('*')
                n_mul_hrr = n_mul_hrr + n
                n = newline.count('/')
                n_mul_hrr = n_mul_hrr + n

        # finally, close the file
        name.close()

def vrr_analyze(workDir,filelist):
    """
    analyze the VRR/VRR contraction part of codes
    """
    global vrr_var_num
    global vrr_vec_nvar
    global n_add_k2   
    global n_mul_k2  
    global n_add_k4 
    global n_mul_k4 
    global n_ignore_vrr 

    # now loop over all of VRR files
    in_k4_part = False
    for iFile in filelist:

        # whether this is the main cpp file?
        isMain = False
        if filelist.index(iFile) == 0:
            isMain = True

        # we only count on VRR files
        if not isMain: 
            if not is_vrr_file(iFile):
                continue

        # open the file
        workFile = workDir + "/" + iFile
        name = open(workFile, "r")
        goto_beginning(name)

        # for main cpp file, in terms of esp; the top loop
        # is over grid point, we need to pass that
        if isMain:
            if re.search(r"esp", iFile) is not None:
                while True:
                    line = name.readline()
                    newline = line.strip()
                    if re.search(r"^for", newline) is not None and re.search(r"nGrids", newline) is not None:
                        break
                    if not line:
                        print "in vrr_analyze, we did not get the the beginning of VRR for esp??\n"
                        sys.exit()

        # for the main cpp file, we need go to the VRR beginning
        if isMain:

            # read in possible VRR result delcare
            while True:
                line = name.readline()
                newline = line.strip()
                if is_comment(line) or line.isspace():
                    continue
                if re.search(r"^for", newline) is not None and re.search(r"{", newline) is not None:
                    break

                # if we reach the end of line, report an error
                if not line:
                    print "in vrr_analyze, we did not get the vrr beginning??\n"
                    sys.exit()

        # now let's begin to parse
        # be careful that we may have K4 part of loop, also
        # in calculating the fmt() part we also have single precision code
        # which we do not want to count
        inSinglePrecisionZone = False
        inCoreVRR = False
        while True:
            line = name.readline()
            newline = line.strip()

            # omit the non-code lines
            if line.isspace():
                continue

            # jump out the end of file
            if not line: break

            # when we stop reading?
            if isMain:
                # whether we are in the core of VRR? when VRR is not in an independent file
                if re.search(r"shell quartet name", newline) is not None:
                    inCoreVRR = True

                # VRR is in independent file
                if re.search(r"hgp_os", newline) is not None and re.search(r"vrr", newline) is not None:
                    inCoreVRR = True

                # if we already in the core of VRR, we can prepare for a stop
                if inCoreVRR:
                    if re.search(r"^}", newline) is not None:
                        break
            else:
                if re.search(r"^}", newline) is not None:
                    break

            # whether this is K4 part?
            if isMain and re.search(r"{", newline) is not None and \
                    re.search(r"^for", newline) is not None:
                in_k4_part = True

            # whether we begin the single precision zone?
            if re.search(r"WITH_SINGLE_PRECISION", newline) is not None \
                            and re.search(r"ifdef", newline) is not None:
                inSinglePrecisionZone = True

            # whether the single precision zone ends?
            if inSinglePrecisionZone:
                if re.search(r"#endif", newline) is not None \
                                or re.search(r"#else", newline) is not None:
                    inSinglePrecisionZone = False

            # if we are in the single precision zone, we omit
            # all of codes
            if inSinglePrecisionZone:
                continue

            # count the number of ignore integrals
            if is_comment(line): 
	        if re.search(r"integrals are omitted", newline) is not None: 
		    tmp = newline.split(" ")
		    pos = tmp.index("totally")
		    n = tmp[pos+1]
		    ignore = int(n)
                    n_ignore_vrr += ignore
            else:

                # whether we declare a vector?
                if re.search(r"DoubleVec", line) is not None:
                    begin = newline.index("(")
                    end   = newline.index(",")
                    num   = newline[begin+1:end]
                    n     = int(num)
                    vrr_vec_nvar = vrr_vec_nvar + n
                elif re.search(r"getNewMemPos", line) is not None:
                    begin = newline.index("(")
                    end   = newline.index(")")
                    num   = newline[begin+1:end]
                    n     = int(num)
                    vrr_vec_nvar = vrr_vec_nvar + n

                # whether this is a new double var?
                if re.search(r"^Double", newline) is not None:
                    vrr_var_num = vrr_var_num + 1

                # counting the number of adding/substracting
                n = newline.count('+')
                if in_k4_part: 
                    n_add_k4 = n_add_k4 +n
                else:
                    n_add_k2 = n_add_k2 +n
                n = newline.count('-')
                if in_k4_part: 
                    n_add_k4 = n_add_k4 +n
                else:
                    n_add_k2 = n_add_k2 +n

                # counting the number of multiplication/dividing
                n = newline.count('*')
                if in_k4_part: 
                    n_mul_k4 = n_mul_k4 +n
                else:
                    n_mul_k2 = n_mul_k2 +n
                n = newline.count('/')
                if in_k4_part: 
                    n_mul_k4 = n_mul_k4 +n
                else:
                    n_mul_k2 = n_mul_k2 +n

        # finally, close the file
        name.close()

def infor_print(fname):
    """
    print out the results
    """
    global hrr_var_num
    global vrr_var_num
    global hrr_vec_nvar
    global vrr_vec_nvar
    global n_add_k2
    global n_mul_k2
    global n_add_k4
    global n_mul_k4
    global n_add_hrr
    global n_mul_hrr
    global n_ignore_vrr 
    global n_ignore_hrr 

    # do we have K4 part?
    has_k4_part = True
    n_total_flops_k4 = n_add_k4 + n_mul_k4 
    n_total_flops_k2 = n_add_k2 + n_mul_k2 
    if n_total_flops_k4 == 0:
        has_k4_part = False

    # do we have hrr part?
    has_hrr_part = True
    n_total_flops_hrr = n_add_hrr + n_mul_hrr 
    if n_total_flops_hrr == 0:
        has_hrr_part = False

    print "*************************************************************************"
    print "for file ", fname
    print "VRR section:"
    print "number of variables     : ", vrr_var_num
    print "number of vector vars   : ", vrr_vec_nvar
    print "No. of ignored integrals: ", n_ignore_vrr 
    print "total K2 FLOPS          : ", n_total_flops_k2
    print "number of + and - in K2 : ", n_add_k2
    print "number of * and / in K2 : ", n_mul_k2
    if has_k4_part:
        print "number of + and - in K4 : ", n_add_k4
        print "number of * and / in K4 : ", n_mul_k4
        print "total K4 FLOPS          : ", n_total_flops_k4
    if has_hrr_part:
        print "HRR/non-RR/derivatives section:"
        print "number of variables     : ", hrr_var_num
        print "number of vector vars   : ", hrr_vec_nvar
        print "No. of ignored integrals: ", n_ignore_hrr 
        print "total K0 FLOPS          : ", n_total_flops_hrr
        print "number of + and -       : ", n_add_hrr
        print "number of * and /       : ", n_mul_hrr
    print "*************************************************************************\n\n"

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

# on the top of dir, we should have 
# three sub-dir whose name like this
projectDir = getWorkPath(topDir,order)

# now the real work
projects = os.listdir(projectDir)
gotIt = False
for iProj in projects:
    workDir = projectDir + "/" + iProj
    if os.path.isdir(workDir) and job == iProj:
        gotIt = True
        break

# now go over the cpp files
if gotIt:
    files = os.listdir(workDir)
    for iFile in files:

        # we only concentrate on cpp file
        if iFile.find("cpp") < 0:
            continue

        # we only do the main cpp file
        if not isMainCPPFile(iFile):
            continue

        # now clear all of variables
        initilize()

        # now generate the file list
        filelist = formFileList(workDir,iFile)

        # perform VRR analyze
        vrr_analyze(workDir,filelist)

        # perform HRR/non-RR/derivatives analyze
        hrr_analyze(workDir,filelist)

        # report what we get
        infor_print(iFile)
else:
    print workDir
    print "We did not find anything matching the work dir"
    sys.exit()


