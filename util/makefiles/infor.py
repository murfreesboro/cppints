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
job_list        = []          # what kind of integral module we would like to include?
maxL            = 5           # maximum L for integral
auxMaxL         = 5           # maximum L for the aux basis sett dimension of integrals
file_size       = 5242880     # this is the default file size that we may keep a caution (5M)
                              # file larger than this size should go with lower compiler flags


#########################################################
#
#  initialize the information
# 
#########################################################
def set_project (fname):

    global makefile_choice
    global job_list
    global maxL
    global auxMaxL
    global file_size 

    # initialize the value to be default
    job_list = [ ]
    makefile_choice = "makefile"
    maxL = 5
    auxMaxL = 5
    file_size = 5242880

    # read in things
    name = open(fname, "r")
    while True:
	line = name.readline()
	if not line: break 

        # skip the empty line and comment lines
	if re.search(r"(?i)#", line) is not None: 
            continue
        if len(line.strip()) == 0:
            continue

        # strp the back space
        l = line.strip()
            
        # now search each key
	if re.search(r"(?i)maxL", l) is not None: 
            tmp = l.split()
            maxL = int(tmp[1])
            if maxL <= 0:
                print line
                print "Illegal value for setting maxL"
		sys.exit()
	elif re.search(r"(?i)auxMaxL", l) is not None: 
            tmp = l.split()
            auxMaxL = int(tmp[1])
            if auxMaxL <= 0:
                print line
                print "Illegal value for setting auxMaxL"
		sys.exit()
	elif re.search(r"(?i)jobs", l) is not None: 
            tmp = l.split()
            for i in range(1,len(tmp)):
                job_name = tmp[i].lower()
                job_list.append(job_name)
	elif re.search(r"(?i)choice", l) is not None: 
            tmp = l.split()
            choice = tmp[1].lower()
            makefile_choice = choice
            if choice != "makefile" and choice != "cmake":
                print "Illegal choice for generating the makefiles/cmakelist.txt"
                sys.exit()
	elif re.search(r"(?i)file_size", l) is not None: 
            tmp  = l.split()
            size = int(tmp[1])
            if size <= 0:
                print line
                print "Illegal value for setting file size"
		sys.exit()
            file_size = size
        else:
            print line
	    print "Illegal line in the infor file, can not be recognized!"
	    sys.exit()
    name.close()

    # final check
    if auxMaxL< maxL:
	print "maxL should be same or larger than the aux max L"
	sys.exit()

#########################################################
#
#  return the value
# 
#########################################################
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

def hasJob(name):
    global job_list
    
    # divide the name from the path
    job = name.lower()
    job_name = os.path.basename(job)
    
    for i in job_list:
        if i == job_name:
            return True
    return False

def largeFile(name):
    global file_size

    # get the size of this file
    size = os.path.getsize(name)
    if size > file_size:
        return True
    else:
        return False


