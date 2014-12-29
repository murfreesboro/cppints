#
#   simple Makefile to build the cppints
#   for debug choice:
#   sqints    : SQINTS_DEBUG
#   rr        : RR_DEBUG
#   rrsqsearch: RRSEARCH_DEBUG
#   all of these debug choice are independent
#
NAME       = cppints
CC		     = g++
CFLAGS     = -Wall -O3
MACRO      = # -DRRSEARCH_DEBUG  #-DRR_DEBUG
INCLUDE    = -Isrc/include 
LIB        = -lboost_filesystem -lboost_system
OBJC       = src/main.o src/shell.o src/basis.o  src/basisutil.o \
				 src/inttype.o src/integral.o  src/shellquartet.o \
				 src/rrbuild.o src/rrsqsearch.o src/rrints.o src/rr.o \
				 src/sqints.o src/infor.o  src/sqintsinfor.o  \
				 src/sqintsprint.o src/codegen.o src/derivinfor.o src/intderiv.o

ALL: $(OBJC) 
	$(CC) -o  $(NAME) $(OBJC) $(LIB)

$(OBJC): %.o:%.cpp
	$(CC) $(INCLUDE) $(CFLAGS) $(MACRO) -c  $< -o $@

.PHONY: clean
clean:
	@cd src; find . -name '*.o' -exec rm -rf {} \;
	@rm -rf cppints
	@rm -rf hgp_os

