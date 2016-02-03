
##########################################################
#          general setting for the makefile              #
#          you need to have the following libs:          #
#          boost, mpfr, gmp                              #
#          Here we assume that boost, mpfr etc.          #
#          are just in default place like /usr/lib       #
##########################################################
NAME                 = test
CC		   	         = g++
CFLAGS       			= -Wall # -O2 #-g 
INCLUDE              = -Iinclude -Ihgp_os/include -I${TBB_ROOT}/include
LIBRARY              = -lm -lmpfr -lgmp -ltbbmalloc -ltbb
MACRO                = -DDEBUG  -DDEBUG_EXPR12 #-DWITH_SINGLE_PRECISION

##########################################################
#                    object files                        #
##########################################################

# these are the functions we wrote/pick up for testing
OBJMAIN    = test.o
OBJMAIN   += ov.o
OBJMAIN   += ki.o
OBJMAIN   += nai.o
OBJMAIN   += eritest.o
OBJMAIN   += expr12test.o
OBJMAIN   += eri.o
OBJMAIN   += tov.o
OBJMAIN   += tki.o
OBJMAIN   += mom.o
OBJMAIN   += norm.o

# hgp module of codes
# main ones
OBJMAIN   += hgp_os/energy/hgp_os_twobodyoverlap.o
OBJMAIN   += hgp_os/energy/hgp_os_kinetic.o
OBJMAIN   += hgp_os/energy/hgp_os_nai.o
OBJMAIN   += hgp_os/energy/hgp_os_eri.o
OBJMAIN   += hgp_os/energy/hgp_os_threebodyoverlap.o
OBJMAIN   += hgp_os/energy/hgp_os_threebodyki.o
OBJMAIN   += hgp_os/energy/hgp_os_mom.o
OBJMAIN   += hgp_os/energy/hgp_os_expr12.o

# real integral codes
OBJONEINT  = $(patsubst %.cpp, %.o,   $(wildcard hgp_os/energy/twobodyoverlap/*.cpp))
OBJONEINT += $(patsubst %.cpp, %.o,   $(wildcard hgp_os/energy/kinetic/*.cpp))
OBJONEINT += $(patsubst %.cpp, %.o,   $(wildcard hgp_os/energy/nai/*.cpp))
OBJONEINT += $(patsubst %.cpp, %.o,   $(wildcard hgp_os/energy/threebodyoverlap/*.cpp))
OBJONEINT += $(patsubst %.cpp, %.o,   $(wildcard hgp_os/energy/threebodyki/*.cpp))
OBJONEINT += $(patsubst %.cpp, %.o,   $(wildcard hgp_os/energy/mom/*.cpp))
OBJTWOINT  = $(patsubst %.cpp, %.o,   $(wildcard hgp_os/energy/eri/*.cpp))
OBJTWOINT += $(patsubst %.cpp, %.o,   $(wildcard hgp_os/energy/expr12/*.cpp))

##########################################################
#                    archive files                       #
##########################################################
ARONEINT   = oneint.a
ARTWOINT   = twoint.a

##########################################################
#                    compiling starts                    #
##########################################################
ALL: $(ARONEINT) $(OBJMAIN) $(ARTWOINT) 
	$(CC) -o  $(NAME) $(OBJMAIN) $(ARONEINT) $(ARTWOINT) $(LIBRARY)

$(ARONEINT): $(OBJONEINT)
	$(AR) -r ${ARONEINT} ${OBJONEINT}

$(ARTWOINT): $(OBJTWOINT)
	$(AR) -r ${ARTWOINT} ${OBJTWOINT}

$(OBJMAIN): %.o:%.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

$(OBJONEINT): %.o:%.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

$(OBJTWOINT): %.o:%.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@


##########################################################
#                        cleaning                        #
##########################################################
.PHONY: clean
clean:
	@find . -name '*.o' -exec rm -rf {} \;
	@find . -name '*.out' -exec rm -rf {} \;
	@find . -name '*~' -exec rm -rf {} \;
	@find . -name '*.a' -exec rm -rf {} \;

