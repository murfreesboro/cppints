
##########################################################
#          general setting for the makefile              #
#          you need to have the following libs:          #
#          boost, mpfr, gmp                              #
#          Here we assume that boost, mpfr etc.          #
#          are just in default place like /usr/lib       #
##########################################################
NAME                 = test
CC		   	         = icpc
CFLAGS       			= -Wall  -O3 #-g 
LOWCFLAGS       		= -Wall  -O0 #-g 
INCLUDE              = -Iinclude -Ihgp_os/include 
LIBRARY              = -lm -lmpfr -lgmp -ltbbmalloc -ltbb
MACRO                = -DDEBUG_EXPR12 -DWITH_SINGLE_PRECISION # -DDEBUG

##########################################################
#                    object files                        #
##########################################################

# for the file eri.cpp, we must do it carefully without -O3
LOWOBJ     = eri.o

# these are the functions we wrote/pick up for testing
OBJMAIN    = test.o
OBJMAIN   += ov.o
OBJMAIN   += ki.o
OBJMAIN   += nai.o
OBJMAIN   += eritest.o
OBJMAIN   += expr12test.o
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
ALL: $(ARONEINT) $(OBJMAIN) $(ARTWOINT) $(LOWOBJ)
	$(CC) -o  $(NAME) $(OBJMAIN) $(LOWOBJ) $(ARONEINT) $(ARTWOINT) $(LIBRARY)

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

$(LOWOBJ): %.o:%.cpp
	$(CC) $(LOWCFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

##########################################################
#                        cleaning                        #
##########################################################
.PHONY: clean
clean:
	@find . -name '*.o' -exec rm -rf {} \;
	@find . -name '*.out' -exec rm -rf {} \;
	@find . -name '*~' -exec rm -rf {} \;
	@find . -name '*.a' -exec rm -rf {} \;

