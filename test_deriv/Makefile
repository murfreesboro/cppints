
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

##########################################################
#               define the deriv order                   #
##########################################################
DERIV_ORDER          = 1
ifeq (${DERIV_ORDER}, 1)
MACRO                = -DORDER1 #-DWITH_SINGLE_PRECISION -DDEBUG  
else
MACRO                = -DORDER2 #-DWITH_SINGLE_PRECISION -DDEBUG
endif

##########################################################
#                    object files                        #
##########################################################
# for the file eri.cpp, we must do it carefully without -O3
LOWOBJ     = eri.o

# these are the functions we wrote/pick up for testing
OBJMAIN    = derivinfor.o
OBJMAIN   += textread.o
OBJMAIN   += ov.o
OBJMAIN   += ki.o
OBJMAIN   += nai.o
OBJMAIN   += eritest.o
OBJMAIN   += norm.o
OBJMAIN   += test.o

# hgp module of codes
# main ones
OBJMAIN1    = hgp_os/first_deriv/hgp_os_twobodyoverlap_d1.o
OBJMAIN1   += hgp_os/first_deriv/hgp_os_kinetic_d1.o
OBJMAIN1   += hgp_os/first_deriv/hgp_os_nai_d1.o
OBJMAIN1   += hgp_os/first_deriv/hgp_os_eri_d1.o
OBJMAIN2    = hgp_os/second_deriv/hgp_os_twobodyoverlap_d2.o
OBJMAIN2   += hgp_os/second_deriv/hgp_os_kinetic_d2.o
OBJMAIN2   += hgp_os/second_deriv/hgp_os_nai_d2.o
OBJMAIN2   += hgp_os/second_deriv/hgp_os_eri_d2.o

# real integral codes
OBJONEINT1  = $(patsubst %.cpp, %.o,   $(wildcard hgp_os/first_deriv/twobodyoverlap/*.cpp))
OBJONEINT1 += $(patsubst %.cpp, %.o,   $(wildcard hgp_os/first_deriv/kinetic/*.cpp))
OBJONEINT1 += $(patsubst %.cpp, %.o,   $(wildcard hgp_os/first_deriv/nai/*.cpp))
OBJTWOINT1  = $(patsubst %.cpp, %.o,   $(wildcard hgp_os/first_deriv/eri/*.cpp))
OBJONEINT2  = $(patsubst %.cpp, %.o,   $(wildcard hgp_os/second_deriv/twobodyoverlap/*.cpp))
OBJONEINT2 += $(patsubst %.cpp, %.o,   $(wildcard hgp_os/second_deriv/kinetic/*.cpp))
OBJONEINT2 += $(patsubst %.cpp, %.o,   $(wildcard hgp_os/second_deriv/nai/*.cpp))
OBJTWOINT2  = $(patsubst %.cpp, %.o,   $(wildcard hgp_os/second_deriv/eri/*.cpp))

##########################################################
#                    archive files                       #
##########################################################
ARONEINT1   = oneint_d1.a
ARTWOINT1   = twoint_d1.a
ARONEINT2   = oneint_d2.a
ARTWOINT2   = twoint_d2.a

##########################################################
#                    compiling starts                    #
##########################################################
ifeq (${DERIV_ORDER}, 1)
ALL: $(ARONEINT1) $(ARTWOINT1) $(OBJMAIN1) $(OBJMAIN) $(LOWOBJ)
	$(CC) -o  $(NAME) $(OBJMAIN) $(OBJMAIN1) $(LOWOBJ) $(ARONEINT1) $(ARTWOINT1) $(LIBRARY)
else
ALL: $(ARONEINT2) $(ARTWOINT2) $(OBJMAIN2) $(OBJMAIN) $(LOWOBJ)
	$(CC) -o  $(NAME) $(OBJMAIN) $(OBJMAIN2) $(LOWOBJ) $(ARONEINT2) $(ARTWOINT2) $(LIBRARY)
endif

$(ARONEINT1): $(OBJONEINT1) 
	$(AR) -r ${ARONEINT1} ${OBJONEINT1}

$(ARTWOINT1): $(OBJTWOINT1)
	$(AR) -r ${ARTWOINT1} ${OBJTWOINT1}

$(ARONEINT2): $(OBJONEINT2) 
	$(AR) -r ${ARONEINT2} ${OBJONEINT2}

$(ARTWOINT2): $(OBJTWOINT2)
	$(AR) -r ${ARTWOINT2} ${OBJTWOINT2}

$(OBJMAIN): %.o:%.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

$(OBJMAIN1): %.o:%.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

$(OBJMAIN2): %.o:%.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

$(OBJONEINT1): %.o:%.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

$(OBJTWOINT1): %.o:%.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

$(OBJONEINT2): %.o:%.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

$(OBJTWOINT2): %.o:%.cpp
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

