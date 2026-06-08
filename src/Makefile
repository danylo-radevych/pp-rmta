# Makefile for rmta.x

include ../../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../../PW/src \
          $(MOD_FLAG)../../../PHonon/PH \
          $(MOD_FLAG)../../../PP/src \
          $(MOD_FLAG)../../../atomic/src \
          $(MOD_FLAG)../../../EPW/src \
          $(MOD_FLAG)../../../EPW/RMTA/src

RMTAOBJS = \
muffin_tin.o \
mh_eta.o \
part_dos.o \
hamann.o \
rmta_lschps_qe.o \
compare.o \
sym_type.o \
mt_var.o \
neighbor.o \
mt_printing.o

EPWOBJS = ../../../EPW/src/libepw.a
LD1OBJS = ../../../atomic/src/libld1.a
PPOBJS = ../../../PP/src/libpp.a
PHOBJS = ../../../PHonon/PH/libph.a
PWOBJS = ../../../PW/src/libpw.a ../../../KS_Solvers/libks_solvers.a ../../../dft-d3/libdftd3qe.a
QEMODS = $(BASEMODS)

all : pw ph pp ld1 librmta.a rmta.x

librmta.a : $(RMTAOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

rmta.x : rmta.o librmta.a $(EPWOBJS) $(LD1OBJS) $(PPOBJS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	rmta.o librmta.a $(EPWOBJS) $(LD1OBJS) $(PPOBJS) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../../bin ; ln -fs ../EPW/RMTA/src/$@ . )


epw :
	cd ../../../ ; make epw

ld1 :
	cd ../../../ ; make ld1

pp :
	cd ../../../ ; make pp

ph :
	cd ../../../ ; make ph

pw :
	cd ../../../ ; make pw


clean :
	- /bin/rm -f *.x *.o *.a *~ *.mod
	- /bin/rm -f ../../../bin/rmta.x librmta.a

include make.depend
