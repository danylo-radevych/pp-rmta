# Makefile for rmta.x

include ../../../make.inc

# location of needed modules and included files (if any)
MODFLAGS = $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../../dft-d3 \
          $(MOD_FLAG)../../../PW/src \
          $(MOD_FLAG)../../../atomic/src \
          $(MOD_FLAG)../../../PP/src \
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
mt_printing.o \
const.o

PPOBJS = ../../../PP/src/libpp.a
LD1OBJS = ../../../atomic/src/libld1.a
PWOBJS = ../../../PW/src/libpw.a ../../../KS_Solvers/libks_solvers.a ../../../dft-d3/libdftd3qe.a
QEMODS = ../../../KS_Solvers/libks_solvers.a ../../../dft-d3/libdftd3qe.a $(BASEMODS)

.PHONY: all pw pp ld1

all : pw librmta.a rmta.x

librmta.a : $(RMTAOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

rmta.x : rmta.o librmta.a $(PPOBJS) $(LD1OBJS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	rmta.o librmta.a $(PPOBJS) $(LD1OBJS) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../../bin ; ln -fs ../EPW/RMTA/src/$@ . )

pp $(PPOBJS) :
	cd ../../../ ; make pp

ld1 $(LD1OBJS) :
	cd ../../../ ; make ld1

pw $(PWOBJS) :
	cd ../../../ ; make pw


clean :
	- /bin/rm -f *.x *.o *.a *~ *.mod
	- /bin/rm -f ../../../bin/rmta.x librmta.a

include make.depend
