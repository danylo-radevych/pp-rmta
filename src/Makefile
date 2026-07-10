# Makefile for rmta.x
# QE_ROOT = $(HOME)/install/rmta-qe

# include ../make.inc
include $(QE_ROOT)/make.inc

# location of needed modules and included files (if any)
MODFLAGS = $(BASEMOD_FLAGS) \
          $(MOD_FLAG)$(QE_ROOT)/dft-d3 \
          $(MOD_FLAG)$(QE_ROOT)/PW/src \
          $(MOD_FLAG)$(QE_ROOT)/atomic/src \
          $(MOD_FLAG)$(QE_ROOT)/PP/src \
          $(MOD_FLAG)./src

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

PPOBJS = $(QE_ROOT)/PP/src/libpp.a
LD1OBJS = $(QE_ROOT)/atomic/src/libld1.a
PWOBJS = $(QE_ROOT)/PW/src/libpw.a $(QE_ROOT)/KS_Solvers/libks_solvers.a $(QE_ROOT)/dft-d3/libdftd3qe.a
QEMODS = $(QE_ROOT)/KS_Solvers/libks_solvers.a $(QE_ROOT)/dft-d3/libdftd3qe.a $(BASEMODS)

.PHONY: all

all : librmta.a rmta.x

librmta.a : $(RMTAOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

rmta.x : rmta.o librmta.a $(PPOBJS) $(LD1OBJS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	rmta.o librmta.a $(PPOBJS) $(LD1OBJS) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../; mkdir -vp bin; cd bin; ln -fs ../src/$@ . )


clean :
	- /bin/rm -f *.x *.o *.a *~ *.mod
	- /bin/rm -f ../bin/rmta.x librmta.a

include make.depend
