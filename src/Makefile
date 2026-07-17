# Makefile for rmta.x

include $(QE_ROOT)/make.inc

# location of needed modules and included files (if any)
MODFLAGS = $(OPTFLAGS) \
          $(BASEMOD_FLAGS) \
          $(MOD_FLAG)$(QE_ROOT)/dft-d3 \
          $(MOD_FLAG)$(QE_ROOT)/PW/src \
          $(MOD_FLAG)$(QE_ROOT)/atomic/src \
          $(MOD_FLAG)$(QE_ROOT)/PP/src \
          $(MOD_FLAG)../src

RMTASRC = \
muffin_tin.f90 \
mh_eta.f90 \
part_dos.f90 \
hamann.f90 \
compare.f90 \
sym_type.f90 \
mt_var.f90 \
neighbor.f90 \
mt_printing.f90 \
const.f90

### object files:
RMTAOBJS = $(addsuffix .o , $(basename $(RMTASRC)))

PPOBJS = $(QE_ROOT)/PP/src/libpp.a
#LD1OBJS = $(QE_ROOT)/atomic/src/libld1.a
LD1OBJS := $(filter-out %/ld1.o, $(wildcard $(QE_ROOT)/atomic/src/*.o))
PWOBJS = $(QE_ROOT)/PW/src/libpw.a $(QE_ROOT)/KS_Solvers/libks_solvers.a $(QE_ROOT)/dft-d3/libdftd3qe.a
QEMODS = $(QE_ROOT)/KS_Solvers/libks_solvers.a $(QE_ROOT)/dft-d3/libdftd3qe.a $(BASEMODS)

.PHONY: all

all : librmta.a rmta.x

librmta.a : $(RMTAOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

rmta.x : rmta.o librmta.a $(PPOBJS) $(LD1OBJS) $(LD1MODS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	rmta.o librmta.a $(PPOBJS) $(LD1OBJS) $(LD1MODS) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../; mkdir -vp bin; cd bin; ln -fs ../src/$@ . )
	- ( cd ../; echo "export QE_ROOT=$(QE_ROOT)" > paths.sh; chmod +x paths.sh )
	- ( cd ../; PPRMTA_ROOT=$$(pwd); echo "export PPRMTA_ROOT=$$PPRMTA_ROOT" >> paths.sh )
# 	@$(SHELL) -c 'PPRMTA_ROOT=$$(pwd)/../; cd ../; echo "export PPRMTA_ROOT=$$PPRMTA_ROOT" >> paths.sh'


clean :
	- /bin/rm -f *.x *.o *.a *~ *.mod
	- /bin/rm -f ../bin/rmta.x librmta.a

depend .depend :
	makedepf90 -free $(RMTASRC) > .depend

include .depend
