# export QE_ROOT = $(HOME)/install/rmta-qe

include make.inc

MAKE = make

.PHONY: all rmta clean

all rmta :
	if test -d src ; then \
	( cd src ; $(MAKE) all || exit 1 ) ; fi

clean :
	rm -vrf bin; (cd src; rm -vrf *.x *.o *.mod *.a *~; cd - || exit 1 )

depend .depend :
	if test -d src ; then \
	( cd src ; $(MAKE) depend || exit 1 ) ; fi
