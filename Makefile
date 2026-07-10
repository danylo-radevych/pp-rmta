export QE_ROOT = $(HOME)/install/rmta-qe

MAKE = make

.PHONY: all rmta clean

all rmta :
	if test -d src ; then \
	( cd src ; $(MAKE) all || exit 1 ) ; fi

clean :
	rm -vrf bin; rm -vrf src/rmta.x src/*.o src/*.mod
