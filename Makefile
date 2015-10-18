include make.of.inc
include make.inc

SOURCE= mesh/wrap.f90 \
        mesh/mesh.f90 \
        Cc.f90

OBJECTS=$(patsubst %.f90, %.o, $(SOURCE))

all: Cc.x

Cc.x: $(OBJECTS)
	$(FC) $(OF_LDFLAGS) *.o -o Cc.x $(OF_LIBS)

%.o: %.f90
	$(FC) -c $^

.PHONY : clean distclean
clean:
	-rm -f *.o *.mod

distclean:
	-rm -f *.o *.mod make.of.inc
	-rm -f Cc
	-rm -f Cc.x
	-rm -rf Make lnInclude
	-rm -f libCc_ofreader*
	-rm -f ./mesh/*.dep
