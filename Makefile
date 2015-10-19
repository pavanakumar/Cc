include make.of.inc
include make.inc

SOURCE= mesh/wrap.f90 \
        mesh/mesh.f90 \
        Cc.f90

MG_SOURCE= ./mg/aratio.c \
	./mg/blas.c \
	./mg/coarsen.c \
	./mg/dfkeysort.c \
	./mg/dkeysort.c \
	./mg/file.c \
	./mg/ifkeysort.c \
	./mg/ifloatsort.c \
	./mg/iintsort.c \
	./mg/ikeysort.c \
	./mg/io.c \
	./mg/kwayfm.c \
	./mg/match.c \
	./mg/memory.c \
	./mg/merge.c \
	./mg/mgridgen.c \
	./mg/refine.c \
	./mg/setup.c \
	./mg/sort.c \
	./mg/util.c

OBJECTS=$(patsubst %.f90, %.o, $(SOURCE))

MG_OBJECTS=$(patsubst %.c, %.o, $(MG_SOURCE))

all: Cc.x

Cc.x: $(OBJECTS) $(MG_OBJECTS)
	$(FC) $(OF_LDFLAGS) *.o -o Cc.x $(OF_LIBS) -lm

%.o: %.f90
	$(FC) -c $^

%.o: %.c
	$(CC) -I./mg -c -fPIC $^

.PHONY : clean distclean
clean:
	-@rm -f *.o *.mod

distclean:
	-@rm -f *.o *.mod
	-@rm -f make.of.inc
	-@rm -f Cc
	-@rm -f Cc.x
	-@rm -rf Make lnInclude
	-@rm -f libCc_ofreader*
	-@rm -f ./mesh/*.dep

