include make.of.inc
include make.inc

SOURCE= mesh/wrap.f90 \
        mesh/polyMeshConstant.f90 \
        mesh/polyMeshMath.f90 \
        mesh/graph.f90 \
        mesh/polyMeshMG.f90 \
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
	./mg/util.c \
        ./mg/sfc.c

TECIO_SOURCE= ./tecio/TranslatedString.cpp \
  ./tecio/alloc.cpp            \
  ./tecio/arrlist.cpp          \
  ./tecio/auxdata.cpp          \
  ./tecio/dataio.cpp           \
  ./tecio/dataio4.cpp          \
  ./tecio/dataset.cpp          \
  ./tecio/dataset0.cpp         \
  ./tecio/datautil.cpp         \
  ./tecio/filestream.cpp       \
  ./tecio/geom2.cpp            \
  ./tecio/q_msg.cpp            \
  ./tecio/q_unicode.cpp        \
  ./tecio/set.cpp              \
  ./tecio/strlist.cpp          \
  ./tecio/strutil.cpp          \
  ./tecio/tassert.cpp          \
  ./tecio/tecxxx.cpp           \
  ./tecio/wrapCc.cpp           

OBJECTS=$(patsubst %.f90, %.o, $(SOURCE))

MG_OBJECTS=$(patsubst %.c, %.o, $(MG_SOURCE))

TECIO_OBJECTS=$(patsubst %.cpp, %.o, $(TECIO_SOURCE))

all: libCc.so

libCc.so: $(OBJECTS) $(MG_OBJECTS) $(TECIO_OBJECTS)
	$(FC) -fPIC -g -shared $(OF_LDFLAGS) $(OBJECTS) $(MG_OBJECTS) $(TECIO_OBJECTS) -o libCc.so $(OF_LIBS) -lm

#Cc.x: libCc.so
#	$(FC) -fPIC -g Cc.f90 $(OF_LDFLAGS) $(OBJECTS) $(MG_OBJECTS) $(TECIO_OBJECTS) -o Cc.x $(OF_LIBS) -lm -lCc

%.o: %.f90
	$(FC) -g -fPIC -fbounds-check -c -o $@ $<

%.o: %.c
	$(CC) -g -I./mg -fPIC -c -o $@ $<

%.o: %.cpp
	$(CXX) -g -I./tecio -fPIC $(TECIO_CFLAGS) -c -o $@ $<


.PHONY : clean distclean
clean:
	-@rm -f *.o *.mod
	-@rm -f mesh/*.o mesh/*.mod
	-@rm -f mg/*.o mg/*.mod
	-@rm -f tecio/*.o

distclean: clean
	-@rm -f make.of.inc
	-@rm -f Cc
	-@rm -f Cc.x
	-@rm -rf libCc.so*
	-@rm -rf mesh/Make mesh/lnInclude
	-@rm -rf libCc_ofreader*
	-@rm -rf Cc.x.*
	-@rm -f ./mesh/*.dep
	-@rm -rf ./lua
	cd luajit-2.0 && $(MAKE) clean
	cd luajit-2.0 && rm -f Makefile.old

