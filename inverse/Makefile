CXX=mpicxx -g -pedantic -Wno-long-long
F77=mpif77 -g
F77_LIBS=
libext = 
incs =   -I.. -I../util -g   -I/usr/local/pkg/hdf/include -DUSE_HDF  -I/usr/local/pkg/mpich2/include -D__MPI__  -DDEBUG_CROSSING -DDEBUG_NODE_PROPAGATE -DDEBUG_TRI_GRID  -DONED -DTWOD -DTHREED  -DHAS_FENV    -D__GD__ -I/usr/include  
libincs =  -L/usr/lib  /lib64/libblas.so /lib64/libgfortran.so.3 -L/usr/local/pkg/hdf/lib  -L/usr/local/pkg/mpich2/lib 
libs = -lgd  /lib64/libblas.so /lib64/libgfortran.so.3  -lmfhdf -ldf -ljpeg -lz  -lmpich -lpthread
build_cpu = x86_64
incsloc = -I. -I../include -I../solver

PETSC = yes
PETSC_ARCH=
PETSC_DIR=/usr/local/pkg/petsc
Petsc_Include=-I/usr/local/pkg/petsc/include
Petsc_Lib=-L/usr/local/pkg/petsc/lib -lpetsc -llapack -lblas -ldl -lm -L/usr/X11R6/lib -lX11

inverse : FT_libs

interface.o : inverse.h interface.cpp
	${CXX} -c interface.cpp $(incsloc) $(incs) ${Petsc_Include}

iexact.o : inverse.h iexact.cpp
	${CXX} -c iexact.cpp $(incsloc) $(incs) ${Petsc_Include}

cim_ellip.o : ../solver/solver.h ../solver/cim_ellip.cpp
	${CXX} -c ../solver/cim_ellip.cpp $(incsloc) $(incs) ${Petsc_Include}

ellip.o : ../solver/solver.h ../solver/ellip.cpp
	${CXX} -c ../solver/ellip.cpp $(incsloc) $(incs) ${Petsc_Include}

solver.o : ../solver/solver.h ../solver/solver.cpp
	${CXX} -c ../solver/solver.cpp $(incsloc) $(incs) ${Petsc_Include}

icartsn.o : inverse.h ../solver/solver.h icartsn.cpp
	${CXX} -c icartsn.cpp $(incsloc) $(incs) ${Petsc_Include}

inverse: inverse.h ../solver/solver.h inverse.cpp interface.o iexact.o icartsn.o ellip.o cim_ellip.o solver.o
	${CXX} -c inverse.cpp $(incsloc) $(incs) ${Petsc_Include}
	$(CXX) -o inverse inverse.o -I../include -L../lib/${build_cpu} interface.o iexact.o icartsn.o ellip.o cim_ellip.o solver.o ${Petsc_Lib} -L../lib/${build_cpu} -lFronTier $(libincs)  $(libs) ${F77_LIBS}

ppe: ppe.cpp
	${CXX} -c ppe.cpp $(incsloc) $(incs)
	$(CXX) -o ppe ppe.o -I../include -L../lib/${build_cpu} -L../lib/${build_cpu} -lFronTier $(libincs)  $(libs) ${F77_LIBS}

-include ../devel-deps.inc

clean:
	rm -rf *.o inverse 
output :
	rm -rf out/*

gif :
	hdf2gif out/intfc.hdf out/intfc.gif

TagsFile:
	ctags *.h *.cpp ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]

tagsfile:
	ctags *.h *.cpp ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]

etagsfile:
	etags *.h *.cpp ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]
