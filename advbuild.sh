#!/bin/bash

HOST=`uname -n`;
NEWONE=""
OLDEONE=""
WITHGAS=0
WITHMPI=0
WITHHDF=0
WITHGD=0
WITHPETSC=0
WITHDEVEL=0
NOMAKE=0

if [[ "$1" == "-h" ]]; then
    echo "$0    High-level build script for FronTier."
    echo "Usage: $0 [-d] [-n] ... [--with-mpich]."
    echo
    echo "    -d            Enable debugging."
    echo "    -n            Just configure. Do not run make."
    echo "    -gas          Enable compilation with gas."
    echo "    -mpi          Enable compilation with mpi."
    echo "    -devel        Enable compilation with devel."
    echo "    -hdf          Enable compilation with hdf."
    echo "    -gd           Enable compilation with gd."
    echo "    -petsc        Enable compilation with petsc."
    echo "    --with-***    Enable other options used in ./configure."
    exit
else
    OPTS="-O3"
    CONF="--with-no-debugging"
    for arg in $* ; do
	if [[ "$arg" == "-d" ]]; then
	    OPTS="-g"
	    COPTS_GCC="-pedantic -Wno-long-long"
	    COPTS_ICC="-Wall"
	    CONF=""
	elif [[ "$arg" == "-n" ]]; then
	    NOMAKE=1
	elif [[ "$arg" == "-gas" ]]; then
	    WITHGAS=1
	elif [[ "$arg" == "-mpi" ]]; then
	    WITHMPI=1
	elif [[ "$arg" == "-hdf" ]]; then
	    WITHHDF=1
	elif [[ "$arg" == "-gd" ]]; then
	    WITHGD=1
	elif [[ "$arg" == "-devel" ]]; then
	    WITHDEVEL=1
	elif [[ "$arg" == "-petsc" ]]; then
	    WITHPETSC=1
        else
	    NEWONE=\ ${OLDONE}\ ${arg}\ 
	    OLDONE=${NEWONE}
	fi
    done
fi
OTHEROPTION=${NEWONE}

# First, run autoconf to generate ./configure
autoconf

# Second, run ./configure. Choose proper compilers based on hostnames
if [[ "${HOST}" == "fenp" ]]; then
    ######################## Blue Gene/P ######################
    echo "Computer is recognized as NYBlue/P."

    export CXX="mpixlcxx_r -DMPICH_IGNORE_CXX_SEEK ${OPTS}"
    export F77="mpixlf77_r ${OPTS}"
    export CC="mpixlc_r ${OPTS}"
    PMAKE="-j8"

    #Set up configure
    if [ ${WITHGAS} == 1 ]; then
	    GASOPTION="--with-gas"
    else
	    GASOPTION=""
    fi
    if [ ${WITHMPI} == 1 ]; then
	    MPIOPTION="--with-mpi"
    else
	    MPIOPTION=""
    fi
    if [ ${WITHPETSC} == 1 ]; then
	    PETSCOPTION="--with-petsc=${PETSC_DIR}"
    else
	    PETSCOPTION=""
    fi
    if [ ${WITHDEVEL} == 1 ]; then
	    DEVELOPTION="--with-devel"
    else
	    DEVELOPTION=""
    fi

    #Configure
    ./configure ${GASOPTION} ${MPIOPTION} --with-extra-libs="-lmpich.cnk" ${CONF} ${OTHEROPTION}
elif [[ "${HOST}" == "fen" ]]; then
    ################################## Blue Gene/L #########################
    echo "Computer is recognized as NYBlue/L."

    export PETSC_ARCH=bgl-ibm-opt
    export PETSC_DIR=/bgl/apps/petsc-2.3.3-p8/petsc-2.3.3-p8

    export HYPRE_DIR=/bgl/apps/petsc-2.3.3-p8/petsc-2.3.3-p8/externalpackages/hypre-2.0.0/bgl-ibm-opt/lib
    export SUPERLU_DIST_DIR=/bgl/apps/petsc-2.3.3-p8/petsc-2.3.3-p8/externalpackages/SuperLU_DIST_2.0-Jan_5_2006/bgl-ibm-opt
    export SUPERLU_DIR=/bgl/apps/petsc-2.3.3-p8/petsc-2.3.3-p8/externalpackages/SuperLU_3.0-Jan_5_2006/bgl-ibm-opt
    export NSS_FILES_DIR=/bgl/BlueLight/ppcfloor/blrts-gnu/powerpc-bgl-blrts-gnu/lib

    export PETSC_LIB="-L${PETSC_DIR}/lib/${PETSC_ARCH} -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lg2c -L/bgl/local/lib -L${HYPRE_DIR} -lHYPRE -L${SUPERLU_DIR} -lsuperlu_3.0 -L${SUPERLU_DIST_DIR} -lsuperlu_dist_2.0  -lc -L${NSS_FILES_DIR} -lnss_files -lnss_dns -lresolv -llapack.rts -lblas.rts -L/usr/lib -ldl -lm"

    export F77_LIBS="-L/opt/ibmcmp/xlsmp/bg/1.7/blrts_lib -L/opt/ibmcmp/xlmass/bg/4.4/blrts_lib -L/opt/ibmcmp/xlf/bg/11.1/blrts_lib -lxlf90 -lxlopt -lxlomp_ser -lxl -lxlfmath -lm -lc -lgcc"

    export CXX="mpixlcxx -DMPICH_IGNORE_CXX_SEEK ${OPTS}"
    export F77="mpixlf77 ${OPTS}"
    export CC="mpixlc ${OPTS}"
    PMAKE="-j8"

    #Set up configure

    if [ ${WITHGAS} == 1 ]; then
	    GASOPTION="--with-gas"
    else
	    GASOPTION=""
    fi
    if [ ${WITHMPI} == 1 ]; then
	    MPIOPTION="--with-mpi"
    else
	    MPIOPTION=""
    fi
    if [ ${WITHPETSC} == 1 ]; then
	    PETSCOPTION="--with-petsc=${PETSC_DIR}"
    else
	    PETSCOPTION=""
    fi
    if [ ${WITHDEVEL} == 1 ]; then
	    DEVELOPTION="--with-devel"
    else
	    DEVELOPTION=""
    fi


    #Configure
    ./configure ${GASOPTION} ${MPIOPTION} --with-extra-libs="-lmpich.rts" ${PETSCOPTION} ${DEVELOPTION} ${CONF} ${OTHEROPTION}
elif [[ "${HOST}" == "sirius" ]]; then
    ############################# Galaxy ####################################
    echo "Computer is recognized as Galaxy."

    export CXX="mpicxx ${OPTS} ${COPTS_GCC}"
    export F77="mpif77 ${OPTS}"
    export CC="mpicc ${OPTS} ${COPTS_GCC}"
    PMAKE="-j2"

    export PETSC_DIR=/usr/local/pkg/petsc-2.3.3/petsc-2.3.3-p11
    export PETSC_ARCH=linux-gnu-c-debug
    export PETSC_LIB="-lblas -llapack -L${PETSC_DIR}/lib/${PETSC_ARCH} -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -ldl -lm -lX11"


    #Set up configure

    if [ ${WITHGAS} == 1 ]; then
	    GASOPTION="--with-gas"
    else
	    GASOPTION=""
    fi
    if [ ${WITHMPI} == 1 ]; then
	    MPIOPTION="--with-openmpi=/usr/local/pkg/openmpi"
    else
	    MPIOPTION=""
    fi
    if [ ${WITHPETSC} == 1 ]; then
	    PETSCOPTION="--with-petsc=${PETSC_DIR}"
    else
	    PETSCOPTION=""
    fi
    if [ ${WITHHDF} == 1 ]; then
	    HDFOPTION="--with-hdf=/usr/local/pkg/HDF4"
    else
	    HDFOPTION=""
    fi
    if [ ${WITHGD} == 1 ]; then
	    GDOPTION="--with-gd=/usr/local/pkg/gd"
    else
	    GDIOPTION=""
    fi
    if [ ${WITHDEVEL} == 1 ]; then
	    DEVELOPTION="--with-devel"
    else
	    DEVELOPTION=""
    fi

    #Configure

    echo The configure command is ./configure ${GASOPTION} ${HDFOPTION} ${MPIOPTION} ${PETSCOPTION} ${GDOPTION} ${DEVELOPTION} ${CONF} ${OTHEROPTION}

    ./configure ${GASOPTION} ${HDFOPTION} ${MPIOPTION} ${PETSCOPTION} ${GDOPTION} ${DEVELOPTION} ${CONF} ${OTHEROPTION}
elif [[ "${HOST}" == "seawulf" ]]; then
    ################################## Seawulf #####################################3
    echo "Computer is recognized as Seawulf."
    if [[ -z "$LD_LIBRARY_PATH" ]]; then
        echo "#### There seems to be a problem in your environment setting."
        echo "#### You need to set LD_LIBRARY_PATH to include /usr/local/pkg/openmpi/lib."
	exit;
    fi

    export CXX="mpicxx ${OPTS} ${COPTS_GCC}"
    export F77="mpif77 ${OPTS}"
    export CC="mpicc ${OPTS} ${COPTS_GCC}"
    PMAKE="-j2"

    export PETSC_DIR=/usr/local/pkg/petsc-2.3.3-p11
    export PETSC_ARCH=linux-gnu-c-debug
    export PETSC_LIB="-lblas -llapack -L${PETSC_DIR}/lib/${PETSC_ARCH} -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -ldl -lm -lX11"

    #Set up configure
    if [ ${WITHGAS} == 1 ]; then
	    GASOPTION="--with-gas"
    else
	    GASOPTION=""
    fi
    if [ ${WITHMPI} == 1 ]; then
	    MPIOPTION="--with-openmpi"
    else
	    MPIOPTION=""
    fi
    if [ ${WITHPETSC} == 1 ]; then
	    PETSCOPTION="--with-petsc=${PETSC_DIR}"
    else
	    PETSCOPTION=""
    fi
    if [ ${WITHHDF} == 1 ]; then
	    HDFOPTION="--with-hdf=/usr/local/pkg/HDF4"
    else
	    HDFOPTION=""
    fi
    if [ ${WITHGD} == 1 ]; then
	    GDOPTION="--with-gd=/usr/local/pkg/gd"
    else
	    GDIOPTION=""
    fi
    if [ ${WITHDEVEL} == 1 ]; then
	    DEVELOPTION="--with-devel"
    else
	    DEVELOPTION=""
    fi

    #Configure
    ./configure ${GASOPTION} ${HDFOPTION} ${MPIOPTION} ${PETSCOPTION} ${DEVELOPTION} ${CONF} ${OTHEROPTION}
elif [[ "${HOST//[0-9]/}" == "honest.ncsa.uiuc.edu" ]]; then
    ################################### NCSA Abe Linux cluster ##########################
    echo "Computer is recognized as NCSA Abe Linux cluster."

    export CXX="mpicxx -DMPICH_IGNORE_CXX_SEEK -Wcheck ${OPTS} ${COPTS_ICC}"
    export F77="mpif77 ${OPTS}"
    export CC="mpicc ${OPTS} ${COPTS_ICC}"
    PMAKE="-j8"

    # Specify Petsc path
    export PETSC_ARCH=abe-intel10-opt
    export PETSC_DIR=/usr/apps/math/petsc/petsc-2.3.3-p7
    export PETSC_LIB="-lblas -llapack -L${PETSC_DIR}/lib/${PETSC_ARCH} -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -ldl -lm -L/usr/X11R6/lib64"

    # add szip library
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/apps/hdf/szip/lib
    export scriptotherlibsinc="${scriptotherlibincs} -L/usr/apps/hdf/szip/lib"
    export scriptotherlibs="${scriptotherlibs} -lsz"

    #Set up configure
    if [ ${WITHGAS} == 1 ]; then
	    GASOPTION="--with-gas"
    else
	    GASOPTION=""
    fi
    if [ ${WITHMPI} == 1 ]; then
	    MPIOPTION="--with-openmpi"
    else
	    MPIOPTION=""
    fi
    if [ ${WITHPETSC} == 1 ]; then
	    PETSCOPTION="--with-petsc=${PETSC_DIR}"
    else
	    PETSCOPTION=""
    fi
    if [ ${WITHHDF} == 1 ]; then
	    HDFOPTION="--with-hdf=/usr/apps/hdf/hdf4/v424"
    else
	    HDFOPTION=""
    fi
    if [ ${WITHGD} == 1 ]; then
	    GDOPTION="--with-gd=/usr/local/pkg/gd"
    else
	    GDIOPTION=""
    fi
    if [ ${WITHDEVEL} == 1 ]; then
	    DEVELOPTION="--with-devel"
    else
	    DEVELOPTION=""
    fi

    #Configure
    ./configure ${GASOPTION} ${HDFOPTION} --with-mpich ${PETSCOPTION} ${DEVELOPTION} ${CONF} ${OTHEROPTION}
else
    ############## Other platforms. Assume g++, OpenMPI, and HDF4 installed in system path #############
    echo "Computer was not recognized. Using generic configure options."

    export CXX="mpicxx ${OPTS} ${COPTS_GCC} -I/usr/include/hdf -I/usr/include/petsc"
    export F77="mpif77 ${OPTS}"
    export CC="mpicc ${OPTS} ${COPTS_GCC} -I/usr/include/hdf -I/usr/include/petsc"
    PMAKE="-j2"

    #Set up configure  ##Please modify for your own platform!
    if [ ${WITHGAS} == 1 ]; then
	    GASOPTION="--with-gas"
    else
	    GASOPTION=""
    fi
    if [ ${WITHMPI} == 1 ]; then
	    MPIOPTION="--with-openmpi"
    else
	    MPIOPTION=""
    fi
    if [ ${WITHPETSC} == 1 ]; then
	    PETSCOPTION="--with-petsc=${PETSC_DIR}"
    else
	    PETSCOPTION=""
    fi
    if [ ${WITHHDF} == 1 ]; then
	    HDFOPTION="--with-hdf"
    else
	    HDFOPTION=""
    fi
    if [ ${WITHGD} == 1 ]; then
	    GDOPTION="--with-gd"
    else
	    GDIOPTION=""
    fi
    if [ ${WITHDEVEL} == 1 ]; then
	    DEVELOPTION="--with-devel"
    else
	    DEVELOPTION=""
    fi

    #Configure
    ./configure ${GASOPTION} ${HDFOPTION} ${DEVELOPTION} ${MPIOPTION} ${CONF} ${OTHEROPTION}
fi

# Third, invoke make
if [ ${NOMAKE} == 0 ]; then
    make ${PMAKE}
fi
