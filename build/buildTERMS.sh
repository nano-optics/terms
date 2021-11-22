# This script is a poor substitute for (c)make...

sysLAPACK=true # links with your system's lapack
debug=false    # enables debug flags 
quad=false     # promotes all real(8) to real(16)

if $quad && $sysLAPACK ; then
    echo "ERROR: Quad precision incompatible with system's LAPACK linking"
    exit 1
fi

if $debug ; then
    flags="-Og -g -fbacktrace -Wall -Wsurprising -fbounds-check -fcheck=all \
    	       -ffpe-trap=overflow,denormal,underflow,invalid"
else
    flags="-O3"
fi

libs="" #    -L/usr/lib/x86_64-linux-gnu -llapack -lblas -lhdf5


if $sysLAPACK ; then
    extra="../src/extLibs/toms644.f"
    flags="$flags -fexternal-blas"
    libs="$libs -llapack -lblas "
    
    
else
    extra=$( ls ../src/extLibs/*.{f,F})
 
fi
  

if $quad ; then
    flags="$flags -freal-8-real-16"
fi
MKLROOT="/usr/lib/x86_64-linux-gnu"

gfortran  -I${MKLROOT}/hdf5/serial/include \
    $extra \
    ../src/eps.f90 \
    ../src/swav.f90 \
    ../src/sphmsv.f90 \
    ../src/miet.f90 \
    ../src/linalg.f90 \
    ../src/HDFfive.f90 \
    ../src/multiscat.f90 \
    ../src/termsProgram.f90 \
     -o terms \
    ${MKLROOT}/hdf5/serial/libhdf5.so \
    ${MKLROOT}/hdf5/serial/libhdf5_fortran.so \
    ${MKLROOT}/hdf5/serial/libhdf5hl_fortran.so \
    ${MKLROOT}/hdf5/serial/libhdf5_hl.so \
    ${MKLROOT}/hdf5/serial/libhdf5_hl_cpp.so \
    ${MKLROOT}/hdf5/serial/libhdf5_cpp.so \
    ${MKLROOT}/libz.so  ${MKLROOT}/libsz.so \
    ${MKLROOT}/libpthread.so  ${MKLROOT}/libm.so \
    ${MKLROOT}/libdl.so \
    $flags $libs \

rm *.mod    
    



