# Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.0.2)


SET (PROJECTNAME "terms")

PROJECT (${PROJECTNAME} Fortran)

ENABLE_LANGUAGE(C) 

FIND_PACKAGE(LAPACK REQUIRED)
FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(HDF5 COMPONENTS Fortran REQUIRED)
# note: can specify a path with HDF5_ROOT to resolve conflicts, 
# e.g. on MacOS after brew install hdf5
# export HDF5_ROOT=/usr/local/Cellar/hdf5/1.12.0_3/

INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})


set(sysLAPACK TRUE) # links with your system's lapack
set(debug FALSE)    # enables debug flags 
set(quad FALSE)     # promotes all real(8) to real(16)

if (quad AND sysLAPACK)
    message("ERROR: Quad precision incompatible with system's LAPACK linking")
    return()
endif()

if (debug)
    SET (CMAKE_Fortran_FLAGS "-Og -g -fbacktrace -Wall -Wsurprising -fbounds-check -fcheck=all \
    	       -ffpe-trap=overflow,denormal,underflow,invalid")
else()
    SET (CMAKE_Fortran_FLAGS "-O3")
endif()



if (sysLAPACK)
     # extra="../src/extLibs/toms644.f"
     # SET (CMAKE_Fortran_FLAGS "-fexternal-blas")
     #libs="" #    -L/usr/lib/x86_64-linux-gnu -llapack -lblas -lhdf5
     #libs="$libs -llapack -lblas "
     string(APPEND CMAKE_Fortran_FLAGS " -fexternal-blas ")
    
    
#else()
     #   extra=$( ls ../src/extLibs/*.{f,F})
 
endif()
  

if (quad)
    SET (CMAKE_Fortran_FLAGS "-freal-8-real-16")
endif()

#MKLROOT="/usr/lib/x86_64-linux-gnu"

ADD_EXECUTABLE(terms         src/eps.f90
			      src/swav.f90		             
                             src/sphmsv.f90
                             src/miet.f90
                             src/linalg.f90
                             src/multiscat.f90
                             src/termsProgram.f90
                             src/HDFfive.f90
                             src/extLibs/toms644.f
			)

target_link_libraries(${PROJECTNAME} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${HDF5_Fortran_LIBRARIES} ${ext_Fortran_LIBRARIES})





