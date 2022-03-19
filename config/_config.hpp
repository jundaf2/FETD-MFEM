// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_CONFIG_HEADER
#define MFEM_CONFIG_HEADER

// MFEM version: integer of the form: (major*100 + minor)*100 + patch.
#define MFEM_VERSION 40301

// MFEM version string of the form "3.3" or "3.3.1".
#define MFEM_VERSION_STRING "4.3.1"

// MFEM version type, see the MFEM_VERSION_TYPE_* constants below.
#define MFEM_VERSION_TYPE ((MFEM_VERSION)%2)

// MFEM version type constants.
#define MFEM_VERSION_TYPE_RELEASE 0
#define MFEM_VERSION_TYPE_DEVELOPMENT 1

// Separate MFEM version numbers for major, minor, and patch.
#define MFEM_VERSION_MAJOR ((MFEM_VERSION)/10000)
#define MFEM_VERSION_MINOR (((MFEM_VERSION)/100)%100)
#define MFEM_VERSION_PATCH ((MFEM_VERSION)%100)

// The absolute path of the MFEM source prefix
#define MFEM_SOURCE_DIR "/home/ubuntu/Desktop/FETD/mfem"

// The absolute path of the MFEM installation prefix
#define MFEM_INSTALL_DIR "/home/ubuntu/Desktop/FETD/mfem/mfem"

// Description of the git commit used to build MFEM.
#define MFEM_GIT_STRING "heads/master-0-gad60dac701b2ab36300a0e8eccd08e51b043f9c3-dirty"

// Build the parallel MFEM library.
// Requires an MPI compiler, and the libraries HYPRE and METIS.
// #define MFEM_USE_MPI

// Enable debug checks in MFEM.
// #define MFEM_DEBUG

// Throw an exception on errors.
// #define MFEM_USE_EXCEPTIONS

// Enable zlib in MFEM.
// #define MFEM_USE_ZLIB

// Enable backtraces for mfem_error through libunwind.
// #define MFEM_USE_LIBUNWIND

// Enable MFEM features that use the METIS library (parallel MFEM).
// #define MFEM_USE_METIS

// Enable this option if linking with METIS version 5 (parallel MFEM).
// #define MFEM_USE_METIS_5

// Use LAPACK routines for various dense linear algebra operations.
// #define MFEM_USE_LAPACK

// Use thread-safe implementation. This comes at the cost of extra memory
// allocation and de-allocation.
// #define MFEM_THREAD_SAFE

// Enable the OpenMP backend.
// #define MFEM_USE_OPENMP

// [Deprecated] Enable experimental OpenMP support. Requires MFEM_THREAD_SAFE.
// #define MFEM_USE_LEGACY_OPENMP

// Internal MFEM option: enable group/batch allocation for some small objects.
// #define MFEM_USE_MEMALLOC

// Which library functions to use in class StopWatch for measuring time.
// For a list of the available options, see INSTALL.
// If not defined, an option is selected automatically.
#define MFEM_TIMER_TYPE 2

// Enable MFEM functionality based on the SUNDIALS libraries.
// #define MFEM_USE_SUNDIALS

// Enable MFEM functionality based on the Mesquite library.
// #define MFEM_USE_MESQUITE

// Enable MFEM functionality based on the SuiteSparse library.
// #define MFEM_USE_SUITESPARSE

// Enable MFEM functionality based on the SuperLU library.
// #define MFEM_USE_SUPERLU
// #define MFEM_USE_SUPERLU5

// Enable MFEM functionality based on the MUMPS library.
// #define MFEM_USE_MUMPS
// #define MFEM_MUMPS_VERSION @MFEM_MUMPS_VERSION@

// Enable MFEM functionality based on the STRUMPACK library.
// #define MFEM_USE_STRUMPACK

// Enable MFEM features based on the Ginkgo library
// #define MFEM_USE_GINKGO

// Enable MFEM functionality based on the AmgX library.
// #define MFEM_USE_AMGX

// Enable secure socket streams based on the GNUTLS library
// #define MFEM_USE_GNUTLS

// Enable Sidre support
// #define MFEM_USE_SIDRE

// Enable the use of SIMD in the high performance templated classes
// #define MFEM_USE_SIMD

// Enable FMS support
// #define MFEM_USE_FMS

// Enable Conduit support
// #define MFEM_USE_CONDUIT

// Enable functionality based on the NetCDF library (reading CUBIT files)
// #define MFEM_USE_NETCDF

// Enable functionality based on the PETSc library
// #define MFEM_USE_PETSC

// Enable functionality based on the SLEPc library
// #define MFEM_USE_SLEPC

// Enable functionality based on the MPFR library.
// #define MFEM_USE_MPFR

// Enable MFEM functionality based on the PUMI library
// #define MFEM_USE_PUMI

// Enable MFEM functionality based on the HIOP library.
// #define MFEM_USE_HIOP

// Enable MFEM functionality based on the GSLIB library
// #define MFEM_USE_GSLIB

// Build the NVIDIA GPU/CUDA-enabled version of the MFEM library.
// Requires a CUDA compiler (nvcc).
// #define MFEM_USE_CUDA

// Build the AMD GPU/HIP-enabled version of the MFEM library.
// Requires a HIP compiler (hipcc).
// #define MFEM_USE_HIP

// Enable functionality based on the RAJA library.
// #define MFEM_USE_RAJA

// Enable functionality based on the OCCA library.
// #define MFEM_USE_OCCA

// Enable functionality based on the libCEED library.
// #define MFEM_USE_CEED

// Enable functionality based on the Caliper library.
// #define MFEM_USE_CALIPER

// Enable functionality based on the Umpire library.
// #define MFEM_USE_UMPIRE

// Enable IO functionality based on the ADIOS2 library.
// #define MFEM_USE_ADIOS2

// Version of HYPRE used for building MFEM.
// #define MFEM_HYPRE_VERSION @MFEM_HYPRE_VERSION@

// Macro defined when PUMI is built with support for the Simmetrix SimModSuite
// library.
// #define MFEM_USE_SIMMETRIX

// Enable interface to the MKL CPardiso library.
// #define MFEM_USE_MKL_CPARDISO

#endif // MFEM_CONFIG_HEADER
