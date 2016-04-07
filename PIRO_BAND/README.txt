PIRO_BAND.  Version 1.0:  Pipelined Rotations for Band reduction.
Copyright (C) 2009-2014, Sivasankaran Rajamanickam, Timothy A. Davis.
PIRO_BAND is also available under other licenses; contact authors for details.
http://www.cise.ufl.edu/research/sparse

This package requires the following other packages:

    * ../SuiteSparse_config Version 4.3.0 or later (a component of SuiteSparse)
    * LAPACK at http://www.netlib.org/lapack/
    * BLAS at http://netlib.org/blas/

This software package reduces a band matrix to bidiagonal form via pipelined
Givens rotations, a precursor to computing the SVD of a band matrix.
See the Doc/ subdirectory for details.

TODO proofread user guide
TODO rename Tcov/mem.h and Tcov/memory.h to Tcov/piro_band_cover_memory.[ch]
TODO add ACM TOMS paper draft to Doc/ (when it's clean), in tech report format

