METIS 5.1.0, with minor modifications by Tim Davis
to incorporate it into SuiteSparse.

This copy of METIS is slightly changed from the original METIS v5.1.0
distribution.  The use of METIS in SuiteSparse is optional, but if used, this
revised version is required.

(1) In metis-5.1.0/include/metis.h, the default integer size has been changed
    from 32 bits to 64 (IDXTYPEWIDTH).  METIS 5.1.0 gives this flexility to the
    user, asking the user to modify this file.  That has been done here, and as
    a result, this file is renamed to SuiteSparse_metis.h.  Getting the
    unmodified libmetis.so in a Linux distro (likely with 32-bit integers)
    combined with a modified metis.h (with 64-bit integers) breaks things
    badly.  So the safest thing is to rename this file as SuiteSparse_metis.h,
    and to rename the compiled library as libsuitesparse_metis.so, to ensure
    the right library is linked.

(3) The files metis-5.1.0/GKlib/GKLib.h and metis-5.1.0/GKlib/memory.c have
    been modified to disable the signal-handing in METIS.  METIS runs out of
    memory, but they break MATLAB (you will get a segfault).  This change is
    essential if METIS is to be used in MATLAB.

(4) The abs and iabs functions in the original metis-5.1.0/libmetis/parmetis.c
    and metis-5.1.0/libmetis/balance.c give compiler warnings when IDXTYPEWIDTH
    is 64, so they have been replaced with a type-agnostic macro, ABS.  This is
    just a compiler warning, so the fix is optional.

(5) Minor formatting changes have been made to avoid compiler warnings of
    misleading indentation (getopt.c, csr.c).

(6) The malloc/calloc/realloc/free functions have been replaced with
    SuiteSparse_config.(malloc/calloc/realloc/free) throughout.  The gkmcore
    feature is disabled since it can conflict with the use of mxMalloc
    in the MATLAB interface to SuiteSparse.  See GKlib/memory.c.

(7) All original files from metis-5.1.0 that have been modified here have been
    placed in ./include/original, ./libmetis/original, or ./GKlib/original

Tim Davis, Oct 31, 2022, Texas A&M University Any changes made by Tim Davis are
released to the original copyright holder, under the original Apache-2.0
license of METIS.

METIS, Copyright 1995-2013, Regents of the University of Minnesota.
Author: George Karypis
SPDX-License-identifier: Apache-2.0

