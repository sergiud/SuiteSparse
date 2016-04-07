TODO write a proper REAMME, like PIRO_BAND

-------------------------------------------------------------------------------
Compiling SkylineSVD
-------------------------------------------------------------------------------

1. "make all" in Source and MATLAB directories will compile the libraries. The
SuiteSparse_config and PIRO_BAND packages are expected in the same level as
this SkylineSVD directory.

2. test_blksky, test_blksky_symbolic tests the numerical R-bidiaonalization
and its symbolic factorization from within MATLAB. The UFget package is
expected to be in the path.

3. blksky_symbolic can be used te get the minimum top row for each column of R.

4. blk_sky will do the numerical R-bidiagonalization. 

The m-files have the usage for both blk_sky and blksky_symbolic. 

5. test_blksky_perf will compare our code against MATLAB's sparse svd. (I need
to add a call to diagonalize the bidiagonal matrix. I haven't done that yet).


