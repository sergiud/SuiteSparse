A set of quick tests for the C-callable PIRO_BAND functions.  To run these
tests, first compile the libpiro_band.a archive by typing 'make' in ../Source
(or 'make' in the top-level PIRO_BAND directory).  Then type 'make' in this
directory.

TODO double-check this output
Sample output:

./test_piro_band
------------- Tests for piro_band_reduce ------- 
------------- Tests for double/real/integer ------- 
m = 10, n= 10, bl=4, bu=4
C - U' ----Norm is 0
A - U * B * V' ----Norm is 4.55022e-16
------------- Tests for double/complex/integer ------- 
m = 10, n= 10, bl=4, bu=4
C - U' ----Norm is 0
A - U * B * V' ----Norm is 6.27228e-16
------------- Tests for double/real/long ------- 
m = 10, n= 10, bl=4, bu=4
C - U' ----Norm is 0
A - U * B * V' ----Norm is 6.12047e-16
------------- Tests for double/complex/long ------- 
m = 10, n= 10, bl=4, bu=4
C - U' ----Norm is 0
A - U * B * V' ----Norm is 5.89239e-16
------------- Tests for float/real/integer ------- 
m = 10, n= 10, bl=4, bu=4
C - U' ----Norm is 0
A - U * B * V' ----Norm is 2.59165e-07
------------- Tests for float/complex/integer ------- 
m = 10, n= 10, bl=4, bu=4
C - U' ----Norm is 0
A - U * B * V' ----Norm is 3.01721e-07
------------- Tests for float/real/long ------- 
m = 10, n= 10, bl=4, bu=4
C - U' ----Norm is 0
A - U * B * V' ----Norm is 3.29668e-07
------------- Tests for float/complex/long ------- 
m = 10, n= 10, bl=4, bu=4
C - U' ----Norm is 0
A - U * B * V' ----Norm is 2.55346e-07
PIRO_BAND:  all tests passed
