function err = test_piro_band_error (m, n, err1, err2, pr)
%TEST_PIRO_BAND_ERROR check and optionally print the error
err = max (err1, err2) ;
if (pr)
    fprintf ('A %3d x %3d: err1 = %15g, err2 = %15g\n',  m, n, err1, err2) ;
elseif (n > 100)
    fprintf ('.') ;
end
if (err > 1e-12)
    error ('test failed!') ;
end

