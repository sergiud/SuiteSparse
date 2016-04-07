%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fail = gp_getFailVector(method, experiment, ew)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Builds the cases where the method fails %%%
%

% DATA_SEGMENT : Store the ids of the offending problems:

metis4Fail_sqsym_noew = [ ] ;
metis4Fail_sqsym_ew = [31, 50, 206, 340, 407, 408, 414, 415, 419, 427, 455, 459, 460, 549, 759, 956, 1337, 1439] ;
metis4Fail_bipart_noew = [ ] ;
metis4Fail_bipart_ew = [ ] ;
metis4Fail_aat_noew = [ ] ;
metis4Fail_aat_ew = [ ] ;

metis5Fail_sqsym_noew = [ ] ;
metis5Fail_sqsym_ew = [ ] ;
metis5Fail_bipart_noew = [ ] ;
metis5Fail_bipart_ew = [ ] ;
metis5Fail_aat_noew = [ ] ;
metis5Fail_aat_ew = [ ] ;

gpFail_sqsym_noew = [ ] ;
gpFail_sqsym_ew = [ ] ;
gpFail_bipart_noew = [ ] ;
gpFail_bipart_ew = [ ] ;
gpFail_aat_noew = [ ] ;
gpFail_aat_ew = [ ] ;

% Configure fail vector based on method and experiment.
switch method
    case 1
        switch experiment
            case 1
                fail = metis4Fail_sqsym_noew ;
                if ew == 1, fail = horzcat(fail, metis4Fail_sqsym_ew); end
            case 2
                fail = metis4Fail_bipart_noew ;
                if ew == 1, fail = horzcat(fail, metis4Fail_bipart_ew) ; end
            case 3
                fail = metis4Fail_aat_noew ;
                if ew == 1, fail = horzcat(fail, metis4Fail_aat_ew) ; end
        end

    case 2
        switch experiment
            case 1
                fail = metis5Fail_sqsym_noew ;
                if ew == 1, fail = horzcat(fail, metis5Fail_sqsym_ew); end
            case 2
                fail = metis5Fail_bipart_noew ;
                if ew == 1, fail = horzcat(fail, metis5Fail_bipart_ew) ; end
            case 3
                fail = metis5Fail_aat_noew ;
                if ew == 1, fail = horzcat(fail, metis5Fail_aat_ew) ; end
        end

    otherwise
        switch experiment
            case 1
                fail = gpFail_sqsym_noew ;
                if ew == 1, fail = horzcat(fail, gpFail_sqsym_ew); end
            case 2
                fail = gpFail_bipart_noew ;
                if ew == 1, fail = horzcat(fail, gpFail_bipart_ew) ; end
            case 3
                fail = gpFail_aat_noew ;
                if ew == 1, fail = horzcat(fail, gpFail_aat_ew) ; end
        end
end



