%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gp_printExecDetails(methodID, experiment, ew)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gp_log('\n') ;
gp_log('===============================================================================') ;
gp_log('\n') ;

gp_log('Method: ') ;
switch methodID
case 1
    gp_log('METIS 4.0 - Default') ;
case 2
    gp_log('METIS 5.0.2 - Default') ;
otherwise
    O = gp_generateOptions(methodID);
    gp_log(O.methodDesc) ;
end
gp_log(' | ') ;

SQSYM  = 1 ;
BIPART = 2 ;
AAT    = 3 ;
gp_log('Experiment: ') ;
switch experiment
    case SQSYM
        gp_log('SQSYM') ;
    case BIPART
        gp_log('BIPART') ;
    case AAT
        gp_log('AAT') ;
    otherwise
        gp_log('Unknown experiment') ;
end
gp_log(' ') ;
switch ew
    case 0
        gp_log('Ignore Edge Weights') ;
    case 1
        gp_log('Consider Edge Weights') ;
end

gp_log('\n') ;
gp_log('===============================================================================') ;
gp_log('\n') ;


