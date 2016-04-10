%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gp_printDetails(k, id, experiment, methodID, ew)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format compact ;

c = clock ;
gp_log(sprintf('%d/%d/%d %d:%d | ', c(3), c(2), c(1), c(4), c(5))) ;
gp_log(sprintf('k: %4d | ', k)) ;
gp_log(sprintf('id: %4d ', id)) ;
gp_log('\n') ;

