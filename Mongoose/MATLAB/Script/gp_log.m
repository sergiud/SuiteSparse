function gp_log(message, echo)

path = '~/gp_logFile.txt' ;
mode = 'a' ;
logfile = fopen(path, mode) ;

% Possibly report to console
do_echo = true;
if nargin == 2, do_echo = echo; end
if do_echo, fprintf(1, message); end

% Definitely report to file
fprintf(logfile, message) ;

fclose(logfile) ;
