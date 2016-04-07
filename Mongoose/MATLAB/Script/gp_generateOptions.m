%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function O = generateOptions(method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If we're just generating a vanilla user params by calling this with no
% arguments then return the default user params.
if nargin == 0, method = 3 ; end

% Configure the UserParams based on the method we want to use.
O = gp_getDefaultOptions() ;

switch method

% FM and QP, one dance
    case 3
        O.methodDesc = 'Mongoose; Default' ;

    case 4
        O.methodDesc = 'Mongoose; G=PA, 1D' ;
        O.guessCutType = 1 ;

    case 5
        O.methodDesc = 'Mongoose; G=GP, 1D' ;
        O.guessCutType = 2 ;

% FM and QP, two dances
    case 6
        O.methodDesc = 'Mongoose; Default, 2D' ;
        O.numDances = 2 ;

    case 7
        O.methodDesc = 'Mongoose; G=PA, 2D' ;
        O.guessCutType = 1 ;
        O.numDances = 2 ;

    case 8
        O.methodDesc = 'Mongoose; G=GP, 2D' ;
        O.guessCutType = 2 ;
        O.numDances = 2 ;

% FM Only, one dance
    case 9
        O.methodDesc = 'Mongoose; Default, noQP' ;
        O.useQPGradProj = false ;

    case 10
        O.methodDesc = 'Mongoose; G=PA, noQP' ;
        O.guessCutType = 1 ;
        O.useQPGradProj = false ;

    case 11 % note that this option doesn't make any sense
        O.methodDesc = 'Mongoose; G=GP, noQP' ;
        O.guessCutType = 2 ;
        O.useQPGradProj = false ;

% QP Only, one dance
    case 12
	O.methodDesc = 'Mongoose; Default, noFM' ;
        O.useFM = false ;

    case 13
	O.methodDesc = 'Mongoose; G=PA, noFM' ;
        O.guessCutType = 1 ;
        O.useFM = false ;

    case 14
	O.methodDesc = 'Mongoose; G=GP, noFM' ;
        O.guessCutType = 2 ;
        O.useFM = false ;

% FM and QP, one dance, PAComm
    case 15
        O.methodDesc = 'Mongoose; PAComm' ;

    case 16
        O.methodDesc = 'Mongoose; G=PA, PAComm' ;
        O.guessCutType = 1 ;

    case 17
        O.methodDesc = 'Mongoose; G=GP, PAComm' ;
        O.guessCutType = 2 ;

    otherwise, O.methodDesc = 'Unknown methodID' ;
        % Use the default user params
end
