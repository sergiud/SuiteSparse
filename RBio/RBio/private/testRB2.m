function testRB2
%testRB2: test the RBio toolbox.  ssget is required.
%
% Example:
%   testRB2
%
% See also ssget, RBread, RBreade, testRB1.

% RBio, Copyright (c) 2009-2022, Timothy A. Davis.  All Rights Reserved.
% SPDX-License-Identifier: GPL-2.0+

Problem = ssget ('Meszaros/farm') ;
% disp (Problem) ;
A = RBread ('farm.rb') ;
if (~isequal (A, Problem.A))
    error ('test failure: farm.rb') ;
end
% fprintf ('mtype: %s\n', RBtype (A)) ;
mtype = RBtype (A) ;
if (any (mtype ~= 'ira'))
    error ('test failure: farm.rb') ;
end

Problem = ssget ('HB/bcsstk01') ;
% disp (Problem) ;
A = RBread ('bcsstk01.rb') ;
if (~isequal (A, Problem.A))
    error ('test failure: bcsstk01.rb') ;
end
% fprintf ('mtype: %s\n', RBtype (A)) ;
mtype = RBtype (A) ;
if (any (mtype ~= 'rsa'))
    error ('test failure: bcsstk01.rb') ;
end

Problem = ssget ('HB/lap_25') ;
% disp (Problem) ;
A = RBread ('lap_25.rb') ;
if (~isequal (A, Problem.A))
    error ('test failure: lap_25.rb') ;
end
A = RBreade ('lap_25.pse') ;
if (~isequal (A, Problem.A))
    error ('test failure: lap_25.pse') ;
end
% fprintf ('mtype: %s\n', RBtype (A)) ;
mtype = RBtype (A) ;
if (any (mtype ~= 'psa'))
    error ('test failure: bcsstk01.rb') ;
end

Problem = ssget ('HB/west0479') ;
% disp (Problem) ;
[A Z] = RBread ('west0479.rb') ;
if (~isequal (A, Problem.A))
    error ('test failure: west0479.rb') ;
end
if (~isequal (Z, Problem.Zeros))
    error ('test failure: west0479.rb') ;
end
[A Z] = RBread ('west0479.rua') ;
if (~isequal (A, Problem.A))
    error ('test failure: west0479.rua') ;
end
if (~isequal (Z, Problem.Zeros))
    error ('test failure: west0479.rua') ;
end
% fprintf ('mtype: %s\n', RBtype (A)) ;
mtype = RBtype (A) ;
if (any (mtype ~= 'rua'))
    error ('test failure: west0479.rua') ;
end

fprintf ('RB test 2: passed\n') ;

