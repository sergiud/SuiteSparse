function cs_install (do_pause)
%CS_INSTALL: compile and install CXSparse for use in MATLAB.
%   Your current working directory must be CXSparse/MATLAB in order to use this
%   function.
%   
%   The directories
%
%       CXSparse/MATLAB/CSparse
%       CXSparse/MATLAB/Demo
%       CXSparse/MATLAB/ssget
%
%   are added to your MATLAB path (see the "pathtool" command to add these to
%   your path permanently, for future MATLAB sessions).
%
%   Next, the MATLAB CXSparse demo program, CXSparse/MATLAB/cs_demo is executed.
%   To run the demo with pauses so you can see the results, use cs_install(1).
%   To run the full MATLAB test programs for CXSparse, run testall in the
%   Test directory.
%
%   Example:
%       cs_install          % install and run demo with no pauses
%       cs_install(1)       % install and run demo with pauses
%
%   See also: cs_demo
%
% CXSparse, Copyright (c) 2006-2022, Timothy A. Davis. All Rights Reserved.
% SPDX-License-Identifier: LGPL-2.1+

fprintf ('Compiling and installing CXSparse\n') ;
if (nargin < 1)
    do_pause = 0 ;
end

if (do_pause)
    input ('Hit enter to continue: ') ;
end
addpath ([pwd '/CSparse']) ;
addpath ([pwd '/Demo']) ;

if (verLessThan ('matlab', '8.4'))
    fprintf ('ssget not installed (MATLAB 8.4 or later required)\n') ;
else
    % install ssget, unless it's already in the path
    try
        % if this fails, then ssget is not yet installed
        index = ssget ;
        fprintf ('ssget already installed:\n') ;
        which ssget
    catch
        index = [ ] ;
    end
    if (isempty (index))
        % ssget is not installed.  Use ./ssget
        fprintf ('Installing ./ssget\n') ;
        try
            addpath ([pwd '/ssget']) ;
        catch me
            disp (me.message) ;
            fprintf ('ssget not installed\n') ;
        end
    end
end

cd ('CSparse') ;
cs_make (1) ;
cd ('../Demo') ;
cs_demo (do_pause)
