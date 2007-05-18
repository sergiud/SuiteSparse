function cs_install (do_pause)
%CS_INSTALL: compile and install CXSparse for use in MATLAB.
%   Your current working directory must be CXSparse/MATLAB in order to use this
%   function.
%   
%   The directories
%
%       CXSparse/MATLAB/CSparse
%       CXSparse/MATLAB/Demo
%       CXSparse/MATLAB/UFget
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
%   Copyright 2006-2007, Timothy A. Davis.
%   http://www.cise.ufl.edu/research/sparse

fprintf ('Compiling and installing CXSparse\n') ;
if (nargin < 1)
    do_pause = 0 ;
end

if (do_pause)
    input ('Hit enter to continue: ') ;
end
addpath ([pwd filesep 'CSparse']) ;
addpath ([pwd filesep 'Demo']) ;

v = getversion ;
if (v >= 7.0)
    addpath ([pwd filesep 'UFget']) ;
else
    fprintf ('UFget not installed (MATLAB 7.0 or later required)\n') ;
end

cd ('CSparse') ;
cs_make (1) ;
cd ('../Demo') ;
cs_demo (do_pause)

%-------------------------------------------------------------------------------
function [v,pc] = getversion
% determine the MATLAB version, and return it as a double.
% only the primary and secondary version numbers are kept.
% MATLAB 7.0.4 becomes 7.0, version 6.5.2 becomes 6.5, etc.
v = version ;
t = find (v == '.') ;
if (length (t) > 1)
    v = v (1:(t(2)-1)) ;
end
v = str2double (v) ;
try
    % ispc does not appear in MATLAB 5.3
    pc = ispc ;
catch
    % if ispc fails, assume we are on a Windows PC if it's not unix
    pc = ~isunix ;
end
