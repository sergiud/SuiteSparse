function kinds = sskinds
%SSKINDS get 'kind' of each problem in the SuiteSparse Matrix Collection.
%
% A problem 'kind' is the domain from which the problem arises.  Examples
% include 'circuit simulation', 'optimization', 'least squares', and 'robotics'.
% See ssget/README.txt for a complete list.
%
% A 'kind' is related to the matrix 'Group'.  The latter is determined by
% person, persons, or organization that created the matrix.  'Group' may also
% refer to another matrix collection that has been included in this collection
% (such as 'HB' for the Harwell/Boeing collection).
%
% kinds = sskinds returns a cell array with one row per matrix.  kinds{id}
% is the same as Problem.kind where Problem=ssget(id).  If an error occurs,
% an empty cell array is returned.
%
% Example:
%
%   kinds = sskinds ;
%   Problem = ssget (42)
%   kinds{42}
%
% results in the following output:
%
%   Problem = 
%        title: 'SYMMETRIC STIFFNEss MATRIX - FRAME WITHIN A SUSPENSION BRIDGE'
%            A: [485x485 double]
%         name: 'HB/bcsstk20'
%           id: 42
%         date: '1984'
%       author: 'J. Lewis'
%           ed: 'I. Duff, R. Grimes, J. Lewis'
%         kind: 'structural problem'
%
%   ans =
%       'structural problem'
%
% To find a list t of the id's all matrices from material science:
%
%   kinds = sskinds ;
%   t = zeros (1, length (kinds)) ;
%   for id = 1:length(kinds)
%       t (id) = ~isempty (strfind (kinds {id}, 'material')) ;
%   end
%   t = find (t) ;
%   for id = t
%       Problem = ssget (id)
%   end
%
% See also ssget.

% ssget, Copyright (c) 2009-2022, Timothy A. Davis. All Rights Reserved.
% SPDX-License-Identifier: BSD-3-clause

params = ssget_defaults ;
statfile = sprintf ('%sfiles/ssstats.csv', params.topdir) ;

try
    f = fopen (statfile, 'r') ;
    nmat = str2double (fgetl (f)) ;     % number of matrices
    s = fgetl (f) ;  %#ok               % date ssstats.csv was created (unused)
    kinds = cell (nmat,1) ;
    for id = 1:nmat
        s = fgetl (f) ;
        t = find (s == ',') ;
        t1 = t (11) ;
        t2 = t (12) ;
        r = s ((t1+1):(t2-1)) ;
        kinds {id} = r ;
    end
catch
    kinds = cell (0,0) ;                % return empty string on error
end

