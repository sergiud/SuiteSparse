function meshdemo(whichdemos)
% MESHDEMO : Demo of mesh partitioning toolkit.
%
% meshdemo(whichdemos);
% The argument is a demo or a list of demos to run.  
% The default is to run them all.
%
% 1:  Various partitioning methods on the Tapir mesh.
% 2:  Dicing a triangular grid into 16 pieces.
% 3:  Vertex separators in graphs.
% 4:  Nested dissection.
%
% John Gilbert, Xerox PARC, August 1994.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.
%
% Modified by Tim Davis, for Matlab 5.1.  July 6, 1998
% Modified by JRG 5-14-01
% Modified by JRG 3 Jul 01 to add Metis to demo 1
% Modified by JRG Feb 02 for Matlab 6

if nargin < 1
    whichdemos = [1 2 4];
end;

clc;
format compact;

disp('          *********************************************')
disp('          *** Mesh Partitioning and Separator Demos ***');
disp('          *********************************************')
disp(' ');
disp('        HELP COPYRIGHT for copyright and licensing notice.');
disp(' ');
disp(' ');
disp(' The file "meshes.mat" contains some sample meshes with coordinates.');
disp(' ');
load meshes;
whos;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(whichdemos==1) % Various partitioning methods on the Tapir mesh.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp('          *********************************************')
disp('          ***      Various Partitioning Methods     *** ');
disp('          *********************************************')
disp(' ');
disp(' ');
    
figure(1);
clf reset;
colordef(1,'black')

disp(' "Tapir" is a test of a no-obtuse-angles mesh generation algorithm');
disp(' due to Bern, Mitchell, and Ruppert.  ');

disp('  ');
disp('gplotg(Tapir,Txy);');
disp('  ');

gplotg(Tapir,Txy);

disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' We''ll try a few different partitioners on the mesh:');

disp(' ');
disp(' First is spectral partitioning, which uses the second eigenvector of');
disp(' the Laplacian matrix of the graph, also known as the "Fiedler vector".');

disp('  ');
disp('specpart(Tapir,Txy);');
specpart(Tapir,Txy);
disp('cutsize(Tapir,ans)');
cutsize(Tapir,ans)
disp('  ');


disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' Second is geometric partitioning, which uses the coordinates to find');
disp(' a separating line or circle in the plane.');

disp('  ');
disp('geopart(Tapir,Txy);');
geopart(Tapir,Txy);
disp('cutsize(Tapir,ans)');
cutsize(Tapir,ans)
disp('  ');


disp(' Hit space to continue ...');
pause;

usechaco = 0;
if usechaco
disp(' ');
disp(' Next is multilevel Kernighan-Lin partitioning, which first partitions');
disp(' a coarsened graph and then improves the split by swapping vertices at');
disp(' several levels of refinement.  (This will only work if you have ');
disp(' Hendrickson and Leland''s "Chaco" package and its Matlab interface.)');

disp('  ');
disp('chaco(Tapir,Txy);');
chaco(Tapir,Txy);
disp('cutsize(Tapir,ans)');
cutsize(Tapir,ans)
disp('  ');

disp(' Hit space to continue ...');
pause;
end % if usechaco

usemetis = 1;
if usemetis
disp(' ');
disp(' Next is a multilevel method from the "Metis" package.');
disp(' (This will only work if you have Metis and its Matlab interface.)');

disp('  ');
disp('metispart(Tapir,Txy);');
metispart(Tapir,Txy);
disp('cutsize(Tapir,ans)');
cutsize(Tapir,ans)
disp('  ');

disp(' Hit space to continue ...');
pause;
end % if usemetis

disp(' ');
disp(' The final method is geometric spectral partitioning, which uses');
disp(' the geometric algorithm but in a spectral coordinate system.');

disp('  ');
disp('gspart(Tapir,Txy);');
gspart(Tapir,Txy);
disp('cutsize(Tapir,ans)');
cutsize(Tapir,ans)
disp('  ');

disp(' The first figure is the spectral coordinate system; ');
disp(' the second is the original mesh coordinates.');
disp('  ');

disp(' Hit space to continue ...');
pause;
disp('  ');

end % whichdemos(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(whichdemos==2) % Dicing a triangular grid into 16 pieces.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp('          *********************************************')
disp('          ***          Recursive Bisection          *** ');
disp('          *********************************************')
disp(' ');

figure(1);
clf reset;
colordef(1,'black')

disp(' ');
disp(' Function "gridt" produces a triangular grid:');
disp(' ');

disp('[A,xy] = gridt(20);');
[A,xy] = gridt(20);

gplotg(A,xy);
nvtx = size(A,1);
nedge = (nnz(A)-nvtx)/2;
xlabel([int2str(nvtx) ' vertices, ' int2str(nedge) ' edges'],'visible','on');

disp(' ');
disp(' (See also grid5, grid7, grid9, grid3d, grid3dt.)');
disp(' ');

disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' Use geometric partitioning to dice the grid into sixteenths.');
disp(' ');

disp('geodice(A,4,xy);');
geodice(A,4,xy);

disp(' ');
disp(' Hit space to continue ...');
pause;
disp('  ');

end % whichdemos(2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(whichdemos==3) % Vertex separators in graphs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp('          *********************************************')
disp('          ***           Vertex Separators           *** ');
disp('          *********************************************')
disp(' ');

figure(1);
clf reset;
colordef(1,'black')

disp(' ');
disp(' A vertex separator is a set of vertices whose removal partitions the graph.');
disp(' We can convert a vertex partition (which is an edge separator) into a');
disp(' vertex separator by finding a minimum cover with bipartitite matching.');
disp(' ');

disp('part = gspart(Eppstein);');

part = gspart(Eppstein);

disp(' ');
disp(' Here''s the vertex partition (edge separator):');
disp(' ');

disp('gplotpart(Eppstein,Exy,part)');

gplotpart(Eppstein,Exy,part)
title('Edge Separator');

disp(' ');
disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' And here''s the vertex separator:');
disp(' ');

disp('sep = vtxsep(Eppstein,part)');

sep = vtxsep(Eppstein,part)

figure(2);
clf reset;
colordef(2,'black')

disp('highlight(Eppstein,Exy,sep);');

highlight(Eppstein,Exy,sep);
title('Vertex Separator');
drawnow;

disp(' ');
disp(' How big is the separator?');
disp(' ');

disp('SeparatorSize = length(sep)');

SeparatorSize = length(sep)
xlabel([int2str(SeparatorSize), ' separator vertices'],'visible','on');

disp(' ');
disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' And how big are the pieces?');
disp(' ');

disp('nonsep = other(sep,Eppstein);');
disp('blocks = components(Eppstein(nonsep,nonsep));');
disp('BlockSizes = hist(blocks,max(blocks))');

nonsep = other(sep,Eppstein);
blocks = components(Eppstein(nonsep,nonsep));
BlockSizes = hist(blocks,max(blocks))

disp(' ');
disp(' Hit space to continue ...');
pause;
disp('  ');

end % whichdemos(3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(whichdemos==4) % Nested dissection.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp('          *********************************************')
disp('          ***           Nested Dissection           *** ');
disp('          *********************************************')
disp(' ');

figure(1);
clf reset;
colordef(1,'black')

disp(' ');
disp(' Nested dissection is a recursively defined sparse matrix permutation');
disp(' based on vertex separators in the graph of a symmetric matrix.');
disp(' Any partitioning algorithm can be used; here we illustrate nested');
disp(' dissection with spectral partitioning.');
disp(' ');
disp(' Here''s the original matrix and its graph...');
disp(' ');

disp('gplotg(Smallmesh,Sxy);');

gplotg(Smallmesh,Sxy);
title('Finite Element Mesh');
figure(2); 
clf reset;
colordef(2,'black')

disp('spy(Smallmesh,3);');

spy(Smallmesh,3);
title('Coefficient Matrix');

disp(' ');
disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' Here''s the result of factoring the matrix with Gaussian elimination.');

disp(' ');
disp('spy(chol(Smallmesh),3);');

spy(chol(Smallmesh),3);
title('Cholesky Factor');

disp(' ');
disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' Now we use nested dissection to reorder the matrix.');

disp(' ');
disp('nd = specnd(Smallmesh);');

nd = specnd(Smallmesh);

disp(' ');
disp(' Permute the rows and columns of the matrix with Matlab''s indexing.');

disp(' ');
disp('spy(Smallmesh(nd,nd),3);');

figure(1)
spy(Smallmesh(nd,nd),3);
title('Coefficient Matrix with Nested Dissection');

disp(' ');
disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' The permuted matrix has a sparser Cholesky factor than the original.');

disp(' ');
disp('spy(chol(Smallmesh(nd,nd)),3);');

figure(1);
spy(chol(Smallmesh(nd,nd)),3);
title('Cholesky Factor with Nested Dissection');

disp(' ');
disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' The elimination tree, which has one vertex per matrix column,');
disp(' shows the recursive structure of the nested dissection ordering.');
disp(' The height of the tree is related to the complexity of factoring');
disp(' the matrix on a parallel computer.');

disp(' ');
disp('etreeplotg(Smallmesh(nd,nd),1);');

figure(2);
etreeplotg(Smallmesh(nd,nd),1);
title('Elimination Tree with Nested Dissection');

disp(' ');
disp(' Hit space to continue ...');
pause;

if exist('symamd')  % symamd is new in Matlab 6, use symmmd on earlier versions
	disp(' ');
	disp(' Another ordering is minimum degree, which often gives higher trees');
	disp(' but lower fill than nested dissection.  Say "help symamd" for details.');
	disp(' Here''s a summary of the complexity on this matrix:');
	md = symamd(Smallmesh);
else
    disp(' ');
	disp(' Another ordering is minimum degree, which often gives higher trees');
	disp(' but lower fill than nested dissection.  Say "help symmmd" for details.');
	disp(' Here''s a summary of the complexity on this matrix:');
	md = symmmd(Smallmesh);
end;

disp(' ');
disp('With no permutation:');
analyze(Smallmesh);

disp(' ');
disp('With nested dissection:');
analyze(Smallmesh(nd,nd));

disp(' ');
disp('With minimum degree:');
analyze(Smallmesh(md,md));

disp(' ');
disp(' Hit space to continue ...');
pause;

end % whichdemos(4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final notes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(' ');
disp(' ');
disp(' ');
disp('  **********************************************************************')
disp('  *** For more information, see meshpart.html or say "help meshpart" ***');
disp('  **********************************************************************')
disp(' ');

format;
