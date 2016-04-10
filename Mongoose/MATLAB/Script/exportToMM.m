function exportToMM(filename, matrix, ew)

format compact;

fid = fopen(filename, 'w') ;

% Output the header.
nnz = size(find(matrix), 1) ;
n = size(matrix, 1) ;
m = size(matrix, 2) ;
fprintf(fid, '%d %d %d\n', n, m, nnz) ;

% Because 'matrix' is symmetric,
% whether we iterate through i first or j first doesn't matter.
for i=1:size(matrix, 1)
    for j=1:size(matrix, 2)
        value = matrix(i,j)+0;
        fprintf(fid, '%d %d %f\n', i, j, value) ;
    end
end

fclose(fid);

gp_log(sprintf('Successfully wrote matrix to \"%s\"\n', filename), false) ;
