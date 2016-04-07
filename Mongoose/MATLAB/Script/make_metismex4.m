cd ~/SuiteSparse/Mongoose/MATLAB/Script/meshpart;
mex -largeArrayDims -I/home/nuri/SuiteSparse/metis-4.0/Lib/ -L/home/nuri/SuiteSparse/metis-4.0/ -lmetis metismex4.c;
