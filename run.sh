# Arg1 is the input file name with the matrix in COO formatand Arg2 is the output file 
time mpirun -np 8 ./main $1 $2
