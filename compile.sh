module load compiler/cuda/7.5/compilervars
module load compiler/gcc/4.9.3/compilervars
module load mpi/mpich/3.1.4/gcc/mpivars
nvcc -I/home/soft/mpich-3.1.4/include/ -L/home/soft/mpich-3.1.4/lib/ -lmpi mat_vec.cu -o main

