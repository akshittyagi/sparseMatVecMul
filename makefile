CC = nvcc
MPIFLAGS = -I/home/soft/mpich-3.1.4/include/ -L/home/soft/mpich-3.1.4/lib/ 
CFLAGS = -lmpi
EXECNAME = main
INFILENAME = input_1024.txt
OUTFILENAME = output.txt
NUMPROCS = 8

all:  parallel

parallel:
	$(CC) $(MPIFLAGS) $(CFLAGS) mat_vec.cu -o $(EXECNAME)

hpcCompile:
	module load compiler/cuda/7.5/compilervars
	module load compiler/gcc/4.9.3/compilervars
	module load mpi/mpich/3.1.4/gcc/mpivars

run: 
	mpirun -np $(NUMPROCS) ./$(EXECNAME) $(INFILENAME) $(OUTFILENAME) 

clean:
	rm -f *~ *.o parallel
	rm $(EXECNAME)