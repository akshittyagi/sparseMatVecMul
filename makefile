CC = nvcc
MPIFLAGS = -I/home/soft/mpich-3.1.4/include/ -L/home/soft/mpich-3.1.4/lib/ 
CFLAGS = -lmpi
EXECNAME = main
INFILENAME = input_1024.txt
OUTFILENAME = output.txt
NUMPROCS = 8

all: parallel

parallel:
	$(CC) $(MPIFLAGS) $(CFLAGS) mat_vec.cu -o $(EXECNAME)

run: parallel
	mpirun -np $(NUMPROCS) ./$(EXECNAME) $(INFILENAME) $(OUTFILENAME) 

clean:
	rm -f *~ *.o parallel
	rm $(EXECNAME)