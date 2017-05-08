# !/bin/sh
# PBS -l select=32:ncpus=2:ngpus=2
### Specify "wallclock time" required for this job, hhh:mm:ss
# PBS -l walltime=01:00:00

# After job starts, must goto working directory. 
# $PBS_O_WORKDIR is the directory from where the job is fired. 
echo "==============================="
echo $PBS_JOBID
cat $PBS_NODEFILE
echo "==============================="
cd $PBS_O_WORKDIR
#job 

ime mpirun -np 2 ./main 1 < input.txt > output.txt

#NOTE
# The job line is an example : users need to change it to suit their applications
# The PBS select statement picks n nodes each having m free processors
# OpenMPI needs more options such as $PBS_NODEFILE
