#include "Util.h"

#define THREADS_PER_BLOCK 1024

int N;
int comm_size, proc_Id;
int dim, num_rows;					//dim of data, number of rows the proc has

vector<int> indices, ptrs;
vector<int> data; 	//sparse matrix ka stuff
vector<int> vec; 					//vector
long long int *own_output, *output;
vector<int> own_mat_indices, own_mat_ptrs;
vector<int> own_mat_data; //own_matrix stuff

string out_file = "Output_";

__device__ long long int multRow(int noOfElems, int *colIndices,int *nonZeroElems,int *vecTOR)
{
	long long int sum = 0;
	for(int j = 0;j<noOfElems;j++)
	{
		long long int num1 = nonZeroElems[j];
		long long int num2 = vecTOR[colIndices[j]];
		sum+=(num1)*(num2);
	}
	return sum;
}

__global__ void multKernel(int *firstElemsRows,int *colIndices,int *nonZeroElems,int numRows,int *vecTOR, long long int *output)
{
	int currRow = blockIdx.x*blockDim.x + threadIdx.x;
	if(currRow<numRows)
	{
		int rowStart = firstElemsRows[currRow];
		int rowEnd = firstElemsRows[currRow+1];
		// if(currRow==numRows-1)
		// {
		// 	long long int lsum = 0;
		// 	for(int i=0;i<rowEnd-rowStart;i++)
		// 		{
		// 			lsum += ((long long int)nonZeroElems[i])*((long long int)vecTOR[colIndices[i]]);
		// 		}//printf("| %d VECTOR %d |",nonZeroElems[i],vecTOR[colIndices[i]]);
		// 	printf("%lld ",lsum);
		// }

		output[currRow] = multRow(rowEnd-rowStart,colIndices+rowStart,nonZeroElems+rowStart,vecTOR);
		
		// if(currRow==numRows-1)
		// printf("%lld\n",output[currRow]);
	}

}

void getInput(char *in_file) {
	if (proc_Id == 0) {
		//get input in proc 0
		ifstream f_in;
		f_in.open(in_file);
		//headers, dim and stuff
		string junk, not_junk;
		string temp;
		int data_item;
		int x, y;
		f_in >> junk >> not_junk;
		f_in >> junk >> dim >> not_junk;
		//actually take in the matrix
		// for (int i = 0; i < dim; i++) {
		// 	bool row_started = false;
		// 	for (int j = 0; j < dim; j++) {
		// 		f_in >> data_item;
		// 		if (data_item != 0) {
		// 			data.push_back(data_item);
		// 			indices.push_back(j);
		// 			if (!row_started) {
		// 				cout << i << endl;
		// 				ptrs.push_back(indices.size() - 1);
		// 				row_started = true;
		// 			}
		// 		}
		// 	}
		// }
		f_in >> temp;
		int xold = -1;
		while (temp[0] != 'B') {
			x = atoi(temp.c_str());
			f_in >> y >> data_item;
			data.push_back(data_item);
			indices.push_back(y);
			if (x != xold) {
				int diff = x - xold - 1;
				while (diff--) {
					ptrs.push_back(indices.size() - 1);
				}
				
				ptrs.push_back(indices.size() - 1);
				xold = x;
			}
			f_in >> temp;
		}
		int endIndex = ptrs.size();
		int differ = dim-endIndex;
		while(differ--)
			ptrs.push_back(indices.size());
		// cout << ptrs.size() << endl;
		//take in vector
		// f_in >> not_junk;
		vec.resize(dim);
		for (int i = 0; i < dim; i++) {
			f_in >> data_item;
			vec[i] = data_item;
		}
		f_in.close();
		//pick up left over rows for proc 0
		num_rows += dim % comm_size;
	}
	//tell everyone about their load
	MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//send the vector
	vec.resize(dim);
	MPI_Bcast(&vec[0], dim, MPI_INT, 0, MPI_COMM_WORLD);

	int chunk = dim/comm_size;

	num_rows += chunk;
	//prepare the output vector
	//send the matrix rows
	
	// if (proc_Id == comm_size - 1) {
	// 	//ready to receive on other procs
	// 	own_mat_ptrs.resize(num_rows + 1);
	// 	cout << proc_Id << " waiting" << endl;
	// 	MPI_IRecv(&own_mat_ptrs[0],own_mat_ptrs.size(),MPI_INT,0,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE,&request11);
	// 	own_mat_data.resize(own_mat_ptrs.back());
	// 	own_mat_indices.resize(own_mat_ptrs.back());
	// 	MPI_IRecv(&own_mat_indices[0],own_mat_indices.size(),MPI_INT,0,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE,&request12);
	// 	MPI_IRecv(&own_mat_data[0],own_mat_data.size(),MPI_INT,0,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE,&request13);
	// 	cout << proc_Id << " received" << endl;
	// }
 	if (proc_Id != 0) {
 		// if (proc_Id == comm_size - 1){
			// own_mat_ptrs.resize(num_rows);
 		// }
 		// else {
 			own_mat_ptrs.resize(num_rows + 1);
 		// }
		// cout << proc_Id << " waiting" << endl;
		MPI_Recv(&own_mat_ptrs[0],own_mat_ptrs.size(),MPI_INT,0,proc_Id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		own_mat_data.resize(own_mat_ptrs.back());
		own_mat_indices.resize(own_mat_ptrs.back());
		// if (proc_Id == comm_size - 1){
			// cout << "PTRS " << own_mat_ptrs.size() << endl;
		// 	// cout << "THIS " << own_mat_ptrs.back() << endl;
			// print(own_mat_ptrs);
		// }
		MPI_Recv(&own_mat_indices[0],own_mat_indices.size(),MPI_INT,0,proc_Id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// if (proc_Id == comm_size - 1){
		// 	cout << proc_Id << " received2" << endl;
		// 	// cout << "THIS " << own_mat_ptrs.back() << endl;
		// }
		MPI_Recv(&own_mat_data[0],own_mat_data.size(),MPI_INT,0,proc_Id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// if (proc_Id == comm_size - 1){
		// 	// cout << proc_Id << " received3" << endl;
		// }
	}
	else {
		int next_index = num_rows;
		int next = ptrs[next_index];
		int temp, old_index;
		vector<int> temp_ptr,temp_indices;
		vector<int> temp_data;
		//own stuff
		own_mat_indices = sub(indices,0,next);
		own_mat_data = sub(data,0,next);
		own_mat_ptrs = sub(ptrs, 0, num_rows);
		own_mat_ptrs.push_back(own_mat_indices.size());		
		//send load to others
		for (int i = 1; i < comm_size - 1; i++) {
			old_index = next_index;
			next_index += chunk;
			temp = next;
			next = ptrs[next_index];
			temp_data = sub(data, temp, next);
			temp_indices = sub(indices, temp, next);
			temp_ptr = sub(ptrs, old_index, next_index);
			mapped_subtract(temp_ptr, temp);
			temp_ptr.push_back(temp_indices.size());
			MPI_Send(&temp_ptr[0],temp_ptr.size(),MPI_INT,i,i, MPI_COMM_WORLD);
			MPI_Send(&temp_indices[0],temp_indices.size(),MPI_INT,i,i, MPI_COMM_WORLD);
			MPI_Send(&temp_data[0],temp_data.size(),MPI_INT,i,i, MPI_COMM_WORLD);
		}
		//final process' load
		temp_data = sub(data,next,data.size());
		temp_indices = sub(indices,next,indices.size());
		temp_ptr = sub(ptrs,next_index,ptrs.size());
		mapped_subtract(temp_ptr, next);
		temp_ptr.push_back(temp_indices.size());
		MPI_Send(&temp_ptr[0],temp_ptr.size(),MPI_INT,comm_size - 1,comm_size - 1, MPI_COMM_WORLD);		
		// cout << "sending to last2 "<< temp_indices.size() << endl;
		MPI_Send(&temp_indices[0],temp_indices.size(),MPI_INT,comm_size - 1,comm_size - 1, MPI_COMM_WORLD);		
		// cout << "sending to last3 "<< temp_data.size() << endl;
		MPI_Send(&temp_data[0],temp_data.size(),MPI_INT,comm_size - 1,comm_size - 1, MPI_COMM_WORLD);
		
	}
}

void getOutput() {
	if (proc_Id == 0) {
		output = new long long int[dim];
		memcpy(output, own_output, num_rows*sizeof(long long int));
		int len, totallen = num_rows;
		for (int i = 1; i < comm_size; i++) {
			MPI_Recv(&len, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(output + totallen, len, MPI_LONG_LONG_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			totallen += len;
		}
	}
	else {
		for (int i = 1; i < comm_size; i++) {
			MPI_Send(&num_rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Send(own_output, num_rows, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);
		}
	}
}

void wrapperForCuda()
{
	int dimension = dim;
	//Get the number of rows being handled by the current process 
	int numRows = num_rows;
		
	//Each partition's rows[i], colsIndices[i] and values[i] 
	int *currPartitionFirstElemsRows;
	currPartitionFirstElemsRows = &own_mat_ptrs[0];

	int *currPartitionColIndices; 
	currPartitionColIndices = &own_mat_indices[0];

	int *currPartitionNonZeroElems;
	currPartitionNonZeroElems = &own_mat_data[0];

	//Commom vector for all processes
	int *vecTOR;
	vecTOR = &vec[0];

	//Device copies for computation
	int *devCurrPartitionFirstElemsRows;
	int *devCurrPartitionColIndices;
	int *devCurrPartitionNonZeroElems;
	int *devVec;
	long long int *devFinalVec;

	int size1 = own_mat_ptrs.size()*sizeof(int);
	int size2 = own_mat_indices.size()*sizeof(int);
	int size3 = own_mat_data.size()*sizeof(int);

	//Current process's computed output
	own_output = (long long int *)malloc(sizeof(long long int)*numRows);
	
	//once
	cudaMalloc((void **)&devFinalVec,numRows*sizeof(long long int));
	cudaMalloc((void **)&devVec,dimension*sizeof(int));
	
	N = numRows;

	cudaMalloc((void **)&devCurrPartitionFirstElemsRows,size1);
	cudaMalloc((void **)&devCurrPartitionColIndices,size2);
	cudaMalloc((void **)&devCurrPartitionNonZeroElems,size3);

	cudaMemcpy(devCurrPartitionFirstElemsRows,currPartitionFirstElemsRows,size1,cudaMemcpyHostToDevice);
	cudaMemcpy(devCurrPartitionColIndices,currPartitionColIndices,size2,cudaMemcpyHostToDevice);
	cudaMemcpy(devCurrPartitionNonZeroElems,currPartitionNonZeroElems,size3,cudaMemcpyHostToDevice);
	cudaMemcpy(devVec,vecTOR,dimension*sizeof(int),cudaMemcpyHostToDevice);
	
	//Tuning for the problem size
	int blocks;
	int thrds;
	if(num_rows<THREADS_PER_BLOCK)
	{
		blocks  = 1;
		thrds = num_rows;
	}
	else
	{
		thrds = THREADS_PER_BLOCK;
		blocks = (num_rows/thrds) + 1;
	}

	multKernel<<<blocks,thrds>>>(devCurrPartitionFirstElemsRows,devCurrPartitionColIndices,devCurrPartitionNonZeroElems,numRows,devVec,devFinalVec);

	cudaMemcpy(own_output,devFinalVec,numRows*sizeof(long long int),cudaMemcpyDeviceToHost);
}

void computeForEachProcess()
{
	wrapperForCuda();
}

void fileWrite(char *name)
{
	if(proc_Id==0)
	{
			ofstream f_out;
			f_out.open(name);
			for (int i = 0; i < dim; i++){
				f_out << output[i] << '\n';
			}
			f_out.close();
	}
}
int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_Id);

	num_rows = 0;
	string in_file = argv[1];	
	getInput(argv[1]);
	
	// cout << "proc" << proc_Id << endl;
	// cout << "indices" << endl;
	// print(own_mat_indices);
	// cout << "data" << endl;
	// print(own_mat_data);
	// cout << "ptrs" << endl;
	// print(own_mat_ptrs);
	// cout << "vec";
	// print(vec);
	// cout << endl;
	// own_output.push_back(proc_Id);
	// getOutput();
	// print(output);
	
		// if (proc_Id == 0) {
		// 	cout << "total data" << endl;
		// 	print(indices);
		// 	print(data);
		// 	print(ptrs);
		// 	print(vec);
		// 	cout << endl;
		// }

	computeForEachProcess();
	getOutput();

	// print(output);

	fileWrite(argv[2]);

	// if(proc_Id==0)
	// 	for(int i=0;i<dim;i++)
	// 		printf("OUTPUT[%d]: %d\n",i,output[i]);

	
	MPI_Finalize();
}