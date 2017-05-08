#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>
#include <cuda.h>

using namespace std;

template <typename T>
void print(vector<T> &vec) {
 for (int e=0;e< vec.size();e++) {
     cout << vec[e] << ' ';
 }
 cout << endl;
}

template <typename T>
vector<T> sub(vector<T> &vec, int s, int e) {
	vector<T> sub;
	sub.resize(e - s);
	for (int i = 0; i < e - s; i++) {
		sub[i] = vec[s + i];
	}
	return sub;
}

template <typename T>
vector<T> sub(vector<T> &vec, int s) {
	vector<T> sub;
	int e = vec.size();
	sub.resize(e - s);
	for (int i = 0; i < e - s; i++) {
		sub[i] = vec[s + i];
	}
	return sub;
}

void mapped_subtract(vector<int> &vec, int offset){
	for (int e =0;e<vec.size();e++) {
		vec[e] -= offset;
	}
}