#include "cuda_runtime.h"
#include <iostream>
#include <chrono>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

/*
double min(double a, double b){
	return(a<b?a:b);
}
double max(double a, double b){
	return(a>b?a:b);
}
*/
__constant__ double accuracy[1];
__constant__ double bounds[4];
__constant__ int stride[1];
struct rectangle{
	double x1;
	double x2;
	double y1;
	double y2;
	double d;
	double mg[4],Mg[4],mm,MM;
	double h,l;
};

__host__ __device__ double g1(double x1,double x2 ,double l1_max){
return(x1*x1 + x2*x2 - l1_max*l1_max);
}
__host__ __device__ double g2(double x1,double x2 ,double l1_min){
return(l1_min*l1_min - x1*x1 - x2*x2);
}

__host__ __device__ double g3(double x1,double x2 ,double l2_max,double l0){
return((x1-l0)*(x1-l0) + x2*x2 - l2_max*l2_max);
}
__host__ __device__ double g4(double x1,double x2 ,double l2_min,double l0){
return(l2_min*l2_min  - (x1-l0)*(x1-l0) - x2*x2);
}

__global__ void kernel(double* device_grid) {
double x1[10],x2[10],y1[10],y2[10];
int i = blockIdx.x * blockDim.x + threadIdx.x;

for(int j = 0;j < stride[0];j++){

x1[j] = bounds[0] + (threadIdx.x*stride[0]+j)*accuracy[0];
x2[j] = bounds[0] + (threadIdx.x*stride[0]+j+1)*accuracy[0];
y1[j] = bounds[2] + (blockIdx.x)*accuracy[0];
y2[j] = bounds[2] + (blockIdx.x+1)*accuracy[0];

double mm,mm1,mm2,mm3,mm4;
double MM,MM1,MM2,MM3,MM4;
double l = 8;
double lmax = 12;
double l0 = 5;
double dia;

mm1 = min(g1(x1[j],0,lmax),g1(x2[j],0,lmax)) + min(g1(lmax,y1[j],lmax),g1(lmax,y2[j],lmax));
MM1 = max(g1(x1[j],0,lmax),g1(x2[j],0,lmax)) + max(g1(lmax,y1[j],lmax),g1(lmax,y2[j],lmax));

mm2 = min(g2(x1[j],0,l),g2(x2[j],0,l)) + min(g2(l,y1[j],l),g2(l,y2[j],l));
MM2 = max(g2(x1[j],0,l),g2(x2[j],0,l)) + max(g2(l,y1[j],l),g2(l,y2[j],l));

mm3 = min(g3(x1[j],lmax,lmax,l0),g3(x2[j],lmax,lmax,l0)) + min(g3(l0,y1[j],lmax,l0),g3(l0,y2[j],lmax,l0));
MM3 = max(g3(x1[j],lmax,lmax,l0),g3(x2[j],lmax,lmax,l0)) + max(g3(l0,y1[j],lmax,l0),g3(l0,y2[j],lmax,l0));

mm4 = min(g4(x1[j],l,l,l0),g4(x2[j],l,l,l0)) + min(g4(l0,y1[j],l,l0),g4(l0,y2[j],l,l0));
MM4 = max(g4(x1[j],l,l,l0),g4(x2[j],l,l,l0)) + max(g4(l0,y1[j],l,l0),g4(l0,y2[j],l,l0));

mm = max(max(mm1,mm2),max(mm3,mm4));
MM = max(max(MM1,MM2),max(MM3,MM4));

dia = sqrt(abs(x2[j]-x1[j])*abs(x2[j]-x1[j]) + abs(y2[j]-y1[j])*abs(y2[j]-y1[j]));
if(MM >= 0){
	if(mm > 0){
		if(((x1[j] <= 0 and x2[j] >= 0)and(y1[j]<=0 and y2[j] >= lmax))or
((x1[j] <= l0 and x2[j] >= l0)and(y1[j]<=0 and y2[j] >= lmax)))
			{
				device_grid[i*stride[0]+j] = x1[j];
			  device_grid[blockDim.x*gridDim.x*stride[0] + i*stride[0]+j] = x2[j];
			  device_grid[blockDim.x*gridDim.x*stride[0]*2 + i*stride[0]+j] = y1[j];
			  device_grid[blockDim.x*gridDim.x*stride[0]*3 + i*stride[0]+j] = y2[j];
			}
		}
	else{
		device_grid[i*stride[0]+j] = x1[j];
	  device_grid[blockDim.x*gridDim.x*stride[0] + i*stride[0]+j] = x2[j];
	  device_grid[blockDim.x*gridDim.x*stride[0]*2 + i*stride[0]+j] = y1[j];
	  device_grid[blockDim.x*gridDim.x*stride[0]*3 + i*stride[0]+j] = y2[j];
	}
	}
	else{
		//
		device_grid[i*stride[0]+j] = x1[j];
		device_grid[blockDim.x*gridDim.x*stride[0] + i*stride[0]+j] = x2[j];
		device_grid[blockDim.x*gridDim.x*stride[0]*2 + i*stride[0]+j] = y1[j];
		device_grid[blockDim.x*gridDim.x*stride[0]*3 + i*stride[0]+j] = y2[j];
	}

}


__syncwarp();

}



int main(){
  int l = 8;
	int l0 = 5;
	double approximation = 0.01;

	double lmax;
	lmax = l*1.5;
  cout<<fixed;
	cout.precision(3);
	//float t1,t2;
	rectangle r1;
	//r1 = {xmin-5, xmax+5, min(ymin,h)-5, ymax+5};
	r1 = {-15,15,0,15};
	r1.d = sqrt(abs(r1.x2-r1.x1)*abs(r1.x2-r1.x1) + abs(r1.y2-r1.y1)*abs(r1.y2-r1.y1));
  double const_bounds[4] = {r1.x1,r1.x2,r1.y1,r1.y2};
	//std::chrono::time_point<std::chrono:: high_resolution_clock> start, end; start = std::chrono::high_resolution_clock::now();
	int n_of_blocks = (r1.y2-r1.y1)/approximation;
	int n_of_threads = (r1.x2-r1.x1)/approximation;
  int offset[1] = {10};

  cout<<"# <<<"<<n_of_blocks<<" , "<<n_of_threads/offset[0]<<">>>\n";
	//double* host_grid = new double[n_of_blocks*n_of_threads*4];
	double* device_grid;

  std::chrono::time_point<std::chrono:: high_resolution_clock> start, end;
	start = std::chrono::high_resolution_clock::now();

	cudaMallocManaged(&device_grid, n_of_blocks * n_of_threads * 4 * sizeof (double));


  cudaMemcpyToSymbol(accuracy, &approximation, sizeof(double));
  cudaMemcpyToSymbol(stride, &offset, sizeof(int));
  cudaMemcpyToSymbol(bounds, &const_bounds, 4*sizeof(double));

	//double N_of_squares = n_of_blocks*n_of_threads;
	//printf("%d\n",blockDim.x*gridDim.x*stride[0]);
  //
	kernel<<<n_of_blocks,n_of_threads/offset[0]>>>(device_grid);
	//cudaMemcpy (host_grid, device_grid, n_of_blocks * n_of_threads * 4 * sizeof (double), cudaMemcpyDeviceToHost);
cudaDeviceSynchronize();

end = std::chrono:: high_resolution_clock::now();
int elapsed_seconds = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
std::time_t end_time = std::chrono::system_clock::to_time_t(end);
cout<< "#. Время выполнения: " << elapsed_seconds << "  microseconds\n";
cout<<"#"<<n_of_blocks*n_of_threads<<"\n";
for(int i = 0; i < n_of_blocks*n_of_threads; i++){
	if((device_grid[i]!=device_grid[i+n_of_blocks*n_of_threads]) and (device_grid[i+n_of_blocks*n_of_threads*2]!=device_grid[i+n_of_blocks*n_of_threads*3])){
	cout<<"["<<device_grid[i]<<":"<<device_grid[i+n_of_blocks*n_of_threads]<<"]:["<<device_grid[i+n_of_blocks*n_of_threads*2]<<":"<<device_grid[i+n_of_blocks*n_of_threads*3]<<"]\n";
}
}

}
