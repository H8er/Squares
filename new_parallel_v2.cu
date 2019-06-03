#include "cuda_runtime.h"
#include <iostream>
#include <chrono>
#include <vector>

using namespace std;

__constant__ double accuracy[1];
__constant__ double bounds[4];
__constant__ int stride[1];
struct rectangle{
	double x1;
	double x2;
	double y1;
	double y2;
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
__forceinline__ __host__ __device__ char interval_eval(double x1,double x2,double y1, double y2, double l, double l0, double lmax){


	double mm,mm1,mm2,mm3,mm4;
	double MM,MM1,MM2,MM3,MM4;


	if(x1<0 and x2>0){
		 x1 = 0;
		 x2 = 0;
	}
	if(x1<l0 and x2>l0){
		 x1 = l0;
		 x2 = l0;
	}


	mm1 = min(g1(x1,0,lmax),g1(x2,0,lmax)) + min(g1(lmax,y1,lmax),g1(lmax,y2,lmax));
	MM1 = max(g1(x1,0,lmax),g1(x2,0,lmax)) + max(g1(lmax,y1,lmax),g1(lmax,y2,lmax));

	mm2 = min(g2(x1,0,l),g2(x2,0,l)) + min(g2(l,y1,l),g2(l,y2,l));
	MM2 = max(g2(x1,0,l),g2(x2,0,l)) + max(g2(l,y1,l),g2(l,y2,l));

	mm3 = min(g3(x1,lmax,lmax,l0),g3(x2,lmax,lmax,l0)) + min(g3(l0,y1,lmax,l0),g3(l0,y2,lmax,l0));
	MM3 = max(g3(x1,lmax,lmax,l0),g3(x2,lmax,lmax,l0)) + max(g3(l0,y1,lmax,l0),g3(l0,y2,lmax,l0));

	mm4 = min(g4(x1,l,l,l0),g4(x2,l,l,l0)) + min(g4(l0,y1,l,l0),g4(l0,y2,l,l0));
	MM4 = max(g4(x1,l,l,l0),g4(x2,l,l,l0)) + max(g4(l0,y1,l,l0),g4(l0,y2,l,l0));

	mm = max(max(mm1,mm2),max(mm3,mm4));
	MM = max(max(MM1,MM2),max(MM3,MM4));
	char a = char(bool(MM <= 0) + bool(mm <= 0));

	return a;
}

__global__ void kernel(char* device_grid) {
double x1[100],x2[100],y1[100],y2[100];
long int i = blockIdx.x * blockDim.x + threadIdx.x;

for(int j = 0;j < stride[0];j++){

	x1[j] = bounds[0] + (threadIdx.x*stride[0]+j)*accuracy[0];
	x2[j] = bounds[0] + (threadIdx.x*stride[0]+j+1)*accuracy[0];
	y1[j] = bounds[2] + (blockIdx.x)*accuracy[0];
	y2[j] = bounds[2] + (blockIdx.x+1)*accuracy[0];

	double l = 8;
	double lmax = 12;
	double l0 = 5;

	device_grid[i*stride[0]+j] = interval_eval(x1[j],x2[j],y1[j],y2[j],l,l0,lmax);
	}

}



int main(){
  int l = 8;
	int l0 = 5;
	double approximation = 0.01;
	double lmax;
	lmax = l*1.5;
  cout<<fixed;
	cout.precision(4);
	rectangle r1;

	long int elapsed_seconds = 0;
	double step = 0.5;
	for(double x_c = -15; x_c < 15;x_c += 2*step){
		int b_s = 0;
		for(double y_c = 0;y_c < 15; y_c += step){
			std::chrono::time_point<std::chrono:: high_resolution_clock> start, end;
			start = std::chrono::high_resolution_clock::now();
			r1.x1 = x_c;
			r1.x2 = x_c + 2*step;
			r1.y1 = y_c;
			r1.y2 = y_c + step;

				if(int(interval_eval(r1.x1,r1.x2,r1.y1,r1.y2,l,l0,lmax))!=0){
					if(int(interval_eval(r1.x1,r1.x2,r1.y1,r1.y2,l,l0,lmax))==2){
						// cout<<"_+_+_+_\n";
						// cout<<"["<<r1.x1<<":"<<r1.x2<<"]:";
						// cout<<"["<<r1.y1<<":"<<r1.y2<<"]\n";
						// cout<<"_+_+_+_\n";
						b_s++;
					}else{
						if(b_s > 0){
							cout<<"_+_+_+_\n";
							cout<<"["<<r1.x1<<":"<<r1.x2<<"]:";
							cout<<"["<<r1.y1-(b_s)*step<<":"<<r1.y2-step<<"]\n";
							cout<<"_+_+_+_\n";
							b_s = 0;
						}
					  double const_bounds[4] = {r1.x1,r1.x2,r1.y1,r1.y2};
						int n_of_blocks = (r1.y2-r1.y1)/approximation;
						int n_of_threads = (r1.x2-r1.x1)/approximation;
					  int offset[1] = {10};
							if(n_of_threads > 1000){
								offset[0] *=10;
							}
					char* device_grid = new char[n_of_blocks * n_of_threads * offset[0]];
					cudaMallocManaged(&device_grid, n_of_blocks * n_of_threads * sizeof (char));
				  cudaMemcpyToSymbol(accuracy, &approximation, sizeof(double));
				  cudaMemcpyToSymbol(stride, &offset, sizeof(int));
				  cudaMemcpyToSymbol(bounds, &const_bounds, 4*sizeof(double));

					kernel<<<n_of_blocks,n_of_threads/offset[0]>>>(device_grid);
					cudaDeviceSynchronize();

					end = std::chrono:: high_resolution_clock::now();
					int t = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
					elapsed_seconds += std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
					std::time_t end_time = std::chrono::system_clock::to_time_t(end);
					//BOUNDARY
					for(int i = 0; i < n_of_blocks*n_of_threads; i++){
						if(int(device_grid[i])==1){
								cout<<"["<<r1.x1+i%(n_of_threads)*approximation<<":"<<r1.x1+(i%(n_of_threads) + 1)*approximation<<"]:";
								cout<<"["<<r1.y1+i/(n_of_threads)*approximation<<":"<<r1.y1+(i/(n_of_threads) + 1)*approximation<<"]\n";
						}
					}
					//internal
					 // for(int i = 0; i < n_of_blocks*n_of_threads; i++){
					 // 	if(int(device_grid[i])==2){
					 // 			cout<<"["<<r1.x1+i%(n_of_threads)*approximation<<":"<<r1.x1+(i%(n_of_threads) + 1)*approximation<<"]:";
					 // 			cout<<"["<<r1.y1+i/(n_of_threads)*approximation<<":"<<r1.y1+(i/(n_of_threads) + 1)*approximation<<"]\n";
					 // 	}
					 // }
					cout<<"_+_+_+_\n";

					int sq = 0;
					int lq = 1;

					for(int j = 0;j < n_of_blocks;j++){
						inner_cycle:
						sq=0;
						for(int i = 0;i < n_of_threads;i++){
							if((int(device_grid[j*n_of_threads+i]) == 2)){
								sq++;
								if(i==n_of_threads-1){
									if(sq!=n_of_threads){
										i++;
										i-=sq;
										cout<<"["<<r1.x1+i*approximation<<":"<<r1.x1+(i + sq)*approximation<<"]:";
										cout<<"["<<r1.y1+j*approximation<<":"<<r1.y1+(j + 1)*approximation<<"]\n";
										i += sq;
										sq = 0;
									}
									else{
										lq++;
										j++;
										goto inner_cycle;
									}
								}
							}
							else{
								if(sq > 0){
										i-=sq;
										cout<<"["<<r1.x1+i*approximation<<":"<<r1.x1+(i + sq)*approximation<<"]:";
										cout<<"["<<r1.y1+j*approximation<<":"<<r1.y1+(j + 1)*approximation<<"]\n";
										i += sq;
										sq = 0;
										if(lq>1){
											cout<<"["<<r1.x1<<":"<<r1.x1+n_of_threads*approximation<<"]:";
											cout<<"["<<r1.y1+(j-lq+1)*approximation<<":"<<r1.y1+(j)*approximation<<"]\n";
											lq=1;
										}
								}
								else{
									if(lq>1){
										cout<<"["<<r1.x1<<":"<<r1.x1+n_of_threads*approximation<<"]:";
										cout<<"["<<r1.y1+(j-lq+1)*approximation<<":"<<r1.y1+(j)*approximation<<"]\n";
										lq=1;
									}
								}
							}
						}
					}
					cout<<"_+_+_+_\n";
					cudaFree(device_grid);
					}
				}
		}
}

cout<< "#. Время выполнения: " << elapsed_seconds << "  microseconds\n";
return 0;
}
