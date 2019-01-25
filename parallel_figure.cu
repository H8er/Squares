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

__global__ void kernel(double xmin,double xmax, double ymin, double ymax, double d,double* global_grid) {
__shared__ double x1[1536];
__shared__ double x2[1536];
__shared__ double y1[1536];
__shared__ double y2[1536];
int i = blockIdx.x * blockDim.x + threadIdx.x;
x1[i] = xmin+(threadIdx.x)*(abs(xmax-xmin)/32.);
x2[i] = xmin+(threadIdx.x+1)*(abs(xmax-xmin)/32.);
y1[i] = ymax-(blockIdx.x+1)*(abs(ymax-ymin)/48.);
y2[i] = ymax-(blockIdx.x)*(abs(ymax-ymin)/48.);
__syncwarp();
double mm,mm1,mm2,mm3,mm4;
double MM,MM1,MM2,MM3,MM4;
double l = 8;
double lmax = 12;
double l0 = 5;
double dia;
d = 1;
mm1 = min(g1(x1[i],0,lmax),g1(x2[i],0,lmax)) + min(g1(lmax,y1[i],lmax),g1(lmax,y2[i],lmax));
MM1 = max(g1(x1[i],0,lmax),g1(x2[i],0,lmax)) + max(g1(lmax,y1[i],lmax),g1(lmax,y2[i],lmax));

mm2 = min(g2(x1[i],0,l),g2(x2[i],0,l)) + min(g2(l,y1[i],l),g2(l,y2[i],l));
MM2 = max(g2(x1[i],0,l),g2(x2[i],0,l)) + max(g2(l,y1[i],l),g2(l,y2[i],l));

mm3 = min(g3(x1[i],lmax,lmax,l0),g3(x2[i],lmax,lmax,l0)) + min(g3(l0,y1[i],lmax,l0),g3(l0,y2[i],lmax,l0));
MM3 = max(g3(x1[i],lmax,lmax,l0),g3(x2[i],lmax,lmax,l0)) + max(g3(l0,y1[i],lmax,l0),g3(l0,y2[i],lmax,l0));

mm4 = min(g4(x1[i],l,l,l0),g4(x2[i],l,l,l0)) + min(g4(l0,y1[i],l,l0),g4(l0,y2[i],l,l0));
MM4 = max(g4(x1[i],l,l,l0),g4(x2[i],l,l,l0)) + max(g4(l0,y1[i],l,l0),g4(l0,y2[i],l,l0));

mm = max(max(mm1,mm2),max(mm3,mm4));
MM = max(max(MM1,MM2),max(MM3,MM4));

dia = sqrt(abs(x2[i]-x1[i])*abs(x2[i]-x1[i]) + abs(y2[i]-y1[i])*abs(y2[i]-y1[i]));
if(MM >= 0){
	if(mm > 0){
		if(((x1[i] <= 0 and x2[i] >= 0)and(y1[i]<=0 and y2[i] >= lmax))or
((x1[i] <= l0 and x2[i] >= l0)and(y1[i]<=0 and y2[i] >= lmax)))
		{////cout<<"aaaaaaaaaaaaa\n";
			if(dia <= d){
				//
				global_grid[i] = x1[i];
				global_grid[i+1536] = x2[i];
				global_grid[i+1536*2] = y1[i];
				global_grid[i+1536*3] = y2[i];
			}
			else{
				//division(VectorOfRectangles);
			}
		}
	}
	else{
		if(dia <= d){
			//
			global_grid[i] = x1[i];
			global_grid[i+1536] = x2[i];
			global_grid[i+1536*2] = y1[i];
			global_grid[i+1536*3] = y2[i];
		}
		else{
			//division(VectorOfRectangles);
		}
	}
	}
	else{
		//
		global_grid[i] = x1[i];
		global_grid[i+1536] = x2[i];
		global_grid[i+1536*2] = y1[i];
		global_grid[i+1536*3] = y2[i];
	}



}






bool restriction_1(double x1,double l1_max, double l2_max, double l0){
	return((x1 >= -l1_max) and (x1 <= l0+l2_max));
}
bool restriction_2(double x2,double l1_max, double l2_max){
	return((x2>=0) and (x2 <= min(l1_max,l2_max)));
}

vector<rectangle> division(vector<rectangle> &VectorOfRectangles){
rectangle temp = VectorOfRectangles.front();
rectangle r1;
rectangle r2;
double l = 1;
double lmax = 4;
double l0 = 5;
temp.h = abs(temp.y2 - temp.y1);
temp.l = abs(temp.x2 - temp.x1);
if(temp.h > temp.l){
	 r1 = {temp.x1,temp.x2,temp.y1,temp.y1+temp.h/2.};
	 r2 = {temp.x1,temp.x2,temp.y1+temp.h/2.,temp.y2};

}
else{
	 r1 = {temp.x1,temp.x1+temp.l/2.,temp.y1,temp.y2};
	 r2 = {temp.x1+temp.l/2.,temp.x2,temp.y1,temp.y2};
}
r1.d = sqrt(abs(r1.x2-r1.x1)*abs(r1.x2-r1.x1) + abs(r1.y2-r1.y1)*abs(r1.y2-r1.y1));
r2.d = sqrt(abs(r2.x2-r2.x1)*abs(r2.x2-r2.x1) + abs(r2.y2-r2.y1)*abs(r2.y2-r2.y1));

r1.mg[0] = min(min(g1(r1.x1,r1.y1,lmax),g1(r1.x1,r1.y2,lmax)),min(g1(r1.x2,r1.y2,lmax),g1(r1.x2,r1.y1,lmax)));
r1.Mg[0] = max(max(g1(r1.x1,r1.y1,lmax),g1(r1.x1,r1.y2,lmax)),max(g1(r1.x2,r1.y2,lmax),g1(r1.x2,r1.y1,lmax)));
r1.mg[1] = min(min(g2(r1.x1,r1.y1,l),g2(r1.x1,r1.y2,l)),min(g2(r1.x2,r1.y2,l),g2(r1.x2,r1.y1,l)));
r1.Mg[1] = max(max(g2(r1.x1,r1.y1,l),g2(r1.x1,r1.y2,l)),max(g2(r1.x2,r1.y2,l),g2(r1.x2,r1.y1,l)));
r1.mm = max(r1.mg[0],r1.mg[1]);
r1.MM = max(r1.Mg[0],r1.Mg[1]);

r2.mg[0] = min(min(g1(r2.x1,r2.y1,lmax),g1(r2.x1,r2.y2,lmax)),min(g1(r2.x2,r2.y2,lmax),g1(r2.x2,r2.y1,lmax)));
r2.Mg[0] = max(max(g1(r2.x1,r2.y1,lmax),g1(r2.x1,r2.y2,lmax)),max(g1(r2.x2,r2.y2,lmax),g1(r2.x2,r2.y1,lmax)));
r2.mg[1] = min(min(g2(r2.x1,r2.y1,l),g2(r2.x1,r2.y2,l)),min(g2(r2.x2,r2.y2,l),g2(r2.x2,r2.y1,l)));
r2.Mg[1] = max(max(g2(r2.x1,r2.y1,l),g2(r2.x1,r2.y2,l)),max(g2(r2.x2,r2.y2,l),g2(r2.x2,r2.y1,l)));
r2.mm = max(r2.mg[0],r2.mg[1]);
r2.MM = max(r2.Mg[0],r2.Mg[1]);

VectorOfRectangles.push_back(r1);
VectorOfRectangles.push_back(r2);

	return VectorOfRectangles;
}

int main(){
	int l = 1;
	int l0 = 5;
	double approximation = 0.1;
	cout.precision(5);
	double lmax;
	lmax = 4;
	float t1,t2;
	rectangle r1;
	//r1 = {xmin-5, xmax+5, min(ymin,h)-5, ymax+5};
	r1 = {-15,15,0,15};
	r1.d = sqrt(abs(r1.x2-r1.x1)*abs(r1.x2-r1.x1) + abs(r1.y2-r1.y1)*abs(r1.y2-r1.y1));

	std::chrono::time_point<std::chrono:: high_resolution_clock> start, end;
			start = std::chrono::high_resolution_clock::now();

	double *host_grid = new double[1536*4];
	double* global_grid;
	cudaMalloc ((void **) &global_grid, 1536*4 * sizeof (double));

	kernel<<<48,32>>>(r1.x1,r1.x2,r1.y1,r1.y2,approximation,global_grid);
	cudaMemcpy (host_grid, global_grid, 1536*4 * sizeof (double), cudaMemcpyDeviceToHost);


	end = std::chrono:: high_resolution_clock::now();

			int elapsed_seconds = std::chrono::duration_cast<std::chrono::microseconds>
															 (end-start).count();
			std::time_t end_time = std::chrono::system_clock::to_time_t(end);

			//std:://cout<< "Время выполнения: " << elapsed_seconds << "  microseconds\n";
	t1 = elapsed_seconds;
start = std::chrono::high_resolution_clock::now();

	vector<rectangle> VectorOfRectangles;
	vector<rectangle> InternalRectangles;
	vector<rectangle> BoundRectangles;

	r1.mg[0] = min(min(g1(r1.x1,r1.y1,lmax),g1(r1.x1,r1.y2,lmax)),min(g1(r1.x2,r1.y2,lmax),g1(r1.x2,r1.y1,lmax)));
  r1.Mg[0] = max(max(g1(r1.x1,r1.y1,lmax),g1(r1.x1,r1.y2,lmax)),max(g1(r1.x2,r1.y2,lmax),g1(r1.x2,r1.y1,lmax)));
  r1.mg[1] = min(min(g2(r1.x1,r1.y1,l),g2(r1.x1,r1.y2,l)),min(g2(r1.x2,r1.y2,l),g2(r1.x2,r1.y1,l)));
  r1.Mg[1] = max(max(g2(r1.x1,r1.y1,l),g2(r1.x1,r1.y2,l)),max(g2(r1.x2,r1.y2,l),g2(r1.x2,r1.y1,l)));
  r1.mm = max(r1.mg[0],r1.mg[1]);
  r1.MM = max(r1.Mg[0],r1.Mg[1]);

	cout<<fixed;
	cout.precision(3);
	VectorOfRectangles.push_back(r1);


  while((VectorOfRectangles.size() > 0)){ //(VectorOfRectangles.size() < 10) and
    if(VectorOfRectangles.front().MM >= 0){
      if(VectorOfRectangles.front().mm > 0){
        if(VectorOfRectangles.front().x1 < 0 and VectorOfRectangles.front().x2 > 0){////cout<<"aaaaaaaaaaaaa\n";
          if(VectorOfRectangles.front().d <= approximation){
            BoundRectangles.push_back(VectorOfRectangles.front());
          }
          else{
            division(VectorOfRectangles);
          }
        }
      }
      else{
        if(VectorOfRectangles.front().d <= approximation){
          BoundRectangles.push_back(VectorOfRectangles.front());
        }
        else{
          division(VectorOfRectangles);
        }
      }
      }
      else{
        InternalRectangles.push_back(VectorOfRectangles.front());
      }
VectorOfRectangles.erase(VectorOfRectangles.begin(),VectorOfRectangles.begin()+1);
}
//cout<<VectorOfRectangles.size()<<" "<<InternalRectangles.size()<<" "<<BoundRectangles.size()<<"\n";

for(int i = 0;i<BoundRectangles.size();i++){
	//cout<<"["<<BoundRectangles[i].x1<<":"<<BoundRectangles[i].x2<<"]:[";	//cout<<BoundRectangles[i].y1<<":"<<BoundRectangles[i].y2<<"]\n";
}

for(int i = 0;i<InternalRectangles.size();i++){
	//cout<<"["<<InternalRectangles[i].x1<<":"<<InternalRectangles[i].x2<<"]:["<<InternalRectangles[i].y1<<":"<<InternalRectangles[i].y2<<"]\n";
}
//cout<<48*1024/sizeof(double)/4;

end = std::chrono:: high_resolution_clock::now();

		elapsed_seconds = std::chrono::duration_cast<std::chrono::microseconds>
														 (end-start).count();
		 end_time = std::chrono::system_clock::to_time_t(end);
t2 = elapsed_seconds;
//std:://cout<< "Время выполнения: " << elapsed_seconds << "  microseconds\n";

for(int i = 0; i<1536;i++){
	if(host_grid[i]!=host_grid[i+1536]){
	cout<<"["<<host_grid[i]<<":"<<host_grid[i+1536]<<"]:["<<host_grid[i+1536*2]<<":"<<host_grid[i+1536*3]<<"]\n";
}
}







}
