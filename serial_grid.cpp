#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <chrono>
using namespace std;


double min(double a, double b){
	return(a<b?a:b);
}
double max(double a, double b){
	return(a>b?a:b);
}

struct rectangle{
	double x1;
	double x2;
	double y1;
	double y2;
	double d;
	double mg[4],Mg[4],mm,MM;
	double h,l;
};

double g1(double x1,double x2 ,double l1_max){
return(x1*x1 + x2*x2 - l1_max*l1_max);
}
double dg1x(double x1){
return(2*x1);
}
double g2(double x1,double x2 ,double l1_min){
return(l1_min*l1_min - x1*x1 - x2*x2);
}
double dg2x(double x1){
return(-2*x1);
}
double g3(double x1,double x2 ,double l2_max,double l0){

return((x1-l0)*(x1-l0) + x2*x2 - l2_max*l2_max);
}
double dg3x(double x1,double l0){
return(2*(x1-l0));
}
double g4(double x1,double x2 ,double l2_min,double l0){

return(l2_min*l2_min  - (x1-l0)*(x1-l0) - x2*x2);
}
double dg4x(double x1,double l0){
return(-2*(x1-l0));
}
bool restriction_1(double x1,double l1_max, double l2_max, double l0){
	return((x1 >= -l1_max) and (x1 <= l0+l2_max));
}
bool restriction_2(double x2,double l1_max, double l2_max){
	return((x2>=0) and (x2 <= min(l1_max,l2_max)));
}



vector<rectangle> division(vector<rectangle> &VectorOfRectangles){
rectangle temp = VectorOfRectangles.front();
rectangle r1,r2;
double l = 8;
double lmax = 12;
double l0 = 5;

//GGGGGGGGGGGGGGGGGGGGGG

r1.mg[0] = min(g1(r1.x1,0,lmax),g1(r1.x2,0,lmax)) + min(g1(lmax,r1.y1,lmax),g1(lmax,r1.y2,lmax));
r1.Mg[0] = max(g1(r1.x1,0,lmax),g1(r1.x2,0,lmax)) + max(g1(lmax,r1.y1,lmax),g1(lmax,r1.y2,lmax));
r1.mg[1] = min(g2(r1.x1,0,l),g2(r1.x2,0,l)) + min(g2(l,r1.y1,l),g2(l,r1.y2,l));
r1.Mg[1] = max(g2(r1.x1,0,l),g2(r1.x2,0,l)) + max(g2(l,r1.y1,l),g2(l,r1.y2,l));
r1.mg[2] = min(g3(r1.x1,lmax,lmax,l0),g3(r1.x2,lmax,lmax,l0)) + min(g3(l0,r1.y1,lmax,l0),g3(l0,r1.y2,lmax,l0));
r1.Mg[2] = max(g3(r1.x1,lmax,lmax,l0),g3(r1.x2,lmax,lmax,l0)) + max(g3(l0,r1.y1,lmax,l0),g3(l0,r1.y2,lmax,l0));
r1.mg[3] = min(g4(r1.x1,l,l,l0),g4(r1.x2,l,l,l0)) + min(g4(l0,r1.y1,l,l0),g4(l0,r1.y2,l,l0));
r1.Mg[3] = max(g4(r1.x1,l,l,l0),g4(r1.x2,l,l,l0)) + max(g4(l0,r1.y1,l,l0),g4(l0,r1.y2,l,l0));
r1.mm = max(max(r1.mg[0],r1.mg[1]),max(r1.mg[2],r1.mg[3]));
r1.MM = max(max(r1.Mg[0],r1.Mg[1]),max(r1.Mg[2],r1.Mg[3]));

r2.mg[0] = min(g1(r2.x1,0,lmax),g1(r2.x2,0,lmax)) + min(g1(lmax,r2.y1,lmax),g1(lmax,r2.y2,lmax));
r2.Mg[0] = max(g1(r2.x1,0,lmax),g1(r2.x2,0,lmax)) + max(g1(lmax,r2.y1,lmax),g1(lmax,r2.y2,lmax));
r2.mg[1] = min(g2(r2.x1,0,l),g2(r2.x2,0,l)) + min(g2(l,r2.y1,l),g2(l,r2.y2,l));
r2.Mg[1] = max(g2(r2.x1,0,l),g2(r2.x2,0,l)) + max(g2(l,r2.y1,l),g2(l,r2.y2,l));
r2.mg[2] = min(g3(r2.x1,lmax,lmax,l0),g3(r2.x2,lmax,lmax,l0)) + min(g3(l0,r2.y1,lmax,l0),g3(l0,r2.y2,lmax,l0));
r2.Mg[2] = max(g3(r2.x1,lmax,lmax,l0),g3(r2.x2,lmax,lmax,l0)) + max(g3(l0,r2.y1,lmax,l0),g3(l0,r2.y2,lmax,l0));
r2.mg[3] = min(g4(r2.x1,l,l,l0),g4(r2.x2,l,l,l0)) + min(g4(l0,r2.y1,l,l0),g4(l0,r2.y2,l,l0));
r2.Mg[3] = max(g4(r2.x1,l,l,l0),g4(r2.x2,l,l,l0)) + max(g4(l0,r2.y1,l,l0),g4(l0,r2.y2,l,l0));
r2.mm = max(max(r2.mg[0],r2.mg[1]),max(r2.mg[2],r2.mg[3]));
r2.MM = max(max(r2.Mg[0],r2.Mg[1]),max(r2.Mg[2],r2.Mg[3]));

//VectorOfRectangles.erase(VectorOfRectangles.begin(),VectorOfRectangles.begin()+1);
VectorOfRectangles.push_back(r1);
VectorOfRectangles.push_back(r2);

	return VectorOfRectangles;
}

bool veccomp(rectangle &a, rectangle &b){
	return(a.d > b.d);
}


vector<rectangle> merge_rect(vector<rectangle> &tomerge){

for(int i = 0; i<tomerge.size();i++){
	for(int j=0;j<tomerge.size();j++){
      if(i!=j){
    		if(tomerge[i].x2 == tomerge[j].x1 and tomerge[i].y2 == tomerge[j].y2 and tomerge[i].y1 == tomerge[j].y1){
          rectangle r1;
          r1 = {tomerge[i].x1,tomerge[j].x2,tomerge[i].y1,tomerge[j].y2};
          tomerge.push_back(r1);
          //cout<<"1\n"<<"["<<r1.x1<<":"<<r1.x2<<"]:["<<r1.y1<<":"<<r1.y2<<"]\n";
                if(i>j){
                  tomerge.erase(tomerge.begin()+i,tomerge.begin()+i+1);
                  tomerge.erase(tomerge.begin()+j,tomerge.begin()+j+1);
                }
                else{
                  tomerge.erase(tomerge.begin()+j,tomerge.begin()+j+1);
                  tomerge.erase(tomerge.begin()+i,tomerge.begin()+i+1);
                }i=0;
        }
        if(tomerge[i].y2 == tomerge[j].y1 and tomerge[i].x2 == tomerge[j].x2 and tomerge[i].x1 == tomerge[j].x1){
          rectangle r1;
          r1 = {tomerge[i].x1,tomerge[j].x2,tomerge[i].y1,tomerge[j].y2};
          tomerge.push_back(r1);
          //cout<<"2\n"<<"["<<r1.x1<<":"<<r1.x2<<"]:["<<r1.y1<<":"<<r1.y2<<"]\n";
                if(i>j){
                  tomerge.erase(tomerge.begin()+i,tomerge.begin()+i+1);
                  tomerge.erase(tomerge.begin()+j,tomerge.begin()+j+1);
                }
                else{
                  tomerge.erase(tomerge.begin()+j,tomerge.begin()+j+1);
                  tomerge.erase(tomerge.begin()+i,tomerge.begin()+i+1);
                }i=0;
        }
      if(tomerge[i].y1 == tomerge[j].y2 and tomerge[i].x2 == tomerge[j].x2 and tomerge[i].x1 == tomerge[j].x1){
        rectangle r1;
        r1 = {tomerge[i].x1,tomerge[j].x2,tomerge[j].y1,tomerge[i].y2};
        tomerge.push_back(r1);
        //cout<<"3\n"<<"["<<r1.x1<<":"<<r1.x2<<"]:["<<r1.y1<<":"<<r1.y2<<"]\n";
              if(i>j){
                tomerge.erase(tomerge.begin()+i,tomerge.begin()+i+1);
                tomerge.erase(tomerge.begin()+j,tomerge.begin()+j+1);
              }
              else{
                tomerge.erase(tomerge.begin()+j,tomerge.begin()+j+1);
                tomerge.erase(tomerge.begin()+i,tomerge.begin()+i+1);
              }i=0;
      }
    if(tomerge[i].x1 == tomerge[j].x2 and tomerge[i].y1 == tomerge[j].y1 and tomerge[i].y2 == tomerge[j].y2){
      rectangle r1;
      r1 = {tomerge[j].x1,tomerge[i].x2,tomerge[i].y1,tomerge[j].y2};
      tomerge.push_back(r1);
      //cout<<"4\n"<<"["<<r1.x1<<":"<<r1.x2<<"]:["<<r1.y1<<":"<<r1.y2<<"]\n";
            if(i>j){
              tomerge.erase(tomerge.begin()+i,tomerge.begin()+i+1);
              tomerge.erase(tomerge.begin()+j,tomerge.begin()+j+1);
            }
            else{
              tomerge.erase(tomerge.begin()+j,tomerge.begin()+j+1);
              tomerge.erase(tomerge.begin()+i,tomerge.begin()+i+1);
            }i=0;

    }

}//i=j
sort(tomerge.begin(),tomerge.end(),veccomp);
//cout<<tomerge.size()<<"\n";
//if(tomerge.size() == 0) break;
}
////cout<<"btw i n j "<<i<<" "<<tomerge.size()<<"\n";
}
return(tomerge);
}



int main(){
	int l = 8;
	int l0 = 5;
	double approximation = 0.1;	
	double lmax;
	lmax = l*1.5;
	double ymax,ymin;
	ymin = sqrt(pow(l,2)-pow(l0/2.,2));
	ymax = sqrt(pow(lmax,2)-pow(l0/2.,2));
	double xmin,xmax;
	double x;
	x = (pow(lmax,2)-pow(l,2)-pow(l0,2))/(2*l0);
	xmin = -x;
	xmax = l0+x;

	double h;
	h = sqrt(pow(l,2)-pow(x,2));
	rectangle r1;
	//r1 = {xmin-5, xmax+5, min(ymin,h)-5, ymax+5};
	r1 = {-15,15,0,15};
	r1.d = sqrt(abs(r1.x2-r1.x1)*abs(r1.x2-r1.x1) + abs(r1.y2-r1.y1)*abs(r1.y2-r1.y1));
	vector<rectangle> VectorOfRectangles;
	vector<rectangle> InternalRectangles;
	vector<rectangle> BoundRectangles;
  rectangle r2 = {r1.x1,r1.x1+approximation,r1.y1,r1.y1+approximation};
  //
	std::chrono::time_point<std::chrono:: high_resolution_clock> start, end;
	start = std::chrono::high_resolution_clock::now();

  while(r2.y2 < r1.y2){
    while(r2.x2 < r1.x2){
      r2.x1 += approximation;
      r2.x2 += approximation;
      r2.mg[0] = min(g1(r2.x1,0,lmax),g1(r2.x2,0,lmax)) + min(g1(lmax,r2.y1,lmax),g1(lmax,r2.y2,lmax));
      r2.Mg[0] = max(g1(r2.x1,0,lmax),g1(r2.x2,0,lmax)) + max(g1(lmax,r2.y1,lmax),g1(lmax,r2.y2,lmax));
      r2.mg[1] = min(g2(r2.x1,0,l),g2(r2.x2,0,l)) + min(g2(l,r2.y1,l),g2(l,r2.y2,l));
      r2.Mg[1] = max(g2(r2.x1,0,l),g2(r2.x2,0,l)) + max(g2(l,r2.y1,l),g2(l,r2.y2,l));
      r2.mg[2] = min(g3(r2.x1,lmax,lmax,l0),g3(r2.x2,lmax,lmax,l0)) + min(g3(l0,r2.y1,lmax,l0),g3(l0,r2.y2,lmax,l0));
      r2.Mg[2] = max(g3(r2.x1,lmax,lmax,l0),g3(r2.x2,lmax,lmax,l0)) + max(g3(l0,r2.y1,lmax,l0),g3(l0,r2.y2,lmax,l0));
      r2.mg[3] = min(g4(r2.x1,l,l,l0),g4(r2.x2,l,l,l0)) + min(g4(l0,r2.y1,l,l0),g4(l0,r2.y2,l,l0));
      r2.Mg[3] = max(g4(r2.x1,l,l,l0),g4(r2.x2,l,l,l0)) + max(g4(l0,r2.y1,l,l0),g4(l0,r2.y2,l,l0));
      r2.mm = max(max(r2.mg[0],r2.mg[1]),max(r2.mg[2],r2.mg[3]));
      r2.MM = max(max(r2.Mg[0],r2.Mg[1]),max(r2.Mg[2],r2.Mg[3]));
    VectorOfRectangles.push_back(r2);
    }
    r2.x1 = r1.x1;
    r2.x2 = r2.x1 + approximation;
    r2.y1 += approximation;
    r2.y2 += approximation;
  }

  while(VectorOfRectangles.size() > 0){ //(VectorOfRectangles.size() < 10) and

  	if(VectorOfRectangles.front().MM >= 0){
  		if(VectorOfRectangles.front().mm > 0){
  			if(((VectorOfRectangles.front().x1 <= 0 and VectorOfRectangles.front().x2 >= 0) and (VectorOfRectangles.front().y1 <= 0 and VectorOfRectangles.front().y2 >= lmax))
  			 or	((VectorOfRectangles.front().x1 <= l0 and VectorOfRectangles.front().x2 >= l0) and (VectorOfRectangles.front().y1 <= 0 and VectorOfRectangles.front().y2 >= lmax))){//cout<<"aaaaaaaaaaaaa\n";

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






	end = std::chrono:: high_resolution_clock::now();
	int elapsed_seconds = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	cout<< "#. Время выполнения: " << elapsed_seconds << "  microseconds\n";


  cout<<fixed;
	cout.precision(3);


  for(int i = 0;i<InternalRectangles.size();i++){
  	cout<<"["<<InternalRectangles[i].x1<<":"<<InternalRectangles[i].x2<<"]:["<<InternalRectangles[i].y1<<":"<<InternalRectangles[i].y2<<"]\n";
  }
  cout<<"_+_+_+_\n";
  for(int i = 0;i<BoundRectangles.size();i++){
  	cout<<"["<<BoundRectangles[i].x1<<":"<<BoundRectangles[i].x2<<"]:["<<BoundRectangles[i].y1<<":"<<BoundRectangles[i].y2<<"]\n";
  }

  for(int i = 0;i<VectorOfRectangles.size();i++){
  //cout<<"["<<VectorOfRectangles[i].x1<<":"<<VectorOfRectangles[i].x2<<"]:["<<VectorOfRectangles[i].y1<<":"<<VectorOfRectangles[i].y2<<"]\n";
  }

}
