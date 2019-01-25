#include <iostream>
#include <math.h>
#include <vector>

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

double g2(double x1,double x2 ,double l1_min){
return(l1_min*l1_min - x1*x1 - x2*x2);
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

r1.mm = min(g1(r1.x1,0,lmax),g1(r1.x2,0,lmax))+min(g1(lmax,r1.y1,lmax),g1(lmax,r1.y2,lmax));
r1.MM = max(g1(r1.x1,0,lmax),g1(r1.x2,0,lmax))+max(g1(lmax,r1.y1,lmax),g1(lmax,r1.y2,lmax));

r2.mm = min(g1(r2.x1,0,lmax),g1(r2.x2,0,lmax))+min(g1(lmax,r2.y1,lmax),g1(lmax,r2.y2,lmax));
r2.MM = max(g1(r2.x1,0,lmax),g1(r2.x2,0,lmax))+max(g1(lmax,r2.y1,lmax),g1(lmax,r2.y2,lmax));


VectorOfRectangles.push_back(r1);
VectorOfRectangles.push_back(r2);

	return VectorOfRectangles;
}

main(){
	int l = 1;
	int l0 = 5;
	double approximation = 0.1;

	double lmax;
	lmax = 4;

	rectangle r1;
	//r1 = {xmin-5, xmax+5, min(ymin,h)-5, ymax+5};
	r1 = {-5,5,-5,5};
	//r1 = {0,0.03,4.93,5};
	r1.d = sqrt(abs(r1.x2-r1.x1)*abs(r1.x2-r1.x1) + abs(r1.y2-r1.y1)*abs(r1.y2-r1.y1));
	vector<rectangle> VectorOfRectangles;
	vector<rectangle> InternalRectangles;
	vector<rectangle> BoundRectangles;

	r1.mm = min(g1(r1.x1,0,lmax),g1(r1.x2,0,lmax))+min(g1(lmax,r1.y1,lmax),g1(lmax,r1.y2,lmax));
  r1.MM = max(g1(r1.x1,0,lmax),g1(r1.x2,0,lmax))+max(g1(lmax,r1.y1,lmax),g1(lmax,r1.y2,lmax));

	VectorOfRectangles.push_back(r1);


  while((VectorOfRectangles.size() > 0)){ //(VectorOfRectangles.size() < 10) and

    if(VectorOfRectangles.front().MM >= 0){
      if(VectorOfRectangles.front().mm > 0){
        if((VectorOfRectangles.front().x1 <= 0 and VectorOfRectangles.front().x2 >= 0) and (VectorOfRectangles.front().y1 <= 0 and VectorOfRectangles.front().y2 >= 0)){//cout<<"aaaaaaaaaaaaa\n";
				//if(VectorOfRectangles.front().y1 <= 0 and VectorOfRectangles.front().y2 >= 0){//cout<<"aaaaaaaaaaaaa\n";

          if(VectorOfRectangles.front().d <= approximation){
            BoundRectangles.push_back(VectorOfRectangles.front());
          }
          else{
            division(VectorOfRectangles);//cout<<"d1\n";
          }
        }
      }
      else{
        if(VectorOfRectangles.front().d <= approximation){
          BoundRectangles.push_back(VectorOfRectangles.front());
        }
        else{
          division(VectorOfRectangles);//cout<<"d2\n";
        }
      }
      }
      else{
        InternalRectangles.push_back(VectorOfRectangles.front());
      }
			VectorOfRectangles.erase(VectorOfRectangles.begin(),VectorOfRectangles.begin()+1);
}
//cout<<VectorOfRectangles.size()<<" "<<InternalRectangles.size()<<" "<<BoundRectangles.size()<<"\n";

for(int i = 0;i<InternalRectangles.size();i++){
	cout<<"["<<InternalRectangles[i].x1<<":"<<InternalRectangles[i].x2<<"]:["<<InternalRectangles[i].y1<<":"<<InternalRectangles[i].y2<<"]\n";
}
cout<<"_+_+_+_\n";
for(int i = 0;i<BoundRectangles.size();i++){
	cout<<"["<<BoundRectangles[i].x1<<":"<<BoundRectangles[i].x2<<"]:["<<BoundRectangles[i].y1<<":"<<BoundRectangles[i].y2<<"]\n";
}
}
