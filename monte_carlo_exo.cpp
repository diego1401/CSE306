#include <random>
#include <stdio.h>      /* printf */

static std :: default_random_engine engine(10) ; // random seed = 10 
static std::uniform_real_distribution<double> uniform(0, 1);

double cube(double x){return x*x*x;}
double square(double x){return x*x;}

void boxMuller(double stdev , double &x, double &y) { 
    double r1 = uniform ( engine ) ;
    double r2 = uniform ( engine ) ;
    x = sqrt(-2 *log(r1))*cos(2 * M_PI*r2)*stdev; 
    y = sqrt(-2 *log(r1))*sin(2 * M_PI*r2)*stdev;
}

double p(double stdev,double x,double y,double z){
    return cube(1./(stdev*sqrt(2 *M_PI))) * exp((square(x) + square(y) + square(z)) / (-2. * square(stdev)));
}

double f(double x,double y,double z){
    return cos(x*y*z);
}

bool in_domain(double x,double low_b=-M_PI_2,double up_b=M_PI_2){
    return (x>low_b) and (x<up_b);
}

int main(){
    double N = 10000;
    double alpha = 0;
    double x1,y1,z1;
    double x2,y2,z2;
    double stdev = 1.;
    for(int i=0;i<N;i+=2){
        boxMuller(stdev,x1,x2);
        boxMuller(stdev,y1,y2);
        boxMuller(stdev,z1,z2);
        // the integration is 0 if outisde of the domain
        if(in_domain(x1) and in_domain(y1) and in_domain(z1)){
            alpha += f(x1,y1,z1)/p(stdev,x1,y1,z1);
        }
        if(in_domain(x2) and in_domain(y2) and in_domain(z2)){
            alpha += f(x2,y2,z2)/p(stdev,x2,y2,z2);
        }
    }
    alpha *= 1/N;
    printf("%f\n",alpha);
    return alpha;

}