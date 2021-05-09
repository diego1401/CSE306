#pragma once
#include <stdio.h>      /* printf */
#include <iostream>
#include <math.h>       /* pow */
#include <stdlib.h>     /* srand, rand */
#include <vector>
#include <random>
#include <string>
#include <list>
// #include <iostream>
// #include <stdio.h>
#include <string.h>
#include <algorithm>
// #include <vector>
static std :: default_random_engine engine(10) ; // random seed = 10 
static std::uniform_real_distribution<double> uniform(0, 1);

//functions to simplify computations
double square(double x){return x*x;}

double pow_5(double x){return x*x*x*x*x;}

class Vector { 
public :
    explicit Vector(double x = 0. , double y = 0. , double z = 0.){ 
    coords[0] = x;
    coords[1] = y;
    coords[2] = z;
    };
    //useful functions
    void print_vector(){
        printf("(%f,%f,%f)\n",coords[0],coords[1],coords[2]);
    }
    int argmax(){
        double max = std::max(coords[0],std::max(coords[1],coords[2]));
        for(int i=0;i<3;i++){
            if(coords[i]==max) return i;
        }
        printf("BRUHHHHHH\n");
        return -1;
    }
    // Vector operations
    Vector& operator+=(const Vector& b) {
    coords[0] += b[0]; coords[1] += b[1]; coords[2] += b[2];
    return *this; 
    }

    Vector& operator-=(const Vector& b) {
    coords[0] -= b[0]; coords[1] -= b[1]; coords[2] -= b[2];
    return *this; 
    }

    double norm(){
        return sqrt((square(coords[0]) + square(coords[1]) + square(coords[2])));
    }
    double norm_squared(){
        return (square(coords[0]) + square(coords[1]) + square(coords[2]));
    }

    void normalize(){
        double norm = this->norm();
        if (norm){
        coords[0] /= norm; coords[1] /= norm; coords[2] /= norm;
        }
    }
    Vector& operator+=(const double b) {
    coords[0] += b; coords[1] += b; coords[2] += b;
    return *this; 
    }

    Vector& operator*=(const double b) {
    coords[0] *= b; coords[1] *= b; coords[2] *= b;
    return *this; 
    }

    Vector& operator-=(const double b) {
    coords[0] -= b; coords[1] -= b; coords[2] -= b;
    return *this; 
    }
    

    const double& operator [] (int i) const { return coords [i] ; } 
    double& operator [](int i) { return coords[i]; }


private :
    double coords[3];
};

//CAREFUL: We need to declare the symmetry
Vector operator*(const Vector& a, const double b){
    return Vector(a[0] * b , a[1] * b , a[2] * b); 
    }

Vector operator*(const double b,const Vector& a){
    return Vector(a[0] * b , a[1] * b , a[2] * b); 
    }

Vector pow(const Vector& a,double exp){
    return Vector(pow(a[0],exp) , pow(a[1],exp)  , pow(a[2],exp));
}

// Vector operations
Vector operator+(const Vector& a, const Vector &b){
    return Vector(a[0] + b[0] , a[1] + b[1] , a[2] + b[2]); 
    }

Vector operator-(const Vector& a, const Vector &b){
    return Vector(a[0] - b[0] , a[1] - b[1] , a[2] - b[2]); 
    }

Vector element_wise_product(const Vector& a, const Vector &b){
    return Vector(a[0] * b[0] , a[1] * b[1] , a[2] * b[2]); 
    }
double dot(const Vector& a, const Vector& b) { 
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

Vector cross(const Vector& a, const Vector& b) { 
    return Vector(a[1] * b[2] - a[2] * b[1],
                  a[2] * b[0] - a[0] * b[2],
                  a[0] * b[1] - a[1] * b[0]);
    }

Vector random_cos(const Vector& N){
    double min = N[0];double argmin = 0;
    for (int i=1;i<3;i++){
        if(abs(N[i])<min){
            min = N[i];
            argmin = i;
        }
    }
    Vector T1;
    if(argmin==0) T1 = Vector(0,N[2],-N[1]);
    else if(argmin==1) T1 = Vector(-N[2],0,N[0]);
    else if(argmin==2) T1 = Vector(-N[1],N[0],0);

    //Makes it a bit faster to comment it but colors change a little bit
    T1.normalize();

    //Use Box-Muller formula
    double r1 = uniform (engine) ;
    double r2 = uniform (engine) ;

    double c = sqrt(1-r2);
    double c2 = 2*M_PI*r1;
    double x = cos(c2) * c;
    double y = sin(c2) * c;
    double z = sqrt(r2);

    Vector T2 = cross(T1,N);
    return x * T1 + y * T2 + z *N;
}

Vector random_pow(const Vector&N, const double& alpha){
    double min = N[0];double argmin = 0;
    for (int i=1;i<3;i++){
        if(abs(N[i])<min){
            min = N[i];
            argmin = i;
        }
    }
    Vector T1;
    if(argmin==0) T1 = Vector(0,N[2],-N[1]);
    else if(argmin==1) T1 = Vector(-N[2],0,N[0]);
    else if(argmin==2) T1 = Vector(-N[1],N[0],0);

    //Makes it a bit faster to comment it but colors change a little bit
    T1.normalize();

    //Use Box-Muller formula
    double r1 = uniform (engine) ;
    double r2 = uniform (engine) ;

    double z = pow(r2,1/(alpha+1));
    double c1 = sqrt(1-z*z);
    double x = cos(2*M_PI*r1) * c1;
    double y = sin(2*M_PI*r1) * c1;

    Vector T2 = cross(T1,N);
    return x * T1 + y * T2 + z *N; // return H
}

// Vector Blinn_PhongBRDF(Vector& p_d, Vector& p_s, 
//         const double& alpha, const Vector& omega_o,const Vector& H,const double& expr){
//     // double diffusive_prob = p_d.norm()/(p_d.norm() + p_s[0]);
//     // double diffusive_prob = 1-p_s[0];
//     double diffusive_prob = p_d.norm() /(p_d + p_s).norm();
//     if(uniform(engine)<diffusive_prob){
//         return p_d*(1/(M_PI*diffusive_prob)) * expr;
//     }
//     else{
//         return p_s * ((alpha + 8)/(alpha+1)) *std::max(0., dot(omega_o,H))* (1/(1-diffusive_prob)) ;
//     }
//     return Vector(0.,0.,0.);
// }

// Vector Blinn_PhongBRDF_specular(const Vector& p_s, const double& alpha,
//                                 const Vector& omega_o,const Vector& H){
//     return p_s* ((alpha + 8)/(alpha+1)) * std::max(0., dot(omega_o,H));
// }
struct Motion{
    double t;
    Vector location;
    Vector speed;
};

struct Intersection
{
    bool intersects;
    Vector incoming_direction;
    Vector P;
    Vector N;
    Vector albedo; // color of the sphere it intersects with
    double length;
    int id;
};
