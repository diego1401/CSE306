#include <stdio.h>      /* printf */
#include <math.h>       /* pow */
#pragma once


class Vector { 
public :
    explicit Vector(double x = 0. , double y = 0. , double z = 0.,double r=0.,double g=0.,double b=0.){ 
    coords[0] = x;
    coords[1] = y;
    coords[2] = z;
    albedo[0] = r;
    albedo[1] = g;
    albedo[2] = b;
    };
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
        return pow(pow(coords[0],2) + pow(coords[1],2) + pow(coords[2],2),0.5);
    }

    Vector normalize(){
        double norm = this->norm();
        coords[0] /= norm; coords[1] /= norm; coords[2] /= norm;
        return *this;
    }
    //Arithmetic
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
    double albedo[3];
};

//Arithmetic operations
Vector operator+(const Vector& a, const double b){
    return Vector(a[0] + b , a[1] + b , a[2] + b); 
    }

Vector operator*(const Vector& a, const double b){
    return Vector(a[0] * b , a[1] * b , a[2] * b); 
    }

Vector operator-(const Vector& a, const double b){
    return Vector(a[0] - b , a[1] - b , a[2] - b); 
    }
// Vector operations
Vector operator+(const Vector& a, const Vector &b){
    return Vector(a[0] + b[0] , a[1] + b[1] , a[2] + b[2]); 
    }

Vector operator-(const Vector& a, const Vector &b){
    return Vector(a[0] - b[0] , a[1] - b[1] , a[2] - b[2]); 
    }

double dot(const Vector& a, const Vector& b) { 
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

Vector cross(const Vector& a, const Vector& b) { 
    return Vector(a[1] * b[2] - a[2] * b[1],
                  a[2] * b[0] - a[0] * b[2],
                  a[0] * b[1] - a[1] * b[0]);
    }


// int main() {
//     Vector myVector1 (1. , 2. , 3.) ; 
//     printf("hola %f \n", myVector1.norm() );
//     myVector1.normalize();
//     printf ("hola %f \n",myVector1[0]);   
//     printf("hola %f \n", myVector1.norm() );
    
// }
