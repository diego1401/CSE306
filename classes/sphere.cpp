#include "vector.cpp"
#include "ray.cpp"
#pragma once

class Sphere {
public:
    explicit Sphere(Vector C, double R= 0.0){ 
    C = C;
    R = R;
    };
    private :
    Vector C;
    double R;
};

// int main() {
//     Vector nulV;
//     Vector myVector1 (1. , 2. , 3.) ; 
//     Sphere Asphere (myVector1,3.0);
    
// }
