#pragma once
#include "vector.cpp"
class Ray {
public:
    explicit Ray(Vector v1, Vector v2,double t){ 
    O = v1; u = v2; u.normalize();time = t;
    };
    Vector getO(){ return O;}
    Vector getu(){ return u;}
    double get_time(){return time;}
private :
    Vector O; Vector u; double time; 
};

