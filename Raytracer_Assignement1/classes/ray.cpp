#pragma once
#include "vector.cpp"
class Ray {
public:
    explicit Ray(Vector v1, Vector v2,double _t){ 
    this->O = v1; this->u = v2; this->u.normalize();this->t = _t;
    };
    Vector O; Vector u; double t; 
};

