#pragma once
#include "ray.cpp"
class Geometry {
public:
    int id; Vector C; Vector albedo; 
    double refract_ind = -100.; 
    bool mirror = false;
    bool light = false; 
    bool diffusive = false;
    bool Phong = false;
    double alpha; 
    Vector p_s;
    Motion motion;
    virtual Intersection intersect(const Ray& r) = 0;
};
