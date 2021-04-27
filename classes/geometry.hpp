#pragma once
#include "ray.cpp"
class Geometry {
public:
    int id; Vector C; Vector albedo; 
    double refract_ind = -100.; 
    bool mirror = false;
    bool light = false; 
    Motion motion;
    virtual Intersection intersect(Ray r) = 0;
};
