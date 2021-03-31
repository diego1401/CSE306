#include "vector.cpp"
#include "sphere.cpp"
#pragma once

class Ray {

public:
    explicit Ray(Vector O, Vector u){ 
    O = O;
    u = u;
    };
    
    Vector getO(){ return O;}

    Vector getu(){ return u;}
private :
    Vector O;
    Vector u;
};

struct point
{
    double x;
    double y;
    double z;
};

struct Intersection
{
    bool intersects;
    point p

};

Intersection intersect(Ray r,Sphere s){
    
}