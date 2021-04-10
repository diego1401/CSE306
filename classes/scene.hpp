#pragma once
#include "sphere.cpp"
#include <vector>


class Scene 
{
public:
    int samples;  
    explicit Scene();
    double get_refr_index_air();
    Intersection intersect(Ray r);
    Vector Lambertian(Vector rho,Intersection inter); 
    Vector get_Q();
    Vector get_S();
    Vector getColor(const Ray& ray, int ray_depth, double refr_index);

private :
    double I;
    double eps;
    double refractive_index_air;  
    
    std::vector<Sphere> s;
    Vector Q;
    Vector S;
    
};