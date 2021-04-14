#pragma once
#include "sphere.cpp"
#include <vector>


class Scene 
{
public:
    int samples;  
    explicit Scene();
    double get_refr_index_air();
    bool Fresnel(Vector omega, Intersection& inter,double n1,double&n2);
    Intersection intersect(Ray r);
    Ray reflect(Vector omega_i,Intersection inter);
    Ray refract(Vector omega_i,Intersection& inter,double n1,double n2);
    Vector Lambertian(Vector rho,Intersection inter); 
    Vector get_Q();
    Vector get_S();
    Vector getColor(const Ray& ray, int ray_depth, double refr_index);

private :
    double I;
    double eps;
    double refractive_index_air;  
    double last_refr_index;
    
    std::vector<Sphere> s;
    Vector Q;
    Vector S;
    
};