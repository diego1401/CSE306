#pragma once
#include "sphere.cpp"
#include <vector>


class Scene 
{
public:
    explicit Scene();
    Sphere get_center_sphere();
    Vector Lambertian(Vector rho,Intersection inter);
    Vector get_Q();
    Vector get_S();
    Intersection intersect(Ray r);
    Vector getColor(const Ray& ray, int ray_depth);

private :
    std::vector<Sphere> s;
    Vector Q;
    Vector S;
    double I;
    double eps;
};