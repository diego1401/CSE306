#pragma once
#include "sphere.cpp"
#include "obj_reader.cpp"

class Scene {
public:
    int samples;  
    double light_source_radius; Vector light;
    double camera_distance; double aperture; Vector Camera;
    double I; double eps;
    double refractive_index_air; double refractive_index_ball; 
    
    std::vector<Geometry*> objects;

    Scene();
    bool Fresnel(Vector omega, Intersection& inter,double& n1,double&n2);
    Intersection intersect(Ray r);
    Ray reflect(Vector omega_i,Intersection& inter,double time);
    Ray refract(Vector omega_i,Intersection& inter,double n1,double& n2,double time);
    Vector direct_light(Vector rho,Intersection& inter,double time); 
    Ray depth_of_field(Vector Q,Vector u,double D,double aperture,double time);
    Vector getColor(const Ray& ray, int ray_depth, double refr_index,double time,bool last_bounce_diffuse);
};