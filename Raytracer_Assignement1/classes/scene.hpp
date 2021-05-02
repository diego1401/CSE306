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
    bool Fresnel(const Vector& omega, Intersection& inter,double& n1,double&n2);
    Intersection intersect(const Ray& r);
    Ray reflect(const Vector& omega_i,Intersection& inter,const double& time);
    Ray refract(const Vector& omega_i,Intersection& inter,double n1,double& n2,const double& time);
    Vector direct_light(const Vector& rho,Intersection& inter,const double& time); 
    Ray depth_of_field(const Vector& Q,const Vector& u,const double& D,const double& aperture,const double& time);
    Vector getColor(const Ray& ray, const int& ray_depth, double refr_index,const double& time,bool last_bounce_diffuse);
};