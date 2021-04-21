#pragma once
#include "sphere.cpp"
#include <vector>


class Scene 
{
public:
    int samples;  
    bool going_in;
    double light_source_radius;
    double last_refractive_index;
    double camera_distance;
    double aperture;
    explicit Scene();
    double get_refr_index_air();
    bool Fresnel(Vector omega, Intersection& inter,double& n1,double&n2);
    Intersection intersect(Ray r);
    Ray reflect(Vector omega_i,Intersection& inter,double time);
    Ray refract(Vector omega_i,Intersection& inter,double n1,double& n2,double time);
    Vector direct_light(Vector rho,Intersection& inter,double time); 
    Ray depth_of_field(Vector Q,Vector u,double D,double aperture,double time);
    // Vector random_point_on_light_sphere();
    Vector get_Q();
    Vector get_S();
    Vector getColor(const Ray& ray, int ray_depth, double refr_index,double time,bool last_bounce_diffuse);

private :
    double I;
    double eps;
    double refractive_index_air; 
    double refractive_index_ball; 
    
    
    
    
    
    std::vector<Sphere> s;
    Vector Q;
    Vector S;
    
};