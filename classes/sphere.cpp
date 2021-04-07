#include "sphere.hpp"


Sphere::Sphere(Vector c, double r,Vector rho,bool flag=false,double n=-1.){ 
        C = c;
        R = r;
        albedo = rho * (1./255.);
        mirror = flag;
        refract_ind = n; // -1 means no refraction
        };

Intersection Sphere::intersect(Ray r){
        Vector u = r.getu();
        u.normalize();
        Vector O = r.getO();
        double t = dot(u,C-O);
        double delta = pow(t,2) - ( pow((O-C).norm(),2) - pow(R,2));

        Intersection inter;
        
        if (delta<0){
            inter.intersects = false;
            return inter;
        }
        inter.P = O;
        if(delta==0){
            inter.length = t;
            inter.P +=t*u;}

        if(delta>0){
            double t1 = t - sqrt(delta);
            double t2 = t + sqrt(delta);
            if(t2<0){
                inter.intersects = false;
                return inter;
                }
            else{
                if (t1>=0){ 
                    inter.length = t1;
                    inter.P  += t1*u;}
                else     { 
                    inter.length = t2;
                    inter.P  += t2*u;};
                }
        }
        inter.intersects = true;
        inter.N = inter.P - C;
        inter.N.normalize();
        inter.albedo = albedo; // initialize albedo
        inter.sphere_id = sphere_id;
        inter.incoming_direction = r.getu();
        return inter;
    };

double Sphere::get_R(){return R;}

Vector Sphere::get_albedo(){return albedo;}

bool Sphere::is_mirror(){return mirror;}

double Sphere::get_refract(){return refract_ind;}

Vector direction(Vector Start, Vector P){
    Vector omega = P - Start;
    omega.normalize();
    return omega;
}
