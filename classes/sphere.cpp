#include "sphere.hpp"


Sphere::Sphere(Vector c, double r,Vector rho,Motion m,bool flag=false,double n=-100.,bool flag2=false){ 
        C = c;
        R = r;
        albedo = rho * (1./255.);
        mirror = flag;
        light = flag2;
        refract_ind = n; // -100 means no refraction
        motion = m;
        };

double Sphere::get_R(){return R;}
bool Sphere::is_mirror(){return mirror;}
bool Sphere::is_light(){return light;}
double Sphere::get_refract(){return refract_ind;}
Vector Sphere::get_albedo(){return albedo;}
Intersection Sphere::intersect(Ray r){
        //treat motion
        Vector moved_C = C+ r.get_time() * motion.speed; //we can do more complex stuff here
        // if(sphere_id==2) {
        //     printf("%f\n",r.get_time());
        // }
        //to put like before just replace by C
        Vector u = r.getu();
        Vector O = r.getO();
        double t = dot(u,moved_C-O);
        double delta = square(t) - ( (O-moved_C).norm_squared() - square(R));

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
        inter.N = inter.P - moved_C;
        inter.N.normalize();
        inter.albedo = albedo; // initialize albedo
        inter.sphere_id = sphere_id;
        inter.incoming_direction = u;
        return inter;
    };
