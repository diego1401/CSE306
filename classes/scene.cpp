#pragma once
#include "sphere.cpp"
#include <vector>


class Scene {
public:
    explicit Scene(){ 
    //camera and light source
    Vector q(0,0,55);
    Vector light(-10,20,40);
    eps = pow(10,-5);
    Q = q;
    S = light;
    I = 8*pow(10,9);
    //All the walls have the same radius
    double R = 940.;
    double dist = 1000;
    //colors
    Vector white(255,255,255);
    Vector red(255,0,0);
    Vector cyan(224,255,255);
    Vector green(0,128,0);
    Vector pink(255,182,193);
    Vector color1(255,230,100);
    //State the center of the walls
    Vector wall_t (0,dist,0); //top wall
    Sphere top_wall(wall_t,R,white);

    Vector wall_bot (0,-dist+40,0); //bottom wall
    Sphere bottom_wall(wall_bot,R,white);

    Vector wall_f (0,0,-dist); //front wall
    Sphere front_wall(wall_f,R,green);

    Vector wall_b (0,0,dist); //back wall
    Sphere back_wall(wall_b,R,pink);

    Vector wall_l (dist,0,0); //left wall
    Sphere left_wall(wall_l,R,color1);

    Vector wall_r (-dist,0,0); //right wall
    Sphere right_wall(wall_r,R,color1);
    // Creating and pushing the first sphere of our interest
    Vector center (0,0,0); //ball
    Sphere center_sphere(center,10,red);

    // We push the walls into our list of spheres
    s.push_back(left_wall);
    s.push_back(right_wall);
    s.push_back(center_sphere);
    s.push_back(top_wall);
    s.push_back(bottom_wall);
    s.push_back(front_wall);
    s.push_back(back_wall);
    };

    Sphere get_center_sphere(){
        return s[0];
    }

    Vector Lambertian(Vector rho,Intersection inter){
        //Preparing parameters
        double d = (S-inter.P).norm();
        Vector omega = S-inter.P;
        omega.normalize();
        // launch ray just above the surface
        Ray visible(inter.P + eps*inter.N,omega);
        Intersection inter_visible = intersect(visible);
        Vector L;
        if (!inter_visible.intersects or (inter_visible.length>d)){
        L = I /(4* pow(M_PI*d,2)) * rho * (std::max(dot(inter.N,omega),0.0));
        }
        
        return L;
    }

    Vector get_Q(){return Q;};

    Vector get_S(){return S;};

    Intersection intersect(Ray r){
        double tmp = 0;
        Intersection inter;
        inter.intersects = false;
        double d = INFINITY;
        for(std::vector<Sphere>::iterator it = s.begin(); it != s.end(); ++it) { 
            Intersection tmp_inter = (*it).intersect(r);
            if(tmp_inter.intersects){
                tmp = tmp_inter.P.norm();
                if(tmp < d){
                    d = tmp;
                    inter = tmp_inter;
                }
            }
        }
        return inter;
    }

private :
    std::vector<Sphere> s;
    Vector Q;
    Vector S;
    double I;
    double eps;
};



// int main() {
// }