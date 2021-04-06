#include "scene.hpp"


Scene::Scene(){ 
    //camera and light source
    Vector q(0,0,55);
    Vector light(-10,20,40);
    eps = pow(10,-6);
    Q = q;
    S = light;
    I = 2*pow(10,10);
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
    
    Sphere top_wall(wall_t,R,white,false);

    Vector wall_bot (0,-dist+40,0); //bottom wall
    Sphere bottom_wall(wall_bot,R,red,false);

    Vector wall_f (0,0,-dist); //front wall
    Sphere front_wall(wall_f,R,green,true);

    Vector wall_b (0,0,dist); //back wall
    Sphere back_wall(wall_b,R,pink,false);

    Vector wall_l (dist,0,0); //left wall
    Sphere left_wall(wall_l,R,color1,false);

    Vector wall_r (-dist,0,0); //right wall
    Sphere right_wall(wall_r,R,pink,false);
    // Creating and pushing the first sphere of our interest
    Vector center (0,0,0); //ball
    Sphere center_sphere(center,10,white,false);

    // We push the walls into our list of spheres
    left_wall.sphere_id = 0;
    s.push_back(left_wall);
    right_wall.sphere_id = 1;
    s.push_back(right_wall);
    center_sphere.sphere_id = 2;
    s.push_back(center_sphere);
    top_wall.sphere_id = 3;
    s.push_back(top_wall);
    bottom_wall.sphere_id = 4;
    s.push_back(bottom_wall);
    front_wall.sphere_id = 5;
    s.push_back(front_wall);
    back_wall.sphere_id = 6;
    s.push_back(back_wall);
    };

Sphere Scene::get_center_sphere(){return s[0];}

Vector Scene::Lambertian(Vector rho,Intersection inter){
        //Preparing parameters
        double d = (S-inter.P).norm();
        Vector omega = S-inter.P;
        omega.normalize();
        //launch ray just above the surface
        Ray visible(inter.P + eps*inter.N,omega);
        Intersection inter_visible = intersect(visible);
        Vector L;
        if (!inter_visible.intersects or (inter_visible.length>d)){
        L = I /(4* pow(M_PI*d,2)) * rho * (std::max(dot(inter.N,omega),0.0));
        }
        return L;
    }

Vector Scene::get_Q(){return Q;};

Vector Scene::get_S(){return S;};

Intersection Scene::intersect(Ray r){
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

Vector Scene::getColor(const Ray& ray, int ray_depth){
    if (ray_depth >= 0){

    Intersection inter = intersect(ray);
    if (inter.intersects){ 
        if (s[inter.sphere_id].is_mirror() ) 
        {
        Vector omega = direction(inter.incoming_direction,inter.P);
        omega.normalize();
        if(inter.N[0] * inter.incoming_direction[0] >0){ // if they are the same sign we change
            //since N should point towards the incoming ray
            omega *= (-1);
        }
        Vector omega_r = omega - 2* dot(omega,inter.N) * inter.N;
        Ray reflected_ray(inter.P + eps*inter.N,omega_r);
        return getColor(reflected_ray,(ray_depth-1)); 
        }
        else{
            return Lambertian(inter.albedo,inter);
        } 
    }
    }
    return Vector(0.,0.,0.); // we terminate when ray_depth is less than 0
};