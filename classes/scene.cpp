#include "scene.hpp"

Scene::Scene(){ 
    //camera and light source
    Vector q(0,0,55);
    Vector light(-10,20,40);

    eps = pow(10,-7);
    Q = q;
    S = light;
    I = 2*pow(10,10);
    samples = 100;
    refractive_index_air = 1.0003;
    //All the walls have the same radius
    double R = 940.;
    double dist = 1000;
    //colors
    Vector white(255,255,255);
    Vector blue(0,0,255);
    Vector red(255,0,0);
    Vector purple(255,0,255);
    Vector green(0,255,0);
    Vector cyan(0,255,255);
    Vector yellow(255,255,0);
    //State the center of the walls
    Vector wall_t (0,dist,0); //top wall
    Sphere top_wall(wall_t,R,red);

    Vector wall_bot (0,-dist,0); //bottom wall
    Sphere bottom_wall(wall_bot,R+50,blue);

    Vector wall_f (0,0,-dist); //front wall
    Sphere front_wall(wall_f,R,yellow);

    Vector wall_b (0,0,dist); //back wall
    Sphere back_wall(wall_b,R,cyan);

    Vector wall_l (dist,0,0); //left wall
    Sphere left_wall(wall_l,R,green);

    Vector wall_r (-dist,0,0); //right wall
    Sphere right_wall(wall_r,R,purple);
    // Creating and pushing the first sphere of our interest
    Vector center (-20,0,0); //ball
    Sphere center_sphere(center,10,white,true);
    // Sphere center_sphere(center,10,white);

    //refractive sphere
    Vector center2 (0,0,0); //ball
    Sphere center_sphere2(center2,10,white,false,1.5);
    // Sphere center_sphere2(center2,10,white);

    // Hollow sphere
    Vector center3 (20,0,0); //ball
    Sphere center_sphere3(center3,10,white,false,1.5);
    Sphere center_sphere4(center3,9.8,white,false,refractive_index_air);
    // Sphere center_sphere3(center3,10,white);
    // Sphere center_sphere4(center3,9.8,white);


    
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
    center_sphere2.sphere_id = 7;
    s.push_back(center_sphere2);
    center_sphere3.sphere_id = 8;
    s.push_back(center_sphere3);
    center_sphere4.sphere_id = 9;
    s.push_back(center_sphere4);
    };

Vector Scene::get_Q(){return Q;};
Vector Scene::get_S(){return S;};
double Scene::get_refr_index_air(){return refractive_index_air;}

Vector Scene::Lambertian(Vector rho,Intersection inter){
        //Preparing parameters
        double d = (S-inter.P).norm();
        Vector omega = S-inter.P; omega.normalize();
        //launch ray just above the surface
        Ray visible(inter.P + eps*inter.N,omega);
        Intersection inter_visible = intersect(visible);
        Vector L;
        if (!inter_visible.intersects or (inter_visible.length>d)){
        L = I /(4* pow(M_PI*d,2)) * rho * (std::max(dot(inter.N,omega),0.0));
        }
        return L;
    }

Intersection Scene::intersect(Ray r){
        double tmp = 0;
        Intersection inter;
        inter.intersects = false;
        double d = INFINITY;
        for(std::vector<Sphere>::iterator it = s.begin(); it != s.end(); ++it) { 
            Intersection tmp_inter = (*it).intersect(r);
            if(tmp_inter.intersects){
                tmp = tmp_inter.length;
                if(tmp < d){
                    d = tmp;
                    inter = tmp_inter;
                }
            }
        }
        return inter;
    }

Vector Scene::getColor(const Ray& ray, int ray_depth,double refr_index){
    if (ray_depth >= 0){
    Intersection inter = intersect(ray);
    if (inter.intersects){ 
        Vector omega = inter.incoming_direction;
        double n2 = s[inter.sphere_id].get_refract();
        
        if (s[inter.sphere_id].is_mirror()) 
        {
        Vector omega_r = omega - 2* dot(omega,inter.N) * inter.N;
        Ray reflected_ray(inter.P + eps*inter.N,omega_r);
        return getColor(reflected_ray,(ray_depth-1),refractive_index_air); 
        }
        else if (n2>0.) // sphere is refractive
        {
            double k0 = pow((refr_index-n2),2)/pow((refr_index+n2),2);
            double R = k0 + (1-k0) * pow(1-abs(dot(inter.N,omega)),5);
            double u = ((double) rand() / (RAND_MAX));
            assert(u<1);
            assert(0<u);
            if (u>R) //refraction
            {
                    //we build omega_T
                if (dot(omega,inter.N)>0){
                    n2 = refractive_index_air;
                    inter.N *= -1.;
                }
                double ratio = refr_index/n2;
                Vector omega_T = ratio * (omega - dot(omega,inter.N)*inter.N);
                //we build omega_N
                double delta = std::max(0.,1 - pow(ratio,2) * (1-pow(dot(omega,inter.N),2)));
                Vector omega_N = -1*inter.N * sqrt(delta);
                //We build omega_t
                Vector omega_t = omega_T + omega_N;
                omega_t.normalize();
                Ray reflected_ray(inter.P - eps*inter.N,omega_t);
                return getColor(reflected_ray,(ray_depth-1),n2); 
            }
            else //reflection
            {
                Vector omega_r = omega - 2* dot(omega,inter.N) * inter.N;
                Ray reflected_ray(inter.P + eps*inter.N,omega_r);
                return getColor(reflected_ray,(ray_depth-1),refractive_index_air); 
            }


        }
        /*
        else if (n2 >0.) //so we have a reflective sphere
        {   
            
            //we build omega_T
            if (dot(omega,inter.N)>0){
                n2 = refractive_index_air;
                inter.N *= -1.;
            }
            double ratio = refr_index/n2;
            Vector omega_T = ratio * (omega - dot(omega,inter.N)*inter.N);
            //we build omega_N
            double delta = std::max(0.,1 - pow(ratio,2) * (1-pow(dot(omega,inter.N),2)));
            Vector omega_N = -1*inter.N * sqrt(delta);
            //We build omega_t
            Vector omega_t = omega_T + omega_N;
            omega_t.normalize();
            Ray reflected_ray(inter.P - eps*inter.N,omega_t);
            return getColor(reflected_ray,(ray_depth-1),n2); 
        }
        */
        else{
            return Lambertian(inter.albedo,inter);
        } 
    }
    }
    return Vector(0.,0.,0.); // we terminate when ray_depth is less than 0
};