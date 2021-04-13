#include "scene.hpp"
double pow_5(double x){return x*x*x*x*x;}

Scene::Scene(){ 
    //camera and light source
    Vector q(0,0,55);
    Vector light(-10,20,40);

    eps = pow(10,-7);
    Q = q;
    S = light;
    I = 2*pow(10,10);
    samples = 1000;
    // samples = 1;
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
    Vector center (-21,0,0); //ball
    Sphere center_sphere(center,10,white,true);
    // Sphere center_sphere(center,10,white);

    //refractive sphere
    Vector center2 (0,0,0); //ball
    Sphere center_sphere2(center2,10,white,false,1.5);
    // Sphere center_sphere2(center2,10,white);

    // Hollow sphere
    Vector center3 (20,0,0); //ball
    Sphere center_sphere3(center3,10,white,false,1.5);
    Sphere center_sphere4(center3,9.5,white,false,refractive_index_air);
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
        L = I / (4* square(M_PI*d)) * rho * (std::max(dot(inter.N,omega),0.0));
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
Ray Scene::reflect(Vector omega_i,Intersection inter){
    Vector omega_r = omega_i - 2* dot(omega_i,inter.N) * inter.N;
    Ray reflected_ray(inter.P + eps*inter.N,omega_r);
    return reflected_ray;
}
Ray Scene::refract(Vector omega,Intersection& inter,double n1,double& n2){
    // We need to pass n2 by reference since it can be subject to change
        //we build omega_T
        // if (inter.sphere_id == 9) inter.N *= -1.;
        if (dot(omega,inter.N)>0){
            n2 = refractive_index_air;
            inter.N *= -1.;
        }
        double ratio = n1/n2;
        double d = dot(omega,inter.N);
        Vector omega_T = ratio * (omega - d*inter.N);
        //we build omega_N
        double delta = 1 - square(ratio) * (1-square(d)) ;
        if (delta <0 ){ // case of total reflection
            return reflect(omega,inter);
        }
        Vector omega_N = -1*inter.N * sqrt(delta);
        //We build omega_t
        Vector omega_t = omega_T + omega_N;
        omega_t.normalize();
        Ray reflected_ray(inter.P - eps*inter.N,omega_t);
        
        return reflected_ray;
}
Vector Scene::getColor(const Ray& ray, int ray_depth,double refr_index){
    if (ray_depth >= 0){
        Intersection inter = intersect(ray);
        if (inter.intersects){ 
            Vector omega = inter.incoming_direction;
            double n2 = s[inter.sphere_id].get_refract();

            if (s[inter.sphere_id].is_mirror()) // we reflect
            {
            return getColor(reflect(omega,inter),(ray_depth-1),refr_index); 
            }
            else if (n2>0.) // sphere is refractive
            {
                double k0 = square(refr_index-n2)/square(refr_index+n2);
                double R = k0 + (1-k0) * pow_5(1-abs(dot(inter.N,omega)));
                double u = ((double) rand() / (RAND_MAX));
                // double u = 1;
                if (u>R) //refraction
                {
                    return getColor(refract(omega,inter,refr_index,n2),(ray_depth-1),n2); 
                }
                else //reflection
                {
                    return getColor(reflect(omega,inter),(ray_depth-1),refr_index); 
                }
            }
            else{
                return Lambertian(inter.albedo,inter);
            } 
        }
    }
    return Vector(0.,0.,0.); // we terminate when ray_depth is less than 0
};