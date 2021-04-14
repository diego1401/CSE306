#include "scene.hpp"
#include <random>


static std :: default_random_engine engine(10) ; // random seed = 10 
static std::uniform_real_distribution<double> uniform(0, 1);
// std::random_device rd;
// std::mt19937 gen(rd());


Vector random_cos(const Vector& N){
    
    Vector T1(N[0],N[1],N[2]);
    double min = INFINITY;double argmin = -1;
    for (int i=0;i<3;i++){
        if(abs(T1[i])<min){
            min = T1[i];
            argmin = i;
        }
    }
    // assert(argmin != -1);
    //force smallest component of N to be 0
    T1[argmin] = 0;
    //swap the other 2
    double swap1 = -1; double swap2 = -1;
    for(int i=0;i<3;i++){
        if (i != argmin){
            if (swap1 == -1){
                swap1 = i;
            }
            else{
                swap2 = i;
            }
        }
    }
    double tmp = T1[swap1];
    double tmp2 = T1[swap2];
    T1[swap1] = -tmp2;
    T1[swap2] = tmp;
    T1.normalize();

    //Use Box-Muller formula
    double r1 = uniform (engine) ;
    double r2 = uniform (engine) ;

    double x = cos(2*M_PI*r1) * sqrt(1-r2);
    double y = sin(2*M_PI*r1) * sqrt(1-r2);
    double z = sqrt(r2);

    Vector T2 = cross(T1,N);
    return x * T1 + y * T2 + z *N;
}

Scene::Scene(){ 
    //camera and light source
    Vector q(0,0,55);
    Vector light(-10,20,40);

    eps = pow(10,-7);
    Q = q;
    S = light;
    I = 2*pow(10,10);
    samples = 1000;
    refractive_index_air = 1.0003;
    last_refr_index = -1; // we have to keep track of the last sphere we visited 
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
    Sphere top_wall(wall_t,R,cyan);

    Vector wall_bot (0,-dist,0); //bottom wall
    Sphere bottom_wall(wall_bot,R+50,blue);

    Vector wall_f (0,0,-dist); //front wall
    Sphere front_wall(wall_f,R,yellow);

    Vector wall_b (0,0,dist); //back wall
    Sphere back_wall(wall_b,R,red);

    Vector wall_l (dist,0,0); //left wall
    Sphere left_wall(wall_l,R,green);

    Vector wall_r (-dist,0,0); //right wall
    Sphere right_wall(wall_r,R,purple);
    // Creating and pushing the first sphere of our interest
    Vector center (-20,0,0); //ball
    Sphere center_sphere(center,10,white,true);

    //refractive sphere
    Vector center2 (0,0,0); //ball
    Sphere center_sphere2(center2,10,white,false,1.5);

    // Hollow sphere
    Vector center3 (20,0,0); //ball
    Sphere center_sphere3(center3,10,white,false,1.5);
    Sphere center_sphere4(center3,9.5,white,false,refractive_index_air);


    
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

double Scene::get_refr_index_air(){return refractive_index_air;}

Vector Scene::get_Q(){return Q;};

Vector Scene::get_S(){return S;};

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

Ray Scene::refract(Vector omega,Intersection& inter,double n1,double n2){
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

bool Scene::Fresnel(Vector omega, Intersection& inter,double n1,double&n2){
    //Returns a boolean that true  -> refractive
    //                       false -> reflective
    //passing inter and n2 as reference because we change them
    if (dot(omega,inter.N)>0){
        n2 = last_refr_index;
        inter.N *= -1.;
    }
    double k0 = square(n1-n2)/square(n1+n2);
    double R = k0 + (1-k0) * pow_5(1-abs(dot(inter.N,omega)));
    double u = uniform (engine) ;
    return u>R;
}
/*
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
                if(Fresnel(omega,inter,refr_index,n2)) //refractive
                {
                    last_refr_index = refr_index;
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
*/
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
                if(Fresnel(omega,inter,refr_index,n2)) //refractive
                {
                    last_refr_index = refr_index;
                    return getColor(refract(omega,inter,refr_index,n2),(ray_depth-1),n2); 
                }
                else //reflection
                {
                    return getColor(reflect(omega,inter),(ray_depth-1),refr_index); 
                }
            }
            else{
                Vector Lo = Lambertian(inter.albedo,inter);
                Ray randomRay(inter.P+eps*inter.N,random_cos(inter.N));
                Lo += element_wise_product(inter.albedo,getColor(randomRay,ray_depth-1,refr_index));
                return Lo;
            } 
        }
    }
    return Vector(0.,0.,0.); // we terminate when ray_depth is less than 0
};



