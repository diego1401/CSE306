#include "scene.hpp"
#include <random>


static std :: default_random_engine engine(10) ; // random seed = 10 
static std::uniform_real_distribution<double> uniform(0, 1);

Vector random_cos(const Vector& N){
    
    Vector T1(N[0],N[1],N[2]);
    double min = INFINITY;double argmin = -1;
    for (int i=0;i<3;i++){
        if(abs(T1[i])<min){
            min = T1[i];
            argmin = i;
        }
    }
    if(argmin==0) T1 = Vector(0,N[2],-N[1]);
    else if(argmin==1) T1 = Vector(-N[2],0,N[0]);
    else if(argmin==2) T1 = Vector(-N[1],N[0],0);
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
    camera_distance = 55;
    // aperture = 0.35;
    Vector q(0,0,camera_distance);
    // Vector light_center(-10,25,-20);
    Vector light_center(-10,20,40);

    eps = pow(10,-7);
    Q = q;
    S = light_center;
    I = 2*pow(10,10);
    samples = 2000;
    light_source_radius = 6;
    refractive_index_air = 1.0003;
    refractive_index_ball = 1.5;
    last_refractive_index = refractive_index_air;
    
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
    //motion
    Motion trivial_motion;
    trivial_motion.location = Vector(0,0,0);
    trivial_motion.speed = Vector(0,0,0);
    trivial_motion.t = 0;

    Vector speed(0,-9.8,0);
    Motion ball2;
    ball2.speed = speed;
    //State the center of the walls
    Vector wall_t (0,dist,0); //top wall
    Sphere top_wall(wall_t,R,red,trivial_motion);

    Vector wall_bot (0,-dist,0); //bottom wall
    Sphere bottom_wall(wall_bot,R+50,blue,trivial_motion);

    Vector wall_f (0,0,-dist); //front wall
    Sphere front_wall(wall_f,R,green,trivial_motion);

    Vector wall_b (0,0,dist); //back wall
    Sphere back_wall(wall_b,R,red,trivial_motion);

    Vector wall_l (dist,0,0); //left wall
    Sphere left_wall(wall_l,R,yellow,trivial_motion);

    Vector wall_r (-dist,0,0); //right wall
    Sphere right_wall(wall_r,R,cyan,trivial_motion);
    // Creating and pushing the first sphere of our interest
    Vector center (10,0,-20); //ball
    Sphere center_sphere(center,10,purple,trivial_motion,true);

    //refractive sphere
    Vector center2 (0,0,0); //ball
    Sphere center_sphere2(center2,10,white,trivial_motion);

    // Hollow sphere
    Vector center3 (-10,10,20); //ball
    Sphere center_sphere3(center3,10,red,ball2,false,refractive_index_ball);
    // Sphere center_sphere4(center3,9.5,purple,false,refractive_index_air);

    //ligth sphere
    Sphere light_sphere(light_center,light_source_radius,white,trivial_motion,false,-100,true);


    
    // We push the walls into our list of spheres
    int counter = 0;
    left_wall.sphere_id = counter; counter++;
    s.push_back(left_wall);
    right_wall.sphere_id = counter; counter++;
    s.push_back(right_wall);
    top_wall.sphere_id = counter; counter++;
    s.push_back(top_wall);
    bottom_wall.sphere_id = counter; counter++;
    s.push_back(bottom_wall);
    front_wall.sphere_id = counter; counter++;
    s.push_back(front_wall);
    back_wall.sphere_id = counter; counter++;
    s.push_back(back_wall);
    light_sphere.sphere_id = counter; counter++;
    s.push_back(light_sphere);
    center_sphere.sphere_id = counter; counter++;
    s.push_back(center_sphere);
    center_sphere2.sphere_id = counter; counter++;
    s.push_back(center_sphere2);
    center_sphere3.sphere_id = counter; counter++;
    s.push_back(center_sphere3);
    // center_sphere4.sphere_id = counter; counter++;
    // s.push_back(center_sphere4);
    
    };

double Scene::get_refr_index_air(){return refractive_index_air;}

Vector Scene::get_Q(){return Q;};

Vector Scene::get_S(){return S;};

Vector Scene::direct_light(Vector rho,Intersection& inter,double time){
        // add direct light
        Vector D = inter.P - get_S(); D.normalize();
        Vector xprime = light_source_radius * random_cos(D) + get_S();
        Vector Nprime = (xprime - get_S()); Nprime.normalize();
        Vector omega_i = (xprime - inter.P); double distance_sq = omega_i.norm_squared(); omega_i.normalize();
        Ray visible(inter.P + eps*inter.N,omega_i,time);
        Intersection inter_visible = intersect(visible);
        double visibility = s[inter_visible.sphere_id].is_light() or !inter.intersects; //if we touch the light
        double pdf = dot(Nprime,D)/ (M_PI * square(light_source_radius));
        Vector L;
        L = rho *(1/M_PI) * visibility  * (I / (4* square(M_PI*light_source_radius)))
                                        * (std::max(dot(inter.N,omega_i),0.))
                                        * (std::max(dot(Nprime,-1*omega_i),0.))
                                        *(1/(distance_sq*pdf));
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

Ray Scene::reflect(Vector omega_i,Intersection& inter,double time){
    Vector omega_r = omega_i - 2* dot(omega_i,inter.N)* inter.N;
    Ray reflected_ray(inter.P + eps*inter.N,omega_r,time);
    return reflected_ray;
}

Ray Scene::refract(Vector omega,Intersection& inter,double n1,double& n2,double time){
        double ratio = n1/n2;
        double d = dot(omega,inter.N);
        Vector omega_T = ratio * (omega - d*inter.N);
        //we build omega_N
        double delta = 1 - square(ratio) * (1-square(d)) ;
        if (delta<0) {
            n2 = n1; // we stay in the same environement
            return reflect(omega,inter,time);//total internal reflection
        }
        last_refractive_index = n1;
        //
        //At this point we are sure we are dealing with a refraction      
        Vector omega_N = -1*inter.N * sqrt(delta);
        //We build omega_t
        Vector omega_t = omega_T + omega_N;
        omega_t.normalize();
        Ray refracted_ray(inter.P - eps*inter.N,omega_t,time);
        return refracted_ray;
}

bool Scene::Fresnel(Vector omega, Intersection& inter,double& n1,double&n2){
    //Returns a boolean such that true  -> refractive
    //                            false -> reflective
    //passing inter and n2 as reference because we change them
    if (dot(omega,inter.N)>0){ //going out of a sphere
        //when going out of a sphere n1 becomes the refractive index
        // of the sphere we are intersecting, and n2 is by construction
        // the last n1.
        //we could improve this later on if we wish to have spheres inside spheres inside spheres..
        if (n2!=refractive_index_air) n2 = refractive_index_air;
        else n2 = refractive_index_ball;
        //////////
        inter.N *= -1.;
        going_in = false;
    }
    else going_in = true;
    double k0 = square(n1-n2)/square(n1+n2);
    double R = k0 + (1-k0) * pow_5(1-abs(dot(inter.N,omega)));
    double u = uniform (engine) ;
    return u<R;
}

Vector Scene::getColor(const Ray& ray, int ray_depth,double refr_index,double time,bool last_bounce_diffuse=false){
    if (ray_depth >= 0){
        Intersection inter = intersect(ray);
        if (inter.intersects){ 
            Vector omega = inter.incoming_direction;
            double n2 = s[inter.sphere_id].get_refract();
            Sphere current_sphere = s[inter.sphere_id];

            if(current_sphere.is_light()){ //hit a light source
                if(last_bounce_diffuse){return Vector(0,0,0);}
                return Vector(1.,1.,1.)* (I/ (4.*square(M_PI*light_source_radius)));
            }
            if (current_sphere.is_mirror()) // we reflect
            {
                return getColor(reflect(omega,inter,time),(ray_depth-1),refr_index,time); 
            }
            else if (n2>0.) // sphere is refractive
            {   
                if(Fresnel(omega,inter,refr_index,n2)){
                    return getColor(reflect(omega,inter,time),(ray_depth-1),refr_index,time); //reflection
                }
                return getColor(refract(omega,inter,refr_index,n2,time),(ray_depth-1),n2,time);
            }
            else{ //diffusive surfaces
                //implementation of soft shadows in direct_light
                Vector Lo = direct_light(inter.albedo,inter,time);
                Ray randomRay(inter.P+eps*inter.N,random_cos(inter.N),time);
                Lo += element_wise_product(inter.albedo,getColor(randomRay,ray_depth-1,refr_index,time,true));
                return Lo;
            } 
        }
    }
    return Vector(0.,0.,0.); // we terminate when ray_depth is less than 0
};

Ray Scene::depth_of_field(Vector Q, Vector u,double D,double aperture,double time){
    // The formulas had to be changed a bit
    Vector P = Q + (D/abs(u[2])) * u;
    double r = uniform (engine)  ; r = sqrt(r)* aperture;
    double theta = uniform (engine) * M_PI*2;
    double x = r* cos(theta); double y = r* sin(theta); double z = Q[2];
    Vector Q_prime(x,y,z);
    Vector u_prime = P - Q_prime; u_prime.normalize();
    return Ray(Q_prime,u_prime,time);
}

