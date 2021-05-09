#include "scene.hpp"
Scene::Scene(){ 
    //camera and light source
    camera_distance = 55;
    Vector q(0,0,camera_distance);
    Vector light_center(-10,20,40);

    eps = pow(10,-7);
    Camera = q; light = light_center; I = 2*pow(10,10);light_source_radius = 5;
    refractive_index_air = 1.0003; refractive_index_ball = 1.5;
    // p_s = 0.1 * Vector(1.,1.,1.);  alpha = 10; 
    samples = 32;
    
    //All the walls have the same radius
    double R = 940.; double dist = 1000;
    //colors
    Vector white(255,255,255);Vector blue(0,0,255);Vector red(255,0,0);
    Vector purple(255,0,255); Vector green(0,255,0);Vector cyan(0,255,255);
    Vector yellow(255,255,0);
    //motions
    Motion trivial_motion,ball2;
    trivial_motion.speed = Vector(0,0,0);
    ball2.speed = Vector(0,-9.8,0);

    //Building the scene
    Vector wall_t (0,dist,0); //top wall
    Geometry* top_wall;
    top_wall = new Sphere(wall_t,R,red,trivial_motion,false,-100.,false);
    
    Vector wall_bot (0,-dist,0); //bottom wall
    Geometry* bottom_wall;
    bottom_wall = new Sphere(wall_bot,R+50,blue,trivial_motion,false,-100.,false);

    Vector wall_f (0,0,-dist); //front wall
    Geometry* front_wall;
    front_wall = new Sphere(wall_f,R,green,trivial_motion,false,-100.,false);

    Vector wall_b (0,0,dist); //back wall
    Geometry* back_wall;
    back_wall= new Sphere(wall_b,R,red,trivial_motion,false,-100.,false);

    Vector wall_l (dist,0,0); //left wall
    Geometry* left_wall; left_wall= new Sphere(wall_l,R,yellow,trivial_motion,false,-100.,false);

    Vector wall_r (-dist,0,0); //right wall
    Geometry* right_wall; right_wall= new Sphere(wall_r,R,cyan,trivial_motion,false,-100.,false);

    Geometry* light_sphere; //ligth sphere
    light_sphere = new Sphere(light_center,light_source_radius,white,trivial_motion,false,-100,true);

    //interesting spheres
    Vector center (10,0,0); Vector center2 (-10,0,20); Vector center3 (20,0,0); Vector center4 (10,0,-20);
    Geometry* center_sphere; 
    center_sphere= new Sphere(center,10,white,trivial_motion,true,-100.,false);
    Geometry* center_sphere2; //refractive sphere
    center_sphere2= new Sphere(center2,10,white,trivial_motion,false,refractive_index_ball,false);
    Geometry*center_sphere3; // Hollow sphere
    center_sphere3 = new Sphere(center3,10,white,trivial_motion,false,refractive_index_ball);
    Geometry*center_sphere4; 
    center_sphere4 = new Sphere(center3,9.5,white,trivial_motion,false,refractive_index_air);
    Geometry* mirror_sphere; //mirror sphere
    mirror_sphere= new Sphere(center4,10,white,trivial_motion,true,-100.,false);

    int counter = 0; //counter in order to id the objects
    //scence

    left_wall->id = counter; counter++;
    objects.push_back(left_wall);
    right_wall->id = counter; counter++;
    objects.push_back(right_wall);
    top_wall->id = counter; counter++;
    objects.push_back(top_wall);
    bottom_wall->id = counter; counter++;
    objects.push_back(bottom_wall);
    front_wall->id = counter; counter++;
    objects.push_back(front_wall);
    back_wall->id = counter; counter++;
    objects.push_back(back_wall);
    light_sphere->id = counter; counter++;
    objects.push_back(light_sphere);

    //objects 

    //Spheres

    // center_sphere->id = counter; counter++;
    // objects.push_back(center_sphere);
    // center_sphere2->id = counter; counter++;
    // objects.push_back(center_sphere2);
    // center_sphere3->id = counter; counter++;
    // objects.push_back(center_sphere3);
    // center_sphere4->id = counter; counter++;
    // objects.push_back(center_sphere4);
    // mirror_sphere->id = counter; counter++;
    // objects.push_back(mirror_sphere);

    //Mr. Cat
    TriangleMesh* mesh;
    mesh = new TriangleMesh;
    //THESE NAMES ARE TO BE CHANGED IN ORDERED TO SWITCH MESHES
    // AS WELL AS UNCOMMENTING THE CORRECT TRANSFORMATION
    const char* file_name = "classes/cat/cat.obj"; 
    const char* texture = "classes/cat/cat_diff.png";
    // const char* file_name = "classes/fox/fox.obj"; 
    // const char* texture = "classes/fox/fox_diff.png";

    mesh->readOBJ(file_name); 
    mesh->Phong = true;
    mesh->p_s = 0.1 * Vector(1.,1.,1.);  mesh->alpha = 100; 
    mesh->id = counter; counter++;
    mesh->motion = trivial_motion;
    printf("number of triangles %d\n",int(mesh->vertices.size()));
    //cat transformation
    for (int i = 0; i < int(mesh->vertices.size()); i++) {
        mesh->vertices[i] = 0.6 * mesh->vertices[i] + Vector(0, -10, 0);
    }

    //fox transformation
    // for (int i = 0; i < int(mesh->vertices.size()); i++) {
    //     mesh->vertices[i] = 0.2 * mesh->vertices[i] + Vector(-10, -10, 20);
    // }

    //push just the mesh
    // cat->id = counter; counter++;
    // objects.push_back(cat); 

    //just one bounding box
    // BoundingBox* box;
    // box = new BoundingBox(cat,0,cat->indices.size());
    // box->id = counter; counter++;
    // objects.push_back(box);
    
    // his cage
    BVH* tree;
    int H = 1024; int W = 512; int c = 1;int channels = 3;
    tree = new BVH(mesh,0,mesh->indices.size(),stbi_load(texture,&H,&W,&c,channels));
    objects.push_back(tree);  
    };

Vector Scene::direct_light(Vector& omega_o, Vector& H,Intersection& inter,const double& time){
        // Spherical lights
        Vector D = inter.P - light; D.normalize();
        Vector xprime = light_source_radius * random_cos(D) + light;
        Vector Nprime = (xprime - light); Nprime.normalize();
        Vector omega_i = (xprime - inter.P); double distance_sq = omega_i.norm_squared(); omega_i.normalize();
        Ray visible(inter.P + eps*inter.N,omega_i,time);
        Intersection inter_visible = intersect(visible);
        double visibility = !inter_visible.intersects or objects[inter_visible.id]->light; //if we touch the light
        double pdf = dot(Nprime,D)/ (M_PI * square(light_source_radius));

        double expr = 1/pdf * (std::max(dot(inter.N,omega_i),0.))
                            * (std::max(dot(Nprime,-1*omega_i),0.));
        Vector L,BDRF;
        if(objects[inter.id]->diffusive){
            BDRF = inter.albedo*(1/M_PI);
            BDRF *= expr;
        } 
        else if(objects[inter.id]->Phong){
            H = random_pow(inter.N,objects[inter.id]->alpha); H.normalize();
            omega_o =  omega_i - 2* dot(omega_i,H)* H; omega_o.normalize();
            BDRF = inter.albedo*(1/(M_PI)) *expr
                 + objects[inter.id]->p_s * ((objects[inter.id]->alpha + 8)/(objects[inter.id]->alpha+1)) *std::max(0., dot(omega_o,H)); 
                 //we already divided by the pdf
            
        }
        L = BDRF * visibility  * (I/ (4*distance_sq* square(M_PI*light_source_radius)));
        return L;
        /*
        //point light source
        Vector omega_i = light - inter.P;
        double d = omega_i.norm();
        omega_i.normalize();
        Vector L;
        double visibility;
        //compute visibility
        Ray visible(inter.P + eps*inter.N,omega_i,time);
        Intersection inter_visible = this->intersect(visible);
        if(!inter_visible.intersects){
         visibility = 1;
        }   
        else if(inter_visible.length > d){
            visibility = 1;
        }
        else visibility = 0;

        L = rho *(1/M_PI) * visibility  * (I / (4* M_PI*square(d)))
                                        * (std::max(dot(inter.N,omega_i),0.));
        return L;
        */
        
    }

Intersection Scene::intersect(const Ray& r){
        Intersection inter;
        inter.intersects = false;
        double d = INFINITY;
        Intersection tmp_inter;
        Geometry* object;
        for(std::vector<Geometry*>::iterator it = objects.begin(); it != objects.end(); ++it) { 
            object = *it;
            tmp_inter = object->intersect(r);
            if(tmp_inter.intersects){
                if(tmp_inter.length < d){
                    inter = tmp_inter; d = inter.length;
                }
            }
        }
        return inter;
    }

Ray Scene::reflect(const Vector& omega_i,Intersection& inter,const double& time){
    Vector omega_r = omega_i - 2* dot(omega_i,inter.N)* inter.N;
    Ray reflected_ray(inter.P + eps*inter.N,omega_r,time);
    return reflected_ray;
}

Ray Scene::refract(const Vector& omega,Intersection& inter,double n1,double& n2,const double& time){
        double ratio = n1/n2;
        double d = dot(omega,inter.N);
        Vector omega_T = ratio * (omega - d*inter.N);
        //we build omega_N
        double delta = 1 - square(ratio) * (1-square(d)) ;
        if (delta<0) {
            n2 = n1; // we stay in the same environement
            return reflect(omega,inter,time);//total internal reflection
        }
        //At this point we are sure we are dealing with a refraction      
        Vector omega_N = -1*inter.N * sqrt(delta);
        //We build omega_t
        Vector omega_t = omega_T + omega_N;
        omega_t.normalize();
        Ray refracted_ray(inter.P - eps*inter.N,omega_t,time);
        return refracted_ray;
}

bool Scene::Fresnel(const Vector& omega, Intersection& inter,double& n1,double&n2){
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
    }
    double k0 = square(n1-n2)/square(n1+n2);
    double R = k0 + (1-k0) * pow_5(1-abs(dot(inter.N,omega)));
    double u = uniform (engine) ;
    return u<R;
    // return false; //for no fresnel
}

Vector Scene::getColor(const Ray& ray, const int& ray_depth,double refr_index,const double& time,bool last_bounce_diffuse=false){
    if (ray_depth >= 0){
        Intersection inter = intersect(ray);
        if (inter.intersects){ 
            Vector omega = inter.incoming_direction;
            double n2 = objects[inter.id]->refract_ind;

            if(objects[inter.id]->light){ //hit a light source
                if(last_bounce_diffuse){return Vector(0,0,0);}
                return Vector(1.,1.,1.)* (I/ (4.*square(M_PI*light_source_radius)));
            }
            if (objects[inter.id]->mirror) //Mirror surfaces
            {
                return getColor(reflect(omega,inter,time),(ray_depth-1),refr_index,time); 
            }
            else if (n2>0.) //Refractive surfaces
            {   
                if(Fresnel(omega,inter,refr_index,n2)){
                    return getColor(reflect(omega,inter,time),(ray_depth-1),refr_index,time); //reflection
                }
                return getColor(refract(omega,inter,refr_index,n2,time),(ray_depth-1),n2,time);
            }
            else if(objects[inter.id]->diffusive){ //diffusive surfaces
                //implementation of soft shadows in direct_light
                Vector omega_o,H;
                Vector Lo = direct_light(omega_o,H,inter,time);
                Ray randomRay(inter.P+eps*inter.N,random_cos(inter.N),time);
                Lo += element_wise_product(inter.albedo,getColor(randomRay,ray_depth-1,refr_index,time,true));
                return Lo;
            }
            else if(objects[inter.id]->Phong){ //Phong surfaces
                Vector omega_o,sampled_H;
                Vector Lo = direct_light(omega_o,sampled_H,inter,time);
                Vector p_d = inter.albedo;
                double diffusive_probability = p_d.norm()/(p_d + objects[inter.id]->p_s).norm();
                if(uniform(engine)<diffusive_probability){ //diffuse
                    Ray randomRay(inter.P+eps*inter.N,random_cos(inter.N),time);
                    Lo += element_wise_product(inter.albedo*(1/diffusive_probability),
                                getColor(randomRay,ray_depth-1,refr_index,time,true));
                }
                else{ //specular
                    omega_o = omega - 2* dot(omega,sampled_H)* sampled_H; omega_o.normalize();
                    Ray randomRay(inter.P+eps*inter.N,omega_o,time);
                    if(dot(randomRay.u,inter.N)<0) return Vector(0.,0.,0.);
                    // Vector BDRF_indirect = Blinn_PhongBRDF_specular(this->p_s,this->alpha,randomRay.u,sampled_H);
                    Vector BDRF_indirect = objects[inter.id]->p_s* ((objects[inter.id]->alpha + 8)/(objects[inter.id]->alpha+1)) 
                                        * std::max(0., dot(omega_o,sampled_H));
                    Lo += element_wise_product(BDRF_indirect* (1/(1-diffusive_probability)), 
                            std::max(0.,dot(inter.N,randomRay.u)) 
                            *getColor(randomRay,ray_depth-1,refr_index,time)); 
                }
                return Lo;
            } 
        }
    }
    return Vector(0.,0.,0.); // we terminate when ray_depth is less than 0
};

Ray Scene::depth_of_field(const Vector& Q, const Vector& u,const double &D,const double& aperture,const double& time){
    // The formulas had to be changed a bit
    Vector P = Q + (D/abs(u[2])) * u;
    double r = uniform (engine)  ; r = sqrt(r)* aperture;
    double theta = uniform (engine) * M_PI*2;
    double x = r* cos(theta); double y = r* sin(theta); double z = Q[2];
    Vector Q_prime(x,y,z);
    Vector u_prime = P - Q_prime; u_prime.normalize();
    return Ray(Q_prime,u_prime,time);
}

