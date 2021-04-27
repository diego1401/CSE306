#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb-master/stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb-master/stb_image.h"

//include my classes
#include "classes/scene.cpp"
#include <chrono>

double spread = 0.5; // or 0.25?
void boxMuller(double stdev , double &x, double &y) { 
    double r1 = uniform ( engine ) ;
    double r2 = uniform ( engine ) ;
    x= sqrt(-2 *log(r1))*cos(2 * M_PI*r2)*stdev *spread; 
    y= sqrt(-2 *log(r1))*sin(2 * M_PI*r2)*stdev *spread;
}

Vector pixel_to_coordinates(Vector Q, int W, int H,double alpha, double x,double y){
    double vx = Q[0] + x + 0.5 - W/2;
    double vy = Q[1] + y + 0.5 - H/2;
    double vz = Q[2] - W/(2*tan(alpha/2));
    return Vector(vx,vy,vz);
}

int main() {
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    int W = 512;
    int H = 512;
    
    // int W = 2560;
    // int H = 1600;
    int max_path_length = 5;
    double alpha = M_PI/3;
    double gamma = 2.2;
    Scene scene;
    Vector Q = scene.Camera;
    double D = scene.camera_distance; //distance of the middle white ball
    scene.aperture = 2.8;
    double aperture = scene.aperture;
    double n_air = scene.refractive_index_air;
    int samples = scene.samples;
    std::vector<unsigned char> image(W*H * 3, 0);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color;
            double x = j; double y = H - 1 - i;
            for (int i = 0;i<samples;i++){
                double x1,y1; double t = uniform(engine);
                boxMuller(1,x1,y1);
                Vector rand_dir = pixel_to_coordinates(Q,W,H,alpha,x+x1,y+y1);
                Vector q = rand_dir - Q;
                q.normalize();
                //Without Dof
                // Ray r(Q,q); color += scene.getColor(r, max_path_length,n_air);
                // Implementation of DoF
                color += scene.getColor(scene.depth_of_field(Q,q,D,aperture,t),max_path_length,n_air,t);
            }

            image[i*W*3+j*3 + 0] = std::min(255,int(pow(color[0]/samples,1./gamma)));
            image[i*W*3+j*3 + 1] = std::min(255,int(pow(color[1]/samples,1./gamma)));
            image[i*W*3+j*3 + 2] = std::min(255,int(pow(color[2]/samples,1./gamma)));
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0); 
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(t2 - t1);
    double time = time_span.count();
    std::cout << "It took me " << int(time/60) << " minutes and " << (int(time)%60) << " seconds"; 
    std::cout << std::endl;
    return 0;
}