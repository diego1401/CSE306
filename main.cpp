#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb-master/stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb-master/stb_image.h"

//include my classes
#include "classes/scene.cpp"
#include <chrono>

Vector pixel_to_coordinates(Vector Q, int W, int H,double alpha, int x,int y){
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
    //parameters
    int max_path_length = 7;
    double alpha = M_PI/3; // pi/4
    double gamma = 2.2;
    Scene scene;
    Vector Q = scene.get_Q();
    Vector S = scene.get_S();
    double n_air = scene.get_refr_index_air();
    int samples = scene.samples;
    //to make the averages
    double i1,i2,i3;
    std::vector<unsigned char> image(W*H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            int x = j;
            int y = H - 1 - i;
            //we compute the coordinates of (x,y)
            Vector coord = pixel_to_coordinates(Q,W,H,alpha,x,y);
            // corresponding ray
            Vector q;
            q = (coord-Q);
            q.normalize();
            Ray r(Q,q);
            i1 = 0;i2 = 0;i3 = 0;
            for (int i = 0;i<samples;i++){
                Vector color = scene.getColor(r, max_path_length,n_air);
                // color = pow(color,1./gamma);
                i1+= color[0];i2+= color[1];i3+= color[2];
            }
            i1 /= samples;i2 /= samples;i3 /= samples;

            image[i*W*3+j*3 + 0] = std::min(255,int(pow(i1,1./gamma)));
            image[i*W*3+j*3 + 1] = std::min(255,int(pow(i2,1./gamma)));
            image[i*W*3+j*3 + 2] = std::min(255,int(pow(i3,1./gamma)));
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0); 
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(t2 - t1);
    std::cout << "It took me " << time_span.count() << " seconds.";
    std::cout << std::endl;
    return 0;
}