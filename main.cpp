#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb-master/stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb-master/stb_image.h"

//include my classes
#include "classes/scene.cpp"

Vector pixel_to_coordinates(Vector Q, int W, int H,double alpha, int x,int y){
    // point q = Q.give_point();
    double vx = Q[0] + x + 0.5 - W/2;
    double vy = Q[1] + y + 0.5 - H/2;
    double vz = Q[2] - W/(2*tan(alpha/2));
    return Vector(vx,vy,vz);
}

int main() {
    int W = 512;
    int H = 512;
    //parameters
    int max_path_length = 5;
    double alpha = M_PI/3; // pi/4
    double gamma = 2.2;
    Scene scene;
    Vector Q = scene.get_Q();
    Vector S = scene.get_S();
    
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
            Vector color = scene.getColor(r, max_path_length);
            color = pow(color,1./gamma);
            image[i*W*3+j*3 + 0] = std::min(255,int(color[0]));
            image[i*W*3+j*3 + 1] = std::min(255,int(color[1]));
            image[i*W*3+j*3 + 2] = std::min(255,int(color[2]));
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    return 0;
}