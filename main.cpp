#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb-master/stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb-master/stb_image.h"

//include my classes
#include "classes/scene.cpp"

Vector pixel_to_coordinates(Vector Q, int W, int H,double alpha, double x,double y){
    point q = Q.give_point();
    double vx = q.x + x + 0.5 - W/2;
    double vy = q.y + y + 0.5 - H/2;
    double vz = q.z - W/(2*tan(alpha/2));
    return Vector(vx,vy,vz);
}

int main() {
    int W = 512;
    int H = 512;
    // double eps = pow(1,-5);
    double alpha = M_PI/3; // pi/4
    Scene scence;
    Vector Q = scence.get_Q();
    Vector S = scence.get_S();
    double gamma = 2.2;
    std::vector<unsigned char> image(W*H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            double x = j;
            double y = H - 1 - i;
            //we compute the coordinates of (x,y)
            Vector coord = pixel_to_coordinates(Q,W,H,alpha,x,y);
            // corresponding ray
            Ray r(Q,coord);
            //compute intersection
            Intersection inter = scence.intersect(r);
            if (inter.intersects){
            inter.L = scence.Lambertian(inter.albedo,inter);
            inter.L = pow(inter.L,1./gamma);
            // inter.L.print_vector();
            image[(i*W + j) * 3 + 0] = int(inter.L[0]);
            image[(i*W + j) * 3 + 1] = int(inter.L[1]);
            image[(i*W + j) * 3 + 2] = int(inter.L[2]);
            // image[(i*W + j) * 3 + 0] = inter.L[0];
            // image[(i*W + j) * 3 + 1] = inter.L[1];
            // image[(i*W + j) * 3 + 2] = inter.L[2];
            }
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    return 0;
}