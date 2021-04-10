#pragma onces
#include "ray.cpp"
#include <vector>
struct Intersection
{
    bool intersects;
    Vector incoming_direction;
    Vector P;
    Vector N;
    Vector albedo; // color of the sphere it intersects with
    double length;
    int sphere_id;
};

class Sphere
{
private:
    Vector C;
    Vector albedo;
    double R;
    double refract_ind;
    bool mirror;
    
public:
    int sphere_id;
    explicit Sphere(Vector c,double r,Vector rho,bool mirror,double n);
    bool is_mirror();
    double get_R();
    double get_refract();
    Intersection intersect(Ray r);
    Vector get_albedo();
};

