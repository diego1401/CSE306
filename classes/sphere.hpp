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
    Vector C; double R;
    Vector albedo; double refract_ind;  bool mirror;  bool light; Motion motion;
    
public:
    int sphere_id;
    explicit Sphere(Vector c,double r,Vector rho,Motion m,bool mirror,double n,bool light);
    bool is_mirror();
    bool is_light();
    double get_R();
    double get_refract();
    Intersection intersect(Ray r);
    Vector get_albedo();
};

