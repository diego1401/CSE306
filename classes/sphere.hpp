#pragma onces
#include "ray.cpp"
#include <vector>
struct Intersection
{
    bool intersects;
    Vector incoming_direction;
    Vector P;
    Vector N;
    Vector L;
    double length;
    Vector albedo; // from the sphere it intersects
    int sphere_id;
};

class Sphere
{
private:
    Vector C;
    double R;
    Vector albedo;
    bool mirror;
    
public:
    int sphere_id;
    explicit Sphere(Vector c,double r,Vector rho,bool mirror);
    Intersection intersect(Ray r);
    bool is_mirror();
    double get_R();
    Vector get_albedo();
};

