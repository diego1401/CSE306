#pragma onces
#include "ray.cpp"
#include <vector>
struct Intersection
{
    bool intersects;
    Vector P;
    Vector N;
    Vector L;
    double length;
    Vector albedo; // from the sphere it intersects
};

class Sphere
{
private:
    Vector C;
    double R;
    Vector albedo;
public:
    explicit Sphere(Vector c,double r,Vector rho);
    Intersection intersect(Ray r);
    double get_R();
    Vector get_albedo();
};

