#pragma onces
#include <vector>
#include "vector.cpp"

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

class Geometry {

public:
    virtual Intersection intersect();

private :
    Vector C;
    Vector albedo;
    double refract_ind;
    bool mirror;
    bool light;
};