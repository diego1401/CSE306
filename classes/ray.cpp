#pragma once
#include "vector.cpp"
//#include "sphere.cpp"


class Ray {

public:
    explicit Ray(Vector v1, Vector v2){ 
    O = v1;
    u = v2;
    };
    
    Vector getO(){ return O;}

    Vector getu(){return u;}
private :
    Vector O;
    Vector u;
};

