#include "svg.cpp"
#include <iostream>   // std::cout
#include <string>     // std::string, std::to_string

Polygon get_bounding_box(Polygon dataset){
    double Bmax_x = 0;
    double Bmax_y = 0;
    // double Bmax_z = INFINITY;
    double Bmin_x = INFINITY;
    double Bmin_y = INFINITY;
    // double Bmin_z = INFINITY;
    for(int i=0;i<dataset.vertices.size();i++){
        Vector curr = dataset.vertices[i];
        Bmax_x = std::max(Bmax_x,curr[0]); Bmax_y = std::max(Bmax_y,curr[1]);
        // Bmax_z = std::max(Bmax_z,curr[2]);
        Bmin_x = std::min(Bmin_x,curr[0]); Bmin_y = std::min(Bmin_y,curr[1]);
        // Bmin_z = std::min(Bmin_z,curr[2]);
    }
    Bmax_x = 1;
    Bmax_y = 1;
    Bmin_x = 0;
    Bmin_y = 0;

    Vector A(Bmin_x,Bmin_y),B(Bmax_x,Bmin_y),C(Bmax_x,Bmax_y),D(Bmin_x,Bmax_y);
    Polygon P({A,B,C,D});
    // A.print_vector();
    // B.print_vector();
    // C.print_vector();
    // D.print_vector();
    return P;
}

Polygon polygon_clipping(Polygon SubjectPolygon, Polygon clipPolygon){
    Polygon* outPolygon = new Polygon();
    for(int i=0;i<clipPolygon.vertices.size();i++){
        outPolygon = new Polygon();
        int curr = i; int prev = (i>0)?(i-1):clipPolygon.vertices.size()-1;
        Vector u = clipPolygon.vertices[prev];
        Vector v = clipPolygon.vertices[curr];

        int n_curr = (clipPolygon.vertices.size()-1>i)? i+1:0 ; int n_prev = i;
        Vector n_edge = clipPolygon.vertices[n_curr] - clipPolygon.vertices[n_prev];
        
        Vector N(v[1]-u[1],u[0]-v[0],0);
        if(dot(n_edge,N)>0){
            N*=-1;
        }

        #define inside1(X) (dot(X-u,N)<=0)
        for(int i=0;i<SubjectPolygon.vertices.size();i++){
            Vector A = SubjectPolygon.vertices[(i>0)?(i-1):SubjectPolygon.vertices.size()-1];
            Vector B = SubjectPolygon.vertices[i]; 
            //Compute the intersection
            
            double t = dot(u-A,N)/dot(B-A,N);
            
            
                Vector P = A + t * (B-A);
                if(inside1(B)){
                    if(!inside1(A)){
                        if(0<=t && t<=1){
                        outPolygon->vertices.push_back(P);
                        }
                    } 
                    outPolygon->vertices.push_back(B);
                }
                else if(inside1(A)) {
                    if(0<=t && t<=1){
                    outPolygon->vertices.push_back(P);
                    }
                }
        }
        SubjectPolygon = *outPolygon;
    }
    return *outPolygon;
}

Polygon clip_bissectors(Polygon SubjectPolygon,int k, Polygon dataset){
    // Polygon SubjectPolygon = get_bounding_box(dataset);
    Vector Pi = dataset.vertices[k];
    Pi.print_vector();
    Polygon* outPolygon = new Polygon();
    for(int j=0;j<dataset.vertices.size();j++){
        if(k==j) continue;
        outPolygon = new Polygon();
        Vector Pj = dataset.vertices[j];
        Vector M = (Pj+Pi)*0.5;
        
        #define inside(X) (dot(X-M,Pj-Pi)<0)
        for(int i=0;i<SubjectPolygon.vertices.size();i++){
            Vector A = SubjectPolygon.vertices[(i>0)?(i-1):SubjectPolygon.vertices.size()-1];
            Vector B = SubjectPolygon.vertices[i]; 
            //Compute the intersection
            
            double t = dot(M-A,Pi-Pj) / dot(B-A,Pi-Pj);
            
            
                Vector P = A + t * (B-A);
                if(inside(B)){
                    if(!inside(A)){
                        if(0<=t && t<=1){
                        outPolygon->vertices.push_back(P);
                        }
                    } 
                    outPolygon->vertices.push_back(B);
                }
                else if(inside(A)) {
                    if(0<=t && t<=1){
                    outPolygon->vertices.push_back(P);
                    }
                }
        }
        SubjectPolygon = *outPolygon;
        save_svg({SubjectPolygon},"tmp" + std::to_string(j)+".svg");
    }
    return *outPolygon;
}
