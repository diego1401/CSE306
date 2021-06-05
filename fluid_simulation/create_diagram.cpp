#include "svg.cpp"
#include <iostream>   // std::cout
#include <string>     // std::string, std::to_string

Polygon get_bounding_box(){
    double Bmax_x = 1;
    double Bmax_y = 1;
    double Bmin_x = 0;
    double Bmin_y = 0;

    Vector A(Bmin_x,Bmin_y),B(Bmax_x,Bmin_y),C(Bmax_x,Bmax_y),D(Bmin_x,Bmax_y);
    Polygon P({A,B,C,D});
    return P;
}

// Polygon polygon_clipping(Polygon SubjectPolygon, Polygon clipPolygon){
//     Polygon* outPolygon = new Polygon();
//     for(int i=0;i<clipPolygon.vertices.size();i++){
//         outPolygon = new Polygon();
//         int curr = i; int prev = (i>0)?(i-1):clipPolygon.vertices.size()-1;
//         Vector u = clipPolygon.vertices[prev];
//         Vector v = clipPolygon.vertices[curr];
//         int n_curr = (clipPolygon.vertices.size()-1>i)? i+1:0 ; int n_prev = i;
//         Vector n_edge = clipPolygon.vertices[n_curr] - clipPolygon.vertices[n_prev];
//         Vector N(v[1]-u[1],u[0]-v[0],0);
//         if(dot(n_edge,N)>0){
//             N*=-1;
//         }
//         #define inside1(X) (dot(X-u,N)<=0)
//         for(int i=0;i<SubjectPolygon.vertices.size();i++){
//             Vector A = SubjectPolygon.vertices[(i>0)?(i-1):SubjectPolygon.vertices.size()-1];
//             Vector B = SubjectPolygon.vertices[i]; 
//             //Compute the intersection
//             double t = dot(u-A,N)/dot(B-A,N);
//                 Vector P = A + t * (B-A);
//                 if(inside1(B)){
//                     if(!inside1(A)){
//                         if(0<=t && t<=1){
//                         outPolygon->vertices.push_back(P);
//                         }
//                     } 
//                     outPolygon->vertices.push_back(B);
//                 }
//                 else if(inside1(A)) {
//                     if(0<=t && t<=1){
//                     outPolygon->vertices.push_back(P);
//                     }
//                 }
//         }
//         SubjectPolygon = *outPolygon;
//     }
//     return *outPolygon;
// }

Polygon PowerCell_of_i(Polygon SubjectPolygon,int i_index, Polygon dataset,const double* weights){
    Vector Pi = dataset.vertices[i_index];
    // Pi.print_vector();
    // Polygon ret_polygon;
    double wj,wi;
    // wi = (weights[i_index]>0)? weights[i_index]: 0;
    wi = weights[i_index];
    for(int j=0;j<dataset.vertices.size();j++){
        if(i_index==j) continue;
        Polygon* outPolygon = new Polygon();
        Vector Pj = dataset.vertices[j];
        Vector M = (Pj+Pi)*0.5;
        wj = weights[j];
        Vector M_prime = M + (wi-wj)/(2*(Pi-Pj).norm_squared())* (Pj-Pi);
        #define inside(X) (dot(X-M_prime,Pj-Pi)<0)
        for(int i=0;i<SubjectPolygon.vertices.size();i++){
            Vector A = SubjectPolygon.vertices[(i>0)?(i-1):SubjectPolygon.vertices.size()-1];
            Vector B = SubjectPolygon.vertices[i]; 
            //Compute the intersection
            double t = dot(M_prime-A,Pi-Pj) / dot(B-A,Pi-Pj);
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
        // ret_polygon = *outPolygon;
        delete outPolygon;
    }
    return SubjectPolygon;
}

std::vector<Polygon> Create_diagram(Polygon dataset,double* weights){
    Polygon SubjectPolygon = get_bounding_box();int n = dataset.vertices.size();
    std::vector<Polygon> Output; Polygon tmp;
    for(int i=0;i<n;i++){
            Output.push_back(tmp);
        }
    #pragma omp parallel for
    for(int i=0;i<n;i++){
        Output[i] = PowerCell_of_i(SubjectPolygon,i,dataset,weights);
    }
    return Output;
}

std::vector<Polygon> Centroidal_Voronoi_Tesselation(Polygon& dataset){
    Polygon SubjectPolygon = get_bounding_box();
    int N = 100;
    int n = dataset.vertices.size();
    double* weights =(double*) malloc(n*sizeof(double)); // have all weights equal to go back
    // Polygon*P = (Polygon*) malloc(n*sizeof(Polygon));
    std::vector<Polygon> P(n);
    //to the previous case
    for (int i=0;i<n;i++){ weights[i] = 1.;}

    std::vector<Polygon> Output(n);
    
    for(int j=0;j<N;j++){
        std::cout << "LLoyd iteration Nº" << j+1 << std::endl;
        
        #pragma omp parallel for
        for(int i=0;i<n;i++){
            P[i] = PowerCell_of_i(SubjectPolygon,i,dataset,weights);
        }
        // std::cout << "filled P" << std::endl;
        // std::vector<Vector> new_P;
        bool is_liquid;
        for(int i=0;i<n;i++){
            // new_P.push_back(P[i].Centroid2d());
            is_liquid = dataset.vertices[i].ret_is_liquid();
            dataset.vertices[i] = P[i].Centroid2d();
            dataset.vertices[i].set_is_liquid(is_liquid);
        }
        // std::cout << "fixed centroids" << std::endl;
        for(int i=0;i<n;i++){
            Output[i] = P[i];
        }
        // Output = P;
        
        }
    free(weights);
    return Output;

}