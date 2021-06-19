#include "../Raytracer_Assignement1/classes/obj_reader.cpp"
#include <set>

std::vector<Vector> Id_Boundary(TriangleMesh Mesh){
    //do stuff
}

double sum_boundary(std::vector<Vector> Boundary){
    double res = 0;
    int n = Boundary.size();
    for(int i=0;i<n;i++){
        res += (Boundary[(i+1)%n]-Boundary[i%n]).norm();
    }
    return res;
}

std::vector<Vector> init_disk(int number_vertices,std::vector<Vector> Boundary){
    int n  = Boundary.size(); double cs = 0;
    std::vector<Vector> vertices_disk(number_vertices);
    double s = sum_boundary(Boundary);
    //first n are boundary nodes
    for(int i=0;i<n;i++){
        double theta = 2 * M_PI * (cs/s);
        vertices_disk[i] = Vector(cos(theta),sin(theta),0);
        cs += (Boundary[(i+1)%n]-Boundary[i%n]).norm();
    }
    //rest are interior points
    for(int i=n;i<number_vertices;i++){
        double theta = 2*M_PI* rand()/RAND_MAX;
        double r = 1 * sqrt(rand()/RAND_MAX);
        vertices_disk[i] = Vector(r*cos(theta),r*sin(theta),0);
    }
    return vertices_disk;
}

std::vector<std::set<Vector>> find_adjecent(TriangleMesh Mesh){
    int number_vertices = Mesh.size();
    std::vector<std::set<Vector>> adjecent(number_vertices);
    for(int i=0;i<number_vertices;i++){
        TriangleIndices triangle = Mesh.indices[i];
        adjecent[triangle.vtxi].insert(Mesh.vertices[triangle.vtxj]);
        adjecent[triangle.vtxi].insert(Mesh.vertices[triangle.vtxk]);

        adjecent[triangle.vtxj].insert(Mesh.vertices[triangle.vtxi]);
        adjecent[triangle.vtxj].insert(Mesh.vertices[triangle.vtxk]);

        adjecent[triangle.vtxk].insert(Mesh.vertices[triangle.vtxj]);
        adjecent[triangle.vtxk].insert(Mesh.vertices[triangle.vtxi]);
    }

    return adjecent;
}

Vector sum(std::set<Vector> s){
    Vector res(0,0,0);
    for (std::set<Vector>::iterator it = s.begin(); it != s.end(); ++it) {
    res += *it;
    }
    return res;
}

std::vector<Vector>  Tutte_embedding(TriangleMesh Mesh,int n_iter = 100){
    std::vector<Vector> Bound_M = Id_Boundary(Mesh); int n = Bound_M.size(); int number_vertices = Mesh.size();
    std::vector<Vector> v = init_disk(Mesh.size(),Bound_M);
    std::vector<std::set<Vector>> adjecent = find_adjecent(Mesh);
    for(int i=0;i<n_iter;i++){
        for(int i=n;i<number_vertices;i++){
            double K = adjecent[i].size();
            v[i] = (1./K) * sum(adjecent[i]);
        }
        // for(int i=0;i<n;i++){
        //     v[i] = v[i];
        // }
    }
    return v;
    
}


int main(){
    return 0;
}