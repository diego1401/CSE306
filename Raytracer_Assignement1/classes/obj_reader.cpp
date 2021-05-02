#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "/Users/Diego/Documents/SEM6/CSE306/Raytracer_Assignement1/stb-master/stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "/Users/Diego/Documents/SEM6/CSE306/Raytracer_Assignement1/stb-master/stb_image.h"
#include "geometry.hpp"
class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class TriangleMesh: public Geometry {
public:
  ~TriangleMesh() {}
	TriangleMesh() {};
	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	
	Vector compute_barycenter(int index){
		TriangleIndices triangle = indices[index];
		Vector A = vertices[triangle.vtxi]; 
		Vector B = vertices[triangle.vtxj]; 
		Vector C = vertices[triangle.vtxk]; 
		return (A+B+C)*(0.3333333);
	}

	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];
		FILE* f;
		f = fopen(obj, "r");
		if(f==NULL){
			printf("error opening the file.\n");
			return;
		}
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;
				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);

	}
	
	Intersection intersect(const Ray& r){
		Intersection inter;
        inter.intersects = false;
        double d = INFINITY;
		Vector O = r.O; Vector u = r.u;
		bool in_triangle;
		Vector A,B,C,e1,e2,N,cros;
		double denominator,denominator_inv,Beta,gamma,alpha,t;

		for(int i = 0;i<indices.size();i++){
			TriangleIndices triangle = this->indices[i];
			A = this->vertices[triangle.vtxi]; 
			B = this->vertices[triangle.vtxj]; 
			C = this->vertices[triangle.vtxk];
			e1 = B-A; e2 = C-A; N = cross(e1,e2);

			denominator = dot(N,u);
			if(denominator){
				cros = cross(A-O,u);
                denominator_inv = 1/denominator;
                Beta = dot(e2,cros)*denominator_inv;
                gamma = -dot(e1,cros)*denominator_inv;
                alpha = 1 - Beta - gamma;
                t = dot(A-O,N)*denominator_inv;
                if(t<d && t>=0){
					in_triangle = ((0<=alpha) && (alpha<=1)) &&
									   ((0<=Beta)  && (Beta<=1))  &&
									   ((0<=gamma) && (gamma<=1));
					if(in_triangle){
						d = t;
						inter.length = t;
						inter.intersects = true;
						inter.P = O + t*u;
						inter.albedo = Vector(1.,1.,1.);
						inter.id = this->id;
						inter.N = N ;
						inter.incoming_direction = u;
					}
                }
			}

        }
        if(inter.intersects) inter.N.normalize();
        return inter;
	}
};

class BoundingBox: public Geometry{
    public:
    Vector Bmin,Bmax,C;
	TriangleMesh* Mesh;
	BoundingBox(){};
	
	bool in_box(const Vector& O){
		//returns wether do vector is in the box or not
		bool condition1 = (O[0]<=this->Bmax[0]) && (O[1]<=this->Bmax[1]) && (O[0]<=this->Bmax[2]);
		bool condition2 = (O[0]>=this->Bmin[0]) && (O[1]>=this->Bmin[1]) && (O[0]>=this->Bmin[2]);
		return condition1&&condition2;
	}

	BoundingBox(TriangleMesh* _Mesh,int start,int end){
		this->Mesh = _Mesh;
		compute_Bounding_Box(_Mesh,start,end);
	}

	void compute_Bounding_Box(TriangleMesh* _Mesh,int start,int end){
		//First we get the biggest and smallest values of x,y,z
		double max_x = 0,max_y = 0,max_z = 0;
		double min_x = INFINITY,min_y = INFINITY,min_z = INFINITY;
		
		for(int i = start;i<end;i++){
			TriangleIndices triangle = _Mesh->indices[i];
			Vector A = _Mesh->vertices[triangle.vtxi]; 
			Vector B = _Mesh->vertices[triangle.vtxj]; 
			Vector C = _Mesh->vertices[triangle.vtxk];
			max_x = std::max({max_x,A[0],B[0],C[0]});
			max_y = std::max({max_y,A[1],B[1],C[1]});
			max_z = std::max({max_z,A[2],B[2],C[2]});

			min_x = std::min({min_x,A[0],B[0],C[0]});
			min_y = std::min({min_y,A[1],B[1],C[1]});
			min_z = std::min({min_z,A[2],B[2],C[2]});
			
		}

		this->Bmax = Vector(max_x,max_y,max_z);
		this->Bmin = Vector(min_x,min_y,min_z);
	}	

	Vector compute_diag(){
		return this->Bmax - this->Bmin;
	}

	Intersection intersect(const Ray& r){
		Vector O = r.O; Vector u = r.u;
		Intersection inter;
		inter.intersects = false;
		//compute tmax
		double tmax_x = (Bmax[0]-O[0])/u[0];
		double tmax_y = (Bmax[1]-O[1])/u[1];
		double tmax_z = (Bmax[2]-O[2])/u[2];
		//compute tmin
		double tmin_x = (Bmin[0]-O[0])/u[0];
		double tmin_y = (Bmin[1]-O[1])/u[1];
		double tmin_z = (Bmin[2]-O[2])/u[2];

		double t0_x,t0_y,t0_z,t1_x,t1_y,t1_z;
		if(tmin_x<tmax_x) {t0_x = tmin_x;t1_x=tmax_x;} else{t0_x = tmax_x;t1_x=tmin_x;}
		if(tmin_y<tmax_y) {t0_y = tmin_y;t1_y=tmax_y;} else{t0_y = tmax_y;t1_y=tmin_y;}
		if(tmin_z<tmax_z) {t0_z = tmin_z;t1_z=tmax_z;} else{t0_z= tmax_z;t1_z=tmin_z;}

		
		//if it intersects with the bounding box we check
		double l = std::max({t0_x,t0_y,t0_z});
		if (std::min({t1_x,t1_y,t1_z})>l){
			//BVH
			inter.intersects = true;
			inter.length = l;
            inter.id = this->id;
			//Just one box
			// inter = this->Mesh->intersect(r);
		}
		return inter;
	}
};

class BVH: public Geometry{
	public:
	TriangleMesh* Mesh;
	BoundingBox* value;
	BVH* lchild; BVH* rchild;
	int start,end,id;
	unsigned char* texture_map; 

	BVH(){};

	BVH(TriangleMesh* _Mesh,int _start,int _end,int _id,unsigned char* _texture_map){
		this->lchild = NULL; this->rchild = NULL;
        this->id = _id; 
		this->Mesh = _Mesh; this->start = _start; this->end = _end;
		this->texture_map = _texture_map;
		this->value = new BoundingBox(this->Mesh,this->start,this->end);
		compute_BVH();
	}

	Intersection intersect_with_mesh(const Ray& r,int start,int end){
		Intersection inter;
        inter.intersects = false;
        double d = INFINITY;
		Vector O = r.O; Vector u = r.u;
		bool in_triangle;
		Vector A,B,C,e1,e2,N,cros,UV;
		double denominator,denominator_inv,Beta,gamma,alpha,t;
		int x,y;
		TriangleIndices triangle;
		for(int i = start;i<end;i++){
			triangle = this->Mesh->indices[i];
			A = this->Mesh->vertices[triangle.vtxi]; 
			B = this->Mesh->vertices[triangle.vtxj]; 
			C = this->Mesh->vertices[triangle.vtxk];
			e1 = B-A; e2 = C-A; N = cross(e1,e2);

			denominator = dot(N,u);
			if(denominator){
				cros = cross(A-O,u);
                denominator_inv = 1/denominator;
                Beta = dot(e2,cros)*denominator_inv;
                gamma = -dot(e1,cros)*denominator_inv;
                alpha = 1 - Beta - gamma;
                t = dot(A-O,N)*denominator_inv;
                if(t<d && t>=0){
					in_triangle = ((0<=alpha) && (alpha<=1)) &&
								  ((0<=Beta)  && (Beta<=1))  &&
								  ((0<=gamma) && (gamma<=1));
					if(in_triangle){
						d = t;
						inter.length = t;
						inter.intersects = true;
						inter.P = O + t*u;
						inter.incoming_direction = u;
						inter.N = alpha * this->Mesh->normals[triangle.ni]
								 + Beta * this->Mesh->normals[triangle.nj]
								 + gamma* this->Mesh->normals[triangle.nk];

						UV= alpha * this->Mesh->uvs[triangle.uvi]
						   + Beta * this->Mesh->uvs[triangle.uvj]
						   + gamma* this->Mesh->uvs[triangle.uvk];
						// printf("%f\n",pow(this->texture_map[int(UV[1]*512*3) + int(UV[0]*3) + 0],2.2));
						x = UV[0]*512; y =  UV[1]*1024;
						x = std::max(std::min(x,512),0);y = 1024 - 1 - std::max(std::min(y,1024),0);
						
						// printf("x=%d,y=%d\n",x,y);
						// UV.print_vector();
						inter.albedo[0] = pow(this->texture_map[int(y*512*3 + x*3) + 0]/255.,2.2);
						inter.albedo[1] = pow(this->texture_map[int(y*512*3 + x*3 )+ 1]/255.,2.2);
						inter.albedo[2] = pow(this->texture_map[int(y*512*3 + x*3) + 2]/255.,2.2);
					}
                }
			}

        }
        if(inter.intersects){
			inter.N.normalize();
			inter.id = this->id;
			// inter.albedo = Vector(1.,1.,1.);
		} 
        return inter;
	}

	Intersection intersect(const Ray& r){
		Intersection inter;
		inter.intersects = false;

		if(!this->value->intersect(r).intersects) return inter;

		std::list<BVH*> nodes_to_visit; nodes_to_visit.push_front(this);
		double best_distance = std :: numeric_limits<double >::max() ;
		Intersection inter_lchild,inter_rchild,inter_with_mesh;
		BVH* curr_node;
		while(!nodes_to_visit.empty()){
			curr_node = nodes_to_visit.back();
			nodes_to_visit.pop_back();
			if(curr_node->lchild){ //if there is a left child we are not in a leaf
			// not this is a full binary tree
				inter_lchild = curr_node->lchild->value->intersect(r);
				inter_rchild = curr_node->rchild->value->intersect(r);
				if(inter_lchild.intersects){
					if(inter_lchild.length<best_distance){
						nodes_to_visit.push_back(curr_node->lchild);
					}
				}
				if(inter_rchild.intersects){
					if(inter_rchild.length<best_distance){
						nodes_to_visit.push_back(curr_node->rchild);
					} 
				}
			}
			else{//we are in a leaf
				inter_with_mesh;
				inter_with_mesh = curr_node->intersect_with_mesh(r,curr_node->start,curr_node->end);
				if(inter_with_mesh.intersects && inter_with_mesh.length < best_distance){
					best_distance = inter_with_mesh.length;
					inter = inter_with_mesh;
					}
			}
		}
		return inter;
	}

	void compute_BVH(){
		Vector diag = this->value->compute_diag();
		Vector mid_diag = this->value->Bmin + 0.5 * diag;
		int longest_axis = diag.argmax();
		int pivot_index = this->start;
		for(int i=this->start;i<this->end;i++){
			Vector barycenter = Mesh->compute_barycenter(i);
			if(barycenter[longest_axis]< mid_diag[longest_axis]){
				std::swap(this->Mesh->indices[i],this->Mesh->indices[pivot_index]);
				pivot_index ++;
			}
		}
		if(pivot_index<=this->start || pivot_index>=this->end-1 || (this->end)- (this->start)<5) return;

		this->lchild = new BVH(this->Mesh,this->start,pivot_index,this->id,this->texture_map);
		this->rchild = new BVH(this->Mesh,pivot_index,this->end,this->id,this->texture_map);
	}
};
