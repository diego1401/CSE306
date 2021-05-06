#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../Raytracer_Assignement1/stb-master/stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "../Raytracer_Assignement1/stb-master/stb_image.h"

//We are going to use vectors
#include "../Raytracer_Assignement1/classes/vector.cpp"
class ImageRetargeting{
public:
    unsigned char* image; int W; int H;
    double* C;
    int* seam; //array of best values of x, values of y are deduced from the index
    ImageRetargeting(unsigned char* _image, int _W, int _H,int N){
        this->image = _image; this->W = _W;this->H = _H;
        C = (double*) malloc(W*H*sizeof(double));
        seam = (int*) malloc(H*sizeof(int));
        for(int i = 0;i<this->H;i++){
            for(int j = 0;j<this->W;j++){
                C[i*this->W + j] = -1.;
            }
        }
        for(int i=0;i<N;i++){
        construct_seam();
        delete_seam();
        }
    }

    double I(int x,int y){
        if(x<0 or x>=this->W) return 0;
        if(y<0 or y>=this->H) return 0;
        return  this->image[y*this->W*3+x*3 + 0] +
                this->image[y*this->W*3+x*3 + 1] +
                this->image[y*this->W*3+x*3 + 2];
    }

    double E(int x,int y){
        if(x<0 or x>=this->W) return 0;
        if(y<0 or y>=this->H) return 0;
        return std::abs(I(x+1,y)-I(x-1,y))
            +  std::abs(I(x,y+1)-I(x,y-1));
    }

    double Cumulated_Energy_Map(int x,int y){
        if(x<0 or x>=this->W) return INFINITY;
        else if(y<0 or y>=this->H) return INFINITY;
        else if(y==0){
            C[y*this->W + x] = E(x,y);
        }
        else if(C[y*this->W + x]==-1){
            C[y*this->W + x] = std::min({Cumulated_Energy_Map(x-1,y-1),
                                         Cumulated_Energy_Map(x,  y-1),
                                         Cumulated_Energy_Map(x+1,y-1)})
                                + E(x,y);
        }
        return C[y*this->W + x];
    }

    int get_start_seam(){
        double min = INFINITY;
        int best_x = 0; double tmp;
        for (int i=0;i<this->W;i++){
            tmp = this->Cumulated_Energy_Map(i,H-1);
            if(tmp<min) best_x = i;min = tmp;
        }
        return best_x;
    }

    void construct_seam(){
        int x = get_start_seam();
        seam[H-1] = x;
        double c1,c2,c3;
        double min;
        for (int y = H-1; y>0; y--)
        {
            c1 = Cumulated_Energy_Map(x-1,y-1);
            c2 = Cumulated_Energy_Map(x,y-1);
            c3 = Cumulated_Energy_Map(x+1,y-1);
            min = std::min({c1,c2,c3});
            if(min==c1) x -= 1;
            else if(min==c3) x +=1;
            seam[y-1] = x;

        }
        
    }

    void delete_seam(){
        for (int j = H-1; j>=0; j--)
        {
            for(int i=seam[j];i<W-1;i++){
                this->image[j*this->W*3+i*3 + 0] = this->image[j*this->W*3+(i+1)*3 + 0];
                this->image[j*this->W*3+i*3 + 1] = this->image[j*this->W*3+(i+1)*3 + 1];
                this->image[j*this->W*3+i*3 + 2] = this->image[j*this->W*3+(i+1)*3 + 2];
            }
        }
        for(int i = 0;i<this->H;i++){
            for(int j = 0;j<this->W;j++){
                C[i*this->W + j] = -1.;
            }
        }
        // this->W--;
    }
};




int main(){
    const char* filename = "/Users/Diego/Documents/SEM6/CSE306/Raytracer_Assignement1/report/images/DoF_image.png";
    //Dimensions, TO BE CHANGED
    int H = 512; int W = 512; int c = 0;int channels = 3;

    unsigned char* image = stbi_load(filename,&H,&W,&c,channels);
    int N = 100;
    ImageRetargeting Map(image,W,H,N);
    //In order to resize the image easily
    W-=N;
    unsigned char* new_image = (unsigned char*) malloc((W)*H*3*sizeof(unsigned char));

    for(int i=0;i<W;i++){
        for(int j=0;j<H;j++){
            new_image[j*W*3+i*3 + 0] = image[j*(W+N)*3+(i)*3 + 0];
            new_image[j*W*3+i*3 + 1] = image[j*(W+N)*3+(i)*3 + 1];
            new_image[j*W*3+i*3 + 2] = image[j*(W+N)*3+(i)*3 + 2];
        }
    }
    // W+=N;
    stbi_write_png("image.png", W, H, channels, new_image, 0); 
}