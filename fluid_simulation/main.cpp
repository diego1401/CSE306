#include "create_diagram.cpp"
#include <stdio.h>
#include "L-BFGS/lbfgs.h"
#include <assert.h>  
#include <random>



class objective_function
{
protected:
    lbfgsfloatval_t *m_x;

public:
    Polygon dataset; double f; Polygon SubjectPolygon; int M,N; //Air and Water Molecules
    std::vector<Vector> velocities;
    double* lambdas;
    double mass_water,mass_air;
    double* final_weights;
    // double eps = 0.004;
    double eps_inv_sq = 1./(0.004*0.004);
    double dt = 0.002;
    double mass = 200;
    Vector g = Vector(0.,-9.8,0.);
    std::vector<Polygon> scene;
    std::string color = "blue";

    void GMS_one_step(){
        // X is dataset[M,...M+N]
        // v is velocities
        //Prev opt W are final_weights
        this->run();
        Vector Fi,Fi_spring;
        for(int i=this->M;i<this->M+this->N;i++){
            Polygon cell = this->scene[i];
            Vector Xi = this->dataset.vertices[i];
            Fi_spring = this->eps_inv_sq * (cell.Centroid2d() - Xi);
            Fi = Fi_spring + mass* g;
            this->dataset.vertices[i] += dt * this->velocities[i];
            this->velocities[i] += this->dt/this->mass * Fi;
        }

    }

    void tranform_weights(const lbfgsfloatval_t *x){
        for(int i=0;i<this->M;i++){
            this->final_weights[i] = x[this->N];
        }
        for(int i=this->M;i<(this->M+this->N);i++){
            this->final_weights[i] = x[i-this->M];
        }
    }

    void color_scene(){
        for(int i=0;i<this->scene.size();i++){
            if(this->dataset.vertices[i].ret_is_liquid()){
                this->scene[i].is_wet = true;
            }
            else{
                this->scene[i].is_wet = false;
            }
        }
    }
    
    objective_function(int _M,int _N,double _f,double sig) : m_x(NULL){
    std::cout << "Initialization is starts!" << std::endl;
    this->f = _f; this->M = _M; this->N = _N; this->SubjectPolygon = get_bounding_box();
    this->lambdas =(double*) malloc((N+1)*sizeof(double));
    //transformed weigths 
    this->final_weights =(double*) malloc((this->N+this->M)*sizeof(double));
    //Initial circle of water
    // double area_water = (double)  this->N / (this->N+this->M);
    double R = 0.3;
    double area_water = M_PI *R*R;
    Vector C(0.5,0.5,0.);  
    
    this->mass_water = area_water/this->N;

    //Initialize air cells
    // double total = 0;
    for(int i =0;i< M;i++){
        Vector tmp((double) rand()/RAND_MAX,(double) rand()/RAND_MAX,0.);
        // this->lambdas[i] = exp(-(tmp-C).norm_squared()/(sig)); total += this->lambdas[i]; 
        this->dataset.vertices.push_back(tmp);
    }
    //Lloyd iterations over the Air cells
    Centroidal_Voronoi_Tesselation(this->dataset);
    
    //Init Water cells
    for(int i =0;i<this->N;i++){
        double r = (double) rand()/RAND_MAX; r = R*sqrt(r);
        double theta = (double) rand()/RAND_MAX * 2*M_PI;
        double x = r * cos(theta); double y = r*sin(theta);
        Vector tmp(x,y,0); tmp+= C; tmp.set_is_liquid(true);
        //Water cells start from the index M
        this->dataset.vertices.push_back(tmp);
        // this->lambdas[i] = this->mass_water;
        //Note we have a correspondance dataset[M] -> lamdas[0], ... ,dataset[M+N-1] -> lamdas[N-1]
        
    }
    this->mass_air = 1. - area_water;
    std::cout << "Initialization is done!" << std::endl;

    //init final_weights
    for(int i=0;i<this->N+this->M;i++){
        this->final_weights[i] = 1.;
    }
    //Normalize
    // for(int i=0;i<this->M;i++){
    //     this->lambdas[i] /= total;
    // }
    
    }

    virtual ~objective_function(){
        if (m_x != NULL) {
            lbfgs_free(m_x);
            m_x = NULL;
        }
        if(this->lambdas!= NULL){
            free(this->lambdas);
         }
        if(this->final_weights!=NULL){
            free(this->final_weights);
        } 
    }

    int run()
    {
        lbfgsfloatval_t fx;
    lbfgsfloatval_t *m_x = lbfgs_malloc(this->N + 1);
    

    if (m_x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
        return 1;
    }

    /* Initialize the variables. */
    for (int i = 0;i < (this->N);i ++) {
        // m_x[i] = 1./((double)N);
        // m_x[i] = this->lambdas[i];
        m_x[i] = this->mass_water;
    }
    m_x[this->N] = this->mass_air;

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);

    param.max_iterations = 300;

    /*
        Start the L-BFGS optimization; this will invoke the callback functions
        evaluate() and progress() when necessary.
    */
    int ret = lbfgs(this->N+1, m_x, &fx, _evaluate, _progress, this, &param);

    /* Report the result. */

    printf("L-BFGS optimization terminated with status code = %d\n", ret);
    std::cout << lbfgs_strerror(ret) << std::endl;
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, m_x[0], m_x[1]);

    this->tranform_weights(m_x);

    std::cout << "Creating Diagram" << std::endl; 
    this->scene = Create_diagram(this->dataset,this->final_weights); std::cout << "Finished" << std::endl;
    this->color_scene(); 
    save_svg(this->scene,"first_frame.svg",this->dataset,this->color);

    // return ret;
    return 0;
    }

protected:
    static lbfgsfloatval_t _evaluate(void *instance,const lbfgsfloatval_t *x,lbfgsfloatval_t *g,const int n,const lbfgsfloatval_t step)
    {
        return reinterpret_cast<objective_function*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(const lbfgsfloatval_t *x,lbfgsfloatval_t *g,const int n,const lbfgsfloatval_t step){
    //Compute G and Delta G
    lbfgsfloatval_t fx = 0.0;
    this->tranform_weights(x);
    lbfgsfloatval_t AreaAir = 0.;
    lbfgsfloatval_t Integral_Air = 0.;
    #pragma omp parallel for
    for(int i=0;i<this->dataset.vertices.size();i++){
        Vector Pi = this->dataset.vertices[i];
        Polygon cell = PowerCell_of_i(this->SubjectPolygon,i,this->dataset,this->final_weights); double Area = cell.Area2D();
        
        lbfgsfloatval_t tmp = 0.0;
        //computing the i-th term of the sum for g(W)
        for(int k=0;k<cell.vertices.size();k++){
            Vector A = cell.vertices[(k>0)?(k-1):cell.vertices.size()-1];
            Vector B = cell.vertices[k]; 
            double xk1=A[0];double xk=B[0];
            double yk1=A[1];double yk=B[1];

            tmp += (xk1*yk -xk*yk1) * (xk1*xk1 + xk1*xk + xk*xk
                                      +yk1*yk1 + yk1*yk + yk*yk
                                    -4.* (Pi[0]*(xk1+xk) + Pi[1]* (yk1 + yk))
                                    + 6.* Pi.norm_squared() );
        }
        tmp = std::abs(tmp)/12.;

        // fx += this->lambdas[i]* x[i] -this->f* x[i] * Area;

        //computing Delta g(W) along i;
        // g[i] = this->f*Area  - this->lambdas[i] ;

        
        //The above is common to both cases,It is different when the weight is involved
        if(i>(this->M-1)){
            int index = i - this->M; //we are dealing with a water molecule
            fx +=  this->f * tmp+ this->mass_water * x[index] - this->f * Area * x[index];
            g[index] = this->f * Area - this->mass_water;
        } 
        else{
            AreaAir += Area;
            Integral_Air += tmp;
        }
    }
    //Contribution Air
    fx += this->mass_air * x[this->N] +this->f * (Integral_Air - x[this->N]  * AreaAir);
    g[this->N] = this->f * AreaAir - this->mass_air;
    
    fx *= -1.;
    return fx;
    }

    static int _progress(void *instance,const lbfgsfloatval_t *x,const lbfgsfloatval_t *g,const lbfgsfloatval_t fx,const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,const lbfgsfloatval_t step,int n,int k,int ls){
        return reinterpret_cast<objective_function*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(const lbfgsfloatval_t *x,const lbfgsfloatval_t *g,const lbfgsfloatval_t fx,const lbfgsfloatval_t xnorm,const lbfgsfloatval_t gnorm,const lbfgsfloatval_t step,
        int n,int k,int ls){
        printf("Iteration %d:\n", k);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
        printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        printf("\n");
        return 0;
    }
};


int main(){
    // Vector e1(0.2,0.2,0),e2(0.8,0.2,0),e3(0.1,0.7,0),e4(0.9,0.6,0);
    // Polygon* P = new Polygon;
    // P->vertices = {e1,e2,e4,e3};
    // Vector f1(0.3,0.3,0),f2(0.8,0.3),f3(0.8,0.7),f4(0.3,0.7);
    // Polygon* clip = new Polygon;
    // clip->vertices = {f1,f2,f3,f4};
    // save_svg({polygon_clipping(*P,*clip)},"image_clipped.svg","red");
    // save_svg({*P,*clip},"image.svg");

    // Vector e1(0.2,0.2,0),e2(0.8,0.2,0),e3(0.1,0.9,0),e4(0.9,0.6,0),
    //        e5(0.55,0.5,0),e6(0.4,0.4,0);
    // Polygon dataset({e1,e2,e3,e4,e5,e6});
    // int N = dataset.vertices.size(); 
    // double weights[N] = {1,1,1,1,1,1}; // have all weights equal to go back
    

    // save_svg({dataset},"polygon.svg");
    // save_svg(Create_diagram(dataset,weights),"image_voronoi.svg",dataset);

    
    int M = 2500;
    int N = 700;
    objective_function obj(M,N,1.,0.02);
    obj.run();
    
    return 0;
}