#include "create_diagram.cpp"
#include <stdio.h>
#include "L-BFGS/lbfgs.h"
#include <assert.h>  
#include <random>
#include <cmath>
#include <chrono>
using namespace std::chrono;


class objective_function
{
protected:
    lbfgsfloatval_t *m_x;

public:
    Polygon dataset; double f; Polygon SubjectPolygon; int M,N; //Air and Water Molecules
    std::vector<Vector> velocities;
    double* lambdas;
    int mode;
    double mass_water,mass_air;
    double* final_weights;
    std::string filename;
    int frame_id = 0;
    // double eps = 0.004;
    double eps_inv_sq = 1./(0.004*0.004);double dt = 0.002; double mass = 200; Vector g = Vector(0.,-9.8,0.);
    std::vector<Polygon> scene;
    std::string color = "blue";

    void GMS_one_step(){
        // X is dataset[M,...M+N]
        // v is velocities
        //Prev opt W are final_weights
        // int ret = this->run();
        Vector Fi,Fi_spring;
        std::cout << "Updating physical quantities" << std::endl;
        this->velocities[this->M+1].print_vector();
        for(int i=this->M;i<this->M+this->N;i++){
            // std::cout << "Cell" << std::endl;
            Polygon cell = this->scene[i];
            // std::cout << "Xi" << std::endl;
            Vector Xi = this->dataset.vertices[i];
            // std::cout << "ok" << std::endl;
            Fi_spring = this->eps_inv_sq * (cell.Centroid2d() - Xi);
            // Fi = Fi_spring + this->mass* this->g;
            Fi = Fi_spring;
            // Fi.print_vector();
            this->dataset.vertices[i] += dt * this->velocities[i];
            this->velocities[i] += (this->dt/this->mass) * Fi;
        }

    }

    void tranform_weights(const lbfgsfloatval_t *x){
        for(int i=0;i<this->M+this->N;i++){
            if(i>(this->M-1)){
                this->final_weights[i] = x[i-this->M];
            }
            else{
                this->final_weights[i] = x[this->N];
            }
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
    
    objective_function(int _mode,std::string _filename,int _M,int _N,double _f,double sig) : m_x(NULL){
    this->mode = _mode;
    if(mode==0){
        Voronoi_Diagram(_filename,_N);
    }
    if(mode==1){
        Power_Diagram(_filename,_N,_f,sig);
    }
    if(mode==2){
        fluid_simulation(_filename,_M,_N,_f,sig);
    }
    
    }

    void fluid_simulation(std::string _filename,int _M,int _N,double _f,double sig)
    {
    std::cout << "Initialization is starts!" << std::endl;
    this->filename = _filename;
    this->f = _f; this->M = _M; this->N = _N; this->SubjectPolygon = get_bounding_box();
    this->lambdas =(double*) malloc((N+1)*sizeof(double));
    //transformed weigths 
    this->final_weights =(double*) malloc((this->N+this->M)*sizeof(double));
    //Initial circle of water
    double R = 0.3; double area_water = M_PI *R*R; Vector C(0.5,0.5,0.);  
    
    this->mass_water = area_water/this->N;
    std::vector<Vector> data_vec_air(this->M);
    Polygon tmp_data(data_vec_air);
    //Initialize air cells
    // double total = 0;
    for(int i =0;i< M;i++){
        double x = (double) rand()/RAND_MAX;
        double y = (double) rand()/RAND_MAX;
        Vector tmp(x,y,0.);
        // this->lambdas[i] = exp(-(tmp-C).norm_squared()/(sig)); total += this->lambdas[i]; 
        // this->dataset.vertices.push_back(tmp);
        tmp_data.vertices[i] = tmp;
    }
    
    //Lloyd iterations over the Air cells
    Centroidal_Voronoi_Tesselation(tmp_data);
    std::vector<Vector> data_vec(this->M+this->N);
    Polygon data(data_vec);
    for(int i=0;i<M;i++){
        data.vertices[i] = tmp_data.vertices[i];
    }
    
    //Init Water cells
    for(int i =this->M;i<this->N+this->M;i++){
        double r = (double) rand()/RAND_MAX; r = R*sqrt(r);
        double theta = (double) rand()/RAND_MAX * 2*M_PI;
        double x = r * cos(theta); double y = r*sin(theta);
        Vector tmp(x,y,0); tmp+= C; tmp.set_is_liquid(true);
        //Water cells start from the index M
        // this->dataset.vertices.push_back(tmp);
        data.vertices[i] = tmp;
        //Note we have a correspondance dataset[M] -> weight[0], ... ,dataset[M+N-1] -> weight[N-1]
        
    }

    this->dataset = data;

    this->mass_air = 1. - area_water;
    std::cout << "Initialization is done!" << std::endl;

    std::vector<Vector> tmp_vel(this->N+this->M);
    //init final_weights init velocities
    for(int i=0;i<this->N+this->M;i++){
        this->final_weights[i] = 1.;
        Vector tmp(0.,0.,0.);
        tmp_vel[i] = tmp;
    }
    this->velocities = tmp_vel;
    
    }

    virtual ~objective_function(){
        if(this->mode==1) free(this->lambdas);
        else{
            free(this->final_weights);
        }
        
    }

    int run(int frames)
    {
    if(mode==1){
    std::cout<< "running" << std::endl;
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *m_x = lbfgs_malloc(this->N);
    

    if (m_x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
        return 1;
    }

    /* Initialize the variables. */
    for (int i = 0;i < (this->N);i ++) {
        m_x[i] = this->lambdas[i];
    }

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);

    param.max_iterations = 3000;
    /*
    Start the L-BFGS optimization; this will invoke the callback functions
    evaluate() and progress() when necessary.
    */
    int ret = lbfgs(this->N, m_x, &fx, _evaluate, _progress, this, &param);

    /* Report the result. */

    printf("L-BFGS optimization terminated with status code = %d\n", ret);
    std::cout << lbfgs_strerror(ret) << std::endl;
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, m_x[0], m_x[1]);

    this->scene = Create_diagram(this->dataset,m_x);
    save_frame(this->scene,this->filename,N);
   
    if (m_x != NULL) {
            lbfgs_free(m_x);
            m_x = NULL;
    }
    std::cout << "Freed" << std::endl;
    return 0;

    }
    else{
    std::cout<< "running" << std::endl;
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *m_x = lbfgs_malloc(this->N + 1);
    

    if (m_x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
        return 1;
    }

    /* Initialize the variables. */
    for (int i = this->M;i < (this->M+this->N);i ++) {
        m_x[i] = this->final_weights[i];
    }
    m_x[this->N] = this->final_weights[0];

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);

    param.max_iterations = 3000;
    for(int i=0;i<frames;i++){
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
        this->scene = Create_diagram(this->dataset,this->final_weights); this->color_scene(); 
        save_frame(this->scene,this->filename,this->M,this->frame_id);
        
        this->frame_id++;
        std::cout << "Frame" << this->frame_id << "saved" << std::endl;
        this->GMS_one_step();
    }
   
    if (m_x != NULL) {
            lbfgs_free(m_x);
            m_x = NULL;
    }
    std::cout << "Freed" << std::endl;
    return 0;
    }
    }

    void Voronoi_Diagram(std::string _filename,int N){
        this->filename = _filename;
        this->final_weights =(double*) malloc((N)*sizeof(double));
        std::vector<Vector> data(N); Polygon tmp_data(data);
        for(int i =0;i< N;i++){
        double x = (double) rand()/RAND_MAX; double y = (double) rand()/RAND_MAX;
        Vector tmp(x,y,0.); tmp_data.vertices[i] = tmp;
        this->final_weights[i] = 1.;
        }
        this->dataset = tmp_data;
        std::vector<Polygon> PowerCells = Create_diagram(this->dataset,this->final_weights);
        this->scene = Create_diagram(this->dataset,this->final_weights);
        save_frame(this->scene,this->filename,N+1);
    }

    void Power_Diagram(std::string _filename,int _N,double _f,double sig){
        std::cout << "Initialization is starts!" << std::endl;
        this->filename = _filename;
        this->f = _f; this->N = _N; this->SubjectPolygon = get_bounding_box();
        this->lambdas =(double*) malloc((this->N)*sizeof(double));
        
        std::vector<Vector> data(this->N);
        Polygon tmp_data(data);
        Vector C(0.5,0.5,0); double total = 0;
        for(int i =0;i< this->N;i++){
            double x = (double) rand()/RAND_MAX;
            double y = (double) rand()/RAND_MAX;
            Vector tmp(x,y,0.);
            this->lambdas[i] = exp(-(tmp-C).norm_squared()/(sig)); total += this->lambdas[i]; 
            tmp_data.vertices[i] = tmp;
        }

        for(int i=0;i<this->N;i++){
            this->lambdas[i] /= total;
        }
        this->dataset = tmp_data;

        std::cout << "Initialization is done!" << std::endl;
        }

protected:
    static lbfgsfloatval_t _evaluate(void *instance,const lbfgsfloatval_t *x,lbfgsfloatval_t *g,const int n,const lbfgsfloatval_t step)
    {
        return reinterpret_cast<objective_function*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(const lbfgsfloatval_t *x,lbfgsfloatval_t *g,const int n,const lbfgsfloatval_t step){
    if(this->mode==1){
    //Compute G and Delta G
    lbfgsfloatval_t fx = 0.0;
    std::vector<Polygon> PowerCells = Create_diagram(this->dataset,(double*)x);
    for(int i=0;i<this->N;i++){
        Vector Pi = this->dataset.vertices[i];
        Polygon cell = PowerCells[i];
        double Area = cell.Area2D();
        
        lbfgsfloatval_t tmp = 0.0;
        //computing the i-th term of the sum for g(W)
        for(int k=0;k<cell.vertices.size();k++){
            Vector A = cell.vertices[(k>0)?(k-1):cell.vertices.size()-1];
            Vector B = cell.vertices[k]; 
            double xk1=A[0];double xk=B[0];
            double yk1=A[1];double yk=B[1];

            tmp += (xk1*yk -xk*yk1) * (xk1*xk1 + xk1*xk + xk*xk
                                      +yk1*yk1 + yk1*yk + yk*yk
                                    - 4.* (Pi[0]*(xk1+xk) + Pi[1]* (yk1 + yk))
                                    + 6.* Pi.norm_squared() );
        }
        tmp = std::abs(tmp)/12.;

        fx += this->lambdas[i]*x[i] + this->f * (tmp - x[i]*Area);
        //computing Delta g(W) along i;
        g[i] = this->f*Area  - this->lambdas[i] ;
    }
    fx *= -1.;
    return fx;
    }
    else{
    //Compute G and Delta G
    lbfgsfloatval_t fx = 0.0;
    this->tranform_weights(x);
    lbfgsfloatval_t AreaAir = 0.;
    lbfgsfloatval_t Integral_Air = 0.;
    std::vector<Polygon> PowerCells = Create_diagram(this->dataset,this->final_weights);
    for(int i=0;i<this->N+this->M;i++){
        Vector Pi = this->dataset.vertices[i];
        // Polygon cell = PowerCell_of_i(this->SubjectPolygon,i,this->dataset,this->final_weights); 
        Polygon cell = PowerCells[i];
        double Area = cell.Area2D();
        
        lbfgsfloatval_t tmp = 0.0;
        //computing the i-th term of the sum for g(W)
        for(int k=0;k<cell.vertices.size();k++){
            Vector A = cell.vertices[(k>0)?(k-1):cell.vertices.size()-1];
            Vector B = cell.vertices[k]; 
            double xk1=A[0];double xk=B[0];
            double yk1=A[1];double yk=B[1];

            tmp += (xk1*yk -xk*yk1) * (xk1*xk1 + xk1*xk + xk*xk
                                      +yk1*yk1 + yk1*yk + yk*yk
                                    - 4.* (Pi[0]*(xk1+xk) + Pi[1]* (yk1 + yk))
                                    + 6.* Pi.norm_squared() );
                                    // + 6.* (Pi[0] *Pi[0] + Pi[1]*Pi[1]));
        }
        tmp = std::abs(tmp);

        // fx += this->lambdas[i]* x[i] -this->f* x[i] * Area;

        //computing Delta g(W) along i;
        // g[i] = this->f*Area  - this->lambdas[i] ;

        
        //The above is common to both cases,It is different when the weight is involved
        if(i>(this->M-1)){
            int index = i - this->M; //we are dealing with a water molecule
            // assert(0<= index && index < this->N && Pi.ret_is_liquid());
            fx +=  this->f  *(tmp/12. - x[index]* Area) + this->mass_water * x[index];
            g[index] = this->f * Area - this->mass_water;
        } 
        else{
            // assert(0<= i && i < this->M && !Pi.ret_is_liquid());
            AreaAir += Area;
            Integral_Air += tmp;
        }
    }
    //Contribution Air
    // std::cout << "Integral Air=" << Integral_Air << std::endl;
    // std::cout << "Area Air=" << AreaAir << std::endl;
    // std::cout << "fx=" << fx << std::endl;
    Integral_Air /= 12.;
    fx += this->mass_air * x[this->N] + this->f * (Integral_Air - x[this->N] * AreaAir);
    g[this->N] = this->f * AreaAir - this->mass_air;
    
    fx *= -1.;
    return fx;
    }
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

    //Voronoi
    // int N = 100;
    // high_resolution_clock::time_point t1 = high_resolution_clock::now();
    // objective_function obj(0,"Voronoi_30k",0,N,0,0);
    // high_resolution_clock::time_point t2 = high_resolution_clock::now();


    //Power Diagram
    int N = 2000;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    objective_function obj(1,"Power_Diagram",0,N,1,0.02);
    obj.run(0);
    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    //Fluid simulation
    // int M = 250; int N = 70;
    // objective_function obj(2,"frames_test/simulation",M,N,1.,0.02);
    // obj.run(1000);


    
    // d = Sequential_Dijkstra_Two_Queue(gra, 0);
    
    duration<double> time_span1 = duration_cast<duration<double> >(t2 - t1);
    std::cout << "Time of exec: " << time_span1.count() << std::endl;

    return 0;
}
