/* File : example.cxx */
#include <fstream>
#include <iostream>
#include "ode_solver.h"
#include "math.h"
#define M_PI 3.14159265358979323846

//std::ofstream cout("/home/vahid/Documents/Complex network/c/output.txt");

/*------------------------------------------------------------*/
void ODE::integrate(const dim1& iAdj) 
{   
    Cij = reshape_2d(iAdj);
    calDegree();

    dim1 y = IC;
    Order1.resize(num_steps);
    Psi1.resize(num_steps);
    Order2.resize(num_steps);
    Psi2.resize(num_steps);
    dim1 r1;
    dim1 r2;
    std::cout.precision(3);
    for (int step = 0; step < num_steps; ++step){
        
        //std::cout<<"step="<<step<<"\n";
        r1 = order_parameter(y);
        //std::cout<<"\n"<<"r1="<<r1[0];
        r2 = order_parameter_k(y);
        Order1[step] = r1[0];
        Psi1[step] = r1[1];
        Order2[step] = r2[0];
        Psi2[step] = r2[1];

        runge_kutta4_integrator(y);
    }
}
/*------------------------------------------------------------*/
void ODE::euler_integrator (dim1 &y )
{
    dim1 f(N);
    f = dydt(y);
    for (int i=0; i<y.size(); i++)
        y[i] += f[i] * dt;
}
/*------------------------------------------------------------*/
void ODE::runge_kutta4_integrator (dim1 &y) 
{
    int n = y.size();
    dim1 k1(n), k2(n), k3(n), k4(n);
    dim1 f(n);
    k1 = dydt(y);
    for (int i=0; i<n; i++)
        f[i]= y[i]+ 0.5 * dt * k1[i];
    k2 = dydt (f);
    for (int i=0; i<n; i++)
        f[i]= y[i]+ 0.5 * dt * k2[i];
    k3 = dydt(f);
    for (int i=0; i<n; i++)
        f[i]= y[i]+ dt * k3[i];
    k4 = dydt(f);
    for (int i=0; i<n; i++)
        y[i] += (k1[i] + 2.0*(k2[i] + k3[i]) + k4[i]) * dt/6.0;
}
/*------------------------------------------------------------*/
dim1 ODE::dydt(const dim1 &x)
{
    double sumx = 0.0;
    dim1 f(N);
    #pragma omp parallel for reduction(+:sumx)
    for (int i=0; i<N; i++)
    {
        sumx = 0;
        for(int j=0; j<N; j++)
        {
            if ((i!=j) && (Cij[i][j]!=0))
                sumx += Cij[i][j] * sin(x[j]-x[i]);
        }
        //f[i] = Omega[i] + K_Over_N * sumx;
        f[i] = Omega[i] + couplingStrength * sumx;
    }
    return f;
}
/*------------------------------------------------------------*/
dim2 ODE::reshape_2d(const std::vector<double>& X1D)
{
    using namespace std;
    dim2 X2D;
    X2D.resize(N);
    for (int i = 0; i < N; i++)
        X2D[i].resize(N);
    for (int i = 0; i < X1D.size(); i++)
    {
        int row = i / N;
        int col = i % N;
        X2D[row][col] = X1D[i];
    }
    return X2D;
}
/*------------------------------------------------------------*/
void ODE::calDegree(){
    Degree.resize(N);
    sum_degree = 0;
    assert (Cij.size()==N && Cij[0].size()==N);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            Degree[i] += Cij[i][j];
            if (Cij[i][j]>1e-8)
                sum_degree ++; 
        }
    if (sum_degree==0){
        std::cout<< "no link in network \n";
        exit(EXIT_FAILURE);
    }
}
/*------------------------------------------------------------*/
void ODE::set_matrices(const dim1& iAdj)
{
    Cij = reshape_2d(iAdj);
    Degree.resize(N);
    sum_degree = 0;
    assert (Cij.size()==N && Cij[0].size()==N);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            Degree[i] += Cij[i][j];
            if (Cij[i][j]>1e-8)
                sum_degree ++; 
        }
    if (sum_degree==0){
        std::cout<< "no link in network \n";
        exit(EXIT_FAILURE);
    }
}
/*------------------------------------------------------------*/
dim1 ODE::order_parameter(const dim1& x)
{
    int    n      = x.size();
    double real_R = 0.;
    double imag_R = 0.;
    for (int i=0; i<n; i++) 
    {
        real_R += cos(x[i]);
        imag_R += sin(x[i]);
    }
    real_R /= (double) n;
    imag_R /= (double) n;
    double r   = sqrt(real_R * real_R + imag_R * imag_R);
    double psi = atan2(imag_R,real_R);
    dim1 result {r, psi};
    return result;
}
/*------------------------------------------------------------*/
dim1 ODE::order_parameter_k(const dim1 &x)
{
    int n = x.size();
    double real_R = 0.;
    double imag_R = 0.;
    for (int i = 0; i < n; i++)
    {
        real_R += Degree[i] * cos(x[i]);
        imag_R += Degree[i] * sin(x[i]);
    }
    real_R /= (double)sum_degree;
    imag_R /= (double)sum_degree;
    double r = sqrt(real_R * real_R + imag_R * imag_R);
    double psi = atan2(imag_R, real_R);
    dim1 result{r, psi};
    return result;
}
/*------------------------------------------------------------*/
dim2 ODE::get_order_parameters()
{
    dim2 result{Order1, Order2};
    return result;
}
/*------------------------------------------------------------*/
// dim2   ODE::get_coordinates() 
// { 
//     return Coordinates;
// }
/*------------------------------------------------------------*/
dim2 kuramoto_correlation(const dim1& x)
{
     /* Calculate Kuramoto correlation*/
    int n = x.size();
    dim2 cor(n,dim1(n));
    
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            cor[i][j] = cos(x[j]-x[i]);
        
    return cor;
}
/*------------------------------------------------------------*/
