/* File : example.h */
#include <vector>
#include <iostream>
#include <algorithm>
#include <assert.h>

typedef std::vector<double> dim1;
typedef std::vector<std::vector<double> > dim2;



class ODE {
  private:
  const int N;
  const double dt;
  double tfinal;
  double K_Over_N;
  int num_steps;
  dim2 Coordinates;
  dim1 Omega;
  int sum_degree;
  double couplingStrength;

public:
  ODE(int iN, double itfinal, 
      double idt, 
      double iCouplingStrength, dim1 iIC, dim1 iOmega) : N(iN), dt(idt) 
  {
    tfinal = itfinal;
    num_steps = int(tfinal/dt);
    IC = iIC;
    Omega = iOmega;
    //K_Over_N = iK / (N+0.0);
    couplingStrength = iCouplingStrength;
  }
  virtual ~ODE() { }
  dim2 Cij;
  dim1 Order1;
  dim1 Order2;
  dim1 Psi1;
  dim1 Psi2;
  dim1 Degree;
  void euler_integrator(dim1 & );
  void runge_kutta4_integrator(dim1 &);

//   dim2 get_coordinates();
  dim2 reshape_2d(const dim1& X1D);

  void integrate(const dim1& iAdj);
  dim1 dydt(const dim1 &x);
  dim1 order_parameter(const dim1& x);
  dim1 order_parameter_k(const dim1 &x);
  dim2 get_order_parameters();
  void calDegree();
  void set_matrices(const dim1& iAdj);

  dim1 IC;
};

dim2 kuramoto_correlation(const dim1& x);
