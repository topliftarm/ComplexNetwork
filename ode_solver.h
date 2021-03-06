/* File : example.h */
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <assert.h>
typedef std::vector<double> dim1;
typedef std::vector<std::vector<double> > dim2;

using namespace std;


class ODE {
  private:
  const int N;
  const double dt;
  double tfinal;
  double K_Over_N;
  int num_steps;
  dim2 Coordinates;
  dim1 Omega;
  int sum_degree, NumbertOfSteps, NumberOfIterations;
  double couplingStrength;
  dim1 createSelfishList(int NumberOfSelfishNodes, dim1 NodesOrder);

public:
  ODE(int iN, double itfinal,
      double idt,
      double iCouplingStrength, dim1 iIC, dim1 iOmega, int _NumbertOfSteps, int _NumberOfIterations) : N(iN), dt(idt)
  {
    tfinal = itfinal;
    num_steps = int(tfinal/dt);
    IC = iIC;
    Omega = iOmega;
    //K_Over_N = iK / (N+0.0);
    couplingStrength = iCouplingStrength;
    NumbertOfSteps = _NumbertOfSteps;
    NumberOfIterations = _NumberOfIterations;

  }
  virtual ~ODE() { }
  dim2 Cij;
  dim1 Order1;
  dim1 Order2;
  dim1 Psi1;
  dim1 Psi2;
  dim1 Degree;
  dim1 AcceptanceRateRewiring;
  void euler_integrator(dim1 & );
  dim1 runge_kutta4_integrator(dim1 y, dim2 CijLocal);
  dim2 rewiring(int indexFocusNode, dim1 nodesOrder, dim2);
  void Print2D(dim2);
  void Print1D(dim1);

  dim1 FinalY;
//   dim2 get_coordinates();
  dim2 reshape_2d(const dim1& X1D);

  double Mean(const dim1&, int);
  dim1 getMeanYPrime();
  dim1 getFinalY();
  dim1 runDynamics(int, dim2, dim1, dim1&, int);
  int calBiEdges(dim2);

  dim1 MeanRinEachIteration;
  void integrate(const dim1& iAdj, bool rewire, string currentPath, int NumberOfSelfishNodes);
  dim1 dydt(const dim1 &x, dim2 CijLocal);
  dim1 order_parameter(const dim1& x);
  dim1 order_parameter_k(const dim1 &x);
  dim2 get_order_parameters();
  dim1 getAcceptanceRewiring();
  dim1 getMeanRinEachIteration();
  void calDegree();
  void set_matrices(const dim1& iAdj);
  dim1 getCij();
  void saveMatrix(string fileName, dim2 data);
  void saveArray(string fileName, dim1 data);

  dim1 IC;
};

dim2 kuramoto_correlation(const dim1& x);
