/* File : example.cxx */
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "ode_solver.h"
#include "math.h"
#include <numeric>
#include <random>
#include <cstdio>
#define M_PI 3.14159265358979323846

//using namespace std;

double ODE::Mean(const dim1& vector, int lastChunk=0){
    double mean;
    if (lastChunk != 0){
        double sum = std::accumulate(vector.end()-lastChunk, vector.end(), 0.0);
        mean = sum / lastChunk ;
    }
    else{
        double sum = std::accumulate(vector.begin(), vector.end(), 0.0);
        mean = sum / vector.size();
    }

    return mean;
}
/*------------------------------------------------------------*/
dim1 ODE::createSelfishList(int NumberOfSelfishNodes, dim1 NodesOrder){
   dim1 selfidhNodes;
   std::random_shuffle(NodesOrder.begin(), NodesOrder.end());
   for(int i=0; i<NumberOfSelfishNodes; i++)
     selfidhNodes.push_back(NodesOrder[i]);

   return selfidhNodes;
}
/*------------------------------------------------------------*/
void ODE::integrate(const dim1& iAdj, bool rewire=false, bool selfish=false, int NumberOfSelfishNodes=0)
{
    //freopen("output.txt","w",stdout);
    //ofstream logfile;
    //logfile.open ("out.txt");
    int TotalRewiring, TotalAcceptedRewiring;
    dim2 NewCij;
    Cij = reshape_2d(iAdj);
    ///calDegree();
    int sumAcceptanceRewirig=0;
    std::cout.precision(10);
    dim1 y, NewY, OldAcceptedY, MeanYPrimeAccepted;
    dim1 RGlobalBeforeRewiring, RGlobalAfterRewiring;
    double averageRbeforRewiring, averageRafterRewiring;
    double OmegaGama, OmegaGamaPrime;
    std::srand(unsigned(std::time(0)));
    int randomNode;
    dim1 nodesOrder, CopyOfNodesOrder;
    for(int i = 0; i < N; i++)
        nodesOrder.push_back(i);
    CopyOfNodesOrder = nodesOrder;
    if(!rewire){
	      y = IC;
        for(int i=0; i<NumberOfIterations; i++){
 	           y = runDynamics(NumbertOfSteps, Cij, y, MeanYPrimeAccepted);
             averageRafterRewiring = Mean(Order1,ceil(NumbertOfSteps/2));
             MeanRinEachIteration.push_back(averageRafterRewiring);
	       }//end for
    }// end !rewire
    else if(!selfish){
          TotalAcceptedRewiring = 0;
	        int TotalSelfishRewiringAccepted = 0;
	        TotalRewiring = 0;
          y = runDynamics(NumbertOfSteps, Cij, IC, MeanYPrimeAccepted);
	        OmegaGama = Mean(MeanYPrimeAccepted, ceil(NumbertOfSteps/2));
          OldAcceptedY = y;
          averageRbeforRewiring = Mean(Order1, ceil(NumbertOfSteps/2));
          MeanRinEachIteration.push_back(averageRbeforRewiring);
	        dim1 selfishNodes;
	        selfishNodes = createSelfishList(NumberOfSelfishNodes, CopyOfNodesOrder);
	        std::cout<<" Number Of Selfish Nodes = "<<NumberOfSelfishNodes<<"\n";
          std::cout<<" Selfish Nodes : \n";
	        //logfile<<" Number Of Selfish Nodes = "<<NumberOfSelfishNodes<<"\n";
	        //logfile<<std::cout<<" Selfish Nodes : \n";
	        Print1D(selfishNodes);
          for(int i=0; i<NumberOfIterations; i++){
              std::cout<<" -- Iteration = "<<i<<"\n";
              //logfile<<" -- Iteration = "<<i<<"\n";
              std::random_shuffle(CopyOfNodesOrder.begin(), CopyOfNodesOrder.end());
              randomNode = CopyOfNodesOrder[0];
              NewCij = rewiring(randomNode, nodesOrder, Cij);
              TotalRewiring++;
              NewY = runDynamics(NumbertOfSteps, NewCij, y, MeanYPrimeAccepted);
              OmegaGamaPrime = Mean(MeanYPrimeAccepted, ceil(NumbertOfSteps/2));
              averageRafterRewiring = Mean(Order1,ceil(NumbertOfSteps/2));
              MeanRinEachIteration.push_back(averageRafterRewiring);
              if(std::find(selfishNodes.begin(), selfishNodes.end(), randomNode)!= selfishNodes.end()){//Selfish check
		              std::cout<<"Checking Selfish Node "<<randomNode<<" ";
	                 //logfile<<"Checking Selfish Node "<<randomNode<<" ";
                 if(abs(Omega[randomNode]-OmegaGamaPrime)<abs(Omega[randomNode]-OmegaGama)){
                     OmegaGama = OmegaGamaPrime;
                     Cij = NewCij;
                     y = NewY;
                     sumAcceptanceRewirig++;
		                 std::cout<<"Accepted ++ \n";
		                 //logfile<<"Accepted ++\n";
		                 TotalSelfishRewiringAccepted++;
                 }else std::cout<<"NotAccepted \n";
              }else{ //NonSelfish check
                    if (averageRbeforRewiring < averageRafterRewiring){
                       Cij = NewCij;
                       averageRbeforRewiring = averageRafterRewiring;
                       sumAcceptanceRewirig++;
                       TotalAcceptedRewiring++;
                       y = NewY;
                       OldAcceptedY = y;
                       OmegaGama = OmegaGamaPrime;
                       }else
                       y = OldAcceptedY;
              }
              if(i % 100 == 0){
                  AcceptanceRateRewiring.push_back(sumAcceptanceRewirig);
                  sumAcceptanceRewirig = 0;
              }
              std::cout<<"TotalSelfishRewiringAccepted = "<<TotalSelfishRewiringAccepted<<"\n";
              //logfile<<"Total Selfish Rewiring Accepted = "<<TotalSelfishRewiringAccepted<<"\n";
          }//end !selfish
      }else if (selfish){//selfish
          std::cout<<" ** in the selfish **"<<"\n";
          //ofstream AcceptanceRewiringFile;
          //AcceptanceRewiringFile.open("AcceptanceRewiring.txt");
          sumAcceptanceRewirig = 0;
          dim1 MeanYPrimeTempAfterRewiring;
          y = runDynamics(NumbertOfSteps, Cij, IC, MeanYPrimeAccepted);
          OmegaGama = Mean(MeanYPrimeAccepted, ceil(NumbertOfSteps/2));
          averageRbeforRewiring = Mean(Order1, ceil(NumbertOfSteps/2));
          MeanRinEachIteration.push_back(averageRbeforRewiring);
          for(int ii=0; ii<NumberOfIterations; ii++){
            std::cout<<"Iteration = "<<ii<<"\n";
            std::random_shuffle(CopyOfNodesOrder.begin(), CopyOfNodesOrder.end());
            randomNode = CopyOfNodesOrder[0];
            NewCij = rewiring(randomNode, nodesOrder, Cij);
            NewY = runDynamics(NumbertOfSteps, NewCij, y, MeanYPrimeTempAfterRewiring);
            OmegaGamaPrime = Mean(MeanYPrimeTempAfterRewiring, ceil(NumbertOfSteps/2));
            averageRafterRewiring = Mean(Order1,ceil(NumbertOfSteps/2));
            MeanRinEachIteration.push_back(averageRafterRewiring);
            if(abs(Omega[randomNode]-OmegaGamaPrime)<abs(Omega[randomNode]-OmegaGama)){
              OmegaGama = OmegaGamaPrime;
              Cij = NewCij;
              y = NewY;
              MeanYPrimeAccepted = MeanYPrimeTempAfterRewiring;
              sumAcceptanceRewirig++;
              //std::cout<<"---- sum Accepted Rewiring = "<<sumAcceptanceRewirig<<"\n";
            }
            if( (ii % 100 == 0) || (ii == NumberOfIterations) ){
                  //std::cout<<"---- %50  AcceptanceRateRewiring = "<<sumAcceptanceRewirig<<"\n";
                  AcceptanceRateRewiring.push_back(sumAcceptanceRewirig);
                  //AcceptanceRewiringFile<<sumAcceptanceRewirig<<"\n";
                  sumAcceptanceRewirig = 0;
            }
          }//for
        }//end selfish
        FinalY = y;
}
/*------------------------------------------------------------*/
dim1 ODE::runDynamics(int _NumbertOfSteps, dim2 _Cij, dim1 _y, dim1 &MeanYPrime){
    dim1 r1, r2;
    double lastPushMeanYPrime;
    //std::cout<<"-- in runDynamics "<<"_NumbertOfSteps = "<<_NumbertOfSteps<<"\n";
    for(int i=0; i<_NumbertOfSteps; i++){
        //std::cout<<"step = "<<i<<"\n";
        r1 = order_parameter(_y);
        //r2 = order_parameter_k(_y);
        //r1 = order_parameter(_y);
        //r2 = order_parameter_k(_y);
        Order1.push_back(r1[0]);
        Psi1.push_back(r1[1]);
        //Order2.push_back(r2[0]);
        //Psi2.push_back(r2[1]);
        //std::cout<<"Before\n";
        //Print2D(_Cij);
        //Print1D(_y);
        _y = runge_kutta4_integrator(_y, _Cij);
        //std::cout<<"After\n";
        //Print1D(_y);
        lastPushMeanYPrime = Mean(dydt(_y, _Cij), 0);
        //std::cout<<"=== Mean Y Prime = "<<lastPushMeanYPrime<<"\n";
        MeanYPrime.push_back(lastPushMeanYPrime);
        //std::cout<<"=== size of vector = "<<MeanYPrime.size()<<"\n";

    }
    return _y;
}
/*------------------------------------------------------------*/
dim2 ODE::rewiring(int indexFocusNode, dim1 nodesOrder, dim2 _Cij){
    std::random_shuffle(nodesOrder.begin(), nodesOrder.end());
    bool needZeroTobeOne, needOneTobeZero;
    needOneTobeZero = 0;
    needZeroTobeOne = 0;
    int j = 0;
    if(indexFocusNode == nodesOrder[j]) j = 1;
    if(_Cij[indexFocusNode][nodesOrder[j]]==1){
        _Cij[indexFocusNode][nodesOrder[j]] = 0;
        //_Cij[nodesOrder[j]][indexFocusNode] = 0;
        needZeroTobeOne = 1;
        //std::cout<<nodesOrder[j]<<"->";
    }else{
           _Cij[indexFocusNode][nodesOrder[j]] = 1;
           //_Cij[nodesOrder[j]][indexFocusNode] = 1;
           needOneTobeZero = 1;
           //std::cout<<nodesOrder[j]<<"->";
    }

    for(int i=j+1; i<N && (needZeroTobeOne || needOneTobeZero); i++){
      if(nodesOrder[i] != indexFocusNode){
        if(_Cij[indexFocusNode][nodesOrder[i]]==1 && needOneTobeZero){
            _Cij[indexFocusNode][nodesOrder[i]] = 0;
            //_Cij[nodesOrder[i]][indexFocusNode] = 0;
            needOneTobeZero = 0;
            //std::cout<<nodesOrder[i]<<"\n";
        }else if(_Cij[indexFocusNode][nodesOrder[i]]==0 && needZeroTobeOne){
                 _Cij[indexFocusNode][nodesOrder[i]] = 1;
                 //_Cij[nodesOrder[i]][indexFocusNode] = 1;
                 needZeroTobeOne = 0;
                 //std::cout<<nodesOrder[i]<<"\n";
        }
      }
    }
    return _Cij;
}
/*------------------------------------------------------------*/
void ODE::euler_integrator (dim1 &y )
{
    dim1 f(N);
    f = dydt(y, Cij);
    for (int i=0; i<y.size(); i++)
        y[i] += f[i] * dt;
}
/*------------------------------------------------------------*/
dim1 ODE::runge_kutta4_integrator (dim1 y, dim2 CijLocal)
{
    int n = y.size();
    dim1 k1(n), k2(n), k3(n), k4(n);
    dim1 f(n);
    k1 = dydt(y, CijLocal);
    for (int i=0; i<n; i++)
        f[i]= y[i]+ 0.5 * dt * k1[i];
    k2 = dydt (f, CijLocal);
    for (int i=0; i<n; i++)
        f[i]= y[i]+ 0.5 * dt * k2[i];
    k3 = dydt(f, CijLocal);
    for (int i=0; i<n; i++)
        f[i]= y[i]+ dt * k3[i];
    k4 = dydt(f, CijLocal);
    for (int i=0; i<n; i++)
        y[i] += (k1[i] + 2.0*(k2[i] + k3[i]) + k4[i]) * dt/6.0;
    return y;
}
/*------------------------------------------------------------*/
dim1 ODE::dydt(const dim1 &x, dim2 CijLocal)
{
    double sumx = 0.0;
    dim1 f(N);
    #pragma omp parallel for reduction(+:sumx)
    for (int i=0; i<N; i++)
    {
        sumx = 0;
        for(int j=0; j<N; j++)
        {
            if ((i!=j) && (CijLocal[i][j]!=0))
                sumx += CijLocal[i][j] * sin(x[j]-x[i]);
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
    int n = x.size();
    double real_R = 0.;
    double imag_R = 0.;
    for(int i=0; i<n; i++)
    {
        real_R += cos(x[i]);
        imag_R += sin(x[i]);
    }
    real_R /= (double) n;
    imag_R /= (double) n;

    double r = sqrt(real_R * real_R + imag_R * imag_R);
    double psi = atan2(imag_R,real_R);

    dim1 result{r, psi};
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
dim1 ODE::getAcceptanceRewiring(){
    return AcceptanceRateRewiring;
}
/*------------------------------------------------------------*/
dim1 ODE::getMeanYPrime(){
    dim1 MeanYPrime;
    return MeanYPrime;
}
/*------------------------------------------------------------*/
dim1 ODE::getCij(){
  dim1 _newCij;
  int i,j;
  for(i=0; i<N; i++)
    for(j=0; j<N; j++)
      _newCij.push_back(Cij[i][j]);

  return _newCij;
}
/*------------------------------------------------------------*/
dim1 ODE::getMeanRinEachIteration(){
  return MeanRinEachIteration;
}
/*------------------------------------------------------------*/
dim1 ODE::getFinalY(){
    return FinalY;
}
/*------------------------------------------------------------*/
void ODE::Print2D(dim2 _Cij){
  int i,j;
  std::cout<<"=========== MATRIX =============\n";
  for(i=0; i<N; i++)
  {
    for(j=0; j<N; j++)
      std::cout<<Cij[i][j]<<" \t";
    std::cout<<"\n";
  }
  std::cout<<"=========== ====== =============\n";
}
/*------------------------------------------------------------*/
void ODE::Print1D(dim1 _y){
  int i;
  //std::cout<<"=========== ARRAY =============\n";
  for(i=0; i<_y.size(); i++)
      std::cout<<_y[i]<<"   \t";
  std::cout<<"\n";
  std::cout<<"=========== ===== =============\n";
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
