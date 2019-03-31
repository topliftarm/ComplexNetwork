/* File : ode_solver.i */
%module ode_solver

%{
#include "ode_solver.h"
%}

%include stl.i
/* instantiate the required template specializations */
namespace std {
    %template(IntVector)     vector<int>;
    %template(DoubleVector)  vector<double>;
    %template(DoubleVector2) vector<vector<double> >;
}
/* Let's just grab the original header file here */
%include "ode_solver.h"
