#ifndef JACOBI_H
#define	JACOBI_H

#include <iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include <armadillo>
#include<time.h>
#include<vector>

using namespace std;
using namespace arma;

void initialize(int,double,mat&,vec&,mat&,int,double);
void jacobirot(mat&,mat&, int, int, int);
void offdiag(mat,int&,int&,int,double&);
vector<double> get_eigenvals(mat,int);
mat get_eigenvecs(mat,mat,int);
#endif /* JACOBI_H */

