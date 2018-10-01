
#include <iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include<time.h>
#include<algorithm>
#include<vector>
#include<armadillo>

#include "jacobi.h"

using namespace arma;

// Performs Jacobi rotation to find eigenvalues and eigenvectors
	void jacobirot(mat& a, mat& v, int p, int q, int n)
	{
		double s, c; // initialization of sin and cos
		if (a(p, q) != 0.0) {
			double t, tau; // initialization of tan and tau, which are only 
			               // used if the non-diagonal element is different than zero
			tau = (a(q,q ) - a(p, p)) / (2 * a(p, q)); // we compute tau

			if (tau >= 0) { // we compute tan
				t = 1.0 / (tau + sqrt(1.0 + tau * tau));
			}
			else {
				t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
			}

			c = 1 / sqrt(1 + t * t); // computation of cos and sin
			s = c * t;
		}
		else {
			c = 1.0; // computation of cos and sin
			s = 0.0;
		}
		double a_pp, a_qq, a_ip, a_iq, v_ip, v_iq; // initialization of the new matrix elements 
												   // and the eigenvector matrix elements
		a_pp = a(p, p);
		a_qq = a(q, q);

		// we compute the new matrix elements
		a(p, p) = c * c*a_pp - 2.0*c*s*a(p, q) + s * s*a_qq;
		a(q, q) = s * s*a_pp + 2.0*c*s*a(p, q) + c * c*a_qq;
		a(p, q) = 0.0;  // hard-coding non-diagonal elements to zero
		a(q, p) = 0.0;  // hard-coding to zero as well
		for (int i = 0; i < n; i++) {
			if (i != p && i != q) {
				a_ip = a(i, p);
				a_iq = a(i, q);
				a(i, p) = c * a_ip - s * a_iq;
				a(p, i) = a(i, p);
				a(i, q) = c * a_iq + s * a_ip;
				a(q, i) = a(i, q);
			}
			//  And finally the new eigenvectors
			v_ip = v(i, p);
			v_iq = v(i, q);

			v(i, p) = c * v_ip - s * v_iq;
			v(i, q) = c * v_iq + s * v_ip;
		}
		return;
	} // end of function jacobirot


//gets only the first three eigenvectors
mat get_eigenvecs(mat a, mat v, int n){
    vector<double>eigenvals=get_eigenvals(a,n);
    mat vecs(3,n);
    for(int i=0;i<3;i++){
        for(int j=0;j<n;j++){
            if(a(j,j)==eigenvals[i]){
                for(int k=0;k<n;k++){
                      vecs(i,k)=v(k,j);
                }
             }
         }
    }
    return vecs;
}

//gets eigenvalues from the diagonalized matrix and puts them in order in a vector
vector<double> get_eigenvals(mat a,int n){
    vector<double>eigen;
    for(int i=0;i<n;i++){
        eigen.push_back(a(i,i));
    }
    sort (eigen.begin(), eigen.begin()+n);
    return eigen;
}

//initialize matrix and vectors
void initialize(int n, double h, mat& a, vec& r, mat& v,int interact,double wr){
    //initialize rho values
    r(0)=h;
    for (int i=1; i<n ;i++){
        r(i)=r(i-1)+h;
    }
    
    //initialize matrix and vector depending on the case we are studying
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            if(i==j && interact==0){ //buckling beam
                a(i,j)=2/(h*h);
                v(i,j)=1;
            }
            else if (i==j && interact==1){ //harmonic oscillator with one electron
				a(i, j) = 2 / (h*h) + r(i)*r(i);
                v(i,j)=1;
            }
			else if (i == j && interact == 2) { //harmonic oscillator with two electrons
				a(i, j) = 2 / (h*h) + wr * wr*r(i)*r(i) + 1 / r(i);
				v(i, j) = 1;
			}
            else if (i==j+1 or i==j-1){
                a(i,j)=-1/(h*h);
            } 
            else{
                a(i,j)=0;
                v(i,j)=0;
            }
        }
    }
}

//find maximum non-diagonal matrix elements
void offdiag(mat a, int& p, int& q, int n, double &max) {
	max = 0.0;
	for (int i = 0; i < n; ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			double aij = fabs(a(i, j));
			if (aij > max)
			{
				max = aij;  p = i; q = j; //"returns" the value and "place" of the max element
			}
		}
	}
}