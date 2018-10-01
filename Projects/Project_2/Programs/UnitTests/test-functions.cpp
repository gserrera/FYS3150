#include "catch.hpp"
#include "jacobi.h"

#define PI 3.14159265
TEST_CASE("Testing max a(i,j)"){
	int n = 3, interact = 0;
    double pmin=0, pmax=10,h = (pmax-pmin)/n;
	n = n - 1;
    mat a = zeros<mat>(n,n);
    mat v = zeros<mat>(n,n);
    vec r(n);
    //initialize matrices and vector
    initialize(n,h,a,r,v,interact,0);
    int p;
    int q;
	double max;
    //find maximum matrix element
	offdiag(a, p, q, n, max);
    REQUIRE(p==0);
    REQUIRE(q==1);
    REQUIRE(max==Approx(0.09));
}

TEST_CASE("Testing eigenvalues"){
    int n=4,interact=0;
    double tol=0.001,wr=0.01, pmin=0, pmax=10,h = (pmax-pmin)/(double(n));
    n=n-1;
    mat a = zeros<mat>(n,n);
    mat v = zeros<mat>(n,n);
    vec r(n);
    //initialize matrices and vector
	initialize(n, h, a, r, v, interact, wr);
	int iterations = 0;
	int p, q;
	double max = 1.0;
	//do jacobi algorithm until convergence
	while (max > tol) {
		offdiag(a, p, q, n, max);
		jacobirot(a, v, p, q, n);
		iterations++;
	}
    //get eigenvalue vector
    vector<double>eigen=get_eigenvals(a,n);
	vector<double>exact(n);
	for (int i = 1; i <= n; i++) {
		exact[i - 1] = 2 / (h*h) - 2 / (h*h)*cos(i*PI / (n + 1));
	}
    REQUIRE(eigen[0]==Approx(exact[0]));
    REQUIRE(eigen[1]==Approx(exact[1]));
    REQUIRE(eigen[2]==Approx(exact[2]));
}
TEST_CASE("Testing eigenvector orthogonality"){
	int n = 4;
    double tol=0.001, pmin=0, pmax=10,h = (pmax-pmin)/(double(n));
    n=n-1;
    mat a = zeros<mat>(n,n);
    mat v = zeros<mat>(n,n);
    vec r(n);
    //initialize matrices and vector
	initialize(n, h, a, r, v, 0, 0);
	int iterations = 0;
	int p, q;
	double max = 1.0;
	//do jacobi algorithm until convergence
	while (max > tol) {
		offdiag(a, p, q, n, max);
		jacobirot(a, v, p, q, n);
		iterations++;
	}
    mat eigenvec=get_eigenvecs(a,v,n);

    //test eigen vector orthogonality
    //dot1=v0*v1=0
    double dot1=eigenvec(0,0)*eigenvec(1,0)+eigenvec(0,1)*eigenvec(1,1)
        +eigenvec(0,2)*eigenvec(1,2);
    //dot2=v0*v0=1
    double dot2=eigenvec(0,0)*eigenvec(0,0)+eigenvec(0,1)*eigenvec(0,1)
        +eigenvec(0,2)*eigenvec(0,2);
    REQUIRE(dot1==Approx(0.000).epsilon(0.01));
    REQUIRE(dot2==Approx(1.000));
}
