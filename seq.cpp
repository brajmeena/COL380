#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;
void LU_Decomp(int);

int main(int argc, char const *argv[])
{
	int n;
	cout<<"Enter n:";
	cin>>n; 
	LU_Decomp(n);
	return 0;
}


void LU_Decomp(int n){	
	vector<int> p(n);
	double *u = (double *)malloc(n*n*sizeof(double));
	double *l = (double *)malloc(n*n*sizeof(double));
	double *a = (double *)malloc(n*n*sizeof(double));
	double *a2 = (double *)malloc(n*n*sizeof(double));
	double *zero = (double *)malloc(n*n*sizeof(double));  //to make p[n][n]
	
    double temp, rand=1000;
	for (int i=0; i<n; i++){
		p[i] = i;
		for (int j=0; j<n; j++){
			temp = (drand48() * rand);
			a[i*n+j] = temp;
			a2[i*n+j] = temp;
			u[i*n+j] = 0;
			zero[i*n+j] = 0;
			if (i==j)	l[i*n+j] = 1;
			else	l[i*n+j] = 0;
		}
	}

	for (int k=0; k<n ; k++){
		double max = 0.0;
		int k2 = k;
		for (int i=k; i<n; i++){
			if (max<abs(a[i*n+k])){
				max = abs(a[i*n+k]);
				k2 = i;
			}
		}
		
		if (max == 0){
			cout<<"ERROR - Singular Matrix"<<endl;
		}
		
		// #### SWAP P[k] and P[K']
		double temp = p[k];
		p[k] = p[k2];
		p[k2] = temp; 
		
		// # swap a(k,:) and a(k',:)
		for (int i=0; i<n; i++){
			temp = a[k*n+i];
			a[k*n+i] = a[k2*n+i];
			a[k2*n+i] = temp;
		}
		
		// swap l(k,1:k-1) and l(k',1:k-1)
		for (int i=0; i<k; i++){
			temp = l[k*n+i];
			l[k*n+i] = l[k2*n+i];
			l[k2*n+i] = temp;
		}

		u[k*n+k] = a[k*n+k];
		
		for (int i=k+1; i<n; i++){
			l[i*n+k] = a[i*n+k]/u[k*n+k];
			u[k*n+i] = a[k*n+i];
		}
		
		for (int i=k+1; i<n; i++){
			for (int j=k+1; j<n; j++){
				a[i*n+j] = a[i*n+j]-(l[i*n+k]*u[k*n+j]);
			}
		}
	}

	for (int i=0; i<n; i++){
		zero[i*n+(p[i])] = 1;   // zero = P[n][n]
	}

	// double PA[n][n];
	double *PA = (double *)malloc(n*n*sizeof(double));
	for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			double sum1=0.0, sum2=0.0;
			for (int k=0; k<n; k++)
			{
				sum1 = sum1 + (zero[i*n+k] * a2[k*n+j]);
				sum2 = sum2 + (l[i*n+k] * u[k*n+j]);
			}
			PA[i*n + j] = sum1;
			a[i*n+j] = sum2;    // a = LU
		}
	}

	// double diff[n][n]; == a2[i][j]
	for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			a2[i*n+j] = (PA[i*n + j] - a[i*n+j]);   //a2=PA-LU
		}
	}

	double ans=0.0;
	for (int j=0; j<n; j++){
		double col = 0;
		for (int i=0; i<n; i++){
			col = col + (a2[i*n+j]*a2[i*n+j]);
		}
		ans = ans + sqrt(col);
	}

	cout<<" ::ANS = "<<ans<<endl;

}

