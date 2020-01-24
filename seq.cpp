#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;

int main(int argc, char const *argv[])
{
	ifstream is("data.txt");
	int n;
	is >> n;

	vector<int> p(n);
	cout<<"Error 0"<<endl;
	// float u[n][n];
	float *u = (float *)malloc(n*n*sizeof(float));
	cout<<"Error 1"<<endl;
	// float l[n][n];
	float *l = (float *)malloc(n*n*sizeof(float));
	cout<<"Error 2"<<endl;
	// float a[n][n];
	float *a = (float *)malloc(n*n*sizeof(float));
	cout<<"Error 3"<<endl;
	// float a2[n][n];
	float *a2 = (float *)malloc(n*n*sizeof(float));
	cout<<"Error 4"<<endl;

    float temp;
	for (int i=0; i<n; i++){
		p[i] = i;
		for (int j=0; j<n; j++){
			is >> temp;
			a[i*n+j] = temp;
			a2[i*n+j] = temp;
			u[i*n+j] = 0;
			// zero[i][j] = 0;
			if (i==j)	l[i*n+j] = 1;
			else	l[i*n+j] = 0;
		}
	}


	for (int k=0; k<n ; k++){
		float max = 0.0;
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
		float temp = p[k];
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
		for (int j=0; j<n; j++){
			a[i*n+j]=0;
		}
	}
	for (int i=0; i<n; i++){
		a[i*n+(p[i])] = 1;
	}

	cout<<"ERROR CHECK"<<endl;
	// float PA[n][n];
	float *PA = (float *)malloc(n*n*sizeof(float));
	cout<<"ERROR PA"<<endl;
	for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			float sum = 0.0;
			for (int k=0; k<n; k++)
			{
				sum = sum + (a[i*n+k] * a2[k*n+j]);
			}
			// PA[i][j] = sum;
			PA[i*n + j] = sum;
		}
	}

	// float LU[n][n]; ==== a[i][j]
	for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			float sum = 0.0;
			for (int k=0; k<n; k++)
			{
				sum = sum + (l[i*n+k] * u[k*n+j]);
			}
			a[i*n+j] = sum;
		}
	}


	// float diff[n][n]; == a2[i][j]
	for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			a2[i*n+j] = (PA[i*n + j] - a[i*n+j]);
		}
	}

	float ans;
	for (int j=0; j<n; j++){
		float col = 0;
		for (int i=0; i<n; i++){
			col = col + (a2[i*n+j]*a2[i*n+j]);
		}
		ans = ans + sqrt(col);
	}

	cout<<" ::ANS = "<<ans<<endl; 
	return 0;
}