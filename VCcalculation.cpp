#include <iostream>
#include <ctime>
#include <algorithm> 
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include <map>
#include <string>
#define t 2.8
#define a 1.42
#define e 1
#define Pi 3.1415
using namespace std;
using namespace Eigen;
std::complex<double> I = 1i; 


using Eigen::MatrixXd;

void saveToCSV(double* matrix, int dim, string filename);
MatrixXcd transformation(int n ,int wD1 ,int wD2 ,double D ,double d );
MatrixXcd hamiltonian(int n ,double k ,double Ey ,int wD1 ,int wD2 ,double D ,double d );
complex<double> vp(double k, int n, MatrixXcd SH, int b);

string IntToString (int h)
{
    ostringstream temp;
    temp<<h;
    return temp.str();
}

int main()
{
    //define your own matrix
    int n = 25;
    double k=0.5;
    double Ey=0.005;
    int b;
    int wD1=1, wD2=1;
    double D=0.1;
    double d=0.05;
    double vc=0;
    double* eigenvalue_arr1;
    eigenvalue_arr1 = (double *)operator new(2*n * sizeof(double));
    double* eigenvalue_arr2;
    eigenvalue_arr2 = (double *)operator new(2*n * sizeof(double));
    double* eigenvalue_arr3;
    eigenvalue_arr3 = (double *)operator new(2*n * sizeof(double));
    complex<double>* eigenvector_arr;
    eigenvector_arr = (complex<double> *)operator new(2*n * sizeof(complex<double>));
    MatrixXcd H(2*n, 2*n), H1(2*n, 2*n), T(2*n, 2*n), TH(2*n, 2*n), SS(2*n, 2*n), Data(200,200);
	clock_t start, stop;
	/*
	for(int i=0;i<=2*n-1;i++)
	{
		for(int m=0;m<=2*n-1;m++)
		{
			SS(i,m)=(m+1)*0.1+(i+1)*0.1*I;
		}
	}
	*/
	start = clock();

	T=transformation(n,wD1,wD2,D,d);
	for(int ee=0;ee<=0;ee++)
	{
		for(int kk=0;kk<=9;kk++)
		{
			H=hamiltonian(n,-Pi/(3*a) + Pi/(600*a)*(kk ),ee*Ey,wD1,wD2,D,d);
			ComplexEigenSolver<MatrixXcd> h(H);
			
		    for (int i = 0; i <= 2*n-1; i++) 
			{
		        eigenvalue_arr1[i] = real(h.eigenvalues()[i]);
		    }
		    sort(eigenvalue_arr1, eigenvalue_arr1 + 2*n);

			H1=hamiltonian(n,-(Pi/(3*a)) + Pi/(600*a)*(kk )+0.00001,ee*Ey,wD1,wD2,D,d);
			ComplexEigenSolver<MatrixXcd> h1(H1);
			
		    for (int i = 0; i <= 2*n-1; i++) 
			{
		        eigenvalue_arr2[i] = real(h1.eigenvalues()[i]);
		    }
		    sort(eigenvalue_arr2, eigenvalue_arr2 + 2*n);

			TH=T*H*(T.transpose());
			SelfAdjointEigenSolver<MatrixXcd> es(TH);
			SS=es.eigenvectors().transpose();
			vc=0;
			for(int b=0;b<=n-1;b++)
			{
				vc=vc+real(vp(-(Pi/(3*a)) + Pi/(600*a)*(kk),n,SS,b))*(eigenvalue_arr2[b]-eigenvalue_arr1[b])/0.00001;
			}
			Data(kk,ee)	= vc;
				
		}
	}
	stop = clock();
	cout << double(stop - start) / CLOCKS_PER_SEC <<endl;
	cout<<Data<<endl;
			/*
		    for (int i = 0; i <= 2*n-1; i++) 
			{
		        eigenvalue_arr3[i] = real(es.eigenvalues()[i]);
		    }
		    cout<<es.eigenvectors()<<endl<<endl;
		    map<double, VectorXcd> index_map;
		    for (int i = 0; i <= 2*n-1; i++) 
			{
		        // use es.eigenvectors().col(0) to get the eigenvector corresponding to the i'th eigenvalue
		        index_map[eigenvalue_arr3[i]] = es.eigenvectors().col(i);
		    }
			sort(eigenvalue_arr3, eigenvalue_arr3 + 2*n);
			
			for(int i=0;i<=2*n-1;i++)
			{
				for(int j=0;j<=2*n-1;j++)
				{
					SS(i,j)=index_map[eigenvalue_arr3[i]][j];
				}
			}
			
			cout<<eigenvalue_arr3[0]<<endl<<endl;
			cout<<SS.row(0)<<endl;
			cout<<(((SS.row(0)).transpose()).adjoint()*((SS.row(0)).transpose())).value()<<endl;
			cout<<eigenvalue_arr3[0]*SS.row(0)<<endl;
			cout<<(TH*(SS.row(0)).transpose()).transpose()<<endl;
			vc=0;
			for(b=0;b<=n-1;b++)
			{

			}
			*/
  //  cout << "H:" << endl << H << endl << endl;
  //  cout << "T:" << endl << T << endl << endl;
  //  TH=T*H*(T.transpose());
  //  cout << "SS:" << endl << SS << endl << endl;
  //  h=vp(k,n,SS,5);
  //  cout << "vpd:" << endl << h << endl << endl;
    
  //  ComplexEigenSolver<MatrixXcd> es(H);
 //   cout << "The eigenvalues of H are:" << endl << es.eigenvalues() << endl << endl;
    // use es.eigenvalues()[i] to get the i'th eigenvalue

 //   cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
    
    
    
    /*
    
    // Now, sort eigenvalue and map eigenvalues with eigenvectors.
    double* eigenvalue_arr;
    eigenvalue_arr = (double *)operator new(2*n * sizeof(double));
    for (int i = 0; i < 2*n-1; i++) {
        eigenvalue_arr[i] = real(es.eigenvalues()[i]);
    }
    
    cout << "Eigenvalue array:" << endl;
    for (int i = 0; i < 2*n; i++) {
        cout << eigenvalue_arr[i] << " ";
    }
    cout << endl << endl;
    
    map<double, VectorXcd> index_map;
    for (int i = 0; i < 2*n-1; i++) {
        // use es.eigenvectors().col(0) to get the eigenvector corresponding to the i'th eigenvalue
        index_map[eigenvalue_arr[i]] = es.eigenvectors().col(i);
    }
    
    sort(eigenvalue_arr, eigenvalue_arr + 2*n);
    cout << "Sorted eigenvalue array:" << endl;
    for (int i = 0; i < 2*n; i++) {
        cout << eigenvalue_arr[i] << " ";
    }
    cout << endl << endl;
    
    // Now eigenvalue_arr is h array of sorted eigenvalues.
    // To find the corresponding eigenvector of eigenvalue_arr[i],
    // use index_map[eigenvalue_arr[i]].
    // And it will return h VectorXcd datatype.
    
    // End of sorting and mapping eigenvalues and eigenvectors.
    
    
    
    
    // ########################################################################
    // To turn VectorXcd into an array, use the following example
    
    VectorXcd v = index_map[eigenvalue_arr[0]];
    
    cout << "eigenvalue: " << endl << eigenvalue_arr[0] << endl;
    cout << "corresponding eigenvector: " << endl;
    
    complex<double>* eigenvector_arr;
    eigenvector_arr = (complex<double> *)operator new(2*n * sizeof(complex<double>));

    for (int i = 0; i < 2*n; i++) {
        // use v[i] to access the i'th element of v, the datatype will be complex<double>
        eigenvector_arr[i] = v[i];
    }
    for (int i = 0; i < 2*n; i++) {
        cout << eigenvector_arr[i] << " ";
    }
    cout << endl;
    // ########################################################################

    
    
    
    
//    cout << "Consider the first eigenvalue, lambda = " << lambda << endl;
//    cout << "If v is the corresponding eigenvector, then lambda * v = " << endl << lambda * v << endl << endl;
//    cout << "... and H * v = " << endl << H.cast<complex<double> >() * v << endl << endl;

    
    // Example for saving double type matrix to csv file.
    double* matrix;
    matrix = (double *)operator new(2*n * 2*n * sizeof(double));
    for(int i = 0; i < 2*n; i++){
        for(int j = 0; j < 2*n; j++){
            matrix[i * 2*n + j] = 1 + i*j;
        }
    }
    
    saveToCSV(matrix, 2*n, "e" + IntToString(43) + "output.csv");
    */
    system ("pause");
}

void saveToCSV(double* matrix, int dim, string filename) {
    ofstream outfile (filename.c_str());
    
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim-1; j++){
            outfile << matrix[i * dim + j] << ", ";
        }
        outfile << matrix[i * dim + dim - 1] << endl;
    }
    outfile.close();
    cout << "Success saving file\n";
}

MatrixXcd transformation(int n ,int wD1 ,int wD2 ,double D ,double d )
{
	MatrixXcd T(2*n, 2*n);
	
	for(int i=0;i<=2*n-1;i++)
	{
		for(int m=0;m<=2*n-1;m++)
		{
			T(i,m)=0;
		}
	}
	
	
	int n1=(n+1)/2; //half of left side atom #
	int n2=(n-1)/2; //half of right side atom #
	double N;
	double d_bond_x=1./sqrt(2*d*d+2*t*t-2*d*sqrt(d*d+t*t))*(d-sqrt(d*d+t*t));
	double d_bond_y=-1./sqrt(2*d*d+2*t*t-2*d*sqrt(d*d+t*t))*t;
	double d_at_bond_x=1./sqrt(2*d*d+2*t*t+2*d*sqrt(d*d+t*t))*(d+sqrt(d*d+t*t));
	double d_at_bond_y=-1./sqrt(2*d*d+2*t*t+2*d*sqrt(d*d+t*t))*t;
	
	double D_bond_x=1./sqrt(2*D*D+2*t*t-2*D*sqrt(D*D+t*t))*(D-sqrt(D*D+t*t));
	double D_bond_y=-1./sqrt(2*D*D+2*t*t-2*D*sqrt(D*D+t*t))*t;
	double D_at_bond_x=1./sqrt(2*D*D+2*t*t+2*D*sqrt(D*D+t*t))*(D+sqrt(D*D+t*t));
	double D_at_bond_y=-1./sqrt(2*D*D+2*t*t+2*D*sqrt(D*D+t*t))*t;

	for( int i=0;i<=2*n-1;i++)
	{
		for( int j=0;j<=2*n-1;j++)
		{
			T(i,j)=0.;
		}
	}
	
	if( n % 2 != 0 )
	{
		if( (((n+1)/2)%2) != 0)
		{

			//5,9,13...
			//bonding
			//left
			//cos
			for(int i=0;i<=(n1+1)/2-1;i++)
			{
				N=0;
				for(int m=0;m<=n2;m++)
				{
					N=N+cos(4*(i+0.5)*Pi/(n+3)*(n2/2-m))*cos(4*(i+0.5)*Pi/(n+3)*(n2/2-m));
				}
				for(int j=0;j<=n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4*(i+0.5)*Pi/(n+3)*(n2/2-j))*d_bond_x;
				}
				for(int j=n1;j<=n1+n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4*(i+0.5)*Pi/(n+3)*(n2/2-(j-n1)))*d_bond_y;
				}
			}
			//5,9,13...
			//bonding
			//left
			//sin
			for(int i=(n1+1)/2;i<=(n1-1);i++)
			{
				N=0;
				for(int m=0;m<=n2;m++)
				{
					N=N+sin(4*(i+1-(n1+1)/2)*Pi/(n+3)*(n2/2-m))*sin(4*(i+1-(n1+1)/2)*Pi/(n+3)*(n2/2-m));
				}
				for(int j=0;j<=n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4*(i+1-(n1+1)/2)*Pi/(n+3)*(n2/2-j))*d_bond_x;
				}
				for(int j=n1;j<=n1+n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4*(i+1-(n1+1)/2)*Pi/(n+3)*(n2/2-(j-n1)))*d_bond_y;
				}
			}
			//5,9,13...
			//bonding
			//Right
			//cos
			for(int i=n1;i<=n1+n2/2-1;i++)
			{
				N=0;
				for(int m=0;m<=n2-1;m++)
				{
					N=N+cos(4*(i-n1+0.5)*Pi/(n+1)*((n2-1.)/2.-m))*cos(4*(i-n1+0.5)*Pi/(n+1.)*((n2-1.)/2.-m));
				}
				for(int j=2*n1;j<=2*n1+n2-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4*(i-n1+0.5)*Pi/(n+1.)*((n2-1.)/2.-j+2*n1))*d_bond_x;
				}
				for(int j=2*n1+n2;j<=2*n-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4*(i-n1+0.5)*Pi/(n+1.)*((n2-1.)/2.-j+2*n1+n2))*d_bond_y;
				}
			}
			//5,9,13...
			//bonding
			//Right
			//sin
			for(int i=n1+n2/2;i<=(n-1);i++)
			{
				N=0;
				for(int m=0;m<=n2-1;m++)
				{
					N=N+sin(4*(i+1-(n1+n2/2.))*Pi/(n+1)*((n2-1.)/2.-m))*sin(4*(i+1-(n1+n2/2.))*Pi/(n+1)*((n2-1.)/2.-m));
				}
				
				for(int j=2*n1;j<=2*n1+n2-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4.*(i+1.-(n1+n2/2.))*Pi/(n+1.)*((n2-1.)/2.-j+2*n1))*d_bond_x;
				}
				for(int j=2*n1+n2;j<=2*n-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4.*(i+1.-(n1+n2/2.))*Pi/(n+1.)*((n2-1.)/2.-j+2*n1+n2))*d_bond_y;
				}
			}
			//5,9,13...
			//anti-bonding
			//left
			//cos
			for(int i=0+n;i<=(n1+1)/2+n-1;i++)
			{
				N=0;
				for(int m=0;m<=n2;m++)
				{
					N=N+cos(4*(i-n+0.5)*Pi/(n+3)*(n2/2-m))*cos(4*(i-n+0.5)*Pi/(n+3)*(n2/2-m));
				}
				for(int j=0;j<=n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4*(i-n+0.5)*Pi/(n+3)*(n2/2-j))*d_at_bond_x;
				}
				for(int j=n1;j<=n1+n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4*(i-n+0.5)*Pi/(n+3)*(n2/2-(j-n1)))*d_at_bond_y;
				}
			}
			//5,9,13...
			//anti-bonding
			//left
			//sin
			for(int i=(n1+1)/2+n;i<=(n1+n-1);i++)
			{
				N=0;
				for(int m=0;m<=n2;m++)
				{
					N=N+sin(4*(i+1-(n1+1)/2-n)*Pi/(n+3)*(n2/2-m))*sin(4*(i+1-(n1+1)/2-n)*Pi/(n+3)*(n2/2-m));
				}
				for(int j=0;j<=n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4*(i+1-(n1+1)/2-n)*Pi/(n+3)*(n2/2-j))*d_at_bond_x;
				}
				for(int j=n1;j<=n1+n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4*(i+1-(n1+1)/2-n)*Pi/(n+3)*(n2/2-(j-n1)))*d_at_bond_y;
				}
			}
			//5,9,13...
			//anti-bonding
			//Right
			//cos
			for(int i=n1+n;i<=n1+n2/2+n-1;i++)
			{
				N=0;
				for(int m=0;m<=n2-1;m++)
				{
					N=N+cos(4*(i-n1-n+0.5)*Pi/(n+1)*((n2-1.)/2.-m))*cos(4*(i-n1-n+0.5)*Pi/(n+1.)*((n2-1.)/2.-m));
				}
				for(int j=2*n1;j<=2*n1+n2-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4*(i-n1-n+0.5)*Pi/(n+1.)*((n2-1.)/2.-j+2*n1))*d_at_bond_x;
				}
				for(int j=2*n1+n2;j<=2*n-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4*(i-n1-n+0.5)*Pi/(n+1.)*((n2-1.)/2.-j+2*n1+n2))*d_at_bond_y;
				}
			}
			//5,9,13...
			//anti-bonding
			//Right
			//sin
			for(int i=n1+n2/2+n;i<=(n+n-1);i++)
			{
				N=0;
				for(int m=0;m<=n2-1;m++)
				{
					N=N+sin(4*(i+1-(n1+n2/2.)-n)*Pi/(n+1)*((n2-1.)/2.-m))*sin(4*(i+1-(n1+n2/2.)-n)*Pi/(n+1)*((n2-1.)/2.-m));
				}
				
				for(int j=2*n1;j<=2*n1+n2-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4.*(i+1.-(n1+n2/2.)-n)*Pi/(n+1.)*((n2-1.)/2.-j+2*n1))*d_at_bond_x;
				}
				for(int j=2*n1+n2;j<=2*n-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4.*(i+1.-(n1+n2/2.)-n)*Pi/(n+1.)*((n2-1.)/2.-j+2*n1+n2))*d_at_bond_y;
				}
			}
			//5,9,13...
			//barrier
			//bonding
			for(int i=0;i<=n-1;i++)
			{
				//Ubarrier1
				for(int j=0;j<=wD1-1;j++)
				{
					T(i,j)=T(i,j)/d_bond_x*D_bond_x;
				}
				//Dbarrier1
				for(int j=n1-wD1;j<=n1-1;j++)
				{
					T(i,j)=T(i,j)/d_bond_x*D_bond_x;
				}	
				//Ubarrier2
				for(int j=0+n1;j<=wD1-1+n1;j++)
				{
					T(i,j)=T(i,j)/d_bond_y*D_bond_y;
				}
				//Dbarrier2
				for(int j=n1-wD1+n1;j<=n1+n1-1;j++)
				{
					T(i,j)=T(i,j)/d_bond_y*D_bond_y;
				}
				//Ubarrier3
				for(int j=0+n1+n1;j<=wD2-1+n1+n1;j++)
				{
					T(i,j)=T(i,j)/d_bond_x*D_bond_x;
				}
				//Dbarrier3
				for(int j=n1-wD2+n1+n2;j<=n1+n1+n2-1;j++)
				{
					T(i,j)=T(i,j)/d_bond_x*D_bond_x;
				}
				//Ubarrier4
				for(int j=0+n1+n1+n2;j<=wD2-1+n1+n1+n2;j++)
				{
					T(i,j)=T(i,j)/d_bond_y*D_bond_y;
				}
				//Dbarrier4
				for(int j=n1-wD2+n1+n2+n2;j<=n1+n1+n2+n2-1;j++)
				{
					T(i,j)=T(i,j)/d_bond_y*D_bond_y;
				}
		    }
		    //5,9,13...
			//barrier
			//antibonding
			for(int i=n;i<=2*n-1;i++)
			{
				//Ubarrier1
				for(int j=0;j<=wD1-1;j++)
				{
					T(i,j)=T(i,j)/d_at_bond_x*D_at_bond_x;
				}
				//Dbarrier1
				for(int j=n1-wD1;j<=n1-1;j++)
				{
					T(i,j)=T(i,j)/d_at_bond_x*D_at_bond_x;
				}	
				//Ubarrier2
				for(int j=0+n1;j<=wD1-1+n1;j++)
				{
					T(i,j)=T(i,j)/d_at_bond_y*D_at_bond_y;
				}
				//Dbarrier2
				for(int j=n1-wD1+n1;j<=n1+n1-1;j++)
				{
					T(i,j)=T(i,j)/d_at_bond_y*D_at_bond_y;
				}
				//Ubarrier3
				for(int j=0+n1+n1;j<=wD2-1+n1+n1;j++)
				{
					T(i,j)=T(i,j)/d_at_bond_x*D_at_bond_x;
				}
				//Dbarrier3
				for(int j=n1-wD2+n1+n2;j<=n1+n1+n2-1;j++)
				{
					T(i,j)=T(i,j)/d_at_bond_x*D_at_bond_x;
				}
				//Ubarrier4
				for(int j=0+n1+n1+n2;j<=wD2-1+n1+n1+n2;j++)
				{
					T(i,j)=T(i,j)/d_at_bond_y*D_at_bond_y;
				}
				//Dbarrier4
				for(int j=n1-wD2+n1+n2+n2;j<=n1+n1+n2+n2-1;j++)
				{

					T(i,j)=T(i,j)/d_at_bond_y*D_at_bond_y;

				}	
			
			}
		}
		else
		{	
			//7,11,15...
			//bonding
			//left
			//sin			
			for(int i=0;i<=n1/2-1;i++)
			{
				N=0;
				for(int m=0;m<=n2;m++)
				{
					N=N+cos(4*(i+0.5)*Pi/(n+3.)*(n2/2.-m))*cos(4.*(i+0.5)*Pi/(n+3)*(n2/2.-m));
					
				}
				
				for(int j=0;j<=n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4*(i+0.5)*Pi/(n+3.)*(n2/2.-j))*d_bond_x;
				}
				for(int j=n1;j<=n1+n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4*(i+0.5)*Pi/(n+3.)*(n2/2.-(j-n1)))*d_bond_y;
				}
			}
			//7,11,15...
			//bonding
			//left
			//sin
			for(int i=n1/2;i<=(n1-1);i++)
			{
				N=0;
				for(int m=0;m<=n2;m++)
				{
					N=N+sin(4*(i+1-(n1)/2.)*Pi/(n+3.)*(n2/2.-m))*sin(4.*(i+1-(n1)/2.)*Pi/(n+3.)*(n2/2.-m));
				}
				for(int j=0;j<=n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4.*(i+1-(n1)/2.)*Pi/(n+3.)*(n2/2.-j))*d_bond_x;
				}
				for(int j=n1;j<=n1+n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4.*(i+1-(n1)/2.)*Pi/(n+3)*(n2/2.-(j-n1)))*d_bond_y;
				}
			}
			//7,11,15...
			//bonding
			//Right
			//cos
			for(int i=n1;i<=n1+n1/2-1;i++)
			{
				N=0;
				for(int m=0;m<=n2-1;m++)
				{
					N=N+cos(4*(i-n1+0.5)*Pi/(n+1)*((n2-1.)/2.-m))*cos(4*(i-n1+0.5)*Pi/(n+1.)*((n2-1.)/2.-m));
				}
				for(int j=2*n1;j<=2*n1+n2-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4.*(i-n1+0.5)*Pi/(n+1.)*((n2-1.)/2.-j+2.*n1))*d_bond_x;
				}
				for(int j=2*n1+n2;j<=2*n-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4.*(i-n1+0.5)*Pi/(n+1.)*((n2-1.)/2.-j+2*n1+n2))*d_bond_y;
				}
			}
			//7,11,15...
			//bonding
			//Right
			//sin
			for(int i=n1+n1/2;i<=(n-1);i++)
			{
				N=0;
				for(int m=0;m<=n2-1;m++)
				{
					N=N+sin(4*(i+1-(n1+n1/2.))*Pi/(n+1)*((n2-1.)/2.-m))*sin(4*(i+1-(n1+n1/2.))*Pi/(n+1)*((n2-1.)/2.-m));
				}
				
				for(int j=2*n1;j<=2*n1+n2-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4.*(i+1.-(n1+n1/2.))*Pi/(n+1.)*((n2-1.)/2.-j+2*n1))*d_bond_x;
				}
				for(int j=2*n1+n2;j<=2*n-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4.*(i+1.-(n1+n1/2.))*Pi/(n+1.)*((n2-1.)/2.-j+2*n1+n2))*d_bond_y;
				}
			}
			//7,11,15...
			//anti-bonding
			//left
			//cos
			for(int i=0+n;i<=(n1)/2+n-1;i++)
			{
				N=0;
				for(int m=0;m<=n2;m++)
				{
					N=N+cos(4*(i-n+0.5)*Pi/(n+3)*(n2/2.-m))*cos(4*(i-n+0.5)*Pi/(n+3.)*(n2/2.-m));
				}
				for(int j=0;j<=n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4*(i-n+0.5)*Pi/(n+3.)*(n2/2.-j))*d_at_bond_x;
				}
				for(int j=n1;j<=n1+n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4*(i-n+0.5)*Pi/(n+3.)*(n2/2.-(j-n1)))*d_at_bond_y;
				}
			}
			//7,11,15...
			//anti-bonding
			//left
			//sin
			for(int i=(n1)/2+n;i<=(n1+n-1);i++)
			{
				N=0;
				for(int m=0;m<=n2;m++)
				{
					N=N+sin(4*(i+1-(n1)/2.-n)*Pi/(n+3.)*(n2/2.-m))*sin(4*(i+1-(n1)/2.-n)*Pi/(n+3.)*(n2/2.-m));
				}
				for(int j=0;j<=n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4.*(i+1-(n1)/2.-n)*Pi/(n+3.)*(n2/2.-j))*d_at_bond_x;
				}
				for(int j=n1;j<=n1+n1-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4.*(i+1-(n1)/2.-n)*Pi/(n+3.)*(n2/2.-(j-n1)))*d_at_bond_y;
				}
			}
			//7,11,15...
			//anti-bonding
			//Right
			//cos
			for(int i=n1+n;i<=n1+n1/2+n-1;i++)
			{
				N=0;
				for(int m=0;m<=n2-1;m++)
				{
					N=N+cos(4.*(i-n1-n+0.5)*Pi/(n+1)*((n2-1.)/2.-m))*cos(4*(i-n1-n+0.5)*Pi/(n+1.)*((n2-1.)/2.-m));
				}
				for(int j=2*n1;j<=2*n1+n2-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4.*(i-n1-n+0.5)*Pi/(n+1.)*((n2-1.)/2.-j+2*n1))*d_at_bond_x;
				}
				for(int j=2*n1+n2;j<=2*n-1;j++)
				{
					T(i,j)=1/(sqrt(N))*cos(4.*(i-n1-n+0.5)*Pi/(n+1.)*((n2-1.)/2.-j+2*n1+n2))*d_at_bond_y;
				}
			}
			//7,11,15...
			//anti-bonding
			//Right
			//sin
			for(int i=n1+n1/2+n;i<=(n+n-1);i++)
			{
				N=0;
				for(int m=0;m<=n2-1;m++)
				{
					N=N+sin(4*(i+1-(n1+n2/2.)-n)*Pi/(n+1)*((n2-1.)/2.-m))*sin(4*(i+1-(n1+n2/2.)-n)*Pi/(n+1)*((n2-1.)/2.-m));
				}
				
				for(int j=2*n1;j<=2*n1+n2-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4.*(i+1.-(n1+n2/2.)-n)*Pi/(n+1.)*((n2-1.)/2.-j+2*n1))*d_at_bond_x;
				}
				for(int j=2*n1+n2;j<=2*n-1;j++)
				{
					T(i,j)=1/(sqrt(N))*sin(4.*(i+1.-(n1+n2/2.)-n)*Pi/(n+1.)*((n2-1.)/2.-j+2*n1+n2))*d_at_bond_y;
				}
			}
			//7,11,15...
			//barrier
			//bonding
			for(int i=0;i<=n-1;i++)
			{
				//Ubarrier1
				for(int j=0;j<=wD1-1;j++)
				{
					T(i,j)=T(i,j)/d_bond_x*D_bond_x;
				}
				//Dbarrier1
				for(int j=n1-wD1;j<=n1-1;j++)
				{
					T(i,j)=T(i,j)/d_bond_x*D_bond_x;
				}	
				//Ubarrier2
				for(int j=0+n1;j<=wD1-1+n1;j++)
				{
					T(i,j)=T(i,j)/d_bond_y*D_bond_y;
				}
				//Dbarrier2
				for(int j=n1-wD1+n1;j<=n1+n1-1;j++)
				{
					T(i,j)=T(i,j)/d_bond_y*D_bond_y;
				}
				//Ubarrier3
				for(int j=0+n1+n1;j<=wD2-1+n1+n1;j++)
				{
					T(i,j)=T(i,j)/d_bond_x*D_bond_x;
				}
				//Dbarrier3
				for(int j=n1-wD2+n1+n2;j<=n1+n1+n2-1;j++)
				{
					T(i,j)=T(i,j)/d_bond_x*D_bond_x;
				}
				//Ubarrier4
				for(int j=0+n1+n1+n2;j<=wD2-1+n1+n1+n2;j++)
				{
					T(i,j)=T(i,j)/d_bond_y*D_bond_y;
				}
				//Dbarrier4
				for(int j=n1-wD2+n1+n2+n2;j<=n1+n1+n2+n2-1;j++)
				{
					T(i,j)=T(i,j)/d_bond_y*D_bond_y;
				}
		    }
		    //7,11,15...
			//barrier
			//antibonding
			for(int i=n;i<=2*n-1;i++)
			{
				//Ubarrier1
				for(int j=0;j<=wD1-1;j++)
				{
					T(i,j)=T(i,j)/d_at_bond_x*D_at_bond_x;
				}
				//Dbarrier1
				for(int j=n1-wD1;j<=n1-1;j++)
				{
					T(i,j)=T(i,j)/d_at_bond_x*D_at_bond_x;
				}	
				//Ubarrier2
				for(int j=0+n1;j<=wD1-1+n1;j++)
				{
					T(i,j)=T(i,j)/d_at_bond_y*D_at_bond_y;
				}
				//Dbarrier2
				for(int j=n1-wD1+n1;j<=n1+n1-1;j++)
				{
					T(i,j)=T(i,j)/d_at_bond_y*D_at_bond_y;
				}
				//Ubarrier3
				for(int j=0+n1+n1;j<=wD2-1+n1+n1;j++)
				{
					T(i,j)=T(i,j)/d_at_bond_x*D_at_bond_x;
				}
				//Dbarrier3
				for(int j=n1-wD2+n1+n2;j<=n1+n1+n2-1;j++)
				{
					T(i,j)=T(i,j)/d_at_bond_x*D_at_bond_x;
				}
				//Ubarrier4
				for(int j=0+n1+n1+n2;j<=wD2-1+n1+n1+n2;j++)
				{
					T(i,j)=T(i,j)/d_at_bond_y*D_at_bond_y;
				}
				//Dbarrier4
				for(int j=n1-wD2+n1+n2+n2;j<=n1+n1+n2+n2-1;j++)
				{
					T(i,j)=T(i,j)/d_at_bond_y*D_at_bond_y;
				}					
			}
		}
	}
	else
	{
		
	}
	return T;
}



MatrixXcd hamiltonian(int n ,double k ,double Ey ,int wD1 ,int wD2 ,double D ,double d )
{
	std::complex<double> I = 1i; 
	MatrixXcd H(2*n, 2*n);
	
	for(int i=0;i<=2*n-1;i++)
	{
		for(int m=0;m<=2*n-1;m++)
		{
			H(i,m)=0;
		}
	}
	
	if( n % 2 != 0)
	{
		//odd part
		//diagonal term and diagonal -t
		int n1=(n+1)/2;
		int n2=(n-1)/2;
		
		for(int i=0; i <= n1-1; i++)
		{
			H(i ,i ) = d;
			H(i, i + n1) = -t;
			H(i + n1, i) = -t;
		}
		for(int i=n1+1-1;i<=n+1-1;i++)
		{
			H(i ,i ) = -d;
		}
		for(int i=n+1+1-1;i<=n+1+n2-1;i++)
		{
			H(i ,i ) = d;
			H(i ,i + n2 ) = -t;
			H(i + n2, i ) = -t;
		}
		for(int i=n+1+n2+1-1;i<=2*n-1;i++)
		{
			H(i ,i ) = -d;
		}
		//Exp term
		for(int i=0, j=2*n-n2+1-1;j<=2*n-1;i++,j++)
		{
			H(i ,j )=-t*exp(-3.*I*k*a);
			H(j ,i )=-t*exp(3.*I*k*a);
			H(i+1 ,j )=-t*exp(-3.*I*k*a);
			H(j ,i+1 )=-t*exp(3.*I*k*a);
		}
		//laddering -t
		for(int i=n1+1-1;i<=n-1;i++)
		{
			H(i ,i + n1) = -t;
			H(i + n1 ,i) = -t;
			H(i+1 ,i + n1) = -t;
			H(i + n1 ,i+1) = -t;
		}
		//Ey and barrier
		if( ((n+1)/2)%2 != 0 )
		{
			//n=5,9,13...
			
			//Ubarrier1
			for(int i=0;i<=wD1-1;i++)
			{	
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-d+D;

			}
			//channel1
			for(int i=wD1+1-1;i<=n1-wD1-1;i++)
			{
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-(i+1-wD1)*sqrt(3.0)*a*e*Ey;
			}
			//Dbarrier1
			for(int i=n1-wD1+1-1;i<=n1-1;i++)
			{
				H(i ,i ) = H(i ,i )-e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-d+D;
			}
			//Ubarrier2
			for(int i=n1+1-1;i<=wD1+n1-1;i++)
			{
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)+d-D;
			}
			//channel2
			for(int i=wD1+1+n1-1;i<=n1-wD1+n1-1;i++)
			{
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-(i+1-n1-wD1)*sqrt(3.0)*a*e*Ey;
			}
			//Dbarrier2
			for(int i=n1-wD1+1+n1-1;i<=n1+n1-1;i++)
			{
				H(i ,i ) = H(i ,i )-e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)+d-D;
			}
			//Ubarrier3
			for(int i=n1+n1+1-1;i<=n1+n1+wD2-1;i++)
			{
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-d+D;
			}
			//channel3
			for(int i=n1+n1+wD2+1-1;i<=n+1+n2-wD2-1;i++)
			{
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-(i+1-n1-n1-wD2+0.5)*sqrt(3.0)*a*e*Ey;
			}
			//Dbarrier3
			for(int i=n+1+n2-wD2+1-1;i<=n+1+n2-1;i++)
			{
				H(i ,i ) = H(i ,i )-e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-d+D;
			}
			//Ubarrier4
			for(int i=n1+n1+1+n2-1;i<=n1+n1+wD2+n2-1;i++)
			{
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)+d-D;
			}
			//channel4
			for(int i=i=n1+n1+wD2+1+n2-1;i<=n+1+n2-wD2+n2-1;i++)
			{
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-(i+1-n1-n1-n2-wD2+0.5)*sqrt(3.0)*a*e*Ey;
			}
			//Dbarrier4
			for(int i=n+1+n2-wD2+1+n2-1;i<=n+1+n2+n2-1;i++)
			{
				H(i ,i ) = H(i ,i )-e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)+d-D;
			}
			
		}
		else
		{
			//n=7,11,15...same as 5,9,13...
			
			//Ubarrier1
			for(int i=0;i<=wD1-1;i++)
			{	
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-d+D;
			}
			//channel1
			for(int i=wD1+1-1;i<=n1-wD1-1;i++)
			{
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-(i+1-wD1)*sqrt(3.0)*a*e*Ey;
			}
			//Dbarrier1
			for(int i=n1-wD1+1-1;i<=n1-1;i++)
			{
				H(i ,i ) = H(i ,i )-e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-d+D;
			}
			
			
			//Ubarrier2
			for(int i=n1+1-1;i<=wD1+n1-1;i++)
			{
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)+d-D;
			}
			//channel2
			for(int i=wD1+1+n1-1;i<=n1-wD1+n1-1;i++)
			{
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-(i+1-n1-wD1)*sqrt(3.0)*a*e*Ey;
			}
			//Dbarrier2
			for(int i=n1-wD1+1+n1-1;i<=n1+n1-1;i++)
			{
				H(i ,i ) = H(i ,i )-e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)+d-D;
			}
			
			
			//Ubarrier3
			for(int i=n1+n1+1-1;i<=n1+n1+wD2-1;i++)
			{
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-d+D;
			}
			//channel3
			for(int i=n1+n1+wD2+1-1;i<=n+1+n2-wD2-1;i++)
			{
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-(i+1-n1-n1-wD2+0.5)*sqrt(3.0)*a*e*Ey;
			}
			//Dbarrier3
			for(int i=n+1+n2-wD2+1-1;i<=n+1+n2-1;i++)
			{
				H(i ,i ) = H(i ,i )-e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-d+D;
			}
			
			
			//Ubarrier4
			for(int i=n1+n1+1+n2-1;i<=n1+n1+wD2+n2-1;i++)
			{
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)+d-D;
			}
			//channel4
			for(int i=i=n1+n1+wD2+1+n2-1;i<=n+1+n2-wD2+n2-1;i++)
			{
				H(i ,i ) = H(i ,i )+e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)-(i+1-n1-n1-n2-wD2+0.5)*sqrt(3.0)*a*e*Ey;
			}
			//Dbarrier4
			for(int i=n+1+n2-wD2+1+n2-1;i<=n+1+n2+n2-1;i++)
			{
				H(i ,i ) = H(i ,i )-e*Ey*(n2*(sqrt(3.0))/2*a-(wD1-1)*sqrt(3.0)*a)+d-D;
			}
		}
	}
	else
	{
		//even part
		//diagonal term and diagonal -t
		int n1=(n)/2;
		int n2=(n)/2;
		
		for(int i=0; i <= n1-1; i++)
		{
			H(i ,i ) = d;
			H(i, i + n1) = -t;
			H(i + n1, i) = -t;
		}
		for(int i=n1+1-1;i<=n-1;i++)
		{
			H(i ,i ) = -d;
		}
		for(int i=n+1+1-1;i<=n+n1-1;i++)
		{
			H(i ,i ) = d;
			H(i ,i + n2 ) = -t;
			H(i + n2, i ) = -t;
		}
		for(int i=n+1+n1-1;i<=2*n-1;i++)
		{
			H(i ,i ) = -d;
		}
		//Exp term
		for(int i=0, j=2*n-n2+1-1;j<=2*n-1;i++,j++)
		{
			H(i ,j )=-t*exp(-3.*I*k*a);
			H(j ,i )=-t*exp(3.*I*k*a);
			H(i+1 ,j )=-t*exp(-3.*I*k*a);
			H(j ,i+1 )=-t*exp(3.*I*k*a);
		}
		//laddering -t
		for(int i=n1+1-1;i<=n-1;i++)
		{
			H(i ,i + n1) = -t;
			H(i + n1 ,i) = -t;
			H(i+1 ,i + n1) = -t;
			H(i + n1 ,i+1) = -t;
		}
	}
    return H;
}

complex<double> vp(double k, int n, MatrixXcd SH, int b)
{
	std::complex<double> I = 1i; 
	complex<double> vpu=0;
	complex<double> vpd=0;
	complex<double> tot1;
	complex<double> tot2;
	
	int s0;
	int j0;
	
	
	
	double N1=0;
	double N2=0;
	
	//vpu
	if(n%2!=0)
	{
		if((n+1)/2%2!=0)
		{
			int nL=(n+1)/2;// # of half left atoms 
			int nr=(n-1)/2;// # of half left atoms
			
			int n1=(nL+1)/2;
			int n2=(nL-1)/2;
			int n3=nr/2;
			
			//normalization const
			for(int m=0;m<=nr;m++)
			{
				N1=N1+cos(2*Pi/(n+3.)*(n3-m))*cos(2*Pi/(n+3.)*(n3-m));
			}
			
			//normalization const
			for(int m=0;m<=nr-1;m++)
			{
				N2=N2+cos(2*Pi/(n+1.)*((nr-1)/2.-m))*cos(2*Pi/(n+1.)*((nr-1)/2.-m));
			}
			
			//bonding left cs mixing 1
			s0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1);	
		
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=0;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(0.5+i+1-1.)*Pi/(n*1.+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2-1;i++)
				{
					tot2=tot2+conj((1.0/I)*SH(b,i)*exp(I*4.*(i+1.-n1)*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
				
				tot1=0;
				for(int i=0;i<=s0-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(-I*4.*(0.5+i+1-1.)*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2-1;i++)
				{
					tot2=tot2+(1.0/I)*SH(b,i)*exp(I*4.*(i+1.-n1)*Pi/(n*1.+3.)*(n3-m*1.));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
						
			}
			
			//bonding left cs mixing 2
			s0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+0.5);
			j0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1);	
		
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=0;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(0.5+i+1-1.)*Pi/(n*1.+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=n1;i<=j0-1;i++)
				{
					tot2=tot2+conj((-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-n1)*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
				
				tot1=0;
				for(int i=0;i<=s0-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(-I*4.*(0.5+i+1-1.)*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				tot2=0;
				for(int i=n1;i<=j0-1;i++)
				{
					tot2=tot2+(-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-n1)*Pi/(n*1.+3.)*(n3-m*1.));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
						
			}
			
			//bonding left cs mixing 3
			s0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+0.5);
			j0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1);	
		
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=n1-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(I*4.*(0.5+i+1-1.)*Pi/(n*1.+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=n1;i<=j0-1;i++)
				{
					tot2=tot2+conj((-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-n1)*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
				
				tot1=0;
				for(int i=s0-1;i<=n1-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(I*4.*(0.5+i+1-1.)*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				tot2=0;
				for(int i=n1;i<=j0-1;i++)
				{
					tot2=tot2+(-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-n1)*Pi/(n*1.+3.)*(n3-m*1.));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
						
			}
			
			//bonding left cs mixing 4
			s0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1);	
		
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=n1-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(I*4.*(0.5+i+1-1.)*Pi/(n*1.+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2-1;i++)
				{
					tot2=tot2+conj((1.0/I)*SH(b,i)*exp(I*4.*(i+1.-n1)*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
				
				tot1=0;
				for(int i=s0-1;i<=n1-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(I*4.*(0.5+i+1-1.)*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2-1;i++)
				{
					tot2=tot2+(1.0/I)*SH(b,i)*exp(I*4.*(i+1.-n1)*Pi/(n*1.+3.)*(n3-m*1.));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
						
			}
			
			//bonding right cs mixing 1
			s0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3);	
		
			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=n1+n2;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(0.5+i-n1-n2+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3-1;i++)
				{
					tot2=tot2+conj((1.0/I)*SH(b,i)*exp(I*4.*(i+1.-n1-n2-n3)*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
				
				tot1=0;
				for(int i=n1+n2;i<=s0-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(-I*4.*(0.5+i-n1-n2+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3-1;i++)
				{
					tot2=tot2+(1.0/I)*SH(b,i)*exp(I*4.*(i+1.-n1-n2-n3)*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
						
			}
			
			//bonding right cs mixing 2
			s0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+0.5);
			j0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3);	
		
			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=n1+n2;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(0.5+i-n1-n2+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				tot2=0;
				for(int i=n1+n2+n3;i<=j0-1;i++)
				{
					tot2=tot2+conj((-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-n1-n2-n3)*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
				
				tot1=0;
				for(int i=n1+n2;i<=s0-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(-I*4.*(0.5+i-n1-n2+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				tot2=0;
				for(int i=n1+n2+n3;i<=j0-1;i++)
				{
					tot2=tot2+(-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-n1-n2-n3)*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
						
			}
			
			//bonding right cs mixing 3
			s0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+0.5);
			j0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3);	
		
			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(I*4.*(0.5+i-n1-n2+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				tot2=0;
				for(int i=n1+n2+n3;i<=j0-1;i++)
				{
					tot2=tot2+conj((-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-n1-n2-n3)*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
				
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(I*4.*(0.5+i-n1-n2+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				tot2=0;
				for(int i=n1+n2+n3;i<=j0-1;i++)
				{
					tot2=tot2+(-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-n1-n2-n3)*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
						
			}
			
			//bonding right cs mixing 4
			s0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3);	
		
			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(I*4.*(0.5+i-n1-n2+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3-1;i++)
				{
					tot2=tot2+conj((1.0/I)*SH(b,i)*exp(I*4.*(i+1.-n1-n2-n3)*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
				
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(I*4.*(0.5+i-n1-n2+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3-1;i++)
				{
					tot2=tot2+(1.0/I)*SH(b,i)*exp(I*4.*(i+1.-n1-n2-n3)*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
			}
			
			//at-bonding left cs mixing 1
			s0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1);	
		
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=n1+n2+n3+n3;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(0.5+i+1-(n1+n2+n3+n3)-1.)*Pi/(n*1.+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1+n2-1;i++)
				{
					tot2=tot2+conj((1.0/I)*SH(b,i)*exp(I*4.*(i+1.-(n1+n2+n3+n3+n1))*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
				
				tot1=0;
				for(int i=n1+n2+n3+n3;i<=s0-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(-I*4.*(0.5+i+1-(n1+n2+n3+n3)-1.)*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1+n2-1;i++)
				{
					tot2=tot2+(1.0/I)*SH(b,i)*exp(I*4.*(i+1.-(n1+n2+n3+n3+n1))*Pi/(n*1.+3.)*(n3-m*1.));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
						
			}
			
			//at-bonding left cs mixing 2
			s0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+0.5);
			j0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1);	
		
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=n1+n2+n3+n3;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(0.5+i+1-(n1+n2+n3+n3)-1.)*Pi/(n*1.+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=n1+n2+n3+n3+n1;i<=j0-1;i++)
				{
					tot2=tot2+conj((-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-(n1+n2+n3+n3+n1))*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
				
				tot1=0;
				for(int i=n1+n2+n3+n3;i<=s0-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(-I*4.*(0.5+i+1-(n1+n2+n3+n3)-1.)*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				tot2=0;
				for(int i=n1+n2+n3+n3+n1;i<=j0-1;i++)
				{
					tot2=tot2+(-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-(n1+n2+n3+n3+n1))*Pi/(n*1.+3.)*(n3-m*1.));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
			}
			
			//at-bonding left cs mixing 3
			s0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+0.5);
			j0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1);	
		
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3+n3+n1-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(I*4.*(0.5+i+1-(n1+n2+n3+n3)-1.)*Pi/(n*1.+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=n1+n2+n3+n3+n1;i<=j0-1;i++)
				{
					tot2=tot2+conj((-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-(n1+n2+n3+n3+n1))*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
				
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3+n3+n1-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(I*4.*(0.5+i+1-(n1+n2+n3+n3)-1.)*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				tot2=0;
				for(int i=n1+n2+n3+n3+n1;i<=j0-1;i++)
				{
					tot2=tot2+(-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-(n1+n2+n3+n3+n1))*Pi/(n*1.+3.)*(n3-m*1.));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
			}
			
			//at-bonding left cs mixing 4
			s0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1);	
		
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3+n3+n1-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(I*4.*(0.5+i+1-(n1+n2+n3+n3)-1.)*Pi/(n*1.+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1+n2-1;i++)
				{
					tot2=tot2+conj((1.0/I)*SH(b,i)*exp(I*4.*(i+1.-(n1+n2+n3+n3+n1))*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
				
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3+n3+n1-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(I*4.*(0.5+i+1-(n1+n2+n3+n3)-1.)*Pi/(n*1.+3.)*(n3-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1+n2-1;i++)
				{
					tot2=tot2+(1.0/I)*SH(b,i)*exp(I*4.*(i+1.-(n1+n2+n3+n3+n1))*Pi/(n*1.+3.)*(n3-m*1.));
				}
				vpu=vpu+(1/(4*N1))*tot1*tot2;	
			}
			
			//at-bonding right cs mixing 1
			s0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+n3);	
		
			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=n1+n2+n3+n3+n1+n2;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(0.5+i-(n1+n2+n3+n3+n1+n2)+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1+n2+n3+n3-1;i++)
				{
					tot2=tot2+conj((1.0/I)*SH(b,i)*exp(I*4.*(i+1.-(n1+n2+n3+n3+n1+n2+n3))*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
				
				tot1=0;
				for(int i=n1+n2+n3+n3+n1+n2;i<=s0-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(-I*4.*(0.5+i-(n1+n2+n3+n3+n1+n2)+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1+n2+n3+n3-1;i++)
				{
					tot2=tot2+(1.0/I)*SH(b,i)*exp(I*4.*(i+1.-(n1+n2+n3+n3+n1+n2+n3))*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
						
			}
			
			//at-bonding right cs mixing 2
			s0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+0.5);
			j0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+n3);	
		
			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=n1+n2+n3+n3+n1+n2;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(0.5+i-(n1+n2+n3+n3+n1+n2)+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				tot2=0;
				for(int i=n1+n2+n3+n3+n1+n2+n3;i<=j0-1;i++)
				{
					tot2=tot2+conj((-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-(n1+n2+n3+n3+n1+n2+n3))*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
				
				tot1=0;
				for(int i=n1+n2+n3+n3+n1+n2;i<=s0-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(-I*4.*(0.5+i-(n1+n2+n3+n3+n1+n2)+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				tot2=0;
				for(int i=n1+n2+n3+n3+n1+n2+n3;i<=j0-1;i++)
				{
					tot2=tot2+(-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-(n1+n2+n3+n3+n1+n2+n3))*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
						
			}
		
			//at-bonding right cs mixing 3
			s0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+0.5);
			j0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+n3);	
		
			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3+n3+n1+n2+n3-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(I*4.*(0.5+i-(n1+n2+n3+n3+n1+n2)+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				tot2=0;
				for(int i=n1+n2+n3+n3+n1+n2+n3;i<=j0-1;i++)
				{
					tot2=tot2+conj((-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-(n1+n2+n3+n3+n1+n2+n3))*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
				
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3+n3+n1+n2+n3-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(I*4.*(0.5+i-(n1+n2+n3+n3+n1+n2)+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				tot2=0;
				for(int i=n1+n2+n3+n3+n1+n2+n3;i<=j0-1;i++)
				{
					tot2=tot2+(-1.0/I)*SH(b,i)*exp(-I*4.*(i+1.-(n1+n2+n3+n3+n1+n2+n3))*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
			}
			
			//at-bonding right cs mixing 4
			s0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+n3);	
		
			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3+n3+n1+n2+n3-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(I*4.*(0.5+i-(n1+n2+n3+n3+n1+n2)+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1+n2+n3+n3-1;i++)
				{
					tot2=tot2+conj((1.0/I)*SH(b,i)*exp(I*4.*(i+1.-(n1+n2+n3+n3+n1+n2+n3))*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
				
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3+n3+n1+n2+n3-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(I*4.*(0.5+i-(n1+n2+n3+n3+n1+n2)+1-1.)*Pi/(n+1.)*((nr-1)/2.-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1+n2+n3+n3-1;i++)
				{
					tot2=tot2+(1.0/I)*SH(b,i)*exp(I*4.*(i+1.-(n1+n2+n3+n3+n1+n2+n3))*Pi/(n+1.)*((nr-1)/2.-m*1.));
				}
				vpu=vpu+(1/(4*N2))*tot1*tot2;	
				
			}
			
			//vpd
			
			vpd=0;
			
			//Bonding Left Cos Left down triangle
			
			s0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+0.5);
			j0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+0.5);

			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=0;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(i+0.5)*Pi/(n+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=0;i<=j0-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(-I*4.*(i+0.5)*Pi/(n+3.)*(n3-m*1.)));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
			}
			
			//Bonding Left Cos Middle top triangle
			
			s0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+0.5);
	
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=n1-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(I*4.*(i+0.5)*Pi/(n+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(I*4.*(i+0.5)*Pi/(n+3.)*(n3-m*1.)));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
			}
			
			//Bonding Left Cos term mixing
			
			s0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+0.5);
	
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=0;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(i+0.5)*Pi/(n+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(I*4.*(i+0.5)*Pi/(n+3.)*(n3-m*1.)));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
				
				tot1=0;
				for(int i=0;i<=s0-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(-I*4.*(i+0.5)*Pi/(n+3.)*(n3-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=n1-1;i++)
				{
					tot2=tot2+SH(b,i)*exp(I*4.*(i+0.5)*Pi/(n+3.)*(n3-m*1.));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
			}
			
			//Bonding Left Sin Left down triangle
			
			s0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1);
			j0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1);
	
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=n1;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(i+1.-n1)*Pi/(n+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=n1;i<=j0-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(-I*4.*(i+1.-n1)*Pi/(n+3.)*(n3-m*1.)));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
			}
			
			//Bonding Left Sin middle top triangle
			
			s0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1);
			j0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1);
	
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=n1+n2-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(i+1.-n1)*Pi/(n+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(-I*4.*(i+1.-n1)*Pi/(n+3.)*(n3-m*1.)));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
			}
			
			//Bonding Left Sin term mixing
			
			s0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1);
			j0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+n1);
	
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=n1;i<=s0-1;i++)
				{
					tot1=tot1+(-1./I)*SH(b,i)*exp(-I*4.*(i+1.-n1)*Pi/(n+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2-1;i++)
				{
					tot2=tot2+conj((1./I)*SH(b,i)*exp(I*4.*(i+1.-n1)*Pi/(n+3.)*(n3-m*1.)));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
				
				tot1=0;
				for(int i=n1;i<=s0-1;i++)
				{
					tot1=tot1+conj((-1./I)*SH(b,i)*exp(-I*4.*(i+1.-n1)*Pi/(n+3.)*(n3-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2-1;i++)
				{
					tot2=tot2+(1./I)*SH(b,i)*exp(I*4.*(i+1.-n1)*Pi/(n+3.)*(n3-m*1.));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
			}
			
			//Bonding right Cos Left down triangle
			
			s0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+0.5);
			j0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+0.5);

			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=n1+n2;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(i-(n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				tot2=0;
				for(int i=n1+n2;i<=j0-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(-I*4.*(i-(n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
			}
			
			//Bonding right Cos Middle top triangle
			
			s0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+0.5);
	
			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(I*4.*(i-(n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(I*4.*(i-(n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
			}
			
			//Bonding right Cos term mixing
			
			s0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+0.5);

			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=n1+n2;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(i-(n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(I*4.*(i-(n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
				
				tot1=0;
				for(int i=n1+n2;i<=s0-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(-I*4.*(i-(n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3-1;i++)
				{
					tot2=tot2+SH(b,i)*exp(I*4.*(i-(n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
			}
			
			//Bonding right Sin Left down triangle
			
			s0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3);
			j0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3);
	
			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=+n1+n2+n3;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(i+1.-(n1 + n2 + n3))*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				tot2=0;
				for(int i=+n1+n2+n3;i<=j0-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(-I*4.*(i+1.-(n1 + n2 + n3))*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
			}
			
			//Bonding right Sin middle top triangle
			
			s0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3);
			j0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3);
	
			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3+n3-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(I*4.*(i+1.-(n1 + n2 + n3))*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(I*4.*(i+1.-(n1 + n2 + n3))*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
			}
			
			//Bonding right Sin term mixing
			
			s0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3);
			j0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3);
	
			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=n1+n2+n3;i<=s0-1;i++)
				{
					tot1=tot1+(-1./I)*SH(b,i)*exp(-I*4.*(i+1.-(n1 + n2 + n3))*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3-1;i++)
				{
					tot2=tot2+conj((1./I)*SH(b,i)*exp(I*4.*(i+1.-(n1 + n2 + n3))*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
				
				tot1=0;
				for(int i=n1+n2+n3;i<=s0-1;i++)
				{
					tot1=tot1+conj((-1./I)*SH(b,i)*exp(-I*4.*(i+1.-(n1 + n2 + n3))*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3-1;i++)
				{
					tot2=tot2+(1./I)*SH(b,i)*exp(I*4.*(i+1.-(n1 + n2 + n3))*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
			}
			
			//at-Bonding Left Cos Left down triangle
			
			s0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+(n1+n2+n3+n3)+0.5);
			j0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+(n1+n2+n3+n3)+0.5);

			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=(n1+n2+n3+n3);i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3)+0.5)*Pi/(n+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=(n1+n2+n3+n3);i<=j0-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3)+0.5)*Pi/(n+3.)*(n3-m*1.)));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
			}
			
			//at-Bonding Left Cos Middle top triangle
			
			s0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+(n1+n2+n3+n3)+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+(n1+n2+n3+n3)+0.5);
	
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=(n1+n2+n3+n3+n1)-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3)+0.5)*Pi/(n+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=(n1+n2+n3+n3+n1)-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3)+0.5)*Pi/(n+3.)*(n3-m*1.)));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
			}
			
			//at-Bonding Left Cos term mixing
			
			s0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+(n1+n2+n3+n3)+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+(n1+n2+n3+n3)+0.5);
	
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=n1+n2+n3+n3;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3)+0.5)*Pi/(n+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3)+0.5)*Pi/(n+3.)*(n3-m*1.)));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
				
				tot1=0;
				for(int i=n1+n2+n3+n3;i<=s0-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3)+0.5)*Pi/(n+3.)*(n3-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1-1;i++)
				{
					tot2=tot2+SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3)+0.5)*Pi/(n+3.)*(n3-m*1.));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
			}

			//anti - Bonding Left Sin Left down triangle
			
			s0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+(n1+n2+n3+n3+n1));
			j0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+(n1+n2+n3+n3+n1));

			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=(n1+n2+n3+n3+n1);i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3+n1)+1.)*Pi/(n+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=(n1+n2+n3+n3+n1);i<=j0-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3+n1)+1.)*Pi/(n+3.)*(n3-m*1.)));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
			}
			
			//anti - Bonding Left Sin middle top triangle
			
			s0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+(n1+n2+n3+n3+n1));
			j0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+(n1+n2+n3+n3+n1));
	
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=(n1+n2+n3+n3+n1+n2)-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3+n1)+1.)*Pi/(n+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=(n1+n2+n3+n3+n1+n2)-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3+n1)+1.)*Pi/(n+3.)*(n3-m*1.)));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
			}
			
			//anti - Bonding Left Sin term mixing
			
			s0=floor(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+(n1+n2+n3+n3+n1));
			j0=ceil(sqrt(3)*fabs(k)*(n+3.)*sqrt(3)*a/(4.*Pi)+(n1+n2+n3+n3+n1));
	
			for(int m=0;m<=nr;m++)
			{
				tot1=0;
				for(int i=(n1+n2+n3+n3+n1);i<=s0-1;i++)
				{
					tot1=tot1+(-1./I)*SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3+n1)+1.)*Pi/(n+3.)*(n3-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=(n1+n2+n3+n3+n1+n2)-1;i++)
				{
					tot2=tot2+conj((1./I)*SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3+n1)+1.)*Pi/(n+3.)*(n3-m*1.)));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
				
				tot1=0;
				for(int i=(n1+n2+n3+n3+n1);i<=s0-1;i++)
				{
					tot1=tot1+conj((-1./I)*SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3+n1)+1.)*Pi/(n+3.)*(n3-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=(n1+n2+n3+n3+n1+n2)-1;i++)
				{
					tot2=tot2+(1./I)*SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3+n1)+1.)*Pi/(n+3.)*(n3-m*1.));
				}
				vpd=vpd+1/(4*N1)*tot1*tot2;
			}

			//anti - Bonding right Cos Left down triangle
			
			s0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+0.5);
			j0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+0.5);

			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=n1+n2+n3+n3+n1+n2;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3+n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				tot2=0;
				for(int i=n1+n2+n3+n3+n1+n2;i<=j0-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3+n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
			}
			
			//anti - Bonding right Cos Middle top triangle
			
			s0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+0.5);
	
			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3+n3+n1+n2+n3-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3+n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1+n2+n3-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3+n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
			}
			
			//anti - Bonding right Cos term mixing
			
			s0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+0.5);
			j0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+0.5);

			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=n1+n2+n3+n3+n1+n2;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3+n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1+n2+n3-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3+n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
				
				tot1=0;
				for(int i=n1+n2+n3+n3+n1+n2;i<=s0-1;i++)
				{
					tot1=tot1+conj(SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3+n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1+n2+n3-1;i++)
				{
					tot2=tot2+SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3+n1+n2)+0.5)*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
			}
		
			//anti - Bonding Right Sin Left down triangle
			
			s0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+n3);
			j0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+n3);

			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=n1+n2+n3+n3+n1+n2+n3;i<=s0-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3+n1+n2+n3)+1.)*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				tot2=0;
				for(int i=n1+n2+n3+n3+n1+n2+n3;i<=j0-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3+n1+n2+n3)+1.)*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
			}
		
			//anti - Bonding Right Sin middle top triangle
			
			s0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+n3);
			j0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+n3);
	
			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=s0-1;i<=n1+n2+n3+n3+n1+n2+n3+n3-1;i++)
				{
					tot1=tot1+SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3+n1+n2+n3)+1.)*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1+n2+n3+n3-1;i++)
				{
					tot2=tot2+conj(SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3+n1+n2+n3)+1.)*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
			}
			
			//anti - Bonding Right Sin term mixing
			
			s0=floor(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+n3);
			j0=ceil(sqrt(3)*fabs(k)*(n+1.)*sqrt(3)*a/(4.*Pi)+n1+n2+n3+n3+n1+n2+n3);

			for(int m=0;m<=nr-1;m++)
			{
				tot1=0;
				for(int i=n1+n2+n3+n3+n1+n2+n3;i<=s0-1;i++)
				{
					tot1=tot1+(-1./I)*SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3+n1+n2+n3)+1.)*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1+n2+n3+n3-1;i++)
				{
					tot2=tot2+conj((1./I)*SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3+n1+n2+n3)+1.)*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
				
				tot1=0;
				for(int i=n1+n2+n3+n3+n1+n2+n3;i<=s0-1;i++)
				{
					tot1=tot1+conj((-1./I)*SH(b,i)*exp(-I*4.*(i-(n1+n2+n3+n3+n1+n2+n3)+1.)*Pi/(n+1.)*((n-3.)/4.-m*1.)));
				}
				tot2=0;
				for(int i=j0-1;i<=n1+n2+n3+n3+n1+n2+n3+n3-1;i++)
				{
					tot2=tot2+(1./I)*SH(b,i)*exp(I*4.*(i-(n1+n2+n3+n3+n1+n2+n3)+1.)*Pi/(n+1.)*((n-3.)/4.-m*1.));
				}
				vpd=vpd+1/(4*N2)*tot1*tot2;
			}
	
		}
		else
			cout<<"under construction"<<endl;
	}
	else
		cout<<"under construction"<<endl;
	

	
		
	
	return vpu/vpd;
}
