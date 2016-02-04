#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "complex.h"
#include <fftw3.h>

using namespace std;

double randNormal(const double, const double);

double f(int x, int y);

int main()
{
  int x0, y0, x;
  int n=384, ii, idx;
  int M = 2*n-1;
  int N = M*M;
  double stackRq;
  double **row, **col, **Rows, **Cols;
  stackRq = ((2*n*n + M*M + n*M)*8.0)/1e3;
  cout << "Requires " << stackRq << " kB in stack " << endl;
  row = new double*[M]; col = new double*[n];
  Rows = new double*[n]; Cols = new double*[n];
  for (int i=0;i<M;i++)
    {
      row[i] = new double[M];
    }
  for (int i=0;i<n;i++)
    {
      Rows[i] = new double[n];
      Cols[i] = new double[n];
      col[i] = new double[M];
    }
	     
  // double Rows[n][n]; double Cols[n][n];
  // double row[M][M];
  // double col[n][M];
 
  fftw_complex Z;//, a;
  double a;
  double map_[n][n];
  // FOR FFT
  fftw_plan p;
  fftw_complex *G, *Gamma, *GammaZ, *F;

  
  // Span n^2*n^2 covariance matrix to gather
  // first rows and columns of each n*n Toeplitz matrix
  x0 = y0 = 0;
  for(int RIdx=0;RIdx<n;RIdx++)
    {
      for(int y=0;y<n;y++)
	{
	  x = RIdx;
	  //1st row : corr between origin (x0,y0) and (x,y). 
	  Rows[RIdx][y] = f(x-x0, y-y0);
	  //1st column : corr between (x,y0) and (x0,y).
	  Cols[RIdx][y] = f(x0-x, y-y0);
	}
    }
  
  //LOOP on the n n*n (small) circulant matrices
  for(int RIdx=0;RIdx<n;RIdx++)
    {
	  // Reconstruct the row
	  //Gather the first n elements (already known)
	  for(int j=0;j<n;j++)
	    {
	      row[RIdx][j] = Rows[RIdx][j];
	      col[RIdx][j] = Cols[RIdx][j];
	    }
	  // Adds the additional elements to make the matrix circulant
	  ii = 0;
	  for(int j=n;j<M;j++)
	    {
	      ii++;
	      row[RIdx][j] = Cols[RIdx][n-ii];
	      col[RIdx][j] = Rows[RIdx][n-ii];
	    }
    }
  //Now pad the row with the two (M) long rows of last circulant blocks.
  // These are transp(C_n),...,transp(C_2)
  // row of transp(C) == column of C
  ii=0;
  for(int RIdx=n;RIdx<M;RIdx++)
    {
      ii++;
      for(int j=0;j<M;j++)
	{
	  row[RIdx][j] = col[n-ii][j];
	}
    }

  // ------- VALIDATED UP TO HERE -------------


  G = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*M*M);
  Gamma = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*M*M);
  p = fftw_plan_dft_2d(M, M, G, Gamma, FFTW_FORWARD, FFTW_ESTIMATE);


  // CREATE VECTOR G PRIOR TO FFT.
  for(int RIdx=0;RIdx<M;RIdx++)
    {
      for (int j=0;j<M;j++)
  	{
	  idx = j + RIdx*M;
	  //idx = RIdx + j*M;

	  G[idx] = row[RIdx][j]/N;
	  // ONE VERIFIES THAT G IS (2*n-1)*(2*n-1) long. OK.
	  // DIVIDES DIRECTLY BY SIZE OF VECTOR (2*n-1)X(2*n-1)
  	}
    }


  // Compute eigen values and free G to save memory

  fftw_execute(p);

  fftw_free(G);
  // Clear plan
  fftw_destroy_plan(p);
  
  //CREATE NEW PLAN FOR NEXT FFT

  GammaZ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  F = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  p = fftw_plan_dft_2d(M, M, GammaZ, F, FFTW_FORWARD, FFTW_ESTIMATE);
  //ALLOCATE GAMMAZ
  
  // COMPUTE DOT PRODUCT BETWENN SQRT(GAMMA)
  // AND GAUSSIAN IMAGINARY RAND. VECTOR Z.
  ii=0;
  for(int i=0;i<N;i++)
    {
      //a = csqrt(Gamma[i]);
      //cout << creal(Gamma[i]) << endl;
      a = creal(Gamma[i]);
      if(a<0){a=0.0;}
      Z = randNormal(0,1) + randNormal(0,1)*I;
      GammaZ[i] = sqrt(a)*Z;
    }

  //CLEAR GAMMA, ALLOCATE F AND DO FFT TO GET F
  fftw_free(Gamma);
  fftw_execute(p);
  // CLEAR GAMMAZ AND FFT PLAN
  fftw_free(GammaZ); fftw_destroy_plan(p);
  

  //EXTRACT SUB BLOCK
  for (int i=0;i<n;i++)
    {
      for (int j=0;j<n;j++)
	{
	  idx = j+M*i;
	  //idx = i + n*j;
	  map_[i][j] = creal(F[idx]);
	}
    }

  //WRITE NOISE MAP ON FILE MATRIX.DAT TO BE READ WITH GNUPLOT (IMAGE MODE)
  ofstream result("result_matrix.dat");
  for(int i=0;i<n;i++)
    {
      for(int j=0;j<n;j++)
	{
	  result << map_[j][i] << " ";
	}
      result << endl;
    }

   result.close();
}

double f(int x, int y)
{
  double h1, h2, a, b, c, e;
  h1 = (double)x;
  h2 = (double)y;
  a = (h1*h1)/(50*50);
  b = (h1*h2)/(50*15);
  c = (h2*h2)/(15*15);

  e = exp(-a-c);
  return (1.0-a-b-c)*e;

  // double d, a;
  // double sigma = 5.0;
  // d = x*x + y*y;
  // a = exp(-d/(2*sigma*sigma));
  // return a;
	  

}

double randNormal(const double mean_, const double sigma_)
{
  /* Return a random number sampled in N(mean_, sigma_).
     Box-Muller method.
  */

  double x1, x2, w;
  do {
    x1 = 2.0 * (rand () / double (RAND_MAX)) - 1.0;
    x2 = 2.0 * (rand () / double (RAND_MAX)) - 1.0;
    w = x1 * x1 + x2 * x2;
  } while (w >= 1.0);

  w = sqrt (-2.0 * log (w)/w);
  const double y1 = x1 * w;
  const double y2 = x2 * w;

  return mean_ + y1 * sigma_;
}
