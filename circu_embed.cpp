#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "complex.h"
//#include <fftw3.h>

using namespace std;

double randNormal(const double, const double);

double f(int x, int y);

int main()
{
  int x0, y0, x;
  int n=3, ii;
  double Rows[n][n]; double Cols[n][n];
  double row[2*n-1][2*n-1]; double col[n][2*n-1];
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
	  for(int j=n;j<2*n-1;j++)
	    {
	      ii++;
	      row[RIdx][j] = Cols[RIdx][n-ii];
	      col[RIdx][j] = Rows[RIdx][n-ii];
	    }
    }
  //Now pad the row with the two (2*n-1) long rows of last circulant blocks.
  // These are transp(C_n),...,transp(C_2)
  // row of transp(C) == column of C
  ii=0;
  for(int RIdx=n;RIdx<2*n-1;RIdx++)
    {
      ii++;
      for(int j=0;j<2*n-1;j++)
	{
	  row[RIdx][j] = col[n-ii][j];
	}
    }
  for(int RIdx=0;RIdx<2*n-1;RIdx++)
    {
      cout << "---------" << endl;
      for (int j=0;j<2*n-1;j++)
  	{
  	  cout << row[RIdx][j] << endl;
  	}
    }

  // Compute eigen values
  ... WIP
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
	  

}
  
