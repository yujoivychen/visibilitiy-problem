#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <sys/time.h>
using namespace std;

//  FUNCTIONS FOR CREATING AND REMOVING MEMORY
double **matrix(int width, int height);
void free_matrix(double **f, int width, int height);
int **imatrix(int width, int height);
void free_matrix(int **f, int width, int height);

void output(ofstream &outfile, double **phi, int width, int height);
double forward_diff(int i, int j, double h, double **phi, int width, int height, int dim);
double backward_diff(int i, int j, double h, double **phi, int width, int height, int dim);
double circle(double x, double y, double center_x, double center_y, double radius);

int main()
{
   
   int width;
   int height;
   int N;

   cout << "How many steps? ";
   cin >> N;
   cout << "N = " << N << endl;


   cout << "Width/Height? ";
   cin >> width;
   height = width;
   cout << "Width = " << width << endl;
   cout << "Height = " << height << endl;


   double **v1 = matrix(width,height);
   double **v2 = matrix(width,height);
   double **phiold_1 = matrix(width,height);
   double **phiold_2 = matrix(width,height);
   double **phiold_3 = matrix(width,height);
   double **phiold_4 = matrix(width,height);
   double **phi_1 = matrix(width,height);
   double **phi_2 = matrix(width,height);
   double **phi_3 = matrix(width,height);
   double **phi_4 = matrix(width,height);
   double **phi = matrix(width,height);
   double **phiold = matrix(width,height);
   double x, y, h, delta_t, temp1, temp2, temp3, temp4, temp_1, temp_2, temp_3, temp_4, temp;
   int i, j, n, dim, r, v, vold;

   ofstream outfile("phi.dat", ios::out | ios::binary);

//Original Shape

   h = (1.0-(-1.0))/(width-1);   //////////////////////////CHANGEABLE (hx and hy)


   for (i = 0; i < width; i++)
      for (j = 0; j < height; j++)
      {
         x = -1.0+h*i;
         y = -1.0+h*j;
         phiold_1[i][j] = circle(x,y,0.2,0.2,0.1);
         phiold_2[i][j] = circle(x,y,-0.2,-0.2,0.1);
         phiold_3[i][j] = circle(x,y,0.2,-0.2,0.1);
         phiold_4[i][j] = circle(x,y,-0.2,0.2,0.1);

         temp1 = min(phiold_1[i][j],phiold_2[i][j]);
         temp2 = min(phiold_3[i][j],phiold_4[i][j]);
         phiold[i][j] = min(temp1, temp2);

      }
   

//Velocity Vectors
   for (i = 0; i < width; i++)
      for (j = 0; j < height; j++)
      {
         x = -1.0+h*i;
         y = -1.0+h*j;
         v1[i][j] = x/(sqrt(x*x+y*y));
         v2[i][j] = y/(sqrt(x*x+y*y));
      }

vold = 0.0;

//set velocity at the origin
v1[100][100]=0;
v2[100][100]=0;

for (i = 0; i < width; i++)
      for (j = 0; j < height; j++){
         v = sqrt(v1[i][j]*v1[i][j] + v2[i][j]*v2[i][j]);
         if (v > vold)
            vold = v;
      }

delta_t = h / (2.0*vold);

//Iterations
for (n = 0; n < N; n++)
   {
      for (i = 0; i < width; i++)
         for (j = 0; j < height; j++)
         {

            if(v1[i][j] <= 0){
               if(v2[i][j] <= 0){
                  temp_1 = v1[i][j]*forward_diff(i, j, h, phiold_1, width, height, 0) + v2[i][j]*forward_diff(i, j, h, phiold_1, width, height, 1);
                  temp_2 = v1[i][j]*forward_diff(i, j, h, phiold_2, width, height, 0) + v2[i][j]*forward_diff(i, j, h, phiold_2, width, height, 1);
                  phi_1[i][j] = phiold_1[i][j] - delta_t*(temp_1);
                  phi_2[i][j] = phiold_2[i][j] - delta_t*(temp_2);
                  temp_3 = v1[i][j]*forward_diff(i, j, h, phiold_3, width, height, 0) + v2[i][j]*forward_diff(i, j, h, phiold_3, width, height, 1);
                  temp_4 = v1[i][j]*forward_diff(i, j, h, phiold_4, width, height, 0) + v2[i][j]*forward_diff(i, j, h, phiold_4, width, height, 1);
                  phi_3[i][j] = phiold_3[i][j] - delta_t*(temp_3);
                  phi_4[i][j] = phiold_4[i][j] - delta_t*(temp_4);
               }
               else{
                  temp_1 = v1[i][j]*forward_diff(i, j, h, phiold_1, width, height, 0) + v2[i][j]*backward_diff(i, j, h, phiold_1, width, height, 1);
                  temp_2 = v1[i][j]*forward_diff(i, j, h, phiold_2, width, height, 0) + v2[i][j]*backward_diff(i, j, h, phiold_2, width, height, 1);
                  phi_1[i][j] = phiold_1[i][j] - delta_t*(temp_1);
                  phi_2[i][j] = phiold_2[i][j] - delta_t*(temp_2);
                  temp_3 = v1[i][j]*forward_diff(i, j, h, phiold_3, width, height, 0) + v2[i][j]*backward_diff(i, j, h, phiold_3, width, height, 1);
                  temp_4 = v1[i][j]*forward_diff(i, j, h, phiold_4, width, height, 0) + v2[i][j]*backward_diff(i, j, h, phiold_4, width, height, 1);
                  phi_3[i][j] = phiold_3[i][j] - delta_t*(temp_3);
                  phi_4[i][j] = phiold_4[i][j] - delta_t*(temp_4);
               }
            }
            else{
               if(v2[i][j] <= 0){
                  temp_1 = v1[i][j]*backward_diff(i, j, h, phiold_1, width, height, 0) + v2[i][j]*forward_diff(i, j, h, phiold_1, width, height, 1);
                  temp_2 = v1[i][j]*backward_diff(i, j, h, phiold_2, width, height, 0) + v2[i][j]*forward_diff(i, j, h, phiold_2, width, height, 1);
                  phi_1[i][j] = phiold_1[i][j] - delta_t*(temp_1);
                  phi_2[i][j] = phiold_2[i][j] - delta_t*(temp_2);
                  temp_3 = v1[i][j]*backward_diff(i, j, h, phiold_3, width, height, 0) + v2[i][j]*forward_diff(i, j, h, phiold_3, width, height, 1);
                  temp_4 = v1[i][j]*backward_diff(i, j, h, phiold_4, width, height, 0) + v2[i][j]*forward_diff(i, j, h, phiold_4, width, height, 1);
                  phi_3[i][j] = phiold_3[i][j] - delta_t*(temp_3);
                  phi_4[i][j] = phiold_4[i][j] - delta_t*(temp_4);
               }
               else{
               /*if(i == 150 && j == 150){
                  cout << phiold_1[i][j] << endl;
                  cout << phiold_1[max(i-1,0)][j] << endl;
                  cout << phiold_1[149][150] << endl;
               }*/
                  temp_1 = v1[i][j]*backward_diff(i, j, h, phiold_1, width, height, 0) + v2[i][j]*backward_diff(i, j, h, phiold_1, width, height, 1);
                  temp_2 = v1[i][j]*backward_diff(i, j, h, phiold_2, width, height, 0) + v2[i][j]*backward_diff(i, j, h, phiold_2, width, height, 1);
                  phi_1[i][j] = phiold_1[i][j] - delta_t*(temp_1);
                  phi_2[i][j] = phiold_2[i][j] - delta_t*(temp_2);
                  temp_3 = v1[i][j]*backward_diff(i, j, h, phiold_3, width, height, 0) + v2[i][j]*backward_diff(i, j, h, phiold_3, width, height, 1);
                  temp_4 = v1[i][j]*backward_diff(i, j, h, phiold_4, width, height, 0) + v2[i][j]*backward_diff(i, j, h, phiold_4, width, height, 1);
                  phi_3[i][j] = phiold_3[i][j] - delta_t*(temp_3);
                  phi_4[i][j] = phiold_4[i][j] - delta_t*(temp_4);
               /*if(i == 150 && j == 150){
                  cout << phiold_1[i][j] << endl;
                  cout << phiold_1[max(i-1,0)][j] << endl;
                  cout << phiold_1[149][150] << endl;
               }*/
               //cout << backward_diff(i, j, h, phiold_2, width, height, 0) << endl;
               }
            }
         
            temp1 = min(phiold_1[i][j],phi_1[i][j]);
            temp2 = min(phiold_2[i][j],phi_2[i][j]);
            temp3 = min(phiold_3[i][j],phi_3[i][j]);
            temp4 = min(phiold_4[i][j],phi_4[i][j]);

         
            temp_1 = min(temp1,temp2);
            temp_2 = min(temp3, temp4);
            temp = min(temp_1,temp_2);
            phi[i][j] = min(temp,phiold[i][j]);
         }

      for (i = 0; i < width; i++)
         for (j = 0; j < height; j++)
         {
            phiold[i][j] = phi[i][j];
            phiold_1[i][j] = phi_1[i][j];
            phiold_2[i][j] = phi_2[i][j];
            phiold_3[i][j] = phi_3[i][j];
            phiold_4[i][j] = phi_4[i][j];
         }
   }


output(outfile,phi,width,height);


   free_matrix(phi,width,height);
   free_matrix(phiold,width,height);
   free_matrix(v1,width,height);
   free_matrix(v2,width,height);
   free_matrix(phiold_1,width,height);
   free_matrix(phiold_2,width,height);
   free_matrix(phiold_3,width,height);
   free_matrix(phiold_4,width,height);
   free_matrix(phi_1,width,height);
   free_matrix(phi_2,width,height);
   free_matrix(phi_3,width,height);
   free_matrix(phi_4,width,height);
}

//**************************************************************************
//FUNCTIONS FOR CREATING AND REMOVING MEMORY
//**************************************************************************

double **matrix(int width, int height)
{
   double **f;
   int i;

   f = new double*[width];
   for (i = 0; i < width; i++)
      f[i] = new double[height];

   return f;
}

void free_matrix(double **f, int width, int height)
{
   int i;

   for (i = 0; i < width; i++)
      delete f[i];
   delete f;
}

int **imatrix(int width, int height)
{
   int **f;
   int i;

   f = new int*[width];
   for (i = 0; i < width; i++)
      f[i] = new int[height];

   return f;
}

void free_matrix(int **f, int width, int height)
{
   int i;

   for (i = 0; i < width; i++)
      delete f[i];
   delete f;
}

void output(ofstream &outfile, double **phi, int width, int height)
{
   int i, j;

   for (j = 0; j < height; j++)
   {
      for (i = 0; i < width; i++)
         outfile << phi[i][j] << " ";
      outfile << endl;
   }
}

double circle(double x, double y, double center_x, double center_y, double r) 
{
   double result;
   result = ((x-center_x)*(x-center_x))+((y-center_y)*(y-center_y))-(r*r);  //////////////////////////CHANGEABLE
   return result;
}

double forward_diff(int i, int j, double h, double **phi, int width, int height, int dim)
{
   double value;

   if (dim == 0) // x-direction
   {
      value = phi[min(i+1,width-1)][j] - phi[i][j];
      value = value/h;
      return value;
   }

   if(dim == 1) // y-direction
   {
      value = phi[i][min(j+1,height-1)] - phi[i][j];
      value = value/h;
      return value;
      
   }
   else{
      exit(1);
   }
}

double backward_diff(int i, int j, double h, double **phi, int width, int height, int dim)
{
   double value;


   if (dim == 0) // x-direction
   {
      value = phi[i][j] - phi[max(i-1,0)][j];
      /*if(i == 150 && j == 150){
         cout << "x " << phi[i][j] << " " << phi[max(i-1,0)][j] << " " << max(i-1,0) << endl;
      }*/
      value = value/h;
      return value;
   }

   if(dim == 1) // y-direction
   {
      value = phi[i][j] - phi[i][max(j-1,0)];
      /*if(i == 150 && j == 150){
         cout << "y " << phi[i][j] << " " << phi[i][max(j-1,0)] << " " << max(j-1,0) << endl;
      }*/
      value = value/h;
      return value;
   }
   else{
      exit(1);
   }
}
