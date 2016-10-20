#include <fstream>
#include <iostream>
#include <stdio.h>
using namespace std;

void print_vector(double*, int);
void conjugate_gradient_method(double **arr, int numr, 
    int numc, double *vb, double *vx_inv);


int main(int argc, char const *argv[])
{
    // matrix A shape
    // Also: 
    // num of unknown parameters, 
    // num of equations in system
    int numrows, numcols; 
    

    ifstream fin("in.dat");
    // Read from file
    // num of unknown parameters, 
    fin >> numcols;
    // num of equations in system
    fin >> numrows;
    fin.close();
    cout << "in.dat was read \n" << flush;
    cout << numrows << "\t" << numcols << endl << flush;


    // A*x == b
    double** a = new double*[numrows];
    for (int i = 0; i < numrows; ++i)
    {
    	a[i] = new double[numcols];
    }

    double x_inv[numcols];
    double b[numrows];



    ifstream faData("AData.dat");
    while ( !faData.eof() ) // reading A, row-major (in file)
    {
        for (int i = 0; i < numrows; ++i)
        {
            for (int j = 0; j < numcols; ++j)
            {
                faData >> a[i][j];
            }
        }
    }
    faData.close();
    cout << "AData.dat was read \n" << flush;
    print_vector(a[0], 10);


    ifstream fbData("bData.dat");
    for (int i = 0; i < numrows; ++i)
    {
        fbData >> b[i];
    }
    fbData.close();
    cout << "bData.dat was read \n" << flush;
    cout << "Vector b[:10]: \n";
    print_vector(b, 10);


    // Vector x should be empty
    for (int i = 0; i < numcols; ++i)
    {
        x_inv[i] = 0;
    }
    cout << "Vector x_inv[:10]: \n";
    print_vector(x_inv, 10);


    // Conjugate Gradient Method (from N.N. Kalitkin)
    conjugate_gradient_method(a, numrows, numcols, b, x_inv);


    ofstream f_res("MyCGM_Results.dat");
    for (int i = 0; i < numcols; ++i)
    {
        f_res << x_inv[i] << " ";
    }
    f_res.close();



}


void print_vector(double *v, int count)
{
    for (int i = 0; i < count; ++i)
    {
        printf("%2.6g \t", v[i]);
    }
    printf("\n");
}


void conjugate_gradient_method(double **arr, int numr, 
    int numc, double *vb, double *vx_inv) 
{
    // algorithm
    cout << arr[0][2] << endl;
}
