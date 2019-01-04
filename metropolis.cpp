#include <iostream>
#include <cstdlib>
#include <time.h>
#include<math.h>
#include<fstream>

using namespace std;

int del(int n, int L) { //////////////////////////////////////////////////////////////
    // takes a cyclic boundary position
    if( n<0 ) return (L - 1);
    else if ( n == L ) return 0;
    else return n;
}

int **randomSpins(int L) { //////////////////////////////////////////////////////////
    int ll = 1;
    //returns random spin configuration of LxL matrix. 
    int **arrpp;

    arrpp = new int*[L];

    for(int i=0; i<L; i++) { 
        arrpp[i] = new int[L];
    }
   for(int h=0; h<L; h++) {
        for(int w=0; w<L; w++) {
            srand( 2*time( NULL ) + ll );
            int k = rand()%2;
            if (k) arrpp[h][w] = 1;    
            else arrpp[h][w] = -1;    
            

                ll = ll*ll + 3*ll + 5;
        }
        cout << "\n";
    }
    return arrpp;
    for(int i = 0; i < L; i++) {
    delete[] arrpp[i];
    }
    delete[] arrpp;
}

double sumsRound2D ( double y ,double z, int m ,int  n ,int** a , int L ) {///////////////////////////////////////////
    return y* (double) a[m][n]*((double) a[m][del(n-1, L)] + a[del(m+1, L)][n] + a[m][del(n+1, L)] + a[del(m-1, L)][n])
    + z*a[m][n]*( a[del(m-1,L)][del(n-1, L)] + a[del(m+1, L)][del(n-1, L)] + a[del(m-1,L)][del(n+1, L)] + a[del(m+1, L)][del(n+1, L)]);
}

double sums2D(double y, double z,int** a, int L) {/////////////////////////////////////////////////////////////////////
    double sum = 0.0;

    for(int i = 0; i < L; i++) {
        for(int j = 0; j < L; j++) {
            sum = sumsRound2D(y,z,i,j,a,L) + sum;
        }
    }
    return sum;
}

    
double weight2D(double oldsum,double newsum){////////////////////////////////////////////////
    return exp(-1.0/2.0*(oldsum - newsum));
}



int** sweep2D(double y, double z,int** a,int L) {////////////////////////////////////////////////////////////////////////////////

    double oldsum = sums2D(y,z,a,L);

    for(int i=0;i < L; i++){
        for(int j=0;j < L; j++){

            a[i][j] = - a[i][j];
            
            double r =(double) rand()/RAND_MAX;
            double newsum = sums2D(y,z,a,L);
            if( weight2D(oldsum,newsum) <=r ) a[i][j] = - a[i][j];
            else oldsum = newsum;

        }
    }

    return a;
    for(int i = 0; i < L; i++) {
    delete[] a[i];
    }
    delete[] a;

}

double energy2D(double y, double z, int** a, double L) { ////////////////////////////////////////////////

  return (double) sums2D(y,z,a,L)/4.0/L/L/y;
}

int** metropolis2D(int iter, double y, double z, int** a, int L){ //////////////////////////////////

    for(int i = 0; i < iter; i++) {

        a = sweep2D(y,z,a,L);

    }

    return a;
    for(int i = 0; i < L; i++) {
    delete[] a[i];
    }
    delete[] a;

}

void convertit(int** arrpp, int L){////////////////////////////////////////////////////////

     for( int i = 0; i < L; i ++ ) {
        for( int j = 0; j < L; j ++ ) {
           if(arrpp[i][j] == 1){
               cout << "+ ";
            }

            else{
                cout << "- ";
                }
        }
        cout << "\n"; 
     }
}



int main ( ) {  ////////////////////////////////////////////////////////////////////////////

    double val = 0;
    double cutoff = 0.05;

    int oldtime = time(NULL);

    int aveiter = 100;
    double* energy = &val;
    int sweeps = 200;

    double y = 0.4;
    double z = 0.0;
    int L = 16;
    int** arrpp = randomSpins(L);


    ofstream myDat2;
    myDat2.open("myDat2_a10_sz23_sw100.dat",ios::app);
    convertit(arrpp,L);

    cout << " " << endl;

while (cutoff <= 1) {
    int innterSTime = time(NULL);
        *energy = 0; 
    for(int i = 0; i < aveiter; i++) {
        arrpp = metropolis2D(sweeps,cutoff,z,arrpp,L);
        *energy = *energy + (double) energy2D(cutoff,z,arrpp,L)/(double) aveiter;

    }

    int innterFTime = time(NULL);

    convertit(arrpp,L);

    cout << "( " << cutoff << " , " << val << " ) " <<  "Time: " << -innterSTime + innterFTime << " seconds" << endl;

    myDat2 << cutoff << "," << val <<"\n";    
    myDat2.close();
    cutoff += 0.05;
}

     cout << "Average energy is: " << *energy << endl;

     int newtime = time(NULL);
     cout<< "Time in seconds: " << newtime - oldtime << endl;
}
