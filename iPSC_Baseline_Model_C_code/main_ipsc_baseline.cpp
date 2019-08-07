// Kernik-Clancy iPSC-CM model
//**********************************************
//Kernik DC, Morotti S, Wu H, Garg P, Duff HJ, Kurokawa J, Jalife J, Wu JC, Grandi E, Clancy CE.
//A computational model of induced pluripotent stem-cell derived cardiomyocytes
//incorporating experimental variability from multiple data sources"
//J Physiol. 2019 Jul 6. doi: 10.1113/JP277724
//**********************************************
//
//Converted to C-code by Mao-Tsuen Jeng
//
//Colleen Clancy Lab @ UC davis
//
//May-21-2019




#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>


#include "integrate_rk2.h"
#include "ipsc_function.h"

int main() {
    using namespace std;
    cout.precision(16);
    
    //% Main File to generate Baseline Model Figures 10-11
    //load ICs_baseline
    int sy0 = 23;
    double y0[sy0], y0n[sy0];
    
    char filename[] = "ICs_baseline.txt";
    FILE *fp;
    fp = fopen( filename, "r" );
    for( int idy = 0; idy < sy0; idy++ ) {
        fscanf( fp, "%lf", &y0[idy] );
        y0n[idy] = y0[idy];
        cout << "Y_init[" << idy << "] = " << y0[idy] << endl;
    }
    fclose( fp );
    
    //load baseline_parameter_inputs
    int sp0 = 87;
    double p0[sp0];
    fp = fopen( "baseline_parameter_inputs.txt", "r" );
    for( int idp = 0; idp < sp0; idp++ ) {
        fscanf( fp, "%lf", &p0[idp] );
        cout << "P_init[" << idp << "] = " << p0[idp] << endl;
    }
    fclose( fp );
    
    double t, dt;
    dt = 1./128;
    t = 0;
    int sc = 30;
    double currents[sc] ;
    for( int idc = 0; idc < sc; idc++ ) {
        currents[idc] = 0;
    }
    
    FILE *fp_currents, *fp_y;
    
    fp_y = fopen( "ys.txt", "w" );
    fp_currents = fopen( "currents.txt", "w" );
    
    for( t = 0; t < 3000; t+=dt ) {
    //%% Run iPSC_function
    //options = odeset('MaxStep',1,'InitialStep',2e-2);
    //run_time=3e3;
    //[Time, values] = ode15s(@ipsc_function,[0, run_time],Y_init, options, baseline_parameter_inputs);
        
        if( fmod( t, 1. ) == 0 ) {
            fprintf( fp_y, "%16.14e\t", t );
            for( int idy = 0; idy < sy0; idy++ ) {
                fprintf( fp_y, "%16.14e\t", y0n[idy] );
            }
            fprintf( fp_y, "\n" );
            
            fprintf( fp_currents, "%16.14e\t", t );
            for( int idc = 0; idc < sc; idc++ ) {
                fprintf( fp_currents, "%16.14e\t", currents[idc] );
            }
            fprintf( fp_currents, "\n" );
        }

        integrate_rk2( ipsc_function, &t, sy0, y0n, dt, p0 , currents );
    
    
        }
    
    fclose( fp_y );
    fclose( fp_currents );
    
}
