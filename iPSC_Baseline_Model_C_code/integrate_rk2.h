/*
 *  integrate_rk12.h
 *  
 *
 *  Created by Mao-Tsuen Jeng on 5/21/19.
 *  Copyright 2019 __Clancy Lab__. All rights reserved.
 *
 *
 *  Use RK2 to compute ODE integration
 *  Integrate with order 2 accuracy 
 *
 */

#ifndef integrate_rk2_H
#define integrate_rk2_H

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
// #include <omp.h>

void integrate_rk2( void (*f)( double , double *, double * , double * , double * ), double * tspan, int sy0, double * y0, double dt, double * pin , double * currents );

// f : function that evaluate ydots ( 1st derivatives ).
// tspan : time interval for integration
// dt : time step size
// y0 : initial value
// t1 : array of saved time.
// sy0 : size of y0

using namespace std;

void integrate_rk2( void (*f)( double , double *, double *, double * , double *  ), double * tspan, int sy0, double * y0, double dt, double * pin , double * currents ){
	
	cout.precision(16);
	
	double tol = 1.0;
	int counter, counter_sy , sy;
	double t, tn, err_y;
	double y1[sy0], y2[sy0], z[sy0], s, sm;
	double dy1[sy0], dy2[sy0];
	double ydot[sy0];
	char runType[] = "ydot";
	int i, j, k;
	//	int st = floor( 0.5 + ( tspan[1] - tspan[0] ) / ( dt * sy ) );
	// cout << st << "\t" << dt << "\t" << sy << endl;
	
	for( i = 0; i < sy0; i++ ){
		y1[i] = y0[i];
	}
	t = tspan[0];
	
	
	//	f( y1,  ydot , Istim , currents, celltype, gendertype, gterm , dt );
	f( t, y1, pin, ydot , currents );
	
	for( i = 0; i < sy0; i++ ) {
		dy1[i] = ydot[i] * dt; // k1
		y2[i] = y1[i] + dy1[i];
	}
	
    t += dt;
	
	//	f( y2, ydot , Istim , currents, celltype, gendertype, gterm , dt );
	f( t, y2, pin, ydot , currents );
	
	for( i = 0; i < sy0; i++ ) {
		dy2[i] = ydot[i] * dt; // k2
		z[i] = y1[i] + 0.5 * dy1[i] + 0.5 * dy2[i];
	}
	
	
	//t += dt;
	// Result using RK2
	sm = 1;
	
	for( i = 0; i < sy0; i++ ) {
		
		if ( z[i] != y2[i] && z[i] == z[i] ) {
			s = pow( tol * dt / ( 2 * fabs( z[i] - y2[i] ) ) , 0.5 );
			// cout << i << "\t" << z[i] << "\t" << y2[i] << "\t" << ydot[i] << endl;
			
			if( sm > s && s == s ) {
				
				sm = s;
				if( 1E-2 > sm ) {
					sm = 1E-2;
					// cout << z[i] << "\t" << y2[i] << endl;
				}
				// cout << s[i] << "\t";
				
			} else if ( s != s ) { // if s is NAN
				// cout << s << endl;
				sm = 1E-2;
			}
			//cout << endl;
		} else if ( z[i] != z[i] ){
			// cout << z[i] << endl;  // z[i] is NAN
			sm = 1E-2;
		}
	}
	if( sm >= 1 ){
		for( i = 0; i < sy0; i++ ) {
			// cout << "dt = " << dt << endl;
			y1[i] = z[i];
			// y1[i] = y2[i];		
		}
	} else {
		
		sy = ceil( 1./ sm );
		//	sy = 1;
		// cout << "Divided by " << sy << " parts.  sm = " << sm  << endl;
		
		dt = dt / sy;
		t = tspan[0];
		for( i = 0; i < sy0; i++ ){
			y1[i] = y0[i];
		}
		
		for( counter_sy = 0; counter_sy < sy; counter_sy ++ ) {
			
			// f( y1,  ydot , Istim , currents, celltype, gendertype, gterm , dt );
			f( t, y1, pin, ydot , currents );
			
			for( i = 0; i < sy0; i++ ) {
				dy1[i] = ydot[i] * dt;
				y2[i] = y1[i] + dy1[i];
			}
			
            t += dt;
			// f( y2, ydot , Istim , currents, celltype, gendertype, gterm , dt );
			f( t, y2, pin, ydot , currents );
			
			for( i = 0; i < sy0; i++ ) {
				dy2[i] = ydot[i] * dt;
				y1[i] = y1[i] + 0.5 * ( dy1[i] + dy2[i] );
			}
			
			// t += dt;
			
		}
		// cout << endl;
	}
	
	for( i = 0; i < sy0; i++ ) {
		y0[i] = y1[i];
		// y1[i] = y2[i];		
	}
}

#endif
