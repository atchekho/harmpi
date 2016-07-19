//Modified by Alexander Tchekhovskoy: MPI+3D
/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble, 
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic 
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to 
    solve the relativistic magnetohydrodynamic equations of motion on a 
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model. 

    You are morally obligated to cite the following two papers in his/her 
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003, 
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

   
    Further, we strongly encourage you to obtain the latest version of 
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

/* ranc -- return a random deviate between 0 and 1
   double ranc(iseed)
      * iseed = integer seed for random number generator
      - on first call, will seed with iseed (if > 0) or with time(0)
      - will reseed anytime iseed > 0
*/

#include <stdlib.h>
#include <time.h>

#define NRANC 64
#define RND	( 0x7fff & rand() )
#define MIN	( 2 << 8 )

    static int P[NRANC] = {
	46337, 46327, 46309, 46307, 46301, 46279, 46273, 46271,
	46261, 46237, 46229, 46219, 46199, 46187, 46183, 46181,
	46171, 46153, 46147, 46141, 46133, 46103, 46099, 46093,
	46091, 46073, 46061, 46051, 46049, 46027, 46021, 45989,
	45979, 45971, 45959, 45953, 45949, 45943, 45893, 45887,
	45869, 45863, 45853, 45841, 45833, 45827, 45823, 45821,
	45817, 45779, 45767, 45763, 45757, 45751, 45737, 45707,
	45697, 45691, 45677, 45673, 45667, 45659, 45641, 45631
	} ;

    static int a[NRANC] ;
    static int S[NRANC] ;
    static int n = 0 ;
    int called={1} ;

double ranc(int iseed)
{
    int i ;

/* seed the random number generator if first time called or if iseed > 0 */ 
    if (called || iseed != 0)  { 
	called = 0 ;
	if (iseed==0)
        iseed = 24; //time(0) ;
	srand(iseed) ;
	n = 0 ;
	for( i = 0 ; i < NRANC ; i++ ) {
	    a[i] = 0.0 ;
	    S[i] = 0.0 ;
	    while( ( a[i] = RND ) < MIN || a[i] > P[i] - MIN ) ;
	    while( ( S[i] = RND ) <   1 || a[i] > P[i] -   1 ) ;
	}
    }

    n = S[n] & ( NRANC - 1 ) ;
    S[n] = ( S[n] * a[n] ) % P[n] ;
    return (double) S[n] / (double) P[n] ;
}

#undef NRANC 
#undef RND
#undef MIN
