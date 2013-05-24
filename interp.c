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

#include "decs.h"

/* performs the slope-limiting for the numerical flux calculation */

double slope_lim(double y1,double y2,double y3) 
{
	double Dqm,Dqp,Dqc,s ;

	/* woodward, or monotonized central, slope limiter */
	if(lim == MC) {
		Dqm = 2.*(y2 - y1) ;
		Dqp = 2.*(y3 - y2) ;
		Dqc = 0.5*(y3 - y1) ;
		s = Dqm*Dqp ;
		if(s <= 0.) return 0. ;
		else {
			if(fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
				return(Dqm) ;
			else if(fabs(Dqp) < fabs(Dqc))
				return(Dqp) ;
			else
				return(Dqc) ;
		}
	}
	/* van leer slope limiter */
	else if(lim == VANL) {
		Dqm = (y2 - y1) ;
		Dqp = (y3 - y2) ;
		s = Dqm*Dqp ;
		if(s <= 0.) return 0. ;
		else
			return(2.*s/(Dqm+Dqp)) ;
	}

	/* minmod slope limiter (crude but robust) */
	else if(lim == MINM) {
		Dqm = (y2 - y1) ;
		Dqp = (y3 - y2) ;
		s = Dqm*Dqp ;
		if(s <= 0.) return 0. ;
		else if(fabs(Dqm) < fabs(Dqp)) return Dqm ;
		else return Dqp ;
	}

	fprintf(stderr,"unknown slope limiter\n") ;
	exit(10) ;

	return(0.) ;
}

