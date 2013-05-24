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
#include "defs.h"

/*****************************************************************/
/*****************************************************************
   main():
   ------

     -- Initializes, time-steps, and concludes the simulation. 
     -- Handles timing of output routines;
     -- Main is main, what more can you say.  

-*****************************************************************/
int main(int argc,char *argv[])
{
	double tdump,timage,tlog ;
	int nfailed = 0 ;
	 
	nstep = 0 ;

	/* Perform Initializations, either directly or via checkpoint */
	system("mkdir -p dumps images");  
	if(!restart_init()) { 
	  init() ;
	} 

	/* do initial diagnostics */
	diag(INIT_OUT) ;

	tdump = t+DTd ;
	timage = t+DTi ;
	tlog = t+DTl ;

	defcon = 1. ;
	while(t < tf) {


		/* step variables forward in time */
		nstroke = 0 ;
		step_ch() ;

		fprintf(stderr,"%10.5g %10.5g %8d %10.5g\n",t,dt,nstep,
				nstroke/(2.*N1*N2)) ;

		/* Handle output frequencies: */
		if(t >= tdump) {
			diag(DUMP_OUT) ;
			tdump += DTd ;
		}
		if(t >= timage) {
			diag(IMAGE_OUT) ;
			timage += DTi ;
		}
		if(t >= tlog) {
			diag(LOG_OUT) ;
			tlog += DTl ;
		}

		/* restart dump */
		nstep++ ;
		if(nstep%DTr == 0) 
			restart_write() ;


		/* deal with failed timestep, though we usually exit upon failure */
		if(failed) {
			restart_init() ;
			failed = 0 ;
			nfailed = nstep ;
			defcon = 0.3 ;
		}
		if(nstep > nfailed + DTr*4.*(1 + 1./defcon)) defcon = 1. ;


	}
	fprintf(stderr,"ns,ts: %d %d\n",nstep,nstep*N1*N2) ;

	/* do final diagnostics */
	diag(FINAL_OUT) ;

	return(0) ;
}


/*****************************************************************/
/*****************************************************************
  set_arrays():
  ----------

       -- sets to zero all arrays, plus performs pointer trick 
          so that grid arrays can legitimately refer to ghost 
          zone quantities at positions  i = -2, -1, N1, N1+1 and 
          j = -2, -1, N2, N2+1

 *****************************************************************/
void set_arrays()
{
	int i,j,k ;

	p     =  (double (*) [N2+4][NPR])(&(  a_p[2][2][0])) ;
	dq    =  (double (*) [N2+4][NPR])(&( a_dq[2][2][0])) ;
	F1    =  (double (*) [N2+4][NPR])(&( a_F1[2][2][0])) ;
	F2    =  (double (*) [N2+4][NPR])(&( a_F2[2][2][0])) ;
	ph    =  (double (*) [N2+4][NPR])(&( a_ph[2][2][0])) ;

	pflag =  (int    (*) [N2+4])(     &( a_pflag[2][2] )) ;

	/* everything must be initialized to zero */
	ZSLOOP(-2,N1+1,-2,N2+1) {
		PLOOP {
			p[i][j][k]   = 0. ;
			ph[i][j][k]  = 0. ;
			dq[i][j][k]  = 0. ;
			F1[i][j][k]  = 0. ;
			F2[i][j][k]  = 0. ;
		}
		pflag[i][j] = 0 ;
	}

	k = 0;
	IMAGELOOP {  
	    failimage[0][k] = failimage[1][k] = failimage[2][k] = failimage[3][k] = failimage[4][k++] = 0 ; 
	}


	/* grid functions */
	conn = (double (*) [N2+4][NDIM][NDIM][NDIM])
		(& ( a_conn[2][2][0][0][0])) ;
	gcon = (double (*) [N2+4][NPG][NDIM][NDIM])
		(& ( a_gcon[2][2][0][0][0])) ;
	gcov = (double (*) [N2+4][NPG][NDIM][NDIM])
		(& ( a_gcov[2][2][0][0][0])) ;
	gdet = (double (*) [N2+4][NPG])
		(& ( a_gdet[2][2][0])) ;

}


/*****************************************************************/
/*****************************************************************
  set_grid():
  ----------

       -- calculates all grid functions that remain constant 
          over time, such as the metric (gcov), inverse metric 
          (gcon), connection coefficients (conn), and sqrt of 
          the metric's determinant (gdet).

 *****************************************************************/
void set_grid()
{
	int i,j,k ;
	double X[NDIM] ;
	struct of_geom geom ;

	/* set up boundaries, steps in coordinate grid */
	set_points() ;
	dV = dx[1]*dx[2] ;

	DLOOPA X[j] = 0. ;

	ZSLOOP(-2,N1+1,-2,N2+1) {
		
		/* zone-centered */
		coord(i,j,CENT,X) ;
		gcov_func(X, gcov[i][j][CENT]) ;
		gdet[i][j][CENT] = gdet_func(gcov[i][j][CENT]) ;
		gcon_func(gcov[i][j][CENT], gcon[i][j][CENT]) ;

		get_geometry(i,j,CENT,&geom) ;
		conn_func(X, &geom, conn[i][j]) ;

		/* corner-centered */
		coord(i,j,CORN,X) ;
		gcov_func(X,gcov[i][j][CORN]) ;
		gdet[i][j][CORN] = gdet_func(gcov[i][j][CORN]) ;

		gcon_func(gcov[i][j][CORN],gcon[i][j][CORN]) ;

		/* r-face-centered */
		coord(i,j,FACE1,X) ;
		gcov_func(X,gcov[i][j][FACE1]) ;
		gdet[i][j][FACE1] = gdet_func(gcov[i][j][FACE1]) ;
		gcon_func(gcov[i][j][FACE1],gcon[i][j][FACE1]) ;

		/* theta-face-centered */
		coord(i,j,FACE2,X) ;
		gcov_func(X,gcov[i][j][FACE2]) ;
		gdet[i][j][FACE2] = gdet_func(gcov[i][j][FACE2]) ;
		gcon_func(gcov[i][j][FACE2],gcon[i][j][FACE2]) ;
	}

	/* done! */
}

