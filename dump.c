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

void dump(FILE *fp)
{
	int i,j,k ;
	double divb ;
	double X[NDIM] ;
	double r,th,vmin,vmax ;
	struct of_geom geom ;
	struct of_state q ;

	/***************************************************************
	  Write header information : 
	***************************************************************/

	fprintf(fp, FMT_DBL_OUT, t        );
	fprintf(fp, FMT_INT_OUT, N1       );
	fprintf(fp, FMT_INT_OUT, N2       );
	fprintf(fp, FMT_DBL_OUT, startx[1]);
	fprintf(fp, FMT_DBL_OUT, startx[2]);
	fprintf(fp, FMT_DBL_OUT, dx[1]    );
	fprintf(fp, FMT_DBL_OUT, dx[2]    );
	fprintf(fp, FMT_DBL_OUT, tf       );
	fprintf(fp, FMT_INT_OUT, nstep    );
	fprintf(fp, FMT_DBL_OUT, a        );
	fprintf(fp, FMT_DBL_OUT, gam      );
	fprintf(fp, FMT_DBL_OUT, cour     );
	fprintf(fp, FMT_DBL_OUT, DTd      );
	fprintf(fp, FMT_DBL_OUT, DTl      );
	fprintf(fp, FMT_DBL_OUT, DTi      );
	fprintf(fp, FMT_INT_OUT, DTr      );
	fprintf(fp, FMT_INT_OUT, dump_cnt );
	fprintf(fp, FMT_INT_OUT, image_cnt);
	fprintf(fp, FMT_INT_OUT, rdump_cnt);
	fprintf(fp, FMT_DBL_OUT, dt       );
	fprintf(fp, FMT_INT_OUT, lim      );
	fprintf(fp, FMT_INT_OUT, failed   );
	fprintf(fp, FMT_DBL_OUT, Rin      );
	fprintf(fp, FMT_DBL_OUT, Rout     );
	fprintf(fp, FMT_DBL_OUT, hslope   );
	fprintf(fp, FMT_DBL_OUT, R0       );

	fprintf(fp,"\n") ;
		
	/***************************************************************
	  Write header information : 
	***************************************************************/

	ZSLOOP(0,N1-1,0,N2-1) {
		coord(i,j,CENT,X) ;
		bl_coord(X,&r,&th) ;

		fprintf(fp, FMT_INT_OUT, i       );
		fprintf(fp, FMT_INT_OUT, j       );
		fprintf(fp, FMT_DBL_OUT, X[1]       );
		fprintf(fp, FMT_DBL_OUT, X[2]       );
		fprintf(fp, FMT_DBL_OUT, r          );
		fprintf(fp, FMT_DBL_OUT, th         );
		PLOOP fprintf(fp, FMT_DBL_OUT, p[i][j][k] );
		
                /* divb flux-ct defn; corner-centered.  Use
		   only interior corners */
		if(i > 0 && j > 0 && i < N1 && j < N2) {
			divb = fabs( 0.5*(
				p[i][j][B1]*gdet[i][j][CENT]
				+ p[i][j-1][B1]*gdet[i][j-1][CENT]
				- p[i-1][j][B1]*gdet[i-1][j][CENT]
				- p[i-1][j-1][B1]*gdet[i-1][j-1][CENT]
				)/dx[1] +
				0.5*(
				p[i][j][B2]*gdet[i][j][CENT]
				+ p[i-1][j][B2]*gdet[i-1][j][CENT]
				- p[i][j-1][B2]*gdet[i][j-1][CENT]
				- p[i-1][j-1][B2]*gdet[i-1][j-1][CENT]
				)/dx[2]) ;
		}
		else divb = 0. ;

		fprintf(fp, FMT_DBL_OUT, divb     );


		if(!failed) {
			get_geometry(i,j,CENT,&geom) ;
			get_state(p[i][j],&geom,&q) ;

			for(k=0;k<NDIM;k++) fprintf(fp,FMT_DBL_OUT,q.ucon[k]) ;
			for(k=0;k<NDIM;k++) fprintf(fp,FMT_DBL_OUT,q.ucov[k]) ;
			for(k=0;k<NDIM;k++) fprintf(fp,FMT_DBL_OUT,q.bcon[k]) ;
			for(k=0;k<NDIM;k++) fprintf(fp,FMT_DBL_OUT,q.bcov[k]) ;

			vchar(p[i][j],&q,&geom,1,&vmax,&vmin) ;
			fprintf(fp, FMT_DBL_OUT, vmin );
			fprintf(fp, FMT_DBL_OUT, vmax );

			vchar(p[i][j],&q,&geom,2,&vmax,&vmin) ;
			fprintf(fp, FMT_DBL_OUT, vmin );
			fprintf(fp, FMT_DBL_OUT, vmax );

			fprintf(fp, FMT_DBL_OUT, geom.g );
		}

		fprintf(fp,"\n") ;
	}
}

