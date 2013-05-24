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

#define FLOOP for(k=0;k<B1;k++)

/* apply floors to density, internal energy */

void fixup(double (* pv)[N2+4][NPR]) 
{
  int i,j ;

  ZLOOP {
    fixup1zone( i, j, pv[i][j] );
  }

}

void fixup1zone( int i, int j, double pv[NPR] ) 
{
  double r,th,X[NDIM],uuscal,rhoscal, rhoflr,uuflr ;
  double f,gamma ;
  double bsq;
  struct of_geom geom ;

  coord(i,j,CENT,X) ;
  bl_coord(X,&r,&th) ;

  rhoscal = pow(r,-1.5) ;
  //rhoscal = pow(r,-2) ;
  uuscal = rhoscal/r ;

  rhoflr = RHOMIN*rhoscal;
  uuflr  = UUMIN*uuscal;

  //compute the square of fluid frame magnetic field (twice magnetic pressure)
  get_geometry(i,j,CENT,&geom) ;
  bsq = bsq_calc(pv,&geom) ;
  
  //tie floors to the local values of magnetic field and internal energy density
  if( rhoflr < bsq / BSQORHOMAX ) rhoflr = bsq / BSQORHOMAX;
  if( uuflr < bsq / BSQOUMAX ) uuflr = bsq / BSQORHOMAX;
  if( rhoflr < pv[UU] / UORHOMAX ) rhoflr = pv[UU] / UORHOMAX;

  if( rhoflr < RHOMINLIMIT ) rhoflr = RHOMINLIMIT;
  if( uuflr  < UUMINLIMIT  ) uuflr  = UUMINLIMIT;

  /* floor on density and internal energy density (momentum *not* conserved) */
  if(pv[RHO] < rhoflr )   pv[RHO] = rhoflr; 
  if(pv[UU]  < uuflr  )   pv[UU]  = uuflr;


  /* limit gamma wrt normal observer */

  if( gamma_calc(pv,&geom,&gamma) ) { 
    /* Treat gamma failure here as "fixable" for fixup_utoprim() */
    pflag[i][j] = -333;
    failimage[3][i+j*N1]++ ;
  }
  else { 
    if(gamma > GAMMAMAX) {
      f = sqrt(
	       (GAMMAMAX*GAMMAMAX - 1.)/
	       (gamma*gamma - 1.)
	       ) ;
      pv[U1] *= f ;	
      pv[U2] *= f ;	
      pv[U3] *= f ;	
    }
  }

  return;
}


/**************************************************************************************
 INTERPOLATION STENCILS:  
 ------------------------
   -- let the stencils be characterized by the following numbering convention:

           1 2 3 
           8 x 4      where x is the point at which we are interpolating 
           7 6 5
*******************************************************************************************/

/* 12345678 */
#define AVG8(pr,i,j,k)  \
        (0.125*(pr[i-1][j+1][k]+pr[i][j+1][k]+pr[i+1][j+1][k]+pr[i+1][j][k]+pr[i+1][j-1][k]+pr[i][j-1][k]+pr[i-1][j-1][k]+pr[i-1][j][k])) 

/* 2468  */
#define AVG4_1(pr,i,j,k) (0.25*(pr[i][j+1][k]+pr[i][j-1][k]+pr[i-1][j][k]+pr[i+1][j][k]))

/* 1357  */
#define AVG4_2(pr,i,j,k) (0.25*(pr[i+1][j+1][k]+pr[i+1][j-1][k]+pr[i-1][j+1][k]+pr[i-1][j-1][k]))

/* + shaped,  Linear interpolation in X1 or X2 directions using only neighbors in these direction */
/* 48  */
#define AVG2_X1(pr,i,j,k) (0.5*(pr[i-1][j  ][k]+pr[i+1][j  ][k]))
/* 26  */
#define AVG2_X2(pr,i,j,k) (0.5*(pr[i  ][j-1][k]+pr[i  ][j+1][k]))

/* x shaped,  Linear interpolation diagonally along both X1 and X2 directions "corner" neighbors */
/* 37  */
#define AVG2_1_X1X2(pr,i,j,k) (0.5*(pr[i-1][j-1][k]+pr[i+1][j+1][k]))
/* 15  */
#define AVG2_2_X1X2(pr,i,j,k) (0.5*(pr[i-1][j+1][k]+pr[i+1][j-1][k]))

/*******************************************************************************************
  fixup_utoprim(): 

    -- figures out (w/ pflag[]) which stencil to use to interpolate bad point from neighbors;

    -- here we use the following numbering scheme for the neighboring cells to i,j:  

                      1  2  3 
                      8  x  4        where "x" is the (i,j) cell or the cell to be interpolated
                      7  6  5

 *******************************************************************************************/

void fixup_utoprim( double (*pv)[N2 + 4][NPR] )  
{
  int i, j, k;
  static int pf[9];

  /* Flip the logic of the pflag[] so that it now indicates which cells are good  */
  ZSLOOP(-2,(N1+1),-2,(N2+1)) { pflag[i][j] = !pflag[i][j] ; } 

  /* Fix the interior points first */
  ZSLOOP(0,(N1-1),0,(N2-1)) { 
    if( pflag[i][j] == 0 ) { 
      pf[1] = pflag[i-1][j+1];   pf[2] = pflag[i][j+1];  pf[3] = pflag[i+1][j+1];
      pf[8] = pflag[i-1][j  ];                           pf[4] = pflag[i+1][j  ];
      pf[7] = pflag[i-1][j-1];   pf[6] = pflag[i][j-1];  pf[5] = pflag[i+1][j-1];

      /* Now the pf's  are true if they represent good points */

//      if(      pf[1]&&pf[2]&&pf[3]&&pf[4]&&pf[5]&&pf[6]&&pf[7]&&pf[8] ){ FLOOP pv[i][j][k] = AVG8(            pv,i,j,k)                   ; }
//      else if(        pf[2]&&       pf[4]&&       pf[6]&&       pf[8] ){ FLOOP pv[i][j][k] = AVG4_1(          pv,i,j,k)                   ; }
//      else if( pf[1]&&       pf[3]&&       pf[5]&&       pf[7]        ){ FLOOP pv[i][j][k] = AVG4_2(          pv,i,j,k)                   ; }
//      else if(               pf[3]&&pf[4]&&              pf[7]&&pf[8] ){ FLOOP pv[i][j][k] = 0.5*(AVG2_1_X1X2(pv,i,j,k)+AVG2_X1(pv,i,j,k)); }
//      else if(        pf[2]&&pf[3]&&              pf[6]&&pf[7]        ){ FLOOP pv[i][j][k] = 0.5*(AVG2_1_X1X2(pv,i,j,k)+AVG2_X2(pv,i,j,k)); }
//      else if( pf[1]&&              pf[4]&&pf[5]&&              pf[8] ){ FLOOP pv[i][j][k] = 0.5*(AVG2_2_X1X2(pv,i,j,k)+AVG2_X1(pv,i,j,k)); }
//      else if( pf[1]&&pf[2]&&              pf[5]&&pf[6]               ){ FLOOP pv[i][j][k] = 0.5*(AVG2_2_X1X2(pv,i,j,k)+AVG2_X2(pv,i,j,k)); }
//      else if(               pf[3]&&                     pf[7]        ){ FLOOP pv[i][j][k] = AVG2_1_X1X2(     pv,i,j,k)                   ; }
//      else if( pf[1]&&                     pf[5]                      ){ FLOOP pv[i][j][k] = AVG2_2_X1X2(     pv,i,j,k)                   ; }
//      else if(        pf[2]&&                     pf[6]               ){ FLOOP pv[i][j][k] = AVG2_X2(         pv,i,j,k)                   ; }
//      else if(                      pf[4]&&                     pf[8] ){ FLOOP pv[i][j][k] = AVG2_X1(         pv,i,j,k)                   ; }

// Old way:
      if(             pf[2]&&       pf[4]&&       pf[6]&&       pf[8] ){ FLOOP pv[i][j][k] = AVG4_1(          pv,i,j,k)                   ; }
      else if( pf[1]&&       pf[3]&&       pf[5]&&       pf[7]        ){ FLOOP pv[i][j][k] = AVG4_2(          pv,i,j,k)                   ; }
      else{ 
//	fflush(stderr);
//	fprintf(stderr,"fixup_utoprim()1: no good stencils1: i j pflag = %d %d : \n %4d %4d %4d \n %4d %4d %4d \n %4d %4d %4d \n\n", 
//		i,j,pflag[i-1][j+1],pflag[i][j+1],pflag[i+1][j+1],pflag[i-1][j],pflag[i][j],pflag[i+1][j],
//		pflag[i-1][j-1],pflag[i][j-1],pflag[i+1][j-1]);
//	fflush(stderr);
	failimage[4][i+j*N1]++ ;
	/* if nothing better to do, then leave densities and B-field unchanged, set v^i = 0 */
        for( k = RHO; k <= UU; k++ ) { pv[i][j][k] = 0.5*( AVG4_1(pv,i,j,k) + AVG4_2(pv,i,j,k) ); }  
	pv[i][j][U1] = pv[i][j][U2] = pv[i][j][U3] = 0.;
      }
      pflag[i][j] = 0;                /* The cell has been fixed so we can use it for interpolation elsewhere */
      fixup1zone( i, j, pv[i][j] ) ;  /* Floor and limit gamma the interpolated value */
    }
  }

  return;
}


#if( DO_FONT_FIX ) 
/***********************************************************************
   set_Katm():

       -- sets the EOS constant used for Font's fix. 

       -- see utoprim_1dfix1.c and utoprim_1dvsq2fix1.c  for more
           information. 

       -- uses the initial floor values of rho/u determined by fixup1zone()

       -- we assume here that Constant X1,r is independent of theta,X2

***********************************************************************/
void set_Katm( void )
{
  int i, j, k, G_type ;
  double prim[NPR], G_tmp;

  j = 1;

  G_type = get_G_ATM( &G_tmp );

  fflush(stdout);
  fprintf(stdout,"G_tmp = %26.20e \n", G_tmp );
  fflush(stdout);

  j = 0;
  for( i = 0 ; i < N1; i++ ) { 

    PLOOP prim[k] = 0.;
    prim[RHO] = prim[UU] = -1.;

    fixup1zone( i, j, prim );
    Katm[i] = (gam - 1.) * prim[UU] / pow( prim[RHO], G_tmp ) ;
    
    fflush(stdout);
    fprintf(stdout,"Katm[%d] = %26.20e \n", i, Katm[i] );
    fflush(stdout);
    
  }

  return;
}
#endif


#undef FLOOP 
