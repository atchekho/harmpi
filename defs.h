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

/*************************************************************************
    GLOBAL ARRAYS SECTION 
*************************************************************************/
double   a_p[N1+4][N2+4][NPR] ;  /* space for primitive vars */
double  a_dq[N1+4][N2+4][NPR] ;  /* slopes */
double  a_F1[N1+4][N2+4][NPR] ;  /* fluxes */
double  a_F2[N1+4][N2+4][NPR] ;  /* fluxes */
double  a_ph[N1+4][N2+4][NPR] ;  /* half-step primitives */
int     a_pflag[N1+4][N2+4];	/* identifies failure points */

/* for debug */
double psave[N1][N2][NPR] ;
double fimage[NIMG][N1*N2];
int    failimage[5][N1*N2];

/* grid functions */
double a_conn[N1+4][N2+4][NDIM][NDIM][NDIM] ;
double a_gcon[N1+4][N2+4][NPG][NDIM][NDIM] ;
double a_gcov[N1+4][N2+4][NPG][NDIM][NDIM] ;
double a_gdet[N1+4][N2+4][NPG] ;

double (*   p)[N2+4][NPR] ;
double (*  dq)[N2+4][NPR] ;
double (*  F1)[N2+4][NPR] ;
double (*  F2)[N2+4][NPR] ;
double (*  ph)[N2+4][NPR] ;
int    (*  pflag)[N2+4];
double (* conn)[N2+4][NDIM][NDIM][NDIM] ;
double (* gcon)[N2+4][NPG][NDIM][NDIM] ;
double (* gcov)[N2+4][NPG][NDIM][NDIM] ;
double (* gdet)[N2+4][NPG] ;


#if(DO_FONT_FIX)
double Katm[N1];
#endif

/*************************************************************************
    GLOBAL VARIABLES SECTION 
*************************************************************************/
/* physics parameters */
double a ;
double gam ;

/* numerical parameters */
double Rin,Rout,hslope,R0,fractheta ;
double cour ;
double dV,dx[NPR],startx[NPR] ;
double dt ;
double t,tf ;
double rcurr,hcurr ;
int istart,istop,jstart,jstop ;
int icurr,jcurr,pcurr,ihere,jhere,phere ;
double dminarg1,dminarg2 ;
int nstep ;
double fval1,fval2;

/* output parameters */
double DTd ;
double DTl ;
double DTi ;
int    DTr ;
int    dump_cnt ;
int    image_cnt ;
int    rdump_cnt ;
int    nstroke ;

/* global flags */
int    failed ;
int    lim ;
double defcon ;

/* diagnostics */
double mdot = 0. ;
double edot = 0. ;
double ldot = 0. ;

/* current local position */
int icurr,jcurr,pcurr ;
