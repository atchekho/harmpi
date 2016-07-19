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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SMALL 	1.e-14 ;

int oN1,oN2,nN1,nN2 ;
double dx1,dx2,rin,rout,dx,dy,hslope ;

int main(int argc, char *argv[])
{
	char **oldimage,**newimage,**cmatrix(int a, int b, int c, int d)  ;
	int i,j,iold,jold,ioldp,joldp ;
	double xmax,ymax,dx2,x1max ;
	double r,th,diold,djold,x1,x2 ;
	void new_coord(int i,int j, double *r, double *th) ;
	void old_coord(double r,double th,double *x1,double *x2) ;


	/* get size of image from arguments */
	if(argc < 8) {
		fprintf(stderr,"Usage: \n" ) ;
	        fprintf(stderr,"\t image_interp oN1 oN2 nN1 nN2 Rin Rout xmax ymax hslope\n") ; 
	        fprintf(stderr,"  where \n");
		fprintf(stderr,"\t oN1     = N1 from simulation \n");
		fprintf(stderr,"\t oN2     = N2 from simulation \n");
		fprintf(stderr,"\t nN1     = # of pixels in x-dir \n");
		fprintf(stderr,"\t nN2     = # of pixels in x-dir \n"); 
		fprintf(stderr,"\t Rin     = Rin used in simulation \n");
		fprintf(stderr,"\t Rout    = Rout used in simulation\n");
		fprintf(stderr,"\t xmax    = interpolated image will span x=[0,xmax] \n");
		fprintf(stderr,"\t ymax    = interpolated image will span y=[-ymax,ymax] \n");
		fprintf(stderr,"\t hslope  = hslope value from run\n");
		fprintf(stderr,"\n\n Example: \n");
		fprintf(stderr," cat im_lrho_0000.r8|awk '(FNR>4) {print}'|image_interp 128 64 128 256 1.321025567223338859 40. 40. 40. 0.3 > im_int.r8 \n");
		fprintf(stderr,"\n\n ");
		exit(0) ;
	}
	sscanf(argv[1],"%d",&oN1) ;
	sscanf(argv[2],"%d",&oN2) ;
	sscanf(argv[3],"%d",&nN1) ;
	sscanf(argv[4],"%d",&nN2) ;
	sscanf(argv[5],"%lf",&rin) ;
	sscanf(argv[6],"%lf",&rout) ;
	sscanf(argv[7],"%lf",&xmax) ;
	sscanf(argv[8],"%lf",&ymax) ;
	sscanf(argv[9],"%lf",&hslope) ;

	x1max = log(rout/rin) ;
	dx1 = x1max/oN1 ;
	dx2 = 1./(double)oN2 ;

	dx = xmax/(double)nN1 ;
	dy = 2.*ymax/(double)nN2 ;

	/* make arrays for images */
	oldimage = cmatrix(0,oN1-1,0,oN2-1) ;
	newimage = cmatrix(0,nN1-1,0,nN2-1) ;

	fprintf(stdout,"RAW\n#   \n%d %d\n%d\n",nN1,nN2,255);

	/* read in old image */
	for(j=oN2-1;j>=0;j--)
	for(i=0;i<oN1;i++) {
		fread(&oldimage[i][j], sizeof(unsigned char), 1, stdin) ;
		/*
		fprintf(stderr,"%d %d %u\n",i,j,oldimage[i][j]) ;
		*/
	}

	/* interpolate to new image */
	for(j=nN2-1;j>=0;j--)
	for(i=0;i<nN1;i++) {
		new_coord(i,j,&r,&th) ;
		old_coord(r,th,&x1,&x2) ;
		/*
		fprintf(stderr,"%d %d %g %g %g %g\n",i,j,r,th,x1,x2) ;
		*/
		if(x1 < 0. || x1 >= x1max || 
		   x2 < 0. || x2 >= 1.) 
			newimage[i][j] = 0 ;

		else {

#if 0
			iold = (int)(x1/dx1 - 0.5) ;
			diold = x1/dx1 - 0.5 - (int)(x1/dx1 - 0.5) ;
			jold = (int)(x2/dx2 - 0.5) ;
			djold = x2/dx2 - 0.5 - (int)(x2/dx2 - 0.5) ;
			ioldp = iold+1 ;
			joldp = jold+1 ;


			/* take care of boundary effects */
			if( diold < 0. ) ioldp = 0 ;
			else if (iold == oN1-1) ioldp = oN1-1 ;
			if( djold < 0. ) joldp = 0 ;
			else if( jold == oN2-1 ) joldp = oN2-1 ;

			/*
			fprintf(stderr,"iold, jold: %d %d %d %d %g %g\n",
					iold,jold,ioldp,joldp,diold,djold) ;
			fprintf(stderr,"old: %u\n",oldimage[iold][jold]) ;
			*/

			newimage[i][j] = (char)(0.5 + 
				(1. - diold)*(1.-djold)*oldimage[iold][jold] +
				(1. - diold)*djold*oldimage[iold][joldp] +
				diold*(1.-djold)*oldimage[ioldp][jold] +
				diold*djold*oldimage[ioldp][joldp]) ;
#endif

			iold = (int)(x1/dx1 - 1.e-20) ;
			jold = (int)(x2/dx2 - 1.e-20) ;

			newimage[i][j] = oldimage[iold][jold] ;
			/*
			fprintf(stderr,"newim:%d %d %d %d %u\n",i,j,iold,jold,newimage[i][j]) ;
			*/
		}
		fwrite(&newimage[i][j], sizeof(unsigned char), 1, stdout) ;
	}
        return(0);
}

void new_coord(int i,int j,double *r,double *th)
{
	double x,y ;

	x = i*dx ;
	y = (j-nN2/2)*dy ;

	*r = sqrt(x*x + y*y) ;
	*th = atan2(x,y) ;	/* deliberately interchanged args */

	/*
	fprintf(stderr,"x,y: %g %g %g %g\n",x,y,*r,*th) ;
	*/

	return ;
}

void old_coord(double r, double th, double *x1, double *x2)
{
	double root_find(double th) ;

	*x1 = log(r/rin) ;
	*x2 = root_find(th) ;

	return ;
}

double root_find(double th)
{
	int i ;
	double X2a,X2b,X2c,tha,thb,thc ;
	double dthdX2, dtheta_func(double y),theta_func(double y) ;

	if(th < M_PI/2.) {
		X2a = 0.-SMALL ;
		X2b = 0.5+SMALL ;
	}
	else {
		X2a = 0.5-SMALL ;
		X2b = 1.+SMALL ;
	}
	tha = theta_func(X2a) ;
	thb = theta_func(X2b) ;

	/* bisect for a bit */
	for(i = 0 ; i < 8 ; i++) {
		X2c = 0.5*(X2a+X2b) ;
		thc = theta_func(X2c) ;

		if((thc - th)*(thb - th) < 0.) X2a = X2c ;
		else X2b = X2c ;

	}

	/* now do a couple of newton-raphson strokes */
	tha = theta_func(X2a) ;
	for(i = 0 ; i < 2 ; i++) {
		dthdX2 = dtheta_func(X2a) ;
		X2a -= (tha - th)/dthdX2 ;
		tha = theta_func(X2a) ;
	}

	return(X2a) ;
}

double theta_func(double x)
{
	return(M_PI*x + 0.5*(1.-hslope)*sin(2.*M_PI*x)) ;
}
double dtheta_func(double x)
{
	return(M_PI*(1. + (1.-hslope)*cos(2.*M_PI*x))) ;
}



char **cmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
        char **m;
	int i ;

        m=(char **)malloc((unsigned) (nrh-nrl+1)*sizeof(char *));
        if (!m) exit(20) ;
        m -= nrl;

        for(i=nrl;i<=nrh;i++) {
                m[i]=(char *)malloc((unsigned) (nch-ncl+1)*sizeof(char));
                if (!m[i]) exit(20) ;
                m[i] -= ncl;
        }
        return m;
}

