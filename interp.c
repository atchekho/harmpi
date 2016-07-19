//Modified by Alexander Tchekhovskoy: MPI+3D

/* M1: no changes needed */

/***

includes all interpolation routines 
linear (MC or other limiter)
parabolic (from collela and woodward)
weno

***/

#include "decs.h"

/* performs the slope-limiting for the numerical flux calculation */

double slope_lim(double y1, double y2, double y3)
{
	double Dqm, Dqp, Dqc, s;

	/* woodward, or monotonized central, slope limiter */
	if (lim == MC) {
		Dqm = 2. * (y2 - y1);
		Dqp = 2. * (y3 - y2);
		Dqc = 0.5 * (y3 - y1);
		s = Dqm * Dqp;
		if (s <= 0.)
			return 0.;
		else {
			if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
				return (Dqm);
			else if (fabs(Dqp) < fabs(Dqc))
				return (Dqp);
			else
				return (Dqc);
		}
	}
	/* van leer slope limiter */
	else if (lim == VANL) {
		Dqm = (y2 - y1);
		Dqp = (y3 - y2);
		s = Dqm * Dqp;
		if (s <= 0.)
			return 0.;
		else
			return (2. * s / (Dqm + Dqp));
	}

	/* minmod slope limiter (crude but robust) */
	else if (lim == MINM) {
		Dqm = (y2 - y1);
		Dqp = (y3 - y2);
		s = Dqm * Dqp;
		if (s <= 0.)
			return 0.;
		else if (fabs(Dqm) < fabs(Dqp))
			return Dqm;
		else
			return Dqp;
	}

	fprintf(stderr, "unknown slope limiter\n");
	exit(10);

	return (0.);
}

void linear_mc(double x1, double x2, double x3, double *lout, double *rout) 
{
	double Dqm,Dqp,Dqc,s;

	Dqm = 2. * (x2 - x1);
	Dqp = 2. * (x3 - x2);
	Dqc = 0.5 * (x3 - x1);

	s = Dqm * Dqp;

	if (s <= 0.)
		s = 0.;
	else {
		if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
			s = Dqm;
		else if (fabs(Dqp) < fabs(Dqc))
			s = Dqp;
		else
			s = Dqc;
	}

	/* reconstruct left, right */
	*lout = x2 - 0.5*s;
	*rout = x2 + 0.5*s;
}

