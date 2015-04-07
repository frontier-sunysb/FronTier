/***************************************************************
FronTier is a set of libraries that implements different types of 
Front Traking algorithms. Front Tracking is a numerical method 
for the solution of partial differential equations whose solutions 
have discontinuities.  


Copyright (C) 1999 by The University at Stony Brook. 
 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

****************************************************************/
/***************************************************************
 *		clmprop.cpp
 * ***********************************************************/

#include <iFluid.h>
#include "climate.h"
/****************************MACRO***************************/
/*Macroscopic equations*/
/****************************MACRO***************************/
MACRO::~MACRO()
{

}

/********construct*************/
MACRO::MACRO(Front &front):T(params->T0),Q(params->qv0),front(&front),params((PARAMS*)front.extra2)
{
	int i;
	IF_PARAMS *iFparams = (IF_PARAMS*)front.extra1;
	P = iFparams->ref_pres;
	Cd = 0.0;
	Wm = -1.0;
	Gamma = 0;
	S = 0.0;
	rho0 = iFparams->rho2;
	dim = front.f_basic->dim;
}

/********computeCondRate****************/
void MACRO::computeCondRate(void)
{
	PARTICLE* particle_array = params->particle_array;
	RECT_GRID gr = front->pp_grid->Global_grid;
	double *L = gr.L;
	double *U = gr.U;
	double a3 = 1.0;
	double K = params->K;
	int i,j;
	int num_drops = params->num_drops;
	Cd = 0;
	for (i = 0; i < dim; i++)
	   a3 *= U[i] - L[i];
	for (i = 0; i < num_drops; i++)
	   Cd += particle_array[i].radius*K*S*particle_array[i].rho;
#if defined(__MPI__)
	pp_gsync();
	pp_global_sum(&Cd,1);
#endif
	Cd *= 4.0*PI/(a3*rho0);
}

/*compute macroscopic temperature*/ 
void MACRO::computeTemp(void)
{
	double Lh = params->Lh;
	double Cp = params->Cp;
	double Source = -Wm*Gamma+Lh/Cp*Cd;
	double dt = front->dt;
	/*First order method*/
	T += dt*Source;
}

/*compute macroscopic vapor mixing ratio*/
void MACRO::computeVapor(void)
{
	double dt = front->dt;
	Q -= dt*Cd*1000.0;
}

/*compute macroscopic supersaturation*/
void MACRO::computeSupersat(void)
{
	double sat_vap_pre, sat_vap_rat;
        double Lh, Rv, Rd, rhoL, Kc, es, D, Cp, ksi;
        double A;
        Lh = params->Lh;
        Rv = params->Rv;
        Rd = params->Rd;
        rhoL = params->rho_l;
        Kc = params->Kc;
        D = params->D;
        ksi = Rd/Rv;
	
	es = 611.2*exp(17.67*(T-273.15)/(T-29.65));
        params->K = 1/((Lh/(Rv*T)-1)*Lh*rhoL/(Kc*T)+rhoL*Rv*T/(D*es));
        printf("Condensation coeffecient = %e\n",params->K);
	sat_vap_pre = 611.2*exp(17.67*(T-273.15)/(T-29.65));
        sat_vap_rat = 621.97 * sat_vap_pre/P;
	S = Q/sat_vap_rat - 1;
	printf("sat_vap_rat = %f\n",sat_vap_rat);
}

/*compute macroscopic Pressure*/
void MACRO::computePressure(void)
{
	double dt = front->dt;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	double *gravity = iFparams->gravity;
	P -= rho0*gravity[dim-1]*Wm*dt;
}

/*solve function*/
void MACRO::solve(void)
{
	computePressure();
	computeSupersat();
	computeCondRate();
	computeTemp();
	computeVapor();
}

/*output function*/
void MACRO::output(void)
{
	static boolean first = YES;
	char fname[100];
	FILE *outfile;
	char *outname = front->out_name;
	if (pp_mynode() != 0)
	    return;
	sprintf(fname,"%s/output",outname);
	if (first == YES)
	{
	    first = NO;
	    outfile = fopen(fname,"w");
	    fprintf(outfile,"Temp    Vapor    Pressure    Supersat    CondRate\n");
	    fclose(outfile);
	}
	outfile = fopen(fname,"a");
	fprintf(outfile,"%f    %f    %f    %f    %f\n",T,Q,P,S,Cd);
	fclose(outfile);
}

void MACRO::initMovieVariables(void)
{
	return;
}

/*API function*/
double MACRO::getSupersat(void){return S;}
double MACRO::getTemp(void){return T;}
double MACRO::getPressure(void){return P;}
double MACRO::getVapor(void){return Q;}
void   MACRO::setW(double new_W){Wm = new_W; return;}


