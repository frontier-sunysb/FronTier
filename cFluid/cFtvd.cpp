/***************************************************************
FronTier is a set of libraries that implements different types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions have 
discontinuities.  

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
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
****************************************************************/

#include <cFluid.h>

static 	double fncomp(double *dd, double s1, double s2, double beta);

#define EPS             0.000001 /* TOLERANCE */
                                 /* Governs convergence of iterations in */
                                 /* Riemann problem and is used to test if */
                                 /* certain fabs(a - b) < EPS */
#define sign(x)         ((x) > 0.0 ? 1 : -1)
#define fndq(d1,d2,d0) (fabs((d0)) > 0.000001 ? ((d2) - (d1)/(d0))/2.0 : (d2))
#define fnga(d1,d2,d0) (fabs((d0)) > 0.000001 ? ((d2) - (d1))/(d0) : 0.0)
#define fnq(a0) (fabs((a0)) > EPS ? fabs((a0)) : 0.5*((a0)*(a0) + EPS*EPS)/EPS)
#define fnl(r1,r2)      (sign((r1))==sign((r2)) ?                       \
                        sign((r1))*std::min(fabs((r1)),fabs((r2))) : 0.0)


void TVD_flux(
	POINTER params, 
	SWEEP *vst, 
	FSWEEP *vflux,
	int n)
{
      double lambda,beta;
      double t1;
      double aa,bb;
      double ga1,ga3,ga4;
      int i,j,k;
      static int max_n = 0;
      static double *dq;
      static double *s,*ga,*ff,**r;
      static double **f,**q,**h,**a,**z,**g,**d,**comp;
      static double *um,*hm,*cm,*vm,*wm,*km,*tm,*tt;
      static double *p;
      SCHEME_PARAMS *scheme_params = (SCHEME_PARAMS*)params;
      double gamma, einf;

      gamma = scheme_params->gamma;
      einf = scheme_params->einf;
      
      if (max_n < n+6)
      {
	  if (max_n != 0)
	     FT_FreeThese(15,um,hm,cm,vm,wm,km,tm,tt,p,h,d,a,z,g,comp);
	  else
	  {
	     FT_VectorMemoryAlloc((POINTER*)&q,5,sizeof(double *));
	     FT_VectorMemoryAlloc((POINTER*)&f,5,sizeof(double *));
	     FT_VectorMemoryAlloc((POINTER*)&dq,5,sizeof(double));
	     FT_VectorMemoryAlloc((POINTER*)&s,5,sizeof(double));
	     FT_VectorMemoryAlloc((POINTER*)&ff,5,sizeof(double));
	     FT_VectorMemoryAlloc((POINTER*)&ga,5,sizeof(double));
	     FT_MatrixMemoryAlloc((POINTER*)&r,5,5,sizeof(double));
	  }
	  max_n = n+6;
	  FT_VectorMemoryAlloc((POINTER*)&um,max_n,sizeof(double));
	  FT_VectorMemoryAlloc((POINTER*)&hm,max_n,sizeof(double));
	  FT_VectorMemoryAlloc((POINTER*)&cm,max_n,sizeof(double));
	  FT_VectorMemoryAlloc((POINTER*)&vm,max_n,sizeof(double));
	  FT_VectorMemoryAlloc((POINTER*)&wm,max_n,sizeof(double));
	  FT_VectorMemoryAlloc((POINTER*)&km,max_n,sizeof(double));
	  FT_VectorMemoryAlloc((POINTER*)&tm,max_n,sizeof(double));
	  FT_VectorMemoryAlloc((POINTER*)&tt,max_n,sizeof(double));
	  FT_VectorMemoryAlloc((POINTER*)&p,max_n,sizeof(double));
	  FT_MatrixMemoryAlloc((POINTER*)&h,5,max_n,sizeof(double));
	  FT_MatrixMemoryAlloc((POINTER*)&d,5,max_n,sizeof(double));
	  FT_MatrixMemoryAlloc((POINTER*)&a,5,max_n,sizeof(double));
	  FT_MatrixMemoryAlloc((POINTER*)&z,5,max_n,sizeof(double));
	  FT_MatrixMemoryAlloc((POINTER*)&g,5,max_n,sizeof(double));
          FT_MatrixMemoryAlloc((POINTER*)&comp,5,max_n,sizeof(double));
      }

      lambda = scheme_params->lambda;
      beta = scheme_params->beta;

      q[0] = vst->dens;
      q[1] = vst->momn[0];
      q[2] = vst->engy;
      q[3] = vst->momn[1];
      q[4] = vst->momn[2];

      f[0] = vflux->dens_flux;
      f[1] = vflux->momn_flux[0];
      f[2] = vflux->engy_flux;
      f[3] = vflux->momn_flux[1];
      f[4] = vflux->momn_flux[2];

      for (i = 0; i < n+6; i++)
      {
      	  p[i] = vst->pres[i];
	  t1 = q[1][i]/q[0][i];
          f[0][i] = q[1][i];
          f[1][i] = q[1][i]*t1 + p[i];
          f[2][i] = (p[i] + q[2][i])*t1;
          f[3][i] = q[3][i]*t1;
          f[4][i] = q[4][i]*t1;
	  tm[i] = sqrt(q[0][i]);
	  tt[i] = tm[i]/q[0][i];
      }

      for (i = 0; i < n+5; i++)
      {
          t1 = tm[i] + tm[i+1];
          um[i] = (tt[i]*q[1][i] + tt[i+1]*q[1][i+1])/t1;
          vm[i] = (tt[i]*q[3][i] + tt[i+1]*q[3][i+1])/t1;
          wm[i] = (tt[i]*q[4][i] + tt[i+1]*q[4][i+1])/t1;
	  km[i] = 0.5*(um[i]*um[i] + vm[i]*vm[i] + wm[i]*wm[i]);
          hm[i] = ((q[2][i] + p[i])/tm[i] + (q[2][i+1] + p[i+1])/tm[i+1])/t1;
	  cm[i] = sqrt((gamma-1.0)*(hm[i] - km[i] + einf));
      }

      for (i = 0; i < n+5; i++)
      {
	  for (j = 0; j < 5; j++)
	  {
	      dq[j] = q[j][i + 1] - q[j][i];
	  }
          d[4][i] = dq[0]*wm[i]-dq[4];
          d[3][i] = dq[0]*vm[i]-dq[3];
	  aa = (gamma-1.0)*((km[i] + einf)*dq[0] + dq[2] - um[i]*dq[1] -
                         vm[i]*dq[3]-wm[i]*dq[4])/cm[i]/cm[i];
          bb = (dq[1] - um[i]*dq[0])/cm[i];
          d[2][i] = (aa + bb)/2.0;
          d[1][i] = dq[0] - aa;
          d[0][i] = (aa - bb)/2.0;
	  tt[i] = (fabs(um[i])-lambda*um[i]*um[i])/2.0;
      }

      for (i = 2; i < n+4; i++){
          comp[1][i] = fncomp(&(d[1][i-2]),tt[i-1],tt[i],beta);
          comp[3][i] = fncomp(&(d[3][i-2]),tt[i-1],tt[i],beta);
          comp[4][i] = fncomp(&(d[4][i-2]),tt[i-1],tt[i],beta);
      }
      for (j = 0; j < 5; j++){
          comp[j][0] = 0.0;
          comp[j][1] = 0.0;
          comp[j][n+5] = 0.0;
          comp[j][n+4] = 0.0;
      }
	
      for (i = 0; i < n+5; i++)
      {
          ga1 = fnga(comp[1][i],comp[1][i+1],d[1][i]);
          ga3 = fnga(comp[3][i],comp[3][i+1],d[3][i]);
          ga4 = fnga(comp[4][i],comp[4][i+1],d[4][i]);
          a[0][i] = um[i] - cm[i];
          a[1][i] = um[i] + ga1;
          a[2][i] = um[i] + cm[i];
          a[3][i] = um[i] + ga3;
          a[4][i] = um[i] + ga4;
	  for (j = 0; j < 5; j++)
	  {
              z[j][i] = (fnq(a[j][i]) - lambda*a[j][i]*a[j][i])*d[j][i]/2.0;
	  }
      }

      for (i = 2; i < n+4; i++)
      {
	  for (j = 0; j < 5; j++)
	  {
              g[j][i] = fnl(z[j][i-1],z[j][i]);
	  }
      }

      for (i = 2; i < n+4; i++)
      {
	  for (j = 0; j < 5; j++)
	  {
              ga[j] = fnga(g[j][i],g[j][i+1],d[j][i]);
               s[j] = g[j][i] + g[j][i + 1] - fnq(a[j][i]+ga[j])*d[j][i];
              ff[j] = f[j][i] + f[j][i + 1];
	  }
          s[1] += comp[1][i]+comp[1][i+1];
          s[3] += comp[3][i]+comp[3][i+1];
          s[4] += comp[4][i]+comp[4][i+1]; 
          
          r[0][0] = 1;
          r[0][1] = 1;
          r[0][2] = 1;
          r[0][3] = 0;
          r[0][4] = 0;  

	  r[1][0] = um[i]-cm[i];
          r[1][1] = um[i];
	  r[1][2] = um[i]+cm[i];
          r[1][3] = 0; 
          r[1][4] = 0;

	  r[2][0] = hm[i] - um[i]*cm[i];
	  r[2][1] = km[i] - einf;
	  r[2][2] = hm[i] + um[i]*cm[i];
	  r[2][3] = -vm[i];
	  r[2][4] = -wm[i];

          r[3][0] = vm[i];
          r[3][1] = vm[i];
          r[3][2] = vm[i];
	  r[3][3] = -1;   
          r[3][4] = 0;

          r[4][0] = wm[i];
          r[4][1] = wm[i];
          r[4][2] = wm[i];
          r[4][3] = 0;
	  r[4][4] = -1;

	  for (j = 0; j < 5; j++)
	  {
              h[j][i + 1] = ff[j];
	      for (k = 0; k < 5; k++)
	      {
		  h[j][i + 1] += s[k]*r[j][k];
	      }
	      h[j][i + 1] *= 0.5;
	  }
      }
      for (i = 3; i < n+3; i++)
      {
	  for (j = 0; j < 5; j++)
	  {
	      f[j][i] = - lambda*(h[j][i + 1] - h[j][i]);
	      if (isnan(f[j][i]))
	      {
		(void) printf("In TVD_flux(): f[%d][%d] = %f\n",j,i,f[j][i]);
		for (k = 0; k < n+6; ++k)
		    printf("q[%d] = %f %f %f %f %f %f\n",k,q[0][k],q[1][k],
				q[2][k],q[3][k],q[4][k],p[k]);
		clean_up(ERROR);
	      }
	  }
      }
      return;
}	/* end TVD_flux */

static 	double fncomp(double *dd, double s1, double s2, double beta)
{
	double eps,s;
	double fcp;
	if (beta == 0.0) return 0.0;
	s1 = s1*(dd[1] - fnl(fnl(dd[0],dd[1]), dd[2]));
	s2 = s2*(dd[2] - fnl(fnl(dd[1],dd[2]), dd[3]));
	eps = 2*fabs(pow(fabs(dd[2]),beta) - pow(fabs(dd[1]),beta))/
		fabs(pow(fabs(dd[1]),beta) + pow(fabs(dd[2]),beta) + 0.00001);
	s = sign(s1);
	fcp = s*std::max(s*fnl(eps*s2, s1),s*fnl(s2,eps*s1));
	return fcp;
}	/* end fncomp */
