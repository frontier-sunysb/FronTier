#include <cFluid.h>

extern double EosPressure(
	STATE *state)
{
	double 		dens = state->dens;
	double 		engy = state->engy;
	double 		*momn = state->momn;
	double 		ke,pres;
	int		i;
	int		dim = state->dim;
	EOS_PARAMS	*eos = state->eos;
	double		gamma = eos->gamma;

	if (dens <= 0.0)
	    return 0.0;
	ke = 0.0;
	for (i = 0; i < dim; ++i)
	    ke += sqr(momn[i]);
	ke *= 0.5/dens;
	pres = (gamma - 1.0)*(engy - ke + dens*eos->einf) - gamma*eos->pinf;
	
	return pres;
}	/* end EosPressure */

extern double EosSoundSpeedSqr(
	STATE *state)
{
	double		pres = state->pres;
	double		dens = state->dens;
	EOS_PARAMS	*eos = state->eos;
	
	return eos->gamma*(pres + eos->pinf)/dens;
}

extern double EosSoundSpeed(
	STATE *state)
{
	return sqrt(EosSoundSpeedSqr(state));
}

extern double EosInternalEnergy(
	STATE *state)
{
	double		pres = state->pres;
	double		dens = state->dens;
	EOS_PARAMS	*eos = state->eos;
	double		gamma = eos->gamma;

	return (pres+gamma*eos->pinf)/(gamma-1) - dens*eos->einf;
}

extern double EosEnergy(
	STATE *state)
{
	int	i,dim = state->dim;
	double	dens = state->dens;
	double	*momn = state->momn;
	double	e;
	
	e = 0.0;
	for (i = 0; i < dim; ++i)
	    e += 0.5*sqr(momn[i])/dens;
	e += EosInternalEnergy(state);

	return e;
}

extern double EosMaxBehindShockPres(
	double M2,		// Mach number sqaure
	STATE *state)
{
	double im2;		// acoustic impedance squared
	double gamma = state->eos->gamma;
	double dens = state->dens;
	double pres = state->pres;
	double c4,c5,p1;

	im2 = M2*gamma*pres*dens;
	c4 = (gamma -1.0)/(gamma + 1.0);
	c5 = 2.0/(gamma + 1.0);
	p1 = c5*im2/dens - c4*pres;
	return p1;
}	/* end EosMaxBehindShockPres */

extern void CovertVstToState(
	STATE		*state,
	SWEEP		*vst,
	EOS_PARAMS	*eos,
	int		ind,
	int		dim)
{
	int	i;

	state->dim = dim;
	state->eos = eos;
	state->dens = vst->dens[ind];
	state->engy = vst->engy[ind];
	for (i = 0; i < dim; ++i)
	    state->momn[i] = vst->momn[i][ind];
	state->pres = EosPressure(state);
}

extern void EosSetTVDParams(
	SCHEME_PARAMS	*scheme_params,
	EOS_PARAMS	*eos)
{
	scheme_params->gamma = eos->gamma;
	scheme_params->einf = eos->einf;
}

