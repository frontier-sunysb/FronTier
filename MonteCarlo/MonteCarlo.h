#ifndef _MONTECARLO_
#define _MONTECARLO_

#include <FronTier.h>
#include <ctime>
#if defined(__GSL__)
#include <gsl/gsl_sf_bessel.h>
#endif /* if defined(__GSL__) */

enum _RANDOM_TYPE {
        ERROR_TYPE = -1,
        GAUSS_NEWTON = 1,
        GAUSS_BM,
        GAUSS_CL,
        EXPONENTIAL,
	POWER,
	MIDDLE,
        CAUCHY,
        STABLE,
	UNIFORM,
	GIG,
	GH
};
typedef enum _RANDOM_TYPE RANDOM_TYPE;

enum _SEED_TYPE {
	FIXED_SEED = 1,
	INPUT_SEED,
	RANDOM_SEED
};
typedef enum _SEED_TYPE SEED_TYPE;

struct _PARAMS {
	/* Random number generator part */
        RANDOM_TYPE rand_type;
	SEED_TYPE seed_type;
	POINTER pdf_params;
	int num_samples;	/* For plotting distribution function */
	double XL;		/* Plotting left end boundary */
	double XU;		/* Plotting right end boundary */
	int nx;			/* Number of intervals */

	/* Stock Monte Carlo simulation part */
	boolean do_monte_carlo;	/* flag to do Monte Carlo Simulation */
	boolean do_option_price;/* flag to do Monte Carlo Simulation */
	double S0;		/* Start value */
	double E;		/* Option strike price */
	double T;		/* End time */
	double mu;		/* Growth rate */
	double sigma;		/* Volatility */
	boolean print_detail;	/* Flag to output individual cases */
	int num_print_sims;	/* Number cases for output */
	int num_steps;		/* Number of time step */
	int num_sims;		/* Number of simulations */
	unsigned short int seeds[3];
};
typedef struct _PARAMS PARAMS;

#endif
