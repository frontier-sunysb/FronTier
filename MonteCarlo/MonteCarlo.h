
#include <FronTier.h>
#include <ctime>

enum _RANDOM_TYPE {
        ERROR_TYPE = -1,
        GAUSS_NEWTON = 1,
        GAUSS_BM,
        GAUSS_CL,
        EXPONENTIAL,
	POWER,
	MIDDLE,
        CAUCHY
};
typedef enum _RANDOM_TYPE RANDOM_TYPE;

enum _SEED_TYPE {
	FIXED_SEED = 1,
	INPUT_SEED,
	RANDOM_SEED
};
typedef enum _SEED_TYPE SEED_TYPE;

struct _GAUSS_PARAMS {
	double sigma;           /* Deviation */
        double mu;              /* Expectation */
};
typedef _GAUSS_PARAMS GAUSS_PARAMS;

struct _EXP_PARAMS {
	double lambda;
};
typedef _EXP_PARAMS EXP_PARAMS;

struct _POWER_PARAMS {
	int power;
};
typedef _POWER_PARAMS POWER_PARAMS;

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
	double S0;		/* Start value */
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

#define                 EPSILON         10e-14

extern double gauss_newton(POINTER,unsigned short int*);
extern double gauss_box_muller(POINTER,unsigned short int*);
extern double gauss_center_limit(POINTER,unsigned short int*);
extern double dist_cauchy(POINTER,unsigned short int*);
extern double dist_exponential(POINTER,unsigned short int*);
extern double dist_power(POINTER,unsigned short int*);
extern double dist_middle(POINTER,unsigned short int*);

