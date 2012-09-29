
#include <FronTier.h>
#include <ctime>

enum _TS_TYPE {
        UNKNOWN_TS_TYPE = -1,
        WHITE_NOISE = 1
};
typedef enum _TS_TYPE TS_TYPE;

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
	UNIFORM
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

struct _UNIFORM_PARAMS {
	double a;
	double b;
};
typedef _UNIFORM_PARAMS UNIFORM_PARAMS;

struct _POWER_PARAMS {
	int power;
};
typedef _POWER_PARAMS POWER_PARAMS;

struct _STABLE_PARAMS {
	double alpha;
	double beta;
	double sigma;
	double mu;
};
typedef _STABLE_PARAMS STABLE_PARAMS;

struct _TIME_PARAMS {
	/* Random number generator part */
        RANDOM_TYPE rand_type;
	SEED_TYPE seed_type;
	POINTER pdf_params;
	double (*random_func)(POINTER,unsigned short int*);

	/* Stock Monte Carlo simulation part */
	TS_TYPE ts_type;	/* Type of time series */
	int T_density;		/* Incident density */
	unsigned short int seeds[3];
};
typedef struct _TIME_PARAMS TIME_PARAMS;

#define                 EPSILON         10e-14

extern double gauss_newton(POINTER,unsigned short int*);
extern double gauss_box_muller(POINTER,unsigned short int*);
extern double gauss_center_limit(POINTER,unsigned short int*);
extern double dist_cauchy(POINTER,unsigned short int*);
extern double dist_exponential(POINTER,unsigned short int*);
extern double dist_power(POINTER,unsigned short int*);
extern double dist_middle(POINTER,unsigned short int*);
extern double dist_uniform(POINTER,unsigned short int*);
extern double dist_stable(POINTER,unsigned short int*);

