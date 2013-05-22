
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

