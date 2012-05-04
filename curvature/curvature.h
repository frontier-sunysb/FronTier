#if defined(__cplusplus)
#   define LINKAGE "C"
#else
#   define LINKAGE 
#endif


#define ct_max(e1,e2) (((e2) > (e1)) ? (e2): (e1))
#define ct_min(e1,e2) (((e2) < (e1)) ? (e2): (e1))
#define ct_set_max(e1,e2,e3)  (e1 = e2, (e1 < e3 ? e1 = e3 : e1))
#define ct_set_min(e1,e2,e3)  (e1 = e2, (e1 > e3 ? e1 = e3 : e1))


extern LINKAGE void  polyfit2d_lhf(
	double normal_out[2], 
	double *curvature_out, 
	double *ngbpts, 
	int ngbpts_dim1, 
	int mid, 
	int deg);


/* Prototype for top-level function. */
extern LINKAGE void  polyfit_lhf3d(
	double nrm_out[3], 
	int *deg_out, 
	double prcurvs_out[2], 
	double maxprdir_out[3], 
	double *xs, 
	int xs_dim1, 
	double *nrm_coor, 
	int nrm_coor_dim1, 
	int degree, 
	int strip);

/* Prototypes for primary functions. */
extern LINKAGE void  eval_curvature_lhf_surf(
	double curvs_out[2], 
	double dir_out[3], 
	double grad[2], 
	double H[2][2]);

extern LINKAGE int eval_vander_bivar(
	double *us, 
	int us_dim1, 
	double *bs, 
	int *bs_dim1, 
	int degree, 
	double *ws, 
	int ws_dim1, 
	int strip);
