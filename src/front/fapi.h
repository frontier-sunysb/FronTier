/*! \file fapi.h
    
    \brief The fapi.h contains the functions used to operate the interface.
 */

/*! \defgroup INITIALIZATION    FronTier Initialization Functions
 *  \defgroup TIME              FronTier Time Control Functions
 *  \defgroup OUTPUT            FronTier Output Functions
 *  \defgroup PROPAGATION       FronTier Propagation Functions
 *  \defgroup GRIDINTFC         FronTier Interface-Grid Functions
 *  \defgroup GEOMETRY          FronTier Interface Geometry Functions
 *  \defgroup OPTIMIZATION      FronTier Interface Optimization Functions
 *  \defgroup BOUNDARY          FronTier Boundary Setting Functions
 *  \defgroup QUERY          	FronTier Interface Query Functions
 *  \defgroup CHECK          	FronTier Interface Check Functions
 *  \defgroup PARALLEL          FronTier Parallel Communication Functions
 *  \defgroup FIELD          	FronTier Field (State) Functions
 *  \defgroup MEMORY          	FronTier Memory Management Functions
 *  \defgroup CREATE          	FronTier Create Hypersurface
 *  \defgroup INFO          	FronTier Print Information
 **/

#include <front/fdecs.h>

                /* Front IMPORTED Function Prototypes*/

#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

/*! \fn void FT_Init(int argc, char **argv, F_BASIC_DATA *f_basic)
    \ingroup INITIALIZATION
    \brief Process the command line arguments, set IO handles, and store
     initialization parameters in f_basic, Including read the number of 
     processors, the dimesension of the problem, partition of the 
     computational domain, restart message, the input file and the output 
     directory.
     -d dim
     -p nx [ny] [nz]
     -i input-name
     -o output-dir-name
     -r restart-dir-name
     -t restart-step
    \param argc @b in	The number of arguments passed by command line
    \param argv @b in	The argument vector passed by command line
    \param f_basic @b out	Structure to store options for initializing the program
 */
   IMPORT  void FT_Init(int argc,
			char **argv, 
			F_BASIC_DATA *f_basic );

/*! \fn void FT_ReadSpaceDomain(char *in_name, F_BASIC_DATA *f_basic)
 *  \ingroup INITIALIZATION
    \brief Read from input file the computational domain information of the 
     problem, including the domain limits, computational grid, and the types 
     of the rectangular boundaries. The information is stored in the structure 
     f_basic.
    \param in_name @b in	The name of input file
    \param f_basic @b out	Structure to store domain information
 */
   IMPORT  void FT_ReadSpaceDomain(char *in_name, 
				    F_BASIC_DATA *f_basic);

/*! \fn void FT_ReadComparisonDomain(char *in_name, F_BASIC_DATA *f_basic)
 *  \ingroup INITIALIZATION
    \brief Read from input file the comparison domain information of the 
     problem, including the domain limits, computational grid, and the types 
     of the rectangular boundaries. The information is stored in the structure 
     f_basic.
    \param in_name @b in	The name of input file
    \param f_basic @b out	Structure to store domain information
 */
   IMPORT  void FT_ReadComparisonDomain(char *in_name, 
				    F_BASIC_DATA *f_basic);

/*! \fn void FT_StartUp(Front *front, F_BASIC_DATA *ft_basic)
 *  \ingroup INITIALIZATION
    \brief Initialize front computational grid, interface and function hooks, 
     default function calls, and if restart, read interface from restart 
     directory. Constructor for Front object.
    \param front @b inout	Pointer to Front.
    \param f_basic @b in	Structure for initialization information.
 */
   IMPORT  void FT_StartUp(Front* front ,
   			     F_BASIC_DATA* f_basic );

/*! \fn void FT_InitDebug(char *inname)
 *  \ingroup INITIALIZATION
    \brief Initialize strings for debugging option of the problem,
     read the input file.
    \param inname @b in	The name of the input file 
 */
   IMPORT  void FT_InitDebug(char *inname );

/*! \fn void FT_InitIntfc(Front *front, 
     				    LEVEL_FUNC_PACK *level_func_pack)
    \ingroup INITIALIZATION
    \brief Initialize the interface interface curves (2D) and surfaces (3D)
     using level function and parameters, or in some cases, point set.
     Install boundary curves (2D) or surfaces (3D). For 2D, the boundary
     is installed by default. For 3D, the boundary is installed on request.
    \param front @b inout	Pointer to Front.
    \param level_func_pack @b in	Structure of level function and parameters, or point set.
 */
   IMPORT  void FT_InitIntfc(Front *front ,
			       LEVEL_FUNC_PACK *level_func_pack );

/*! \fn void FT_ClipIntfcToSubdomain(Front *front) 
    \ingroup INITIALIZATION
    \brief Clip the Initial interface of the front to a parallel subdomain.
    \param front @b inout	Pointer to Front.
 */
   IMPORT  void FT_ClipIntfcToSubdomain(Front *front);

/*! \fn void FT_InitFrontVeloFunc(Front *front,
                   VELO_FUNC_PACK *velo_func_pack)
    \ingroup INITIALIZATION
    \brief Initialize front velocity function for front point propagation.
     The velocity function use point and other related structures as input,
     must also supply parameters needed for the velocity function.
    \param front @b inout	Pointer to Front.
    \param velo_func_pack @b in	Structure containing velocity function and parameters.
 */
   IMPORT  void    FT_InitFrontVeloFunc(Front *front ,
 				 VELO_FUNC_PACK *velo_func_pack );


/*! \fn void FT_ReadTimeControl(char *in_name, Front *front)
 *  \ingroup INITIALIZATION
    \brief Read the time domain and control information from input file,
     including maximum time, maximum step, restart print interval, 
     movie frame output interval, CFL factor, and redistribution step interval.
    \param in_name @b in	The name of input file
    \param front @b inout	Pointer to Front for computation
 */
   IMPORT  void FT_ReadTimeControl(char *in_name ,
	   				 Front *front );

/*! \fn void FT_ResetTime(Front *front)
 *  \ingroup INITIALIZATION
    \brief Reset the time to 0.0, time step to 0, print interval index to 0,
     movie interval index to 0.
    \param front @b inout	Pointer to Front.
 */
   IMPORT  void FT_ResetTime(Front *front );

/*! \fn POINTER *FT_CreateLevelHyperSurfs(RECT_GRID *rgr, INTERFACE *intfc, int neg_comp, int pos_comp, double (*func)(POINTER,double*), POINTER   func_params, int w_type, int *num_hs)
 *  \ingroup INITIALIZATION
    \brief This function creates a set of hypersurfaces (curves in 2D and 
     surfaces in 3D) using a level function (provided by the caller). The
     function return the handle for the array of hyper surfaces as (POINTER*).
    \param rgr @b in Pointer to the rectangular grid in which the hyper surfaces are constructed.
    \param intfc @b in Pointer to the interface structure.
    \param neg_comp @b in Region index of the negative side of the level surfaces.
    \param pos_comp @b in Region index of the positive side of the level surfaces.
    \param func @b in Pointer level set function, level surface is the zero set.
    \param func_params @b in Pointer level set function parameters.
    \param w_type @b in Wave type of the hyper surfaces.
    \param num_hs @b out Address of number of segments of the hyper surfaces.
 */
   IMPORT  POINTER *FT_CreateLevelHyperSurfs(
			RECT_GRID *rgr,
        		INTERFACE *intfc,
        		int neg_comp,
        		int pos_comp,
        		double    (*func)(POINTER,double*),
        		POINTER   func_params,
        		int       w_type,
        		int       *num_hs);

/*! \fn void FT_Propagate(Front *front)
 *  \ingroup PROPAGATION
    \brief Propagate the Front for one time step. The process includes  
     advancing the front forward by one step to new positions, redistributing
     new interface mesh and resolving physical and topological bifurcations.
     The interface in the front is replaced by a new and propagated one and
     the old one is freed.
    \param front Pointer to Front to be propagated
 */

   IMPORT  void FT_Propagate(Front *front);

/*! \fn void FT_RedistMesh(Front *front)
 *  \ingroup OPTIMIZATION
    \brief This is an independent call for redistribution and optimization
     of the interface mesh. A parallel communication of front should be called
     following this call to ensure that the redistribution is consistent
     globally. 
    \param front @b inout Pointer to Front to be redistributed
 */
   IMPORT  void FT_RedistMesh(Front *front);

/*! \fn void FT_OptimizeMesh(Front *front, SCALED_REDIST_PARAMS params)
 *  \ingroup OPTIMIZATION
    \brief This is an independent call for redistribution and optimization
     of the interface mesh. A parallel communication of front should be called
     following this call to ensure that the redistribution is consistent
     globally. The function will redistribute all curves and surfaces (3D).
     It will recursively do either 10 times or when nothing to be done.
    \param front @b inout Pointer to Front to be redistributed
    \param params @b in structure of scaled redistribution parameters.
 */
   IMPORT  void FT_OptimizeMesh(Front *front, SCALED_REDIST_PARAMS params);

/*! \fn void FT_OptimizeSurfMesh(Front *front, SURFACE *surf, SCALED_REDIST_PARAMS params)
 *  \ingroup OPTIMIZATION
    \brief This is an independent call for redistribution and optimization
     of a surf mesh. No parallel communication is called after redistribution.
     Return YES (nothing_done == YES) if no more triangle to be redistributed.
    \param front @b inout Pointer to Front to be redistributed
    \param surf @b inout Pointer to the surface to be redistributed
    \param params @b in structure of scaled redistribution parameters.
 */
   IMPORT  boolean FT_OptimizeSurfMesh(Front *front, SURFACE *surf, SCALED_REDIST_PARAMS params);

/*! \fn void FT_OptimizeCurveMesh(Front *front, CURVE *curve, SCALED_REDIST_PARAMS params)
 *  \ingroup OPTIMIZATION
    \brief This is an independent call for redistribution and optimization
     of a curve mesh. No parallel communication is called after redistribution.
     Return YES (nothing_done == YES) if no more bond to be redistributed.
    \param front @b inout Pointer to Front to be redistributed
    \param curve @b inout Pointer to the curve to be redistributed
    \param params @b in structure of scaled redistribution parameters.
 */
   IMPORT  boolean FT_OptimizeCurveMesh(Front *front, CURVE *curve, SCALED_REDIST_PARAMS params);

/*! \fn void FT_SetCurveSpacing(Front *front, double scaled_spacing)
 *  \ingroup OPTIMIZATION
    \brief This function set the optimal spacing for curve redistribution.
     The scaled spacing is in the unit of regular grid spacing h.
    \param front @b inout Pointer to Front to be redistributed
    \param scaled_spacing @b in value of front spacing to be set.
 */

   IMPORT  void FT_SetCurveSpacing(Front *front, double scaled_spacing);

/*! \fn void FT_OptimizeCurveMeshWithEqualBonds(Front *front, CURVE *curve)
 *  \ingroup OPTIMIZATION
    \brief This is an independent call for redistribution and optimization
     of curve. A parallel communication of front should be called
     following this call to ensure that the redistribution is consistent
     globally. 
    \param front @b inout Pointer to Front of the curve to be redistributed
    \param curve @b inout Pointer to curve to be redistributed
 */
   IMPORT  void FT_OptimizeCurveMeshWithEqualBonds(Front *front, CURVE *curve);

/*! \fn void FT_SetTimeStep(Front *front)
 *  \ingroup TIME
    \brief Calculate front->dt for next time step using the recorded maximum
     speed from previous time step, it step is reduced by a CFL factor to
     to ensure numerical stability.
     \param front @b inout	Pointer to the Front.
 */
   IMPORT  void FT_SetTimeStep(Front *front );

/*! \fn void FT_SetOutputCounter(Front *front)
 *  \ingroup INITIALIZATION
    \brief This function is used in restart to set the printing index and
     movie output index to appropiate number according to the restart time.
    \param front @b inout	Pointer to Front.
 */
   IMPORT  void FT_SetOutputCounter(Front *front );

/*! \fn void FT_TimeControlFilter(Front *front)
 *  \ingroup TIME
    \brief To further reduce time step if the printing or movie output time
     is smaller than the time after next time step. Increment priting
     or/and movie output index if either of both of them are met.
    \param front @b inout	Pointer to Front.
 */
   IMPORT  void FT_TimeControlFilter(Front *front );

/*! \fn boolean FT_IsSaveTime(Front *front)
 *  \ingroup TIME
    \brief Signals that time is reached for printing restart files. returns
     YES to indicate that restart file should be generated. This function
     does not write the output.
    \param front @b in	Pointer to Front.
 */

   IMPORT  boolean FT_IsSaveTime(Front *front );

/*! \fn boolean FT_IsDrawTime(Front *front)
 *  \ingroup TIME
    \brief Signals that time is reached for output of a movie frame. returns
     YES to indicate that a movie frame should be generated. This function
     does not write the movie frame.
    \param front @b in	Pointer to Front.
 */

   IMPORT  boolean FT_IsDrawTime(Front *front );

/*! \fn boolean FT_TimeLimitReached(Front *front)
 *  \ingroup TIME
    \brief Signals that time is reached for termination of the run. Return
     YES either maximum time has been reached or maximum time step has been 
     reached.
    \param front @b in	Pointer to Front.
 */

   IMPORT  boolean FT_TimeLimitReached(Front *front );

/*! \fn void FT_RecordMaxFrontSpeed(int dir, double speed, POINTER state, double *coords, Front *front)
 *  \ingroup TIME
    \brief This function compare and record maximum speed of the front.
     The recorded speed will be used to determine the time step which
     must satisfy the CFL condition.
    \param dir @b in	Direction of the speed, e. g. 0(x), 1(y), or 2(z).
    \param speed @b in	The speed to be compared and recorded (if bigger than max).
    \param state @b in	Pointer to the state, optional, can be set to NULL.
    \param coords @b in	Coordinates of the point where the speed occurs.
    \param front @b inout	Pointer to Front.
 */
   IMPORT  void FT_RecordMaxFrontSpeed(
				int dir ,
				double speed ,
				POINTER state ,
				double *coords ,
   				Front* front );

/*! \fn boolean FT_AddTimeStepToCounter(Front *front)
 *  \ingroup PROPAGATION
    \brief Add front->dt to front->time after propagation.
    \param front @b inout	Pointer to Front.
 */

   IMPORT  void FT_AddTimeStepToCounter(Front *front );

/*! \fn void FT_Save(Front *front)
 *  \ingroup OUTPUT
    \brief Output front geometric data to the directory of out_name.
     The data can be used for restart of the run. 
    \param front @b in	Pointer to Front.
 */

   IMPORT  void FT_Save(Front *front);

/*! \fn void FT_Draw(Front *front)
 *  \ingroup OUTPUT
    \brief Output a movie frame, currently includes GD, hdf, vtk formats.
    \param front @b in	Pointer to Front.
 */

   IMPORT  void FT_Draw(Front *front);

/*! \fn void FT_XgraphSampleLine(char *dirname,char *varname,boolean data_in_domain,int size,double *x,double *var)
    \ingroup OUTPUT
    \brief  This function output variable sample along a grid line as
     xgraph data file. It considers parallelization of subdomains.
    \param dirname @b in Name of the directory for the output file.
    \param varname @b in Name of the variable for the output file.
    \param data_in_domain @b in Yes if the subdomain contains data.
    \param size @b in Size of the data in the subdomain.
    \param x @b in Horizontal axis containing coordinate of the sample line.
    \param var @b in Vertical axis containing variable data.
 */
   IMPORT  void FT_XgraphSampleLine(char *dirname,
			char *varname,
			boolean data_in_domain,
			int size,
			double *x,
			double *var);

/*! \fn void FT_MakeGridIntfc(Front *front)
 *  \ingroup GRIDINTFC
    \brief Make a duplicate interface whose topological grid is the
     expanded dual grid, currently with buffer of 4h for PERIODIC and
     SUBDOMAIN boundary, 1h for other boundaries. Install crossings
     of the interface and the expanded dual grid and store them in the
     interface table. These crossings can be used to interact with
     various PDE solvers.
    \param front @b inout	Pointer to Front.
 */

   IMPORT  void FT_MakeGridIntfc(Front *front );

/*! \fn void FT_FreeGridIntfc(Front *front)
 *  \ingroup GRIDINTFC
    \brief Delete and free space of grid crossing interface made by
     the function FT_MakeGridIntfc().
    \param front @b inout	Pointer to Front.
 */

   IMPORT  void FT_FreeGridIntfc(Front *front );

/*! \fn void FT_MakeCompGridIntfc(Front *front)
 *  \ingroup GRIDINTFC
    \brief Make a duplicate interface whose topological grid is the
     expanded comp grid, currently with buffer of 4h for PERIODIC and
     SUBDOMAIN boundary, 1h for other boundaries. Install crossings
     of the interface and the expanded dual grid and store them in the
     interface table. These crossings can be used to interact with
     various PDE solvers.
    \param front @b inout	Pointer to Front.
 */

   IMPORT  void FT_MakeCompGridIntfc(Front *front );

/*! \fn void FT_FreeGridIntfc(Front *front)
 *  \ingroup GRIDINTFC
    \brief Delete and free space of grid crossing interface made by
     the function FT_MakeCompGridIntfc().
    \param front @b inout	Pointer to Front.
 */

   IMPORT  void FT_FreeCompGridIntfc(Front *front );

/*! \fn void FT_FreeOldGridIntfc(Front *front)
 *  *  \ingroup GRIDINTFC
 *      \brief Delete and free space of grid crossing interface of
 *           front->old_grid_intfc
 *               \param front @b inout       Pointer to Front.
 *                */

   IMPORT  void FT_FreeOldGridIntfc(Front *front );

/*! \fn void FT_FreeFront(Front *front)
 *  \ingroup GRIDINTFC
    \brief Delete and free space occupied by the front including grid_intfc 
     if still there, and interf.
    \param front @b inout	Pointer to Front.
 */

   IMPORT  void FT_FreeFront(Front *front );

/*! \fn void FT_FreeMainIntfc(Front *front)
 *  \ingroup GRIDINTFC
    \brief Delete and free space occupied by the primary interface of the front 
    \param front @b inout	Pointer to Front.
 */

   IMPORT  void FT_FreeMainIntfc(Front *front );

/*! \fn boolean FT_NormalAtGridCrossing(Front *front, int *icoords, GRID_DIRECTION dir, int comp, double *nor, HYPER_SURF **hs, double *crx_coords)
 *  \ingroup GRIDINTFC
    \brief Standing at grid icoords, looking to the direction dir, this
     function looks for the nearest interface cross on the grid line segment.
     The function returns YES if the crossing exists, in such case, the
     crossing coordinates are copied to crx_coords, the corresponding
     hyper surface (curce in 2D and surface in 3D) is assigned to hs,
     and the normal vector to the side of comp. If no crossing exists, 
     the function return NO;
    \param front @b in	Pointer to Front.
    \param icoords @b in	Grid point coordinate indices.
    \param dir @b in	Direction to which the crossing is to be found.
    \param comp @b in	Component (domain index) of the grid point at icoord.
    \param nor @b out	normal vector at the crossing to the side of comp.
    \param hs @b out	Crossing hyper surface (curve in 2D and surface in 3D).
    \param crx_coords @b out	Crossing coordinates.
 */

   IMPORT  boolean FT_NormalAtGridCrossing(Front *front ,
   				int *icoords ,
   				GRID_DIRECTION  dir ,
   				int  comp ,
   				double *nor ,
   				HYPER_SURF **hs ,
   				double *crx_coords );

/*! \fn void FT_ComputeVolumeFraction(Front *front, int num_phases, COMPONENT *compsCELL_PART *cell_part)
 *  \ingroup GRIDINTFC
    \brief  This function compute fraction volume of each grid cell of
    a conservative variable on topological grid. The conservative
    variable is marked by the solution component as an input variable.
    \param front @b in	Pointer to Front.
    \param num_phases @b in number of phases in the domain..
    \param comps @b in component index of each phase..
    \param cell_part @b inout Cell partition structure.
 */

   IMPORT  void FT_ComputeVolumeFraction(Front *front ,
   				int num_phases,
				COMPONENT *comps,
				CELL_PART *cell_part);

/*! \fn boolean FT_CoordsAtGridCrossing(Front *front, int *icoords, GRID_DIRECTION dir, int comp, HYPER_SURF **hs, double *crx_coords)
 *  \ingroup GRIDINTFC
    \brief Standing at grid icoords, looking to the direction dir, this
     function looks for the nearest interface cross on the grid line segment.
     The function returns YES if the crossing exists, in such case, the
     crossing coordinates are copied to crx_coords, the corresponding
     hyper surface (curce in 2D and surface in 3D) is assigned to hs,
     If no crossing exists, the function return NO;
    \param front @b in	Pointer to Front.
    \param icoords @b in	Grid point coordinate indices.
    \param dir @b in	Direction to which the crossing is to be found.
    \param comp @b in	Component (domain index) of the grid point at icoord.
    \param hs @b out	Crossing hyper surface (curve in 2D and surface in 3D).
    \param crx_coords @b out	Crossing coordinates.
 */

   IMPORT  boolean FT_CoordsAtGridCrossing(Front *front ,
				INTERFACE *grid_intfc,
   				int *icoords,
   				GRID_DIRECTION  dir,
   				int  comp,
   				HYPER_SURF **hs,
   				double *crx_coords);

/*! \fn boolean FT_StateStructAtGridCrossing(Front *front, int *icoords, GRID_DIRECTION dir, int comp, POINTER *state, HYPER_SURF **hs, double *crx_coords)
 *  \ingroup GRIDINTFC
    \brief Standing at grid icoords, looking to the direction dir, this
     function looks for the nearest interface cross on the grid line segment.
     The function returns YES if the crossing exists, in such case, the
     crossing coordinates are copied to crx_coords, the corresponding
     hyper surface (curce in 2D and surface in 3D) is assigned to hs,
     and the state on the side of comp is assigned to the state pointer.
     If no crossing exists, the function return NO;
    \param front @b in	Pointer to Front.
    \param icoords @b in	Grid point coordinate indices.
    \param dir @b in	Direction to which the crossing is to be found.
    \param comp @b in	Component (domain index) of the grid point at icoord.
    \param state @b out	State at the crossing on the side of comp.
    \param hs @b out	Crossing hyper surface (curve in 2D and surface in 3D).
    \param crx_coords @b out	Crossing coordinates.
 */

   IMPORT  boolean FT_StateStructAtGridCrossing(Front *front ,
				INTERFACE *grid_intfc,
   				int *icoords,
   				GRID_DIRECTION  dir,
   				int  comp,
   				POINTER *state,
   				HYPER_SURF **hs,
   				double *crx_coords);

/*! \fn boolean FT_StateVarAtGridCrossing(Front *front, INTERFACE *grid_intfc, int *icoords, GRID_DIRECTION dir, int comp, double (*state_func)(Locstate), double *ans, double *crx_coords)
 *  \ingroup GRIDINTFC
    \brief This function performs the same way as the function
     FT_StateStructAtGridCrossing() except it assigns a specific state
     variable instead of the whole state structure. Since front does not
     know the state structure of application, an application specific
     state retrieving function state_func() must be supplied.
    \param front @b in	Pointer to Front.
    \param grid_intfc @b in	Pointer to the grid_intfc.
    \param icoords @b in	Grid point coordinate indices.
    \param dir @b in	Direction to which the crossing is to be found.
    \param comp @b in	Component (domain index) of the grid point at icoord.
    \param state_func @b in	State function for the requested state variable.

    \param ans @b 
    \param crx_coords @b out	Crossing coordinates.
 */

   IMPORT  boolean FT_StateVarAtGridCrossing(Front *front ,
   				int *icoords ,
   				GRID_DIRECTION  dir ,
   				int  comp ,
				double (*state_func)(Locstate) ,
				double *ans ,
   				double *crx_coords );

/*! \fn boolean FT_StateVarAtPrevStepCrossing(Front *front, INTERFACE *grid_intfc, int *icoords, GRID_DIRECTION dir, int comp, double (*state_func)(Locstate), double *ans, double *crx_coords)
 *  \ingroup GRIDINTFC
    \brief This function performs the same way as the function
     FT_StateVarAtGridCrossing() except it checks the grid_intfc 
     of previous time step. 
    \param front @b in	Pointer to Front.
    \param grid_intfc @b in	Pointer to the grid_intfc.
    \param icoords @b in	Grid point coordinate indices.
    \param dir @b in	Direction to which the crossing is to be found.
    \param comp @b in	Component (domain index) of the grid point at icoord.
    \param state_func @b in	State function for the requested state variable.

    \param ans @b 
    \param crx_coords @b out	Crossing coordinates.
 */

   IMPORT  boolean FT_StateVarAtPrevStepCrossing(Front *front ,
   				int *icoords ,
   				GRID_DIRECTION  dir ,
   				int  comp ,
				double (*state_func)(Locstate) ,
				double *ans ,
   				double *crx_coords );

/*! \fn HYPER_SURF *FT_HyperSurfAtGridCrossing(Front *front, int *icoords, GRID_DIRECTION dir,int wave_type)
 *  \ingroup GRIDINTFC
    \brief Sitting at icoords and look to the direction dir, this function 
     detects the nearest hyper surface (curve in 2D and surface in 3D) 
     on the grid segment. Return pointer to hyper surface if there is
     one, return NULL if no crossing hyper surface is found.
    \param front @b in	Pointer to Front.
    \param icoords @b in	Grid point coordinate indices.
    \param dir @b in	Direction to which the crossing is to be found.
    \param wave_type @b in wave type of the hyper surface to search, if ANY_WAVE_TYPE, will return any hyper surface at the crossing.
 */

   IMPORT  HYPER_SURF *FT_HyperSurfAtGridCrossing(Front *front ,
   				int *icoords ,
   				GRID_DIRECTION  dir,
				int wave_type);

/*! \fn void FT_GetCrossingHistory(Front *front, int *icoords, GRID_DIRECTION dir,double *crx_coords_new,double *crx_coords_old)
 *  \ingroup GRIDINTFC
    \brief Sitting at icoords and look to the direction dir, this function 
     find the crossing history, the current and past crossing potitions
    \param front @b in	Pointer to Front.
    \param icoords @b in	Grid point coordinate indices.
    \param dir @b in	Direction to which the crossing is to be found.
    \param crx_coords_new @b out current crossing coordinates.
    \param crx_coords_old @b out previous crossing coordinates.
 */

   IMPORT  void FT_GetCrossingHistory(Front *front ,
   				int *icoords ,
   				GRID_DIRECTION  dir,
				double *crx_coords_new,
				double *crx_coords_old);

/*! \fn boolean FT_IntrpStateVarAtCoords(Front *front, int comp, double *coords, double *var_array, double (*state_func)(POINTER), double *ans, double *default_ans)
 *  \ingroup GRIDINTFC
    \brief Interpolate a state variable at a space point with coords. If 
     comp == NO_COMP, it interpolates with no regard of interface. Otherwise
     it will interpolate in the subdomain of comp. The state_func() is needed
     to tell the function how to retrieve the variable from the interface
     state. The interpolated variable is assigned in ans. Return YES if
     the interpolation is successful.
    \param front @b in	Pointer to Front.
    \param comp @b in	Component in which the state should be interpolated.
    \param var_array @b in	Array of the variable on the expanded dual grid.
    \param state_func() @b in	Function to retrieve the variable from the interface state pointer.
    \param ans @b out	Address of the interpolated variable.
    \param default_ans @b in	Address of default solution, if NULL, the function will look for solution at nearest interface point.
 */

   IMPORT  boolean FT_IntrpStateVarAtCoords(Front *front ,
   				int comp , 
				double *coords , 
				double *var_array , 
				double (*state_func)(POINTER) , 
				double *ans,
				double *default_ans);

/*! \fn boolean FT_CompGridIntrpStateVarAtCoords(Front *front, int comp, double *coords, double *var_array, double (*state_func)(POINTER), double *ans, double *default_ans)
 *  \ingroup GRIDINTFC
    \brief Interpolate a state variable at a space point with coords 
     on computational grid. If 
     comp == NO_COMP, it interpolates with no regard of interface. Otherwise
     it will interpolate in the subdomain of comp. The state_func() is needed
     to tell the function how to retrieve the variable from the interface
     state. The interpolated variable is assigned in ans. Return YES if
     the interpolation is successful.
    \param front @b in	Pointer to Front.
    \param comp @b in	Component in which the state should be interpolated.
    \param var_array @b in	Array of the variable on the expanded dual grid.
    \param state_func() @b in	Function to retrieve the variable from the interface state pointer.
    \param ans @b out	Address of the interpolated variable.
    \param default_ans @b in	Address of default solution, if NULL, the function will look for solution at nearest interface point.
 */

   IMPORT  boolean FT_CompGridIntrpStateVarAtCoords(Front *front ,
   				int comp , 
				double *coords , 
				double *var_array , 
				double (*state_func)(POINTER) , 
				double *ans,
				double *default_ans);

/*! \fn boolean FT_NearestRectGridVarInRange(Front *front, int comp, double *coords, double *var_array, int range, double *ans)
 *  \ingroup GRIDINTFC
    \brief Find the state variable on rectangular grid point which
     has the same component as the input and is nearest to the input
     coordinate. Return YES if such point is found, and NO if no such point 
     is found. In the latter case, the value of the ans is set to zero.
    \param front @b in	Pointer to Front.
    \param comp @b in	Component in which the state should be interpolated.
    \param var_array @b in	Array of the variable on the expanded dual grid.
    \param range @b in	Rnage of search in number of grid cells.
    \param ans @b out	Address of the interpolated variable.
 */

   IMPORT  boolean FT_NearestRectGridVarInRange(Front *front ,
   				int comp , 
				double *coords , 
				double *var_array , 
				int range,
				double *ans);

/*! \fn boolean FT_FindNearestIntfcPointInRange(Front *front, int comp, double *coords, double *intfc_point, double *t, HYPER_SURF_ELEMENT **hse, HYPER_SURF **hs, int range)
 *  \ingroup GRIDINTFC
    \brief Given a space coordinate coords, this function tries to find the
     nearest point on the interface within range, together with its associated 
     hyper surface element (bond in 2D and triangle in 3D) and hyper surface 
     (curve in 2D and surface in 3D). If no hyper surface is within range, the
     function returns NO;
    \param front @b in	Pointer to Front.
    \param comp @b in	Component index of the space coordinates.
    \param coords @b in	Coordinates of the space point.
    \param bdry @b in	Whether boundary hypersurface should be used.
    \param intfc_point @b out	Coordinates of the nearest point on the interface.
    \param t @b out	Interpolation factors.
    \param hse @b out	Residing hyper surface element (bond in 2D and tri in 3D).
    \param hs @b out	Residing hyper surface (curve in 2D and surface in 3D).
    \param range @b in	Range for the search in unit of grid spacing.
 */

   IMPORT  boolean FT_FindNearestIntfcPointInRange(Front *front ,
   					int comp ,
   					double *coords ,
					USE_BOUNDARIES bdry,
   					double *intfc_point ,
   					double *t ,
   					HYPER_SURF_ELEMENT **hse ,
   					HYPER_SURF **hs ,
					int range );

/*! \fn void FT_NormalAtPoint(POINT *p, Front *front, double *normal, int comp)
 *  \ingroup GEOMETRY
    \brief Get the normal vector at the point p.
    \param p @b in	Pointer to a valid point on interface of the front.
    \param front @b in	Pointer to Front.
    \param normal @b out	The normal vector.
    \param comp @b in	The component of the side which the normal is pointing to.
 */

   IMPORT  void FT_NormalAtPoint(POINT *p ,
   				Front *front ,
   				double *normal ,
				int comp );

/*! \fn void FT_CurvatureAtPoint(POINT *p, Front *front, double *curvature)
 *  \ingroup GEOMETRY
    \brief Get the curvature at the point p.
    \param p @b in	Pointer to a valid point on interface of the front.
    \param front @b in	Pointer to Front.
    \param curvature @b out	The address where the curvature is to be assigned.
 */

   IMPORT  void FT_CurvatureAtPoint(POINT *p ,
   				Front *front ,
   				double *curvature );


/*! \fn double FT_GridSizeInDir(double *dir, Front *front)
 *  \ingroup GEOMETRY
    \brief The grid size in direction dir is (dir dot h). If h are the same
     in all directions (square mesh), it returns the h.
    \param dir @b in	The unit vector of the direction.
    \param front @b in	Pointer to Front.
 */

   IMPORT  double FT_GridSizeInDir(double *dir ,
   				Front *front );

/*! \fn Nor_stencil *FT_CreateNormalStencil(Front *front, POINT *p, int comp, int num_pts)
 *  \ingroup GEOMETRY
    \brief This function create a normal stencil at the interface point p
     with size (number of points) num_pts. The normal stencil in in the
     ambient with component comp.
    \param front @b in	Pointer to Front.
    \param p @b in	Pointer to the point.
    \param comp @b in	Index of the region of the stencil.
    \param num_pts @b in	Number of point in the normal stencel.
 */

   IMPORT  Nor_stencil *FT_CreateNormalStencil(Front *front ,
                                POINT *p ,
                                int comp ,
                                int num_pts );

/*! \fn void FT_ReflectPointThroughBdry(Front *front, HYPER_SURF *hs, double *coords, int comp, double *coords_bdry, double *coords_ref, double *nor)
 *  \ingroup GEOMETRY
    \brief Given the coordinates coords, this function find the reflected
     coordinates coordsrefl through the hypersurface hs, it also provide
     the normal vector nor at the reflection. Return NO if conditions not
     satisfied.
    \param front @b in	Pointer to Front.
    \param hs @b in  Pointer to the hypersurface (curve in 2D, surface in 3D).
    \param comp @b in  Component of the ambient.
    \param coords @b in	Coordinates to be reflected.
    \param coords_bdry @b out boundary point of reflection.
    \param coords_ref @b out Coordinates after reflection.
    \param nor @b out Normal vector at reflection.
 */

   IMPORT  boolean FT_ReflectPointThroughBdry(
				Front *front ,
				HYPER_SURF *hs ,
				double *coords ,
				int comp ,
				double *coords_bdry ,
				double *coords_ref ,
   				double *normal );

/*! \fn void FT_SetDirichletBoundary(Front *front, void (*state_func)(double*,HYPER_SURF*,Front*,POINTER,POINTER), const char *state_func_name, POINTER state_func_params, POINTER state, HYPER_SURF *hs)
 *  \ingroup BOUNDARY
    \brief This function sets state of a Dirichlet boundary as a hyper surface.
     It provide two methods to set the boundary: (a) set the boundary with
     a constant state, in this case, the address of the state as a pointer
     must be supplied as a pointer while state_func, state_func_name and
     state_func_params must be set to NULL; (b) set the boundary state through
     a state function, in this case, state must be set to NULL while
     state_func, state_func_name and state_func_params must be supplied.
     Currently only flow-through and time-dependent functions are available.
    \param front @b inout	Pointer to Front.
    \param state_func @b in	Function pointer where state is evaluated in option (b).
    \param state_func_name @b in	Function name as a char string, for option (b).
    \param state_func_params @b in	Associated function parameters, for option (b).
    \param state @b in	The address of the constant state, for option (a).
    \param hs @b inout	The pointer to the boundary as a hyper surface structure.
 */

   IMPORT  void FT_SetDirichletBoundary(Front *front ,
   				void (*state_func)(double*,HYPER_SURF*,Front*,POINTER,POINTER) ,
				const char *state_func_name ,
				POINTER state_func_params ,
				POINTER state ,
				HYPER_SURF *hs );

/*! \fn void FT_InsertDirichletBoundary(Front *front, void (*state_func)(double*,HYPER_SURF*,Front*,POINTER,POINTER), const char *state_func_name, POINTER state_func_params, POINTER state, HYPER_SURF *hs, int index)
 *  \ingroup BOUNDARY
    \brief This function insert state to a Dirichlet boundary.
     The six (four for 2D) sides of the boundary memory has been allocated, 
     this function fills in the content of the boundary state.
     It provide two methods to set the boundary: (a) set the boundary with
     a constant state, in this case, the address of the state as a pointer
     must be supplied as a pointer while state_func, state_func_name and
     state_func_params must be set to NULL; (b) set the boundary state through
     a state function, in this case, state must be set to NULL while
     state_func, state_func_name and state_func_params must be supplied.
     Currently only flow-through and time-dependent functions are available.
    \param front @b inout	Pointer to Front.
    \param state_func @b in	Function pointer where state is evaluated in option (b).
    \param state_func_name @b in	Function name as a char string, for option (b).
    \param state_func_params @b in	Associated function parameters, for option (b).
    \param state @b in	The address of the constant state, for option (a).
    \param hs @b inout	The pointer to the boundary as a hyper surface structure.
    \param index @b The boundary surface index.
 */

   IMPORT  void FT_InsertDirichletBoundary(Front *front ,
   				void (*state_func)(double*,HYPER_SURF*,Front*,POINTER,POINTER) ,
				const char *state_func_name ,
				POINTER state_func_params ,
				POINTER state ,
				HYPER_SURF *hs ,
				int index);

/*! \fn void FT_MixedBoundaryHypSurfs(INTERFACE *intfc, int idir, int nb, int w_type, int *num_hs)
 *  \ingroup BOUNDARY
    \brief This function returns a set of hyper surfaces (curves in 2D and 
     surfaces in 3D) whose wave type matches the input w_type. If the 
     boundary on idir and side nb is not MIXED_TYPE_BOUNDARY, the function
     returns NULL. Otherwise it will return an array of hyper surfaces.
     The function also assign num_hs, the total number of hyper surfaces
     in the returned set.
    \param intfc @b in	POINTER to the input interface.
    \param idir @b in Direction 0 for x-dir, 1 for y-dir, 2 for z-dir.
    \param nb @b in Side, 0 for lower side, 1 for upper side.
    \param w_type @b in Wave type of the hyper surfaces looking for.
    \param num_hs @b out Number of hyper surfaces in the matching set.
 */

   IMPORT  HYPER_SURF **FT_MixedBoundaryHypSurfs(
   			INTERFACE *intfc,
			int idir,
			int nb,
			int w_type,
			int *num_hs);

/*! \fn void FT_PromptSetMixedTypeBoundary2d(char *in_name, Front *front)
 *  \ingroup BOUNDARY
    \brief This function sets the node types and curve wave types of nodes
     and curves at the MIXED_TYPE_BOUNDARY side(s). This function is only for
     2D. It will prompt and set nodes in ascending order and then curves
     between the node pairs.
    \param in_name @b in	Character string of input file name to open.
    \param front @b inout	Pointer to Front.
 */

   IMPORT  void FT_PromptSetMixedTypeBoundary2d(
   			char *in_name,
   			Front *front);

/*! \fn int FT_Dimension()
 *  \ingroup QUERY
    \brief This function returns the spatial dimension of current run. */

   IMPORT  int FT_Dimension();

/*! \fn int FT_RectBoundaryType(Front *front, int dir, int side)
 *  \ingroup QUERY
    \brief This function returns rectangular boundary type in direction
     dir and side side.
    \param front @b in	Pointer the front.
    \param dir @b in Direction of the boundary.
    \param side @b in Side of the rectangular domain, 0 for lower side and 1 for upper side.
 */

   IMPORT  int FT_RectBoundaryType(Front *front ,
				int dir,
				int side);

/*! \fn boolean FT_FrontContainWaveType(Front *front, int w_type)
 *  \ingroup QUERY
    \brief This function tells if front contains hyper-surface of certain
     type. It returns YES if it does, NO if it does not contain.
    \param front @b in	Pointer the front.
    \param w_type @b in Type of wave for the hyper-surface.
 */

   IMPORT  boolean FT_FrontContainWaveType(Front *front , int w_type);

/*! \fn HYPER_SURF *FT_RectBoundaryHypSurf(INTERFACE *intfc, int wave_type, int dir, int side)
 *  \ingroup QUERY
    \brief This function looks for a boundary hyper surface (curve in 2D and
     surface in 3D) at the rectangular domain border with the given wave type
     of the hyper surface in direction dir and side side. Returns NULL if
     no match is found and return pointer to the hyper surface if found.
    \param intfc @b in	Pointer to an interface of the front.
    \param wave_type @b in Wave type of the hyper surface requested.
    \param dir @b in Direction of the boundary, 0 for x-dir, 1 for y-dir and 2 for z-dir.
    \param side @b in Side of the rectangular domain, 0 for lower side and 1 for upper side.
 */

   IMPORT  HYPER_SURF *FT_RectBoundaryHypSurf(INTERFACE *intfc ,
				int wave_type,
				int dir,
				int side);

/*! \fn HYPER_SURF **FT_InteriorHypSurfs(INTERFACE *intfc, int wave_type, int *num_hs)
 *  \ingroup QUERY
    \brief This function looks for interior hyper surfaces (curve in 2D and
     surface in 3D) with the given wave type. Returns NULL if no match is 
     found and return pointer to an array of hyper surfaces if found. This
     function allocates the memory for the array of pointers to hyper surfaces.
    \param intfc @b in	Pointer to an interface of the front.
    \param wave_type @b in Wave type of the hyper surfaces requested.
    \param num_hs @b out Address to interger of number of hyper surfs found.
 */

   IMPORT  HYPER_SURF **FT_InteriorHypSurfs(INTERFACE *intfc ,
				int wave_type,
				int *num_hs);

/*! \fn void FT_ParallelExchIntfcBuffer(Front *front)
 *  \ingroup PARALLEL
    \brief This is a root level parallel communication function for front
     interface geometry and states. It will cut the old buffer parts of the
     interface and patch it with new parts received from other subdomains
     or periodically shifted sides. This is a synchronous function and must
     be called synchronously by every processor.
    \param front @b inout	Pointer to Front.
 */
   IMPORT  void FT_ParallelExchIntfcBuffer(Front* front );

/*! \fn void FT_ParallelExchGridArrayBuffer(double *grid_array, Front *front, int *symmetry)
 *  \ingroup PARALLEL
    \brief This is a parallel communication function for a double array
     on the expanded dual grid of the grid_intfc in front. It will cut 
     the old buffer parts of the array and patch it with new buffer parts 
     received from other subdomains or periodically shifted sides. This is a 
     synchronous function and must be called synchronously by every processor.
    \param grid_array @b inout A double array of variable on expanded duel grid.
    \param front @b in	Pointer to Front.
    \param symmetry @b in Flag of symmetry in each direction (in reflection).
 */
   IMPORT  void FT_ParallelExchGridArrayBuffer(double *grid_array ,
   				Front* front ,
				int* symmetry);

/*! \fn void FT_ParallelExchGridStructArrayBuffer(POINTER pstruct, Front *front, int psize)
 *  \ingroup PARALLEL
    \brief This is a parallel communication function for a struct array
     on the expanded dual grid of the grid_intfc in front. It will cut 
     the old buffer parts of the array and patch it with new buffer parts 
     received from other subdomains or periodically shifted sides. This is a 
     synchronous function and must be called synchronously by every processor.
    \param pstruct @b inout pointer to the struct array on expanded duel grid.
    \param front @b in	Pointer to Front.
    \param psize @b in size of data struct in byte.
 */
   IMPORT  void FT_ParallelExchGridStructArrayBuffer(POINTER pstruct ,
   				Front* front ,
				int psize);


/*! \fn void FT_ParallelExchGridIntArrayBuffer(int *iarray, Front *front)
 *  \ingroup PARALLEL
    \brief This is a parallel communication function for a integer array
     on the expanded dual grid of the grid_intfc in front. It will cut 
     the old buffer parts of the array and patch it with new buffer parts 
     received from other subdomains or periodically shifted sides. This is a 
     synchronous function and must be called synchronously by every processor.
    \param iarray @b inout A double array of variable on expanded duel grid.
    \param front @b in	Pointer to Front.
 */
   IMPORT  void FT_ParallelExchGridIntArrayBuffer(int *iarray ,
   				Front* front);

/*! \fn void FT_ParallelExchGridVectorArrayBuffer(double **vec_grid_array, Front *front)
 *  \ingroup PARALLEL
    \brief This is a parallel communication function for a vector double array
     on the expanded dual grid of the grid_intfc in front. It will cut 
     the old buffer parts of the array and patch it with new buffer parts 
     received from other subdomains or periodically shifted sides. This is a 
     synchronous function and must be called synchronously by every processor.
     Reflection of vector at the RLECTION_BOUNDARY is considered.
    \param vec_grid_array @b inout A vector double array of variable on expanded duel grid.
    \param front @b in	Pointer to Front.
 */
   IMPORT  void FT_ParallelExchGridVectorArrayBuffer(
				double **vec_grid_array ,
   				Front* front);

/*! \fn void FT_ParallelExchCompGridArrayBuffer(double *grid_array, Front *front, int *symmetry)
 *  \ingroup PARALLEL
    \brief This is a parallel communication function for a double array
     on the expanded comp grid of the grid_intfc in front. It will cut 
     the old buffer parts of the array and patch it with new buffer parts 
     received from other subdomains or periodically shifted sides. This is a 
     synchronous function and must be called synchronously by every processor.
    \param grid_array @b inout A double array of variable on expanded comp grid.
    \param front @b in Pointer to Front.
    \param symmetry @b in Flag of symmetry in each direction (in reflection).
 */
   IMPORT  void FT_ParallelExchCompGridArrayBuffer(double *grid_array ,
   				Front* front ,
				int* symmetry);

/*! \fn void FT_ParallelExchCellIndex(Front *front, int *lbuf, int *ubuf, POINTER ijk_to_I)
 *  \ingroup PARALLEL
    \brief This is a parallel communication function for the cell index
     on the expanded dual grid of the grid_intfc in front. The cell index
     translate the nD (n=2,3) icoordinates to a one dimensional index
     sequence. The indices are parallely globalized.
    \param front @b in	Pointer to Front.
    \param lbuf @b in size of buffer on the lower side.
    \param ubuf @b in size of buffer on the upper side.
    \param ijk_to_I @b inout Pointer to array of indices on the expanded dual grid (ij_to_I for 2D).
 */
   IMPORT  void FT_ParallelExchCellIndex(Front* front ,
				int *lbuf,
				int *ubuf,
				POINTER ijk_to_I);

/*! \fn void FT_ParallelExchCompGridCellIndex(Front *front, int *lbuf, int *ubuf, POINTER ijk_to_I)
 *  \ingroup PARALLEL
    \brief This is a parallel communication function for the dual cell index
     on the expanded comp grid of the comp_grid_intfc in front. The cell index
     translate the nD (n=2,3) icoordinates to a one dimensional index
     sequence. The indices are parallely globalized.
    \param front @b in	Pointer to Front.
    \param lbuf @b in size of buffer on the lower side.
    \param ubuf @b in size of buffer on the upper side.
    \param ijk_to_I @b inout Pointer to array of indices on the expanded dual grid (ij_to_I for 2D).
 */
   IMPORT  void FT_ParallelExchCompGridCellIndex(Front* front ,
				int *lbuf,
				int *ubuf,
				POINTER ijk_to_I);

/*! \fn INTERFACE *FT_CollectHypersurfFromSubdomains(Front *front, int *owner, int w_type)
 *  \ingroup PARALLEL
    \brief This function collect pieces of surface given wave type and 
     patch them to the surface in the owner subdomain. The owner subdomain
     is given by its rectangular subdomain index owner[] = (ix,iy,iz).
     The collected surfaces are packed in a interface structure.
     The function return the interface if the collection is succsessful.
     The function return NULL if the collection is not succsessful.
    \param front @b inout	Pointer to Front.
    \param owner @b in	Owner subdomain coordinates.
    \param w_type @b in	Wave type of hypersurface to be collected.
 */
   IMPORT  INTERFACE *FT_CollectHypersurfFromSubdomains(Front* front,
				int *owner,
				int w_type);

/*! \fn void FT_GetStatesAtPoint(POINT *p,  HYPER_SURF_ELEMENT *hse, HYPER_SURF *hs,  POINTER *sl, POINTER *sr)
 *  \ingroup FIELD
    \brief This function retrieves the left and right states at a point.
     Since a point can be shared by different entities, the associated
     hyper surface element (bond in 2D and tri in 3D) and hyper surface
     (curve in 2D and surface in 3D) must be provided as input.
    \param p @b in Pointer to point.
    \param hse @b in Pointer to the associated hyper surface element.
    \param hs @b in Pointer to the associated hyper surface.
    \param sl @b out Address for the left state.
    \param sr @b out Address for the right state.
 */
   IMPORT  void FT_GetStatesAtPoint(POINT *p,
				HYPER_SURF_ELEMENT *hse,
				HYPER_SURF *hs,
				POINTER *sl,
				POINTER *sr);

/*! \fn void FT_ScalarMemoryAlloc(POINTER *a, int size)
 *  \ingroup MEMORY
    \brief This function allocate the memory for a scalar.
    \param a @b out Address of the scalar pointer.
    \param size @b in size of the scalar entity.
 */
   IMPORT  void FT_ScalarMemoryAlloc(POINTER* a,
				int size);

/*! \fn void FT_VectorMemoryAlloc(POINTER *a, int n1, int size)
 *  \ingroup MEMORY
    \brief This function allocate the memory for a vector.
    \param a @b out Address of the vector pointer.
    \param n1 @b in Dimension of the vector.
    \param size @b in size of the vector entity.
 */
   IMPORT  void FT_VectorMemoryAlloc(POINTER* a,
				int n1,
				int size);

/*! \fn void FT_MatrixMemoryAlloc(POINTER *a, int n1, int n2, int size)
 *  \ingroup MEMORY
    \brief This function allocate the memory for a matrix.
    \param a @b out Address of the matrix pointer.
    \param n1 @b in First dimension of the matrix.
    \param n2 @b in Second dimension of the matrix.
    \param size @b in size of the matrix entity.
 */
   IMPORT  void FT_MatrixMemoryAlloc(POINTER* a,
				int n1,
				int n2,
				int size);

/*! \fn void FT_TriArrayMemoryAlloc(POINTER *a, int n1, int n2, int n3, int size)
 *  \ingroup MEMORY
    \brief This function allocate the memory for a tri-array.
    \param a @b out Address of the tri-array pointer.
    \param n1 @b in First dimension of the tri-array.
    \param n2 @b in Second dimension of the tri-array.
    \param n3 @b in Third dimension of the tri-array.
    \param size @b in size of the tri-array entity.
 */
   IMPORT  void FT_TriArrayMemoryAlloc(POINTER* a,
				int n1,
				int n2,
				int n3,
				int size);

/*! \fn void FT_QuadArrayMemoryAlloc(POINTER *a, int n1, int n2, int n3, int n4, int size)
 *  \ingroup MEMORY
    \brief This function allocate the memory for a quad-array.
    \param a @b out Address of the quad-array pointer.
    \param n1 @b in First dimension of the quad-array.
    \param n2 @b in Second dimension of the quad-array.
    \param n3 @b in Third dimension of the quad-array.
    \param n4 @b in Fourth dimension of the quad-array.
    \param size @b in size of the quad-array entity.
 */
   IMPORT  void FT_QuadArrayMemoryAlloc(POINTER* a,
				int n1,
				int n2,
				int n3,
				int n4,
				int size);

/*! \fn void FT_QuinArrayMemoryAlloc(POINTER *a, int n1, int n2, int n3, int n4, int n5, int size)
 *  \ingroup MEMORY
    \brief This function allocate the memory for a quin-array.
    \param a @b out Address of the quin-array pointer.
    \param n1 @b in First dimension of the quin-array.
    \param n2 @b in Second dimension of the quin-array.
    \param n3 @b in Third dimension of the quin-array.
    \param n4 @b in Fourth dimension of the quin-array.
    \param n5 @b in Fifth dimension of the quin-array.
    \param size @b in size of the quin-array entity.
 */
   IMPORT  void FT_QuinArrayMemoryAlloc(POINTER* a,
				int n1,
				int n2,
				int n3,
				int n4,
				int n5,
				int size);

/*! \fn void FT_SexArrayMemoryAlloc(POINTER *a, int n1, int n2, int n3, int n4, int n5, int n6, int size)
 *  \ingroup MEMORY
    \brief This function allocate the memory for a sex-array.
    \param a @b out Address of the sex-array pointer.
    \param n1 @b in First dimension of the sex-array.
    \param n2 @b in Second dimension of the sex-array.
    \param n3 @b in Third dimension of the sex-array.
    \param n4 @b in Fourth dimension of the sex-array.
    \param n5 @b in Fifth dimension of the sex-array.
    \param n6 @b in Sixth dimension of the sex-array.
    \param size @b in size of the sex-array entity.
 */
   IMPORT  void FT_SexArrayMemoryAlloc(POINTER* a,
				int n1,
				int n2,
				int n3,
				int n4,
				int n5,
				int n6,
				int size);

/*! \fn void FT_SetGlobalIndex(Front *front)
 *  \ingroup INITIALIZATION
    \brief This function set global index for all point, curves, surfaces 
     in the interface.
    \param n @b in Number of items whose memory are to be freed.
    \param ... @b Pointers to the addresses of these items.
 */
   IMPORT  void FT_SetGlobalIndex(Front *front);

/*! \fn void FT_SetSurfGlobalIndex(Front *front)
 *  \ingroup INITIALIZATION
    \brief This function set global index for all surfaces in the interface.
    \param n @b in Number of items whose memory are to be freed.
    \param ... @b Pointers to the addresses of these items.
 */
   IMPORT  void FT_SetSurfGlobalIndex(Front *front);

/*! \fn void FT_SetCurveGlobalIndex(Front *front)
 *  \ingroup INITIALIZATION
    \brief This function set global index for all curves in the interface.
    \param n @b in Number of items whose memory are to be freed.
    \param ... @b Pointers to the addresses of these items.
 */
   IMPORT  void FT_SetCurveGlobalIndex(Front *front);

/*! \fn void FT_FreeThese(int n, ...)
 *  \ingroup MEMORY
    \brief This function free memory of items allocated by FT_...MemoryAlloc()
     functions. The number of arguments is flexible, but needs to equal to
     the input integer n.
    \param front @b inout front in which the interfac global index is to be set.
 */
   IMPORT  void FT_FreeThese(int n,...);

/*TMP will move when mature*/
IMPORT  boolean FT_StateStructAtGridCrossing2(Front *front ,
                                int *icoords ,
                                GRID_DIRECTION  dir ,
                                int  comp ,
                                POINTER *state ,
                                HYPER_SURF **hs ,
                                HYPER_SURF_ELEMENT **hse ,
                                double *crx_coords );

/*! \fn double FT_ComputeTotalVolumeFraction(Front *front, COMPONENT comp_of_vol)
 *  \ingroup GRIDINTFC
    \brief This function compute the total volume fraction associated the
     comp_of_vol based on crossing in front->grid_intfc. It returns the
     total volume fraction as a double precision value.
    \param front @b pointer to the front.
    \param comp_of_vol @b component index of the volume fraction.
 */
   IMPORT  double FT_ComputeTotalVolumeFraction(
				Front *front,
				COMPONENT comp_of_vol);

/*! \fn double FT_ComputeGridVolumeFraction(Front *front, COMPONENT comp_of_vol,
 POINTER *grid_vol_frac)
 *  \ingroup GRIDINTFC
    \brief This function compute the volume fraction on the expanded dual grid
     associated the comp_of_vol based on crossing in front->grid_intfc. 
     It passes the address of the grid volume fraction to the pointer
     grid_vol_frac. In 2D, it is a double**, in 3D it is double***.
    \param front @b pointer to the front.
    \param comp_of_vol @b component index of the volume fraction.
    \param grid_vol_frac @b ananomous pointer pointing to the addess of the volume fraction on the grid.
 */
   IMPORT  void FT_ComputeGridVolumeFraction(
				Front *front,
				COMPONENT comp_of_vol,
				POINTER *grid_vol_frac);

/*! \fn double FT_CurveSegLengthConstr(CURVE *c, BOND *bs, BOND *be, int nb, double seg_length, REDISTRIBUTION_DIRECTION dir)
 *  \ingroup OPTIMIZATION
    \brief This function set curve to constrained length (seg_length), starting
     from the bond bs and end at bond be, with total of nb bonds. The cutting
     direction is the input dir.
    \param curve @b pointer to the curve.
    \param bs @b pointer to the start bond.
    \param be @b pointer to the end bond.
    \param nb @b number of bonds.
    \param seg_length @b total length of the curve.
    \param dir @b direction of the redistribution (FORWARD_REDISTRIBUTION/BACKWARD_REDISTRIBUTION.
 */
   IMPORT  void FT_CurveSegLengthConstr(
				CURVE *c, 
				BOND *bs, 
				BOND *be, 
				int nb, 
				double seg_length, 
				REDISTRIBUTION_DIRECTION dir);

/*! \fn CURVE *FT_MakeNodeArrayCurve(Front *front, int num_nodes, double **node_array, COMPONENT neg_comp, COMPONENT pos_comp, boolean is_closed_curve, scale_factor, int w_type)
 *  \ingroup INSERT
    \brief This function inserts a curve into the front with given
     array of points, if is_closed_curve is true, the curve is closed.
    \param front @b inout Pointer to the front in which curve is inserted.
    \param num_nodes @b in number of nodal points of the curve.
    \param node_array @b in an array of node coordinates.
    \param neg_comp @b in index for negative side of the curve (inner side).
    \param pos_comp @b in index for positive side of the curve (outer side).
    \param is_closed_curve @b in indicate whether the curve is closed.
    \param scale_factor @b in length factor of bond to grid spacing.
    \param w_type @b in wave type of the curve.
 */

   IMPORT  CURVE *FT_MakeNodeArrayCurve(Front *front,int num_nodes,double **node_array,COMPONENT neg_comp,COMPONENT pos_comp,boolean is_closed_curve,double scale_factor, int w_type);

/*! \fn CURVE *FT_MakePointArrayCurve(Front *front, int num_points, double **point_array, COMPONENT neg_comp, COMPONENT pos_comp, boolean is_closed_curve)
 *  \ingroup INSERT
    \brief This function inserts a curve into the front with given
     array of points, if is_closed_curve is true, the curve is closed.
    \param front @b inout Pointer to the front in which curve is inserted.
    \param num_points @b in number of points in the given array.
    \param point_array @b in an array of point coordinates.
    \param neg_comp @b in index for negative side of the curve (inner side).
    \param pos_comp @b in index for positive side of the curve (outer side).
    \param is_closed_curve @b in indicate whether the curve is closed.
    \param w_type @b in wave type of the curve.
 */

   IMPORT  CURVE *FT_MakePointArrayCurve(Front *front,int num_points,double **point_array,COMPONENT neg_comp,COMPONENT pos_comp,boolean is_closed_curve,int w_type);

/*! \fn void FT_MakeEllipticSurf(Front *front, double *center, double *radius, COMPONENT neg_comp, COMPONENT pos_comp, int w_type,int refinement_level,SURFACE **surf)
 *  \ingroup INSERT
    \brief This function inserts an elliptic surface into the front with given
     information of center, radii, components, and wave type.
    
    \param front @b inout Pointer to the front in which surface is inserted.
    \param center @b in center of the ellipsoid.
    \param radius @b in radii of the ellipsoid.
    \param neg_comp @b in index for negative side of the surface (inner side).
    \param pos_comp @b in index for positive side of the surface (outer side).
    \param w_type @b int wave type of the surface.
    \param refinement_level @b int refinement level of the surface.
    \param surf @b out surface made by this function.
 */

   IMPORT  void FT_MakeEllipticSurf(Front *front,double *center,double *radius,COMPONENT neg_comp,COMPONENT pos_comp,int w_type,int refinement_level,SURFACE **surf);

/*! \fn void FT_MakeSphericalSurf(Front *front, double *center, double radius, COMPONENT neg_comp, COMPONENT pos_comp, int w_type,int refinement_level,SURFACE **surf)
 *  \ingroup INSERT
    \brief This function inserts an spherical surface into the front with given
     information of center, radii, components, and wave type.
    
    \param front @b inout Pointer to the front in which surface is inserted.
    \param center @b in center of the sphere.
    \param radius @b in radius of the sphere.
    \param neg_comp @b in index for negative side of the surface (inner side).
    \param pos_comp @b in index for positive side of the surface (outer side).
    \param w_type @b int wave type of the surface.
    \param refinement_level @b int refinement level of the surface.
    \param surf @b out surface made by this function.
 */

   IMPORT  void FT_MakeSphericalSurf(Front *front,double *center,double radius,COMPONENT neg_comp,COMPONENT pos_comp,int w_type,int refinement_level,SURFACE **surf);

/*! \fn void FT_MakeDumbBellSurf(Front *front, double x0, double x1,double y0,double z0,double R,double r,COMPONENT neg_comp, COMPONENT pos_comp, int w_type,SURFACE **surf)
 *  \ingroup INSERT
    \brief This function inserts a dumbbell surface into the front with given
     information of its parameters, components, and wave type.
    
    \param front @b inout Pointer to the front in which surface is inserted.
    \param x0 @b in x-coordinate of the left center.
    \param x0 @b in x-coordinate of the right center.
    \param y0 @b in y-coordinate of the axis.
    \param z0 @b in z-coordinate of the axis.
    \param R @b in radius of the two end spheres.
    \param r @b in radius cylinder connecting the two spheres.
    \param neg_comp @b in index for negative side of the surface (inner side).
    \param pos_comp @b in index for positive side of the surface (outer side).
    \param w_type @b in wave type of the surface.
    \param surf @b out surface made by this function.
 */

   IMPORT  void FT_MakeDumbBellSurf(Front *front,double x0,double x1,double y0,double z0,double R,double r,COMPONENT neg_comp,COMPONENT pos_comp,int w_type,SURFACE **surf);

/*! \fn void FT_MakeProjectileSurf(Front *front, double *center, double R,double r,double h,COMPONENT neg_comp, COMPONENT pos_comp, int w_type,SURFACE **surf)
 *  \ingroup INSERT
    \brief This function inserts a projectile surface into the front with given
     information of its parameters, components, and wave type.
    
    \param front @b inout Pointer to the front in which surface is inserted.
    \param center @b in center-coordinate of the projectile.
    \param R @b in cylindrical radius of the projectile.
    \param r @b in head height of the projectile.
    \param h @b in butt height of the projectile.
    \param neg_comp @b in index for negative side of the surface (inner side).
    \param pos_comp @b in index for positive side of the surface (outer side).
    \param w_type @b in wave type of the surface.
    \param surf @b out surface made by this function.
 */

   IMPORT  void FT_MakeProjectileSurf(Front *front,double *center,double R,double r,double h,COMPONENT neg_comp,COMPONENT pos_comp,int w_type,SURFACE **surf);

/*! \fn CURVE *FT_MakeParametricCurve(Front *front, COMPONENT neg_comp, COMPONENT pos_comp, int w_type, boolean (*func)(POINTER,double,double*), POINTER func_params, int refinement_level, boolean is_closed)
 *  \ingroup INSERT
    \brief This function inserts a parametric curve into the front with given
     information of its function, parameters, components, and wave type.
    \param front @b inout Pointer to the front in which curve is inserted.
    \param neg_comp @b in index for negative side of the curve (left side).
    \param pos_comp @b in index for positive side of the curve (right side).
    \param w_type @b in wave type of the curve.
    \param func @b in parametric function for the curve.
    \param func_params @b in anonymous pointer to parameters for the function.
    \param refinement_level @b in refinement level on grid spacing.
    \param is_closed @b in whether the curve is closed.
 */

   IMPORT  CURVE *FT_MakeParametricCurve(Front *front, COMPONENT neg_comp, COMPONENT pos_comp, int w_type, boolean (*func)(POINTER,double,double*), POINTER func_params, int refinement_level, boolean is_closed);


/*! \fn void FT_RotateSurface(SURFACE *surf,double *center,double phi,double theta)
 *  \ingroup INSERT
    \brief This function rorate surface with azimuthal angle theta and 
     polar angle phi about the given center.
    
    \param surf @b inout Pointer to the surface to be rotated.
    \param center @b in center-coordinate for the rotation.
    \param phi @b in polar angle of the rotation.
    \param theta @b in azimuthal angle of the rotation.
 */

   IMPORT  void FT_RotateSurface(SURFACE *surf,double *center,double phi,double theta);

/*! \fn void FT_MakeCuboidSurf(Front *front,double *center,double *edge,COMPONENT neg_comp,COMPONENT pos_comp,int w_type,SURFACE **surf)
 *  \ingroup INSERT
    \brief This function inserts a cuboid surface into the front with given
     information of its parameters, components, and wave type.    
    \param front @b inout Pointer to the front in which surface is inserted.
    \param center @b in center of the cuboid.
    \param edge @b in edge of the cuboid.
    \param neg_comp @b in index for negative side of the surface (inner side).
    \param pos_comp @b in index for positive side of the surface (outer side).
    \param w_type @b in wave type of the surface.
    \param surf @b out surface made by this function
*/

    IMPORT void FT_MakeCuboidSurf(Front *front,double *center,double *edge,COMPONENT neg_comp,COMPONENT pos_comp,int w_type,SURFACE **surf);

/*! \fn void FT_MakeCylinderSurf(Front *front,double *center,double radius, double height, COMPONENT neg_comp,COMPONENT pos_comp,int w_type,SURFACE **surf)
 *  \ingroup INSERT
    \brief This function inserts a cylinder surface into the front with given
     information of its parameters, components, and wave type.
    \param front @b inout Pointer to the front in which surface is inserted.
    \param center @b in center of the cylinder.
    \param edge @b in radius of the cylinder.
    \param height @b in height of the cylinder.
    \param neg_comp @b in index for negative side of the surface (inner side).
    \param pos_comp @b in index for positive side of the surface (outer side).
    \param w_type @b in wave type of the surface.
    \param surf @b out surface made by this function
*/

    IMPORT void FT_MakeCylinderSurf(Front *front,double *center,double radius, double height, COMPONENT neg_comp,COMPONENT pos_comp,int w_type,SURFACE **surf);

/*! \fn void FT_MakeConeSurf(Front *front,double *center,double slope, double height, COMPONENT neg_comp,COMPONENT pos_comp,int w_type,SURFACE **surf)
 *  \ingroup INSERT
    \brief This function inserts a cone surface into the front with given
     information of its parameters, components, and wave type.
    \param front @b inout Pointer to the front in which surface is inserted.
    \param center @b in vertex of the cone.
    \param slope @b in slope of the cone.
    \param height @b in height of the cone.
    \param neg_comp @b in index for negative side of the surface (inner side).
    \param pos_comp @b in index for positive side of the surface (outer side).
    \param w_type @b in wave type of the surface.
    \param surf @b out surface made by this function
*/

    IMPORT void FT_MakeConeSurf(Front *front,double *center,double slope, double height, COMPONENT neg_comp,COMPONENT pos_comp,int w_type,SURFACE **surf);

/*! \fn void FT_MakeTetrahedronSurf(Front *front,double *center,double radius, COMPONENT neg_comp,COMPONENT pos_comp,int w_type,SURFACE **surf)
 *  \ingroup INSERT
    \brief This function inserts a tetrahedron surface into the front with given
     information of its parameters, components, and wave type.
    \param front @b inout Pointer to the front in which surface is inserted.
    \param center @b in center of the tetrahedron.
    \param radius @b in circumcircle of the tetrahedron.
    \param neg_comp @b in index for negative side of the surface (inner side).
    \param pos_comp @b in index for positive side of the surface (outer side).
    \param w_type @b in wave type of the surface.
    \param surf @b out surface made by this function
*/
  
    IMPORT void FT_MakeTetrahedronSurf(Front *front,double *center,double radius,COMPONENT neg_comp,COMPONENT pos_comp,int w_type,SURFACE **surf);

/*! \fn void FT_MakePlaneSurf(Front *front, double *plane_nor,double *plane_pt,boolean reset_bdry_comp,COMPONENT neg_comp, COMPONENT pos_comp, int w_type,SURFACE **surf)
 *  \ingroup INSERT
    \brief This function inserts a plane surface into the front with given
     information of its parameters, components, and wave type.
    
    \param front @b inout Pointer to the front in which surface is inserted.
    \param plane_nor @b in normal vector of the plane.
    \param plane_pt @b in coordinates of a point on the plane.
    \param reset_bdry_comp @b in if YES, reset boundary component.
    \param neg_comp @b in index for negative side of the surface (inner side).
    \param pos_comp @b in index for positive side of the surface (outer side).
    \param w_type @b in wave type of the surface.
    \param surf @b out surface made by this function.
 */

   IMPORT  void FT_MakePlaneSurf(Front *front,double *plane_nor,double *plane_pt,boolean reset_bdry_comp,COMPONENT neg_comp,COMPONENT pos_comp,int w_type,SURFACE **surf);

/*! \fn void FT_InstallSurfEdge(SURFACE *surf,int hsbdry_type)
 *  \ingroup INSERT
    \brief This function install a curve at the surface boundary. The
     boundary of the surface should contain triangles with NULL side.
    
    \param surf @b inout Pointer to the surface where edge to be installed.
    \param hsbdry_type @b in Boundary curve type.
 */

   IMPORT  void FT_InstallSurfEdge(SURFACE *surf,int hsbdry_type);

/*! \fn void FT_CutSurfBdry(SURFACE *surf,boolean constr_func(POINTER,double*),POINTER func_params,double **insert_coords,int num_pts,int insert_idir)
 *  \ingroup INSERT
    \brief This function cut surface at the boundary defined by the constrain
     function. It first insert significant corner points to make sure it 
     cuts the surface accurately.
    
    \param surf @b inout Pointer to the surface to be cut.
    \param constr_func @b in constrain function defining the boundary.
    \param func_params @b in anonymous pointer of function parameters.
    \param insert_coords @b in coordinates of points to be inserted before cut.
    \param num_pts @b in number of points to be inserted before cut.
    \param insert_idir @b in insert point along this direction.
 */

   IMPORT  void FT_CutSurfBdry(SURFACE *surf,boolean constr_func(POINTER,double*),POINTER func_params,double **insert_coords,int num_pts,int insert_idir);

/*! \fn void FT_MakeEllipticCurve(Front *front, double *center, double *radius, COMPONENT neg_comp, COMPONENT pos_comp, int w_type,int refinement_level,CURVE **curve)
 *  \ingroup INSERT
    \brief This function inserts an elliptic curve into the front with given
     information of center, radii, components, and wave type.
    
    \param front @b inout Pointer to the front in which curve is inserted.
    \param center @b in center of the ellipse.
    \param radius @b in radii of the ellipse.
    \param neg_comp @b in index for negative side of the curve (inner side).
    \param pos_comp @b in index for positive side of the curve (outer side).
    \param w_type @b int wave type of the curve.
    \param refinement_level @b int level of refinement.
    \param curve @b out curve made by this function.
 */

   IMPORT  void FT_MakeEllipticCurve(Front *front,double *center,double *radius,COMPONENT neg_comp,COMPONENT pos_comp,int w_type,int refinement_level,CURVE **curve);

/*! \fn void FT_PrintWaveType(int w_type)
 *  \ingroup INFO
    \brief This function print wave type as a string.
    \param w_type @b in Wave type as enumerate number.
 */

   IMPORT  void FT_PrintWaveType(int w_type);

/*! \fn void FT_PrintBoundaryType(int dir, int side)
 *  \ingroup INFO
    \brief This function print boundary type as a string.
    \param dir @b in Direction of the boundary.
    \param side @b in Side of the boundary.
 */

   IMPORT  void FT_PrintBoundaryType(int dir,int side);

/*! \fn int FT_BoundaryType(int dir, int side)
 *  \ingroup QUERY
    \brief This function return boundary type as an enumerated.
    \param dir @b in Direction of the boundary.
    \param side @b in Side of the boundary.
 */

   IMPORT  int FT_BoundaryType(int dir,int side);

/*! \fn double *FT_GridIntfcTopL(Front*)
 *  \ingroup GRIDINTFC
    \brief This function return lower bounds of grid domain.
    \param front @b in Pointer to front.
 */

   IMPORT  double *FT_GridIntfcTopL(Front*);

/*! \fn double *FT_GridIntfcTopU(Front*)
 *  \ingroup GRIDINTFC
    \brief This function return upper bounds of grid domain.
    \param front @b in Pointer to front.
 */

   IMPORT  double *FT_GridIntfcTopU(Front*);

/*! \fn double *FT_GridIntfcToph(Front*)
 *  \ingroup GRIDINTFC
    \brief This function return grid spacing of grid domain.
    \param front @b in Pointer to front.
 */

   IMPORT  double *FT_GridIntfcToph(Front*);

/*! \fn COMPONENT *FT_GridIntfcTopComp(Front*)
 *  \ingroup GRIDINTFC
    \brief This function return components of grid domain.
    \param front @b in Pointer to front.
 */

   IMPORT  COMPONENT *FT_GridIntfcTopComp(Front*);

/*! \fn int *FT_GridIntfcTopGmax(Front*)
 *  \ingroup GRIDINTFC
    \brief This function return mesh sizes of grid domain.
    \param front @b in Pointer to front.
 */

   IMPORT  COMPONENT *FT_GridIntfcTopGmax(Front*);

/*! \fn RECT_GRID *FT_GridIntfcTopGrid(Front*)
 *  \ingroup GRIDINTFC
    \brief This function return grid structure of grid domain.
    \param front @b in Pointer to front.
 */

   IMPORT  RECT_GRID *FT_GridIntfcTopGrid(Front*);

/*! \fn boolean FT_NextTopGridIcoordsInDir(Front*,int *icoords,GRID_DIRECTION dir,int *icoords_next)
 *  \ingroup GRIDINTFC
    \brief This function find the next coordinate index while standing at the
    input icoords and looking to the direction dir. It returns NO if the
    the input icoords is already at the edge of the domain.
    \param front @b in Pointer to front.
    \param icoords @b in input coordinate index.
    \param dir @b in direction to which of the next grid point.
    \param icoords_next @b out the next coordinate index.
 */

   IMPORT  boolean FT_NextTopGridIcoordsInDir(
   	Front*,
	int *icoords,
	GRID_DIRECTION dir,
	int *icoords_next);

/*! \fn boolean FT_AdjTopGridIcoords(Front*,int *icoords,GRID_DIRECTION *dir,int *icoords_next)
 *  \ingroup GRIDINTFC
    \brief This function find the adjacent coordinate index while standing at 
    the input icoords and looking to all directions dir. It returns NO if the
    the input icoords is already at the edge of the domain. The number of
    directions depends on the dimension. 
    \param front @b in Pointer to front.
    \param icoords @b in input coordinate index.
    \param dir @b in directions to which of adjacent next grid point.
    \param icoords_next @b out the next coordinate index.
 */

   IMPORT  boolean FT_AdjTopGridIcoords(
   	Front*,
	int *icoords,
	GRID_DIRECTION *dir,
	int *icoords_next);

/*! \fn boolean FT_NearestGridIcoordsWithComp(Front*,int *icoords,COMPONENT comp,COMPONENT comp, int *gmax, int *icoords_next)
 *  \ingroup GRIDINTFC
    \brief This function find the nearest coordinate index while standing at 
    the input icoords and looking to all directions dir. It returns NO if the
    the input icoords is already at the edge of the domain. 
    \param front @b in Pointer to front.
    \param icoords @b in input coordinate index.
    \param comp @b in component of the cell to be searched.
    \param comp_map @b in the coomponent map of the grid.
    \param gmax @b in upper index bounds in all directions.
    \param icoords_next @b out the next coordinate index.
 */

   IMPORT  boolean FT_NearestGridIcoordsWithComp(
	Front *front,
	int *icoords,
        COMPONENT comp,
        COMPONENT *comp_map,
        int *gmax,
        int *icoords_next);


/*! \fn void FT_AddHdfMovieVariable(
 * 	Front *front,
        boolean preset_bound,
        boolean untracked,
        COMPONENT obst_comp,
        const char *var_name,
        double *var_field,
        double (*getStateFunc)(POINTER),
        double max_var,
        double min_var)
 *  \ingroup INITIALIZATION
    \brief Initialize a variable and information for hdf movie output.
    \param front @b in	Pointer to the front  
    \param preset_bound @b in Flag whether to set bounds of variable
    \param untracked @b in Flag whether to untrack front  
    \param obst_comp @b in The obstacle component
    \param var_name @b in Names of the variable
    \param idir @b in Normal direction of the plane of the variable (in 3D)
    \param var_field @b in Field values of the variable
    \param getStateFunc @b in Function to get variable from interface state
    \param max_var @b in  Ceilling of the variable (if use preset_bound)
    \param min_var @b in  Floor of the variable (if use preset_bound)
 */
   IMPORT  void FT_AddHdfMovieVariable(
	Front *front,
        boolean preset_bound,
        boolean untracked,
        COMPONENT obst_comp,
        const char *var_name,
	int idir,
        double *var_field,
        double (*getStateFunc)(POINTER),
        double max_var,
        double min_var);

/*! \fn void FT_AddVtkVectorMovieVariable(
 * 	Front *front,
        const char *var_name,
        double **var_field)
 *  \ingroup INITIALIZATION
    \brief Initialize a variable and information for vtk vector movie output.
    \param front @b in	Pointer to the front  
    \param var_name @b in Names of the variable
    \param var_field @b in Field values of the variable
 */
   IMPORT  void FT_AddVtkVectorMovieVariable(
	Front *front,
        const char *var_name,
        double **var_field);

/*! \fn void FT_AddVtkScalarMovieVariable(
 * 	Front *front,
        const char *var_name,
        double *var_field)
 *  \ingroup INITIALIZATION
    \brief Initialize a variable and information for vtk scalar movie output.
    \param front @b in	Pointer to the front  
    \param var_name @b in Names of the variable
    \param var_field @b in Field values of the variable
 */
   IMPORT  void FT_AddVtkScalarMovieVariable(
	Front *front,
        const char *var_name,
        double *var_field);

/*! \fn void FT_AddVtkIntfcMovieVariable(
 * 	Front *front,
        const char *var_name)
 *  \ingroup INITIALIZATION
    \brief Initialize the variable name for vtk interface movie output.
    \param front @b in	Pointer to the front  
    \param var_name @b in Names of the variable
 */
   IMPORT  void FT_AddVtkIntfcMovieVariable(
	Front *front,
        const char *var_name);

/*! \fn void FT_ResetDomainAndGrid(Front *front, 
   	Front *front,
   	double *L,
   	double *U,
   	int *gmax)
    \ingroup INITIALIZATION
    \brief Reset computational grid of front to new domain defined by
     input L,U, and gmax (mesh size). This function is used for modification
     of restart run. Caution: parallelization has not been considered.
    \param front @b inout	Pointer to Front.
    \param L @b in New limit of lower boundary.
    \param U @b in New limit of upper boundary.
    \param gmax @b in New mesh size.
 */
   IMPORT  void FT_ResetDomainAndGrid(Front *front,
	double *L,
        double *U,
        int *gmax);

/*! \fn boolean FT_CoordsInSubdomain(Front *front, double *coords)
 *  \ingroup PARALLEL
    \brief This function test if the coordinate is inside my subdomain.
     The function returns YES if it is in, otherwise it returns NO.
    \param front @b in	Pointer to Front.
    \param coords @b in	Testing coordinate.
 */
   IMPORT  boolean FT_CoordsInSubdomain(Front* front,
				double *coords);

/*! \fn void FT_PrintTimeStamp(Front *front)
 *  \ingroup OUTPUT
    \brief Output time information including current time, step and predicted
     next time step size.
    \param front @b in	Pointer to Front.
 */

   IMPORT  void FT_PrintTimeStamp(Front *front);

/*! \fn void FT_MakeCrossCylinderSurf(Front *front,double *center1,double * center2, double radius1, double radius2,  double height1, double height2, COMPONENT neg_comp,COMPONENT pos_comp,int w_type,SURFACE **surf)
 *  \ingroup INSERT
    \brief This function inserts a surface in which a cylinder go through another cylinder into the front with given
     information of its parameters, components, and wave type.
    \param front @b inout Pointer to the front in which surface is inserted.
    \param center1 @b in center of one cylinder.
    \param center2 @b in center of another cylinder.
    \param radius1 @b in radius of one cylinder.
    \param radius2 @b in radius of another cylinder.
    \param height1 @b in height of one cylinder.
    \param height2 @b in height of another cylinder.
    \param neg_comp @b in index for negative side of the surface (inner side).
    \param pos_comp @b in index for positive side of the surface (outer side).
    \param w_type @b in wave type of the surface.
    \param surf @b out surface made by this function
*/
    IMPORT void FT_MakeCrossCylinderSurf(Front *front,
				double *center1,
				double *center2,
				double radius1,
				double radius2,
				double height1,
				double height2, 
				COMPONENT neg_comp,
				COMPONENT pos_comp,
				int w_type,
				SURFACE **surf);

/*! \fn void FT_MakeBowlSurf(Front *front,double *center,double radius1,double radius2,double radius3,double height1,double height2,COMPONENT neg_comp,COMPONENT pos_comp,int w_type,SURFACE **surf)
 *  \ingroup INSERT
    \brief This function inserts a bowl surface into the front with given
     information of its parameters, components, and wave type.
    \param front @b inout Pointer to the front in which surface is inserted.
    \param center @b in center of the sphere.
    \param radius1 @b in radius of the outer sphere.
    \param radius2 @b in radius of the cylinder.
    \param radius3 @b in radius of the inner sphere.
    \param height1 @b in height of body of the bowl.
    \param height2 @b in height of bottom of the bowl.
    \param neg_comp @b in index for negative side of the surface (inner side).
    \param pos_comp @b in index for positive side of the surface (outer side).
    \param w_type @b in wave type of the surface.
    \param surf @b out surface made by this function
*/
    IMPORT void FT_MakeBowlSurf(Front *front,
				double *center,
				double radius1,
				double radius2,
				double radius3,
				double height1,
				double height2, 
				COMPONENT neg_comp,
				COMPONENT pos_comp,
				int w_type,
				SURFACE **surf);

/*! \fn void FT_MakePlatformSurf(Front *front,double *center,double radius, double height, double slope, COMPONENT neg_comp,COMPONENT pos_comp,int w_type,SURFACE **surf)
 *  \ingroup INSERT
    \bridf This function inserts a platform surface derived by cutting the head of the cone curface into the front with given information of its parameters, components, and wave type.
    \param front @b inout Pointer to the front in which surface is inserted.
    \param center @b in center of the bottom of the platform.
    \param radius @b in radius of the bottom of the platform. 
    \param height @b in height of the platform.
    \param slope @b in slope of the platform.
    \param neg_comp @b in index for negative side of the surface (inner side).
    \param pos_comp @b in index for positive side of the surface (outer side).
    \param w_type @b in wave type of the surface.
    \param surf @b out surface made by this function
*/

    IMPORT void FT_MakePlatformSurf(Front *front,
			double *center,
			double radius, 
			double height, 
			double slope, 
			COMPONENT neg_comp,
			COMPONENT pos_comp,
			int w_type,
			SURFACE **surf);

/*! \fn void FT_MakeStellatedoctahedronSurf(Front *front, double *center,double edge,COMPONENT neg_comp, COMPONENT pos_comp, int w_type,SURFACE **surf)
 *  \ingroup INSERT
    \brief This function inserts a cristal surface into the front with given
     information of its parameters, components, and wave type.
    \param front @b inout Pointer to the front in which surface is inserted.
    \param center @b in center of the upper tetrahedron.
    \param edge @b in edge of the tetrahedron.
    \param neg_comp @b in index for negative side of the surface (inner side).
    \param pos_comp @b in index for positive side of the surface (outer side).
    \param w_type @b in wave type of the surface.
    \param surf @b out surface made by this function
*/

    IMPORT void FT_MakeStellatedOctahedronSurf(Front *front,
			double *center,
			double edge,
			COMPONENT neg_comp,
			COMPONENT pos_comp,
			int w_type,
			SURFACE **surf);

/*! \fn void FT_InitSurfVeloFunc(SURFACE *surf,const char* vfunc_name,POINTER vparams,int (*)vfunc(POINTER,struct _Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*));
    \ingroup INITIALIZATION
    \brief Initialize surface velocity function for propagation.
     The velocity function use point and other related structures as input,
     must also supply parameters needed for the velocity function.
    \param surf @b inout Pointer to the surface.
    \param vparams @b in Pointer to velocity function parameters.
    \param vfunc @b in Pointer to velocity function.
 */
   IMPORT  void    FT_InitSurfVeloFunc(
			SURFACE *surf,
			const char *vfunc_name,
			POINTER vparams,
			int (*vfunc)(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*));

/*! \fn void FT_InitCurveVeloFunc(CURVE *curve,const char* vfunc_name,POINTER vparams,int (*)vfunc(POINTER,struct _Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*));
    \ingroup INITIALIZATION
    \brief Initialize curve velocity function for propagation.
     The velocity function use point and other related structures as input,
     must also supply parameters needed for the velocity function.
    \param curve @b inout Pointer to the curve.
    \param vparams @b in Pointer to velocity function parameters.
    \param vfunc @b in Pointer to velocity function.
 */
   IMPORT  void    FT_InitCurveVeloFunc(
			CURVE *curve,
			const char *vfunc_name,
			POINTER vparams,
			int (*vfunc)(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*));

/*! \fn void FT_InitNodeVeloFunc(NODE *node,const char* vfunc_name,POINTER vparams,int (*)vfunc(POINTER,struct _Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*));
    \ingroup INITIALIZATION
    \brief Initialize node velocity function for propagation.
     The velocity function use point and other related structures as input,
     must also supply parameters needed for the velocity function.
    \param node @b inout Pointer to the node.
    \param vparams @b in Pointer to velocity function parameters.
    \param vfunc @b in Pointer to velocity function.
 */
   IMPORT  void    FT_InitNodeVeloFunc(
			NODE *node,
			const char *vfunc_name,
			POINTER vparams,
			int (*vfunc)(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*));


/*! \fn boolean FT_CheckSurfCompConsistency(Front *front, SURFACE *surf)
 *  \ingroup CHECK
    \brief This function check if the outer (positive) component of a 
     surface is consistent with the ambient component.
    \param front @b in Pointer to the front.
    \param surf @b in Pointer to surface to be checked.
 */

   IMPORT  boolean FT_CheckSurfCompConsistency(
			Front *front,
			SURFACE *surf);

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif
