#ifndef SSSHEADER_H
#define SSSHEADER_H

/***** External libraries ****************************************************/
#ifdef __cplusplus
extern "C" {
#endif
  /*
#define _XOPEN_SOURCE 600
#define _GNU_SOURCE
  */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <sys/types.h>
#include <ctype.h>
#include <dirent.h>
#include <sys/timeb.h>
#include <signal.h>
#include <unistd.h>

#ifndef NO_GFX_MODE
#include <gd.h>
#endif

/*****************************************************************************/

/***** Parallel calculation **************************************************/

#ifdef MPI

/* MPI library header */

#include <mpi.h>

/* mpirun executable path */

#define MPIRUN_PATH "mpirun"

#endif

#ifdef OPEN_MP

/* OpenMP library header */

#include <omp.h>

#endif

/* MPI mode: 1 = divide source size, 2 = divide number of active cycles */

#define MPI_MODE1
#define OLD_IFP

/*****************************************************************************/

/***** Code name and version *************************************************/

#define CODE_NAME     "Serpent"
#define CODE_VERSION  "2.1.26"
#define CODE_DATE     "April 21, 2016"
#define CODE_AUTHOR   "serpent@vtt.fi"

/*****************************************************************************/

/***** Constants *************************************************************/

/* Zero and infinity */

#define INFTY       1e+37
#define ZERO        1e-37

/* Natural constants */

#define M_NEUTRON   1.0086649670000E+00  /* neutron mass in amu              */
#define N_AVOGADRO  6.0220434469282E-01  /* Avogadro's constant in 1E-24     */
#define SHAKE       1E-8                 /* ancient time unit used at LANL   */
#define KELVIN      8.6173E-11           /* MeV per kelvin                   */
#define E_RESTMASS  0.5109989            /* Electron rest mass in MeV        */
#define MEV         1.6021765314E-13     /* J per MeV                        */
#define ETOV2       1.9131290731E+18     /* square speed (cm2/s2) per energy */
#define BARN        1E-24                /* barn in cm2                      */
#define FS_CONST    7.2973525698E-03     /* Fine-structure constant (~1/137) */
#define GRAYH       5.7678382E-07        /* MeV/s to Gy/h                    */

#ifdef wanhat

#define SPD_C       29979245800.0        /* Speed of light in cm/s           */
#define NEUTRON_E0  939.668854490302     /* Neutron rest mass in MeV         */

#else
#define SPD_C       29979250000.0        /* Speed of light in cm/s           */
#define NEUTRON_E0  9.3958000000000E+02  /* Neutron rest mass in MeV         */

#endif

/* U-235 fission Q-value and energy deposition per fission */

#define U235_FISSQ 1.9372E+02
#define U235_FISSE 202.27*MEV

/* Math constants */

#define LOG2    0.69314718055995
#define PI      3.14159265358979
#define SQRT2   1.41421356237309
#define SQRT3   1.73205080756888
#define SIN30   0.50000000000000
#define COS30   0.86602540378444
#define TAN30   0.57735026918963
#define SIN45   0.70710678118655
#define SIN60   0.86602540378444
#define COS60   0.50000000000000
#define SQRTPI  1.77245385090552

/* Random number stride between neutrons = 2^STRIDE (Fast implementations */
/* for STRIDE values 12-19 are implemented for use in serial mode. (TVi)) */

#define STRIDE 16

/* Extrapolation length for boundary distances */

#define EXTRAP_L 1E-6

/* very negative null pointer for data arrays */

#define NULLPTR -1E+6

/* yes and no */

#define YES 1
#define NO  0

/* Maximum array sizes */

#define MAX_STR                  256
#define MAX_FP_NUCLIDES          1500
#define MAX_ISOTOPES             1000
#define MAX_INPUT_PARAMS         1000000
#define MAX_CELL_SURFACES        10000
#define MAX_LATTICE_ITEMS        100000
#define MAX_SURFACE_PARAMS       300
#define MAX_EGRID_NE             100000000
#define MAX_PRECURSOR_GROUPS     8
#define MAX_GEOMETRY_LEVELS      10000
#define MAX_EXT_K_GEN            5
#define MAX_GENERATIONS          1000000000

/* Tracking errors */

#define TRACK_ERR_INF_LOOP     1
#define TRACK_ERR_CELL_SEARCH  2
#define TRACK_ERR_LATTICE      3
#define TRACK_ERR_NO_MATERIAL  4
#define TRACK_ERR_OUTSIDE      5

/* Maximum allowed number of OpenMP threads */

#define MAX_OMP_THREADS 1000

/* Limiting values for parameters (used for sanity checks only) */

#define MAX_XS       1E+12  /* Maximum microscopic cross section  */
#define E_CHECK_MAX  220.0  /* Maximum energy */

/* Energy grid types */

#define EG_INTERP_MODE_LIN 1
#define EG_INTERP_MODE_LOG 2

/* Run mode options */

#define MODE_NORMAL         0
#define MODE_REPLAY         2
#define MODE_MPI            4

/* Simulation modes */

#define SIMULATION_MODE_CRIT    1
#define SIMULATION_MODE_SRC     2
#define SIMULATION_MODE_DYN     3
#define SIMULATION_MODE_DELDYN  4

/* UFS mode */

#define UFS_MODE_NONE  0
#define UFS_MODE_COL   1
#define UFS_MODE_FLUX  2
#define UFS_MODE_FISS  3

/* K-eff iteration modes */

#define ITER_MODE_NONE    0
#define ITER_MODE_ALBEDO  1

/* Burnup modes (solution method for Bateman's equations) */

#define BUMODE_TTA  1
#define BUMODE_CRAM 2

/* Special run modes in coefficient calculation */

#define SPECIAL_COEF_MODE_HIS_ONLY 1
#define SPECIAL_COEF_MODE_COE_ONLY 2

/* Normalization in burnup mode */

#define BURN_NORM_ALL       1
#define BURN_NORM_BURN      2
#define BURN_NORM_NOT_BURN  3

/* Radioactive decay source sampling modes */

#define RAD_SRC_MODE_ANALOG    1
#define RAD_SRC_MODE_IMPLICIT  2

/* Array types */

#define DATA_ARRAY     1
#define RES1_ARRAY     2
#define ACE_ARRAY      3
#define PRIVA_ARRAY    4
#define BUF_ARRAY      5
#define RES2_ARRAY     6

/* Solid 3D geometry types */

#define SOLID_GEO_TYPE_OF   1
#define SOLID_GEO_TYPE_STL  2
#define SOLID_GEO_TYPE_IFC  3

/* STL geometry cell search mode */

#define STL_SEARCH_MODE_NONE  0
#define STL_SEARCH_MODE_FAST  1
#define STL_SEARCH_MODE_SAFE  2

/* STL fail flags */

#define STL_RAY_TEST_FAIL_PARA   -1000
#define STL_RAY_TEST_FAIL_EDGE   -2000
#define STL_RAY_TEST_FAIL_MESH   -3000
#define STL_RAY_TEST_FAIL_LOOP   -4000
#define STL_RAY_TEST_FAIL_STUCK  -5000
#define STL_FACET_OVERLAP        -6000

/* XS data types */

#define XS_TYPE_SAB         3
/*
#define XS_TYPE_NONE        0
#define XS_TYPE_CONTINUOUS  1
#define XS_TYPE_DOSIMETRY   2
#define XS_TYPE_DECAY       4
*/

/* THERM interpolation types */

#define THERM_INTERP_NONE      0
#define THERM_INTERP_STOCHMIX  1
#define THERM_INTERP_MAKXSF    2
#define THERM_INTERP_OTF       3

/* Calculation of analog reaction rates */

#define ARR_MODE_NONE  0
#define ARR_MODE_BALA  1
#define ARR_MODE_ALL   2

/* Nuclide types */

#define NUCLIDE_TYPE_NONE        0
#define NUCLIDE_TYPE_TRANSPORT   1
#define NUCLIDE_TYPE_DOSIMETRY   2
#define NUCLIDE_TYPE_SAB         3
#define NUCLIDE_TYPE_DECAY       4
#define NUCLIDE_TYPE_PHOTON      5
#define NUCLIDE_TYPE_DBRC        6
#define NUCLIDE_TYPE_TRANSMUXS   7

/* Nuclide flags */

#define NUCLIDE_FLAG_INITIAL              1
#define NUCLIDE_FLAG_AP                   2
#define NUCLIDE_FLAG_DP                   4
#define NUCLIDE_FLAG_FP                   8
#define NUCLIDE_FLAG_BP                  16
#define NUCLIDE_FLAG_TRANSPORT_DATA      32
#define NUCLIDE_FLAG_DOSIMETRY_DATA      64
#define NUCLIDE_FLAG_DECAY_DATA         128
#define NUCLIDE_FLAG_NFY_DATA           256
#define NUCLIDE_FLAG_SFY_DATA           512
#define NUCLIDE_FLAG_BRA_DATA          1024
#define NUCLIDE_FLAG_URES_AVAIL        2048
#define NUCLIDE_FLAG_URES_USED         4096
#define NUCLIDE_FLAG_FISSILE           8192
#define NUCLIDE_FLAG_SAB_DATA         16384
#define NUCLIDE_FLAG_SRC              32768
#define NUCLIDE_FLAG_PHOTON_DATA      65536
#define NUCLIDE_FLAG_DBRC            131072
#define NUCLIDE_FLAG_DEP             262144
#define NUCLIDE_FLAG_DELNU_PREC      524288
#define NUCLIDE_FLAG_TRANSMU_DATA   1048576
#define NUCLIDE_FLAG_NEW_DAUGHTERS  2097152
#define NUCLIDE_FLAG_TMS            4194304

/* Branching data types */

#define BRA_TYPE_NONE  0
#define BRA_TYPE_FIX   1
#define BRA_TYPE_ENE   2

/* Reaction types */

#define REACTION_TYPE_SUM         1
#define REACTION_TYPE_PARTIAL     2
#define REACTION_TYPE_SPECIAL     3
#define REACTION_TYPE_DECAY       4
#define REACTION_TYPE_TRA_BRANCH  5
#define REACTION_TYPE_DEC_BRANCH  6

/* Macroscopic reaction MT's (also used with detector response functions) */

#define MT_MACRO_TOTXS           -1
#define MT_MACRO_ABSXS           -2
#define MT_MACRO_ELAXS           -3
#define MT_MACRO_FISSXS          -6
#define MT_MACRO_HEATXS          -4
#define MT_MACRO_PHOTXS          -5
#define MT_MACRO_NSF             -7
#define MT_MACRO_FISSE           -8
#define MT_MACRO_MAJORANT        -9
#define MT_MACRO_RECOILE        -10
#define MT_SOURCE_RATE          -11

#define MT_NEUTRON_DENSITY      -15
#define MT_MACRO_INLPRODXS      -16

#define MT_MACRO_TOTPHOTXS      -25
#define MT_MACRO_HEATPHOTXS     -26
#define MT_PHOTON_PULSE_HEIGHT  -27

#define MT_MACRO_TMP_MAJORANTXS -30

#define MT_USER_DEFINED        -100
#define MT_PHOTON_DOSE         -200

/* Pre-defined photon attenuation coefficients -201 ... -248 */

/* Detector output conversion */

#define DET_CONVERT_NONDE  0
#define DET_CONVERT_GDOSE  1

/* Energy grid types */

#define GRID_TYPE_LIN  1
#define GRID_TYPE_LOG  2

/* User energy grid types */

#define EG_TYPE_ARB     1
#define EG_TYPE_UNI_E   2
#define EG_TYPE_UNI_L   3
#define EG_TYPE_PREDEF  4

/* User time binning types */

#define TB_TYPE_ARB       1
#define TB_TYPE_UNI_T     2
#define TB_TYPE_UNI_LOGT  3

/* Angular distribution types */

#define ANG_TYPE_EQUIBIN 1
#define ANG_TYPE_TABULAR 2

/* General options */

#define OPT_USED               1
#define OPT_EXISTS             2
#define OPT_BURN_MAT           4
#define OPT_FISSILE_MAT        8
#define OPT_PHYSICAL_MAT      16
#define OPT_REPLACED_MAT      32
#define OPT_INCLUDE_MAJORANT  64

/* Sort modes */

#define SORT_MODE_ASCEND         1
#define SORT_MODE_DESCEND        2
#define SORT_MODE_ASCEND_PRIVA   3
#define SORT_MODE_DESCEND_PRIVA  4

/* Fission yield types (used to identify fission reactions in   */
/* transmu chains, these must be different from all ENDF MT's). */

#define FISSION_YIELD_TYPE_NFY  -1
#define FISSION_YIELD_TYPE_SFY  -2

/* Data sizes in bytes */

#define KILO        1024.0
#define MEGA     1048576.0
#define GIGA  1073741824.0

/* Input parameter types */

#define PTYPE_LOGICAL 1
#define PTYPE_REAL    2
#define PTYPE_INT     3

/* Plotter modes */

#define PLOT_MODE_YZ 1
#define PLOT_MODE_XZ 2
#define PLOT_MODE_XY 3

/* Track plot tail colours */

#define TRACK_PLOT_NCOL 16

/* Boundary conditions */

#define BC_BLACK       1
#define BC_REFLECTIVE  2
#define BC_PERIODIC    3

/* Transformation types */

#define TRANS_TYPE_UNI   1
#define TRANS_TYPE_SURF  2
#define TRANS_TYPE_FILL  3

/* Discontinuity factor region types */

#define ADF_REG_TYPE_FUEL 1
#define ADF_REG_TYPE_REF  2

/* Relative width of ADF mid and corner regions */

#define ADF_MID_WIDTH  0.1
#define ADF_CORN_WIDTH 0.1

/* Surface list operators */

#define SURF_OP_NOT    -1
#define SURF_OP_AND    -2
#define SURF_OP_OR     -3
#define SURF_OP_LEFT   -4
#define SURF_OP_RIGHT  -5

/* Timers */

#define TOT_TIMERS                19

#define TIMER_TRANSPORT            1
#define TIMER_TRANSPORT_ACTIVE     2
#define TIMER_TRANSPORT_TOTAL      3
#define TIMER_TRANSPORT_CYCLE      4
#define TIMER_BURNUP               5
#define TIMER_BURNUP_TOTAL         6
#define TIMER_PROCESS              7
#define TIMER_PROCESS_TOTAL        8
#define TIMER_BATEMAN              9
#define TIMER_BATEMAN_TOTAL       10
#define TIMER_RUNTIME             11
#define TIMER_INIT                12
#define TIMER_INIT_TOTAL          13
#define TIMER_VOLUME_CALC         14
#define TIMER_OMP_PARA            15
#define TIMER_MPI_OVERHEAD        16
#define TIMER_MPI_OVERHEAD_TOTAL  17
#define TIMER_FINIX               18
#define TIMER_MISC                19

/* Geometry errors */

#define GEOM_ERROR_NO_CELL         -1
#define GEOM_ERROR_MULTIPLE_CELLS  -2
/*

#define GEOM_ERROR_POINTER_ERROR   -3
*/

/* MPI methods */

#define MPI_METH_BC  1
#define MPI_METH_RED 2

/* Mesh types */

#define MESH_TYPE_CARTESIAN    1
#define MESH_TYPE_CYLINDRICAL  2
#define MESH_TYPE_SPHERICAL    3
#define MESH_TYPE_HEXX         4
#define MESH_TYPE_HEXY         5
#define MESH_TYPE_ORTHOGONAL   6
#define MESH_TYPE_ADAPTIVE     7

/* Mesh content */

#define MESH_CONTENT_RES   1
#define MESH_CONTENT_DAT   2
#define MESH_CONTENT_PTR   3

/* Mesh content data types */

#define MESH_CONTENT_DATA_TET  1
#define MESH_CONTENT_DATA_STL  2

/* Predictor and corrector steps */

#define PREDICTOR_STEP  0
#define CORRECTOR_STEP  1

/* Mesh plot types (numbers 1-3 are reserved for default) */

#define MPL_TYPE_FLUXPOW    4
#define MPL_TYPE_COLPT      5
#define MPL_TYPE_COLWGT     6
#define MPL_TYPE_GAMMAHEAT  7
#define MPL_TYPE_DET        8
#define MPL_TYPE_DET_IMP    9
#define MPL_TYPE_FLUXTEMP  10
#define MPL_TYPE_DENSITY   11
#define MPL_TYPE_DT_NEFF   12
#define MPL_TYPE_DT_GEFF   13

/* Color palettes for mesh plots */

#define PALETTE_HOT             1
#define PALETTE_COLD            2
#define PALETTE_HOTCOLD         3
#define PALETTE_JET             4
#define PALETTE_BW              5
#define PALETTE_HSV             6
#define PALETTE_SPRING          7
#define PALETTE_SUMMER          8
#define PALETTE_AUTUMN          9
#define PALETTE_WINTER         10
#define PALETTE_GREEN_PURPLE   11
#define PALETTE_PURPLE_ORANGE  12
#define PALETTE_BLUE_RED       13

/* Color scales */

#define COLOR_SCALE_LIN  1
#define COLOR_SCALE_LOG  2

/* Surface current detector type */

#define SURF_DET_TYPE_SURF 1
#define SURF_DET_TYPE_UNIV 2

/* TMS modes (järjestyksellä on väliä) */

#define TMS_MODE_NONE 0
#define TMS_MODE_MG   1
#define TMS_MODE_CE   2

/* Entropy calculation mode */

#define ENTROPY_CALC_NONE        0
#define ENTROPY_CALC_SRC         1
#define ENTROPY_CALC_ALL         2

/* Memory allocation */

#define MEM_ALLOC    1
#define MEM_REALLOC  2
#define MEM_FREE     3
#define MEM_ALLOW    4
#define MEM_DENY     5

/* Multi-physics interface types */


#define IFC_TYPE_PT_AVG     1
#define IFC_TYPE_REG_MESH   2
#define IFC_TYPE_FUNC       3
#define IFC_TYPE_TET_MESH   4
#define IFC_TYPE_FUEP       5
#define IFC_TYPE_FPIP       6
#define IFC_TYPE_OPENFOAM   7
#define IFC_TYPE_OF_MAT     8
#define IFC_TYPE_OF_SOLID   9

/* OpenFOAM file types */

#define OF_FILE_POINTS    1
#define OF_FILE_FACES     2
#define OF_FILE_OWNER     3
#define OF_FILE_NEIGHBOUR 4
#define OF_FILE_TEMP      5
#define OF_FILE_DENSITY   6
#define OF_FILE_MATERIAL  7
#define OF_FILE_MAP       8

/* Material-wise delta-tracking modes */

#define DT_MAT_BLOCK  1
#define DT_MAT_FORCE  2

/* Wielandt iteration modes */

#define WIELANDT_MODE_NONE   0
#define WIELANDT_MODE_FIX_K  1
#define WIELANDT_MODE_FIX_P  2

/* Source files */

#define SRC_FILE_BUF_SIZE 1000000

#define SRC_FILE_TYPE_SERPENT1       1
#define SRC_FILE_TYPE_S1_RENORM      2
#define SRC_FILE_TYPE_FUSION_PLASMA  3
#define SRC_FILE_TYPE_WGT_RENORM     4

/* Plot only mode */

#define STOP_AFTER_PLOT_NONE    0
#define STOP_AFTER_PLOT_GEOM    1
#define STOP_AFTER_PLOT_TRACKS  2

/* Restart options */

#define RESTART_CHECK     1
#define RESTART_OVERRIDE  2
#define RESTART_REPLACE   3

/* Material divisor flags (keksi noille paremmat nimet) */

#define MAT_DIV_TYPE_NONE    0
#define MAT_DIV_TYPE_PARENT  1
#define MAT_DIV_TYPE_S1      2
#define MAT_DIV_TYPE_NEW     3

/* Materials in burnup output */

#define BURN_OUT_MAT_DIV     1
#define BURN_OUT_MAT_PARENT  2
#define BURN_OUT_MAT_BOTH    3

/* Tracking modes */

#define TRACK_MODE_DT 1
#define TRACK_MODE_ST 2

/* Super-imposed detector types */

#define SUPERDET_TYPE_CURRENT 1
#define SUPERDET_TYPE_TLEFLUX 2

/* Detector type flags */

#define DETECTOR_TYPE_CUMU     -1
#define DETECTOR_TYPE_UNI_E    -2
#define DETECTOR_TYPE_UNI_L    -3

#define DETECTOR_TYPE_MULTI     2
#define DETECTOR_TYPE_DIVI      3

/* Track types */

#define TRACK_END_STRT   0
#define TRACK_END_SURF   1
#define TRACK_END_LEAK   2
#define TRACK_END_COLL   3
#define TRACK_END_VIRT   4
#define TRACK_END_TCUT   5
#define TRACK_END_WWIN   6

/* Additional types */

#define TRACK_END_CAPT   7
#define TRACK_END_FISS   8
#define TRACK_END_SCAT   9
#define TRACK_END_ECUT  10
#define TRACK_END_WCUT  11
#define TRACK_END_BC    12

/* Last track point for plotter */

#define TRACK_PLOT_LAST 1000

/* Fission matrix types */

#define FISSION_MATRIX_TYPE_MAT  1
#define FISSION_MATRIX_TYPE_UNI  2
#define FISSION_MATRIX_TYPE_LVL  3
#define FISSION_MATRIX_TYPE_XYZ  4

/* Event types (0-11 are from TRACK_END...) */

#define EVENT_TYPE_FISS  13

/* Event record flags */

#define RECORD_EVENT_PLOTTER     1
#define RECORD_EVENT_IFP         2
#define RECORD_EVENT_IMPORTANCE  4

/* Photon production in neutron reactions */

#define PHOTON_PROD_ANA  1
#define PHOTON_PROD_IMP  2

/* Photon constants */

#define PHOTON_NSS_MAX  40

/* Signaling modes for coupled calculation */

#define SIG_MODE_NONE   0
#define SIG_MODE_POSIX  1
#define SIG_MODE_FILE   2

/* Face indices for mesh cells */

#define MESH_CELL_FACE_LEFT    0
#define MESH_CELL_FACE_RIGHT   1
#define MESH_CELL_FACE_FRONT   2
#define MESH_CELL_FACE_BACK    3
#define MESH_CELL_FACE_BOTTOM  4
#define MESH_CELL_FACE_TOP     5

/* Different types of mesh cells */

#define MESH_CELL_TYPE_TET     0
#define MESH_CELL_TYPE_PYRAMID 1
#define MESH_CELL_TYPE_PRISM   2
#define MESH_CELL_TYPE_HEX     3
#define MESH_CELL_TYPE_POLY    4

/* Majorant extra xs types */

#define MAJORANT_EXTRA_FP_POISON_ITER   1
#define MAJORANT_EXTRA_MIXTURE_ITER     2

/* Different modes for tracking precursor concentrations */

#define PREC_MODE_NONE  0
#define PREC_MODE_MESH  1
#define PREC_MODE_POINT 2

/*****************************************************************************/

/***** Structures ************************************************************/

/* Complex number */

typedef struct {
  double re;
  double im;
} complex;


/* Matrix in compressed sparse column format*/

struct ccsMatrix{

  long m;            /* rivit */
  long n;            /* sarakket */
  long nnz;          /* nollasta eroavien lkm */

  long *colptr;      /* osoittimet 1. nollasta eroavaan joka sarakkeessa */
                     /* pituus on (n+1) */
  long *rowind;      /* rivi-indeksit */
  long *colind;      /* sarake-ineksit */
  long *rowptr;
  long *next;
  complex *values;   /* nollasta eroavien arvot */

};

/* Data structure to store nuclide data in depletion files */

struct depnuc {
  long ZAI;
  long Z;
  double AW;
  double lambda;
  double dh;
  double sf;
  double gI;
  double ingtox;
  double inhtox;
};

/*****************************************************************************/

/***** Function prototypes ***************************************************/

void AddBranching(long);

void AddBuf(double, double, long, long, long, ...);

void AddBuf1D(double, double, long, long, long);

void AddChains(long, long, long);

void AddItem(long, long);

void AddMesh(long, double, double, double, double, long);

void AddMeshIdx(long, double, long, long, long, long);

long AddNuclide(char *, long, char *, double, long, long);

void AddPrivateRes(long, double, long);

double *AddPts(double *, long *, const double *, long);

void AddSabData();

void AddSearchMesh(long, long, double, double, double, double, double, double);

void AddSortItem(long, long, long, long, long);

void AddStableNuclides();

void AddStat(double, long, ...);

long AddSTLPoint(long ***, long, long, long, double, double, double);

void AddValuePair(long, double, double, long);

void AdjustEnergyGrid(long, long, const double *);

void AdjustSabData(long);

void AllocInterfaceStat();

void AllocMacroXS();

void AllocMicroXS();

void AllocParticleStack(long, long);

void AllocPrecDet();

long AllocPrivateData(long, long);

void AllocValuePair(long);

void AllocStatHistory(long);

void Alpha(double, double, double *);

double AlphaXS(double);

void ApplyGCSymmetries(long);

void ARESOutput();

long AtoI(char *, char *, char *, long);

double AtoF(char *, char *,char *, long);

double AtomicRelaxation(long, long, long, long, double, double, double, double,
			double, long);

void AverageTransmuXS(long, double, double, long);

void AziRot(double, double *, double *, double *, long);

double B1FluxCorr(long gcu, double E);

long B1Flux(long, long, const double *, const double *, const double *,
	    const double *, const double *, const double *, const double *,
	    double *, double *);

void B1Solver();

void BanksToStore();

long BoundaryConditions(long *, double *, double *, double *, double *,
			double *, double *, double *, long);

void BroadCrossSection(long, long, long, double, double, double, double *,
		       long *, long, long);

void BsfunN(long, double *, double *, complex **, complex **, complex **);

double BufMean(long, ...);

double BufN(long, ...);

double BufVal(long, ...);

double BufWgt(long, ...);

void BurnMatCompositions();

void BurnMaterials(long, long);

long BurnMatrixSize(long);

void BurnupCycle();

complex c_add (complex, complex);

complex c_con (complex);

complex c_div (complex, complex);

complex c_mul (complex, complex);

double  c_norm(complex);

complex c_sub (complex, complex);

complex c_sqrt(complex);

complex c_exp(complex);

void CacheXS();

void CalcMicroGroupXS();

void CalculateActivities();

void CalculateBytes();

void CalculateDTMajorants();

void CalculateEntropies();

void CalculateMasses();

void CalculateMGXS();

void CalculateRelAlpha();

void CalculateRelPopSize();

void CalculateTetCenter(long, long, long, long);

void CalculateTransmuXS(long, long);

void CalculateUresMajorants(long);

void ccsMatrixColSpace(struct ccsMatrix *, long , long );

void ccsMatrixCopy(struct ccsMatrix *, struct ccsMatrix *);

void ccsMatrixIsort(struct ccsMatrix *);

struct ccsMatrix *ccsMatrixNew(long, long, long);

void ccsMatrixFree(struct ccsMatrix  *);

void ccsMatrixPrint(struct ccsMatrix *);

void CellCount(long, long, long, long);

void CellVolumes();

void CheckCoefCalc();

void CheckDuplicates();

void CheckNuclideData();

void CheckPolyhedMesh(long);

void CheckReaListSum(long, long, double, long, long);

void CheckUnused();

void ClearBuf();

void ClearInterfaceStat();

void ClearMicroGroupXS();

void ClearPrivateData(long);

void ClearPrivateRes();

void ClearRelTransmuXS();

void ClearStat(long);

void ClearTransmuXS();

void CloseList(long);

void CoefCycle();

void CoefOutput();

void CombineActinides();

void CombineFissionYields();

long CompareStr(long, long);

void ComplexRea(long, long, double *, double, double, double, double *,
		double *, double *, double, double *, double, double *, long);

void ColDet(long, long, double, double, double, double, double, double, double,
	    double, double, double, double, long);

void CollectBuf();

void CollectBurnData();

void CollectDet();

void CollectDynData();

void CollectParallelData();

void CollectPrecDet();

void CollectResults();

void CollectVRMeshData();

long Collision(long, long, double, double, double, double *, double *,
	       double *, double *, double *, double, long);

void ComptonScattering(long, long, long, double *, double, double, double,
		       double *, double *, double *, double, double, long);

void ContribDet(long, long, long, long, long, double, double, double, long);

void CoordExpans(long, double *, double *, double *, double, long);

void CoordTrans(long, double *, double *, double *, double *, double *,
		double *, long);

void CountDynSrc();

void CreateGeometry();

long CreateMesh(long, long, long, long, long, long, const double *, long);

long CreateUniverse(long, char *, long);

void CSplineConstruct(const double *, const double *, long, double, double,
                      double **, long, double *);

double CSplineIntegrate(double **, double *, double *, long, double, double);

double CSplineIntegrate0(double *, double *, long, double, double, double,
                         double);

void CSplineInterpolate(double **, double *, long, double *, double *, long);

void CSplineInterpolate0(double *, double *, long, double, double, double *,
                         double *, long, long);

double CylDis(double, double, double, double, double);

void DecayMeshPrecDet();

void DecayPointPrecDet();

double DensityFactor(long, double, double, double, double, long);

void DepletionPolyFit(long, long);

long DetBin(long, long, double, double, double, double, double, long);

void DetectorOutput();

long DetIdx(long, long, long, long, long, long, long, long, long, long);

double DetResponse(long, long, long, double, double, long);

void DFPos(long, double, double, double, long*, long *, long *, long *,
	   long *, long *, double *, double *, double *);

double *dfSol(long, long, double **, double *, complex **, complex **,
	      complex *);

void DFSolver(long, long, long, const double *, double, double **,
	      double **, complex **, complex **, complex *, complex *);

int Die(char *, ...);

void DiffCoefED(long, double, double, double, double, double,double,
		double, double, double, long);

void Disperse();

void Disperse2();

void DistributeMaterialData();

void DivideBurnMat();

void DivideMeshCell(long, long, long [8], long (*)[8], long (*)[6],
		    long [6], long [6], long *, long *, long *, long);

void DividePolyhedCell(long, long, long *);

void DividePolyhedFace(long, long, long, long, long, long *,
		       long *, long, long*);

void DivideZone(long, long *, long, long);

double DopMicroXS(long, long, double, double *, double, long);

void DopplerBroad();

double DTMajorant(long, double, long);

long DuplicateItem(long);

long DuplicateParticle(long, long);

void ElasticScattering(long, long, double *, double *, double *, double *,
		       long);

void Element(char *, char *, char *);

double ENDFColF(long, char *);

long ENDFColI(long, char *);

double ENDFInterp(long, double, double, double, double, double);

void ENDFNewLine(char *, char *, FILE *);

void Error(long, ...);

void EstimateRuntime();

void ExpandPrivateArrays();

long EventFromBank(long);

void EventToBank(long);

double FillSTLMesh(long, long, double, double, double);

void FinalizeMPI();

long FindTetCell(long, double, double, double, long);

void FindInterfaceRegions(long, long, long, long, double, double, double, long);

long FindLatticeRegion(long, long, double *, double *, double *, long *, long);

void FindMaterialPointers();

long FindNestRegion(long, long, double, double, double, long);

long FindNuclideData(char *, long, char *, double, long, long);

long FindPBRegion(long, long, double *, double *, double *, long *, long *,
		  long);

void FindRowIndexes(struct ccsMatrix *);

long FindUniverseCell(long, double, double, double, long *, long);

long FindSTLSolid(long, double, double, double, double, double, double,
		  long, long);

long *FindXSLimits(long, long, double, double, double, double *, long *);

long FinixPtrFromFpe(long);

long FirstItem(long);

long FissMtxIndex(long, long);

void FissMtxOutput();

void Fission(long, long, long, double *, double, double, double, double,
	     double *, double *, double *, double, double *, long);

void FixHexMesh(long, long);

void FixPolyhedMesh(long);

void FlushBank();

void FlushPrecSource();

void FormTransmuPaths(long, long, double, double, long, long);

void FreeMem();

long FromBank(long);

long FromSrc(long);

long FromStack(long, long);

long FromStore(long, long);

long FromTrkBank(long, long);

long FromQue(long);

long GaussianSubst(struct ccsMatrix *, complex *, complex *, complex *);

void GetBankedPrecursors();

void GetBurnIDs();

double *GetImportantPts(long, long *, long *, long *);

long GetLatticeIndexes(double, double, double, double, double, double, long *,
		       long *, long *, long);

long GetLineNumber(char *, long);

char **GetParams(char *, char *, long *, long *, long, long, char *);

double GetTemp(long, long);

char *GetText(long);

void GeometryPlotter(long);

double GridFactor(long, double, long);

long GridSearch(long, double);

void hessFactorization(long, complex **, complex **, complex **);

void HexNewTetFace(long, long, long, long (*)[6], long [6], long [6], long,
		   long, long, long, long, long *, long *,
		   long *, long, long);

void HexRotateCell(long *, long (*)[8]);

void HexRotateFace(long *, long, long);

double HisMean(long, long, ...);

double HisRelErr(long, long, ...);

double HisVal(long, long, ...);

void HomoFlux(long);

long ICMIdx(long, double, double, double, double *, double *, double *,
	    double *, double *, double *, double *, double *, double *);

char *IdxStr(long, long);

void IFCPoint(long, double *, double *, long);

long InCell(long, double, double, double, long, long);

void InelasticScattering(long, double *, double *, double *, double *, long);

void InitData();

void InitHistories();

void InitMPI(int, char **);

void InitOMP();

void InitPrecDet();

void InitPrecDetSource();

void InitSignal();

long InSuperCell(long, long, double, double, double, long);

long InterpolateData(const double *, double *, long, const double *,
		     const double *, long, long, long *, long *);

void InterpolateNubar(double *, long);

long InterpolateSab(long);

void IntersectionList(long, long *, long);

long InTetCell(long, long, double, double, double, long, long);

void InvokeBranch(long);

void IsotopeFractions(long);

long IsotoZAI(char *);

void IsotropicDirection(double *, double *, double *, long);

void IterateCC();

void IterateExternal();

void IterateKeff();

void KleinNishina(double, double *, double *, long);

long LastItem(long);

void Leak(long, double, double, double, double, double, double, double,
	  double, long);

void LevelScattering(long, double *, double *, double *, double *, long);

long LIFOListSize(long);

void LinkGCUMaterials(long, long);

void LinkReactions();

void LinkSabData();

void LUdecomposition(long, complex **, complex *, complex *);

double MacroUresCorr(long, double, double, long);

double MacroXS(long, double, long);

double MajorantXS(long, double, long);

double *MakeArray(double, double, long, long);

struct ccsMatrix *MakeBurnMatrix(long, long);

struct ccsMatrix *MakeBurnMatrixMSR(long, double **, long);

void MakeDepletionZones(long, long, long, long, long, long *);

long MakeEnergyGrid(long, long, long, long, const double *, long);

void MakePalette(long *, long *, long *, long, long);

void MakeRing(long);

void MaterialBurnup(long, double *, double *, double, double, long, long);

void MaterialTotals();

void MaterialVolumes();

void MatlabOutput();

void matProduct(long, long, long, complex **, complex **, complex **);

long MatPtr(long, long);

double *MatrixExponential(struct ccsMatrix *, double *, double);

void MaxSurfDimensions(long, double *, double *, double *, double *,
			  double *, double *);

double MaxwellEnergy(double, long);

double Mean(long, ...);

void *Mem(long, ...);

long MemCount();

void MeshCellConnectDiags(long [8], long (*)[8], long);

void MeshCellFromCGNS(long, long, long *, long *, long *, long (*)[6],
		      long (*)[4], long *, long *, long);

void MeshCellGetFace(long [8], long *, long, long);

void MeshCellGetFacePos(long *, long, long);

void MeshCellIndirection(long*, long);

void MeshCellRotateLists(long [8], long (*)[4], long *, long *, long (*)[6], long);

long MeshCellType(long, long);

double MeshCellVol(long, double, double, double);

long MeshIndex(long, double, double, double);

void MeshPlotter();

long MeshPtr(long, double, double, double);

double MeshTot(long);

double MeshVal(long, double, double, double);

double MGXS(long, double, long, long);

void MicroCalc();

double MicroMajorantXS(long, double, long);

double MicroXS(long, double, long);

double MinXS(long, double, long);

void MORAOutput();

long MoveDT(long, double, double, long *, double *, double *, double *,
	    double*, double *, double *, double *, double *, double, long);

void MoveItemFirst(long);

void MoveItemRight(long);

long MoveST(long, double, double, long *, double *, double *, double *,
	    double *, double *, double, double, double, long);

void MoveStore();

void MPITransfer(double *, double *, long, long, long);

long MyParallelMat(long, long);

double NearestBoundary(long);

double NearestPBSurf(long, double, double, double, double, double, double,
		     long);

double NearestMeshBoundary(long, double, double, double, double, double,
			   double, long *);

double NearestSTLSurf(long, double, double, double, double, double, double,
		      long);

double NearestUMSHSurf(long, double, double, double, double, double, double,
		       long);

void NestVolumes();

long NewItem(long, long);

long NewLIFOItem(long, long);

void NewReaList(long, long);

long NewStat(char *, long, ...);

long NextReaction(long, long *, double *, double *, double *, long);

long NextWord(char *, char *);

double NormCoef(long);

void NormalizeCompositions();

void NormalizeCritSrc();

long NormalizeDynSrc();

void NormalizePrecDet();

void Note(long, ...);

double Nubar(long, double, long);

long NumericGauss(struct ccsMatrix *, complex *, complex, complex *);

void Nxn(long, long, double *, double, double, double, double *, double *,
	 double *, double, double *, double, double *, long);

void OmpTester();

FILE *OpenDataFile(long, char *);

double OTFSabXS(long, double, double, long);

void OTFSabScattering(long, double *, double *, double *, double *, long);

void OverrideIDs();

void PairProduction(long, long, long, double, double, double, double, double,
		    double, double, double, double, long);

void parlett(long, complex *, complex **, complex **, complex **);

long ParseCommandLine(int, char **);


void Photoelectric(long, long, long, double, double, double, double, double,
		   double, double, double, double, long);

double PhotonMacroXS(long, double, long);

double PhotonMicroXS(long, double, long);

void PhotonProd(long, double, double, double, double, double, double,
		double, double, double, long);

#ifndef NO_GFX_MODE

long PlotTracks(long, gdImagePtr, long, long, double, double, double, double,
		double, double, long, long);

#endif

void PoisonEq();

double PoisonXS(long, double, long, long);

double PolarAngle(double, double);

long PolyPInF(long, long *, long);

long PolySameFace(long *, long * , long);

double PotCorr(long, double, double);

void PreallocMem(long, long);

void PrecDet(long, long, double, double, double, double, double, double,
	     double, double, long);

void PrecursorPopControl();

void PrepareCCIter();

void PrepareTransportCycle();

void PrintCoeVals(FILE *, long);

void PrintCompositions(long);

void PrintCoreDistr();

void PrintCycleOutput();

void PrintDepMatrix(long, struct ccsMatrix *, double t, double *, double *,
		    long id);

void PrintDepOutput();

void PrintDepVals(FILE *, char *, struct depnuc *, long, long, double **,
		  double *, double *);

void PrintGammaSpectra();

void PrintGeometryData();

void PrintHistoryOutput();

void PrintInterfaceOutput();

void PrintMaterialData();

void PrintMeshCell(long [8], long);

void PrintMVar(FILE *, long);

void PrintNuclideData(long, long);

void PrintPBData();

void PrintPrecDet();

void PrintProgress(long, long);

void PrintReactionLists();

void PrintTitle();

void PrintTMSDiagnostics();

void PrintValues(FILE *, char *, long, long, long, long, long, long);

void ProcessBC();

void ProcessBraData(long);

void ProcessBurnMat();

void ProcessBurnupEGroups();

void ProcessCells();

void ProcessCellMesh();

void ProcessCPD();

void ProcessComplementCells();

void ProcessCompton(long, long);

void ProcessDecayData(long);

void ProcessDecaySrc();

void ProcessDepHis();

void ProcessDetectors();

void ProcessDivisors();

void ProcessEDistributions(long, long);

void ProcessEntropy();

void ProcessEvents();

void ProcessFissionYields(long);

void ProcessFissMtx();

void ProcessGC();

void ProcessICM();

void ProcessIFCFB(long, long);

void ProcessIFCFunc(long, long);

void ProcessIFCPtAvg(long, long);

void ProcessIFCRegMesh(long, long);

void ProcessIFCTetMesh(long, long);

void ProcessInterface();

void ProcessInventory();

void ProcessLattices();

void ProcessMaterials();

void ProcessMeshPlots();

void ProcessMixture(long, long);

void ProcessMSR();

void ProcessMuDistributions(long);

void ProcessNests();

void ProcessNubarData(long);

void ProcessNuclides();

void ProcessPairProduction(long, long);

void ProcessPBGeometry();

void ProcessPhotoelectric(long, long);

void ProcessPhotonAtt();

void ProcessPhotonProd(long);

void ProcessPhotonRea(long);

void ProcessPoisons();

void ProcessPrecDet();

void ProcessRayleigh(long, long);

void ProcessReactionLists();

void ProcessRelaxation();

void ProcessReprocessors();

void ProcessSources();

void ProcessStats();

void ProcessSTLGeometry();

void ProcessSymmetries();

void ProcessTmpData();

void ProcessTimeBins();

void ProcessTransformations();

void ProcessTTB();

void ProcessUMSHGeometry();

void ProcessUresData(long);

void ProcessUserEGrids();

void ProcessVR();

void ProcessXSData();

void PulseDet(long, long, double, double, double, double, double, long);

void PutCompositions();

void PutMeshIdx(long, double, long, long, long);

void PutPoisonConc();

long PutText(char *);

void QRfactorization(long, complex **, complex **, complex **);

long RadGammaSrc(long, long, double *, double *, double, long);

double Rand64(unsigned long *);

double RandF(long);

void RayleighScattering(long, double, double *, double *, double *, long);

void ReactionCount();

void ReactionCutoff();

long ReactionTargetZAI(long);

char *ReactionMT(long);

void ReadACEFile(long);

void ReadBRAFile();

void ReadDecayFile();

void ReadDirectoryFile();

void ReadFissionYields();

long ReadIFCBins(long, FILE *, double, long);

void ReadIFCFB(long, long);

void ReadIFCFBLims(FILE *, long, long);

void ReadIFCFunc(long, long);

void ReadIFCOFMesh(long, long);

void ReadIFCPtAvg(long, long);

void ReadIFCRegMesh(long, long);

void ReadIFCTetMesh(long, long);

long ReadInfix(long, long *, long *);

void ReadInput(char *);

void ReadInterface();

void ReadMaterialComp();

double ReadMesh(long, long, long, long);

long ReadMeshPtr(long, long, long, long);

void ReadOFBatches(long);

char *ReadOFData(FILE *, long);

void ReadOFHeader(FILE *, long *, long *, long *);

void ReadPBGeometry();

void ReadPhotonData();

void ReadPlasmaSrc(long);

void ReadRestartFile(long);

void ReadSourceFile(long, double *, double *, double *, double *, double *,
		    double *, double *, double *, double *);

void ReadSTLGeometry();

char *ReadTextFile(char *);

void ReadUMSHGeometry();

long ReallocMem(long, long);

double ReaMulti(long, long, double, long);

void RecoilDet(long, double, double, double, double, double, double, double,
	       double, double, double, long);

long ReDistributeQues();

void ReDistributeStacks();

void ReduceBuffer();

void ReducePrivateRes();

void RefreshInventory();

unsigned long ReInitRNG(long);

void RelaxInterfacePower();

void RelaxTransmuXS();

double RelErr(long, ...);

long RemoveFlaggedItems(long, long, long, long);

void RemoveItem(long);

void RemoveVoidCells();

void Rendezvous(long *, double *);

void ReplaceItem(long, long);

void ReplacePhotonData();

void Reprocess(long);

void ResetOption(long, long);

void ResetPoisonConc();

void ResetTimer(long);

void ResizeDynSrc();

void ResizeFissionSrc();

double ResponseFunction(double, long);

void RetrieveComposition(long, long);

void RIACycle();

void RROutput();

void SabScattering(long, double *, double *, double *, double *, long);

void SampleENDFLaw(long, long, double, double *, double *, long);

void SampleDelnu();

void SampleMeshDelnu(long, long, long);

double SampleMu(long, long, double, long);

long SampleNu(double, long);

void SamplePlasmaSrc(long, double *, double *, double *, double *, double *,
		     double *, double *, double *, double *, long);

void SamplePointDelnu(long, long, long);

long SamplePrecursorGroup(long, double, long);

double SamplePTable(long, double, long);

long SampleReaction(long, long, double, double, long);

long SampleSrcPoint(long, long, long);

void schurFactorization(long, complex **, complex **, complex **);

void Score(long, long, double, double, double, double, double, double, double,
	   double, double, double, double, double, long);

void ScoreAdjoint(long, long, long, double, double, double, double, double,
		  double, double, double, double, double, double, long);

void ScoreAlb(long, double, double, double, double, double, double, double,
	      double, double, long);

void ScoreCapture(long, long, double, long);

void ScoreCPD(double, double, double, long);

void ScoreDF(double, double, double, double, double, double, double, double,
	     double, long);

void ScoreFission(long, long, double, double, double, long, double,
		  double, double, double, long, long, long);

void ScoreGC(double, double, double, double, double, double, double, double,
	     long, long, double, double, double, double, long);

void ScoreICMTrk(long, double, double, double, double, double, double, double,
		 double, double, long);

void ScoreICMCol(double, double, double, double, double, double, long, long,
		 double, double, double, double, double, double, long);

void ScoreInterfaceFlux(double, double, double, double, double, double, double,
			long);

void ScoreInterfacePower(double, double, double, double, double, double, long);

void ScoreMesh(long, long, double, double, double, double, double, double,
	       double, double, double, long);


void ScorePinPower(long, long, double, double, double, long);

void ScorePoison(double, long, double, double, double, long);

void ScorePB(double, double, long);

void ScoreScattering(long, long, double, double, double, double, double, long);

void ScoreSurf(long, double *, double *, double *, double, double, double,
	       double, double, double, double, double, double, long);

void ScoreTimeConstants(double, double, long, long, long);

void ScoreTransmuXS(double, long, double, double, long);

void ScoreUFS(double, long, double, double, double, double, double, double,
	      long);

void ScoreWWDCurr(double, double, double, double, double, double, long, long);

long SearchArray(const double *, double, long);

long SeekList(long, long, double, long);

long SeekListStr(long, long, char *);

long SetCoefCalc(long);

void SetDepStepSize(long, long);

void SetDirectPointers(long);

void SetFissE();

void SetNormalization(long);

void SetOption(long, long);

void SetOptimization();

void SetPathLevels(long, long);

void SetPrecursorGroups();

void SetSTLMeshPointers(long);

void ShareInputData();

void ShuntingYard(long, long *, long);

void SignalExternal(int);

void SignalHandler(int);

double **MSRReaList(double **, long, long);

void SortAll();

void SortArray(double *, long);

void SortList(long, long, long);

double Speed(long, double);

void SrcDet(long, long, double, double, double, double, double, double,
	    double, double, double, long);

void StartTimer(long);

long StatBin(long, ...);

void StatTests();

void StdComp(char *, char *);

double StdDev(long, ...);

double STLFacetDistance(long, double, double, double, double, double, double,
			long);

void STLMatFinder();

long STLRayTest(long, long, double, double, double, double, double, double,
		long, long);

void StopAtBoundary (long *, double *, double *, double *, double *, double,
		     double, double, long);

long StopAtWWBound(long, double *, double *, double *, double, double,
		   double, double, double *, double *, long *, long);

void StopCCIter();

void StopCI();

void StopTimer(long);

void StoreComposition(long, double, double);

void StoreHistoryPoint(long, long, long, double, double, double, double,
		       double, double, double, double, double, double, long);

void StoreTransmuXS(long, long, long, long);

void StoreValuePair(long, double, double, long);

void SumDivCompositions();

double SumPrivateData(long);

double SumPrivateRes(long);

void SumTotXS(long);

void SuperDet(long, double, double, double, double, double, double, double,
	      double, double, double, long);

void SurfaceNormal(long, double, double, double, double *, double *, double *,
		   long);

double SurfaceDistance(long, const double *, long, long, double, double,
		       double, double, double, double, long);

void SurfaceSrc(long, long, double *, double *, double *, double *, double *,
		double *, long);

double SurfaceVol(long);

void SwapItems(long, long);

void SwapUniverses(long, long);

struct ccsMatrix *SymbolicLU(struct ccsMatrix *);

double SymmetryBoundary(long, double, double, double, double, double, double);

void SystemStat();

void TargetVelocity(long, double, double *, double *, double *, double, double,
		    double, double, long);

long TestASCIIFile(char *);

void TestDOSFile(char *);

double TestParam(char *, char *, long, char *, long, ...);

void TestSTLGeometry();

void TestSTLSolids(long, long, long);

long TestSurface(long, double, double, double, long, long);

double TestValuePair(long, double, long);

void TestXS();

void TetPutBoundingBox(long, long, long[4]);

double TetraVol(long, long);

double *ThinGrid(double *, long *, double);

long TimeCutoff(long, long, long *, double *, double *, double *,
		double *, double, double, double, double, double, double,
		double, long, long);

char *TimeIntervalStr(double);

char *TimeStamp();

char *TimeStr(long);

double TimerCPUVal(long);

double TimerVal(long);

void TmpMajorants();

void ToBank(long, long);

void ToQue(long, long);

double TorusDis(double, double, double, double, double, double, double,
		double, double);

void ToStack(long, long);

void ToStore(long, long, long);

double TotXS(long, long, double, long);

void Tracking(long);

void TrackingError(long, double, long, long, long);

long TrackMode(long, long, double, double, double, long, long);

void TransportCycle();

complex trapz(long, double *, complex *);

double Truncate(double, long);

double *TTA(struct ccsMatrix *, double *, double);

double TTAChain(int, double, double *, double *);

void TTALoop (long, double, long, double *, double, double, double, double,
	      double *, struct ccsMatrix *,long);

double TTB(long, long, double, double, double, double, double, double, double,
	 double, double, long, long);

double UFSFactor(double, double, double, long);

void UnionizeGrid();

void UniSym(long, double *, double *, double *, double *, double *, double *);

void UniverseBoundaries();

void UpdateCIStop(long, double*, long);

double UresDiluMicroXS(long, double, long);

double UresFactor(long, double, long);

void UserIFC(long, double *, double *, double, double, double, double,
	     long, const double *);

void UserSrc(long, double *, double *, double *, double *,  double *, double *,
	     double *, double *, double *, long);

void UserSurf(long, long, const double *, double *, double *, double *,
	      double *, double *, double *, long *, double *, double, double,
	      double, double, double, double);

double ValuePairIdx(long, long);

double ValuePairVal(long, long);

double vectorNorm(long, complex *);

long VrTester(long, double, double, double, double, double, double, double,
	      double *, long);

void VolumesMC();

void Warn(char *, ...);

long WeightWindow(long, long, double, double, double, double, double, double,
		  double, double *, double, long, long);

long WhereAmI(double, double, double, double, double, double, long);

double *WorkArray(long, long, long, long);

void WriteCIMomFluxes();

void WriteDynSrc();

void WriteTetMeshtoGeo();

void WriteDepFile();

void WriteICMData();

void WriteSourceFile(long, double, double, double, double, double, double,
		     double, double, double, double, long);

void WriteUMSHtoSTL();

double WWDis(double, double, double, double, double, double);

double WWImportance(double, double, double, double);

void XSPlotter();

char *ZAItoIso(long, long);

double ZDis(double, double, double);

void ZoneCount(long, long, long);

/*****************************************************************************/

/***** Function prototypes for FINIX coupling ********************************/

#ifdef FINIX

/* Function prototypes */

void CollectFinix();
void CreateFinixIFC();
void DistributeFinix();
void FreeFinix();
void IterateFinix();
void PrintFinix();
void ProcessFinix();
void ReadFinixIFC();
void RunFinix(long, long);
void UpdateFinixIFC();
void UpdateFinixPower(long, long);
void WriteFinixIFC();
void WriteFinixInputFile(long, long);

#else

/* Replace by dummy definitions */

#define CollectFinix();
#define CreateFinixIFC();
#define DistributeFinix();
#define IterateFinix();
#define PrintFinix();
#define ProcessFinix();
#define ReadFinixIFC();
#define RunFinix(a);
#define UpdateFinixIFC();
#define UpdateFinixPower(a,b)
#define WriteFinixIFC(a);
#define WriteFinixInputFile(a,b);

#endif

/*****************************************************************************/

/***** Function prototypes for internal coupling *****************************/

#ifdef SerpentINT

/* Function prototypes */

void InitInternal();
void InternalIFCPoint();
void IterateInternal();
void ScoreInternaPower();

#else

/* Replace by dummy definitions */

#define InitInternal();
#define InternalIFCPoint();
#define IterateInternal();
#define ScoreInternalPower();

#endif

/*****************************************************************************/

/***** Functions replaced by macros in non-OpenMP mode ***********************/

#ifdef OPEN_MP

#define OMP_THREAD_NUM omp_get_thread_num()

#ifdef DEBUG

void AddPrivateData(long, double, long);
double GetPrivateData(long, long);
void PutPrivateData(long, double, long);

#else

#define AddPrivateData(a,b,c)(PRIVA[a + c*(long)RDB[DATA_REAL_PRIVA_SIZE]] += (double)b)
#define GetPrivateData(a,b)(PRIVA[a + b*(long)RDB[DATA_REAL_PRIVA_SIZE]])
#define PutPrivateData(a,b,c)(PRIVA[a + c*(long)RDB[DATA_REAL_PRIVA_SIZE]] = (double)b)

#endif

#else

#define OMP_THREAD_NUM 0

#ifdef DEBUG

void AddPrivateData(long, double, long);
double GetPrivateData(long, long);
double *PrivateDataPtr(long, long);
void PutPrivateData(long, double, long);

#else

#define AddPrivateData(a,b,c)(PRIVA[a] += (double)b)
#define GetPrivateData(a,b)(PRIVA[a])
#define PrivateDataPtr(a,b)(&PRIVA[a])
#define PutPrivateData(a,b,c)(PRIVA[a] = (double)b)

#endif

#endif


/*****************************************************************************/

/***** Functions replaced by macros in non-debug mode ************************/

#ifdef DEBUG

void CheckPointer(char *, char *, long, long);

void CheckValue(char *, char *, char *, double, double, double);

long ListPtr(long, long);

long NextItem(long);

long PrevItem(long);

double GetPrivateRes(long);

#else

#define CheckPointer(a, b, c, d)

#define CheckValue(a, b, c, d, e, f)

#define ListPtr(a, b)((long)RDB[(long)RDB[a + LIST_PTR_DIRECT] + b + 1])

#define NextItem(a)((long)RDB[a + LIST_PTR_NEXT])

#define PrevItem(a)((long)RDB[a + LIST_PTR_PREV])

#define GetPrivateRes(a)(RES2[a])

#endif

/*****************************************************************************/

/***** Additional macros *****************************************************/

#define ListCommon(a)((long)RDB[a + LIST_PTR_COMMON])
#define ListSize(a)((long)RDB[ListCommon(a) + LIST_COMMON_N_ITEMS])
#define ListRoot(a)((long)RDB[ListCommon(a) + LIST_COMMON_PTR_ROOT])
#define OMPPtr(a, b)((long)RDB[a] + b)
#define LorentzFactor(a)(a/NEUTRON_E0 + 1.0)

/*****************************************************************************/

/***** Global arrays and variables *******************************************/
/*                                                                           */
/* Serpent 2 stores data in the following arrays:                            */
/*                                                                           */
/* WDB    Writable main data block that contains nuclear data, geometry,     */
/*        pointters and calculation parameters. Must be OpenMP protected     */
/*        when written during threaded routines.                             */
/*                                                                           */
/* RDB    Pointer to WDB array, but type-cast to constant double. This is to */
/*        induce a compiler warning when trying to write in the main data    */
/*        at the wrong place.                                                */
/*                                                                           */
/* PRV    Data block for storing OpenMP private values. Divided into         */
/*        segments and can be accessed by special routines without atomic or */
/*        critical pragmas.                                                  */
/*                                                                           */
/* BUF    Buffer for storing cycle/batch wise data for statistics. Divided   */
/*        into segments or accessed with atomic pragmas by special routines. */
/*                                                                           */
/* RES1   First results array, used for storing statistics. Not accessed by  */
/*        OpenMP threads.                                                    */
/*                                                                           */
/* RES2   Second results array, containing large tables of integral data for */
/*        which the statistics is not needed. Divided into segments or       */
/*        accessed with atomic pragmas by special routines.                  */
/*                                                                           */
/* ASCII  Data array for storing text strings.                               */
/*                                                                           */
/* ACE    Data array for storing cross sections and nuclide data before      */
/*        processing. Freed before transport cycle.                          */
/*                                                                           */
/*****************************************************************************/

/* Arrays */
/* Added extern to compile on Mac, LMK 6/2016 */

extern double *ACE;
extern double *WDB;
extern double *PRIVA;
extern double *BUF;
extern double *RES1;
extern double *RES2;

extern const double *RDB;

extern char *ASCII;

/* This size must be larger than cache line width */

#define RNG_SZ 100

extern unsigned long *SEED;
extern unsigned long *SEED0;

/*****************************************************************************/

/* Output pointers */

extern FILE *err;
extern FILE *out;

/* Number of mpi tasks and id */

extern int mpitasks;
extern int mpiid;

/* Random number seed */

extern unsigned long parent_seed;

/* Collision counter */

/* Timers */
/* Added tag for struct, LMK 6/2016 */
struct Timer {
  double t0;
  double t;
  double cpu_t0;
  double cpu_t;
  int on;
};

extern struct Timer timer[TOT_TIMERS + 1];

/*****************************************************************************/
#ifdef __cplusplus
} // closing curly bracket
#endif

#endif
