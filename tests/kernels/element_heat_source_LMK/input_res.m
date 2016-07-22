
% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.26' ;
COMPILE_DATE              (idx, [1: 20])  = 'Apr 21 2016 10:47:14' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
INPUT_FILE_NAME           (idx, [1:  5])  = 'input' ;
WORKING_DIRECTORY         (idx, [1: 77])  = '/Users/kerblm/projects/SerpentCouplerII/tests/kernels/element_heat_source_LMK' ;
HOSTNAME                  (idx, [1: 15])  = 'inl604169.local' ;
CPU_TYPE                  (idx, [1:  7])  = 'Unknown' ;
START_DATE                (idx, [1: 24])  = 'Fri Jul 22 12:40:20 2016' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Fri Jul 22 12:40:23 2016' ;

% Run parameters:

POP                       (idx, 1)        = 4000 ;
CYCLES                    (idx, 1)        = 200 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1342234 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 4 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:   4]) = [  9.95971E-01  1.00978E+00  9.96395E-01  9.97849E-01  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 59])  = '/Users/kerblm/Desktop/c757mnyws00/xsdata/sss_endfb7u.xsdata' ;
DECAY_DATA_FILE_PATH      (idx, [1:  3])  = 'N/A' ;
SFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
NFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  8.09949E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.23521E-02 0.00839  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.87648E-01 0.00010  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  3.85387E-01 0.00344  0.00000E+00 0.0E+00 ];
IFC_COL_EFF               (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  9.84162E-01 0.00043  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_SAMPLING_EFF          (idx, [1:   2]) = [  8.06100E-01 0.00477 ];
TOT_COL_EFF               (idx, [1:   4]) = [  3.86437E-01 0.00308  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.65554E+00 0.00262  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_FAIL_STAT             (idx, [1:   6]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  1.09727E+01 0.00578  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  1.01792E+01 0.00626  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.57620E+01 0.01058  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  8.11696E-01 0.00152  0.00000E+00 0.0E+00 ];

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 30 ;
SOURCE_POPULATION         (idx, 1)        = 120582 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  4.01940E+03 0.01027 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  4.01940E+03 0.01027 ];
SIMULATION_COMPLETED      (idx, 1)        = 0 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.51541E-01 ;
RUNNING_TIME              (idx, 1)        =  5.76667E-02 ;
INIT_TIME                 (idx, [1:  2])  = [  1.86500E-02  1.86500E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  1.16667E-04  1.16667E-04 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.88833E-02  0.00000E+00  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  1.98011E-01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 2.62788 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.43180E+00 0.00492 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  6.79480E-01 ;

% Memory usage:

ALLOC_MEMSIZE             (idx, 1)        = 503.74;
MEMSIZE                   (idx, 1)        = 314.69;
XS_MEMSIZE                (idx, 1)        = 13.24;
MAT_MEMSIZE               (idx, 1)        = 1.54;
RES_MEMSIZE               (idx, 1)        = 13.21;
MISC_MEMSIZE              (idx, 1)        = 286.71;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 189.05;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 3 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 33210 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 1 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 4 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 4 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 71 ;
TOT_TRANSMU_REA           (idx, 1)        = 4 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 0 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 0 ;
TMS_MODE                  (idx, 1)        = 2 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  0.00000E+00 ;
TOT_DECAY_HEAT            (idx, 1)        =  0.00000E+00 ;
TOT_SF_RATE               (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  0.00000E+00 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  0.00000E+00 ;
INHALATION_TOXICITY       (idx, 1)        =  0.00000E+00 ;
INGESTION_TOXICITY        (idx, 1)        =  0.00000E+00 ;
SR90_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
TE132_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
I131_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
I132_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
CS134_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
CS137_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
TOT_PHOTON_SRC            (idx, 1)        =  0.00000E+00 ;

% Normaliation coefficient:

NORM_COEF                 (idx, [1:   4]) = [  1.26772E+13 0.00701  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
U235_FISS                 (idx, [1:   4]) = [  7.76589E+15 0.00730  1.00000E+00 0.0E+00 ];
U235_CAPT                 (idx, [1:   4]) = [  1.52572E+15 0.01160  5.63495E-01 0.01006 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  2.50000E+05 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  1.95986E+16 0.00039 ];
TOT_FISSRATE              (idx, [1:   2]) = [  7.71433E+15 0.0E+00 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  2.68484E+15 0.00585 ];
TOT_ABSRATE               (idx, [1:   2]) = [  1.03992E+16 0.00151 ];
TOT_SRCRATE               (idx, [1:   2]) = [  5.07087E+16 0.00701 ];
TOT_FLUX                  (idx, [1:   2]) = [  8.29363E+17 0.00574 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  4.02638E+16 0.00802 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  5.06630E+16 0.00646 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  5.08144E+17 0.00728 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.54055E+00 0.00039 ];
FISSE                     (idx, [1:   2]) = [  2.02270E+02 2.7E-09 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  3.88764E-01 0.00796  3.85470E-01 0.00754  4.07308E-03 0.07334 ];
IMP_KEFF                  (idx, [1:   2]) = [  3.87565E-01 0.00628 ];
COL_KEFF                  (idx, [1:   2]) = [  3.87035E-01 0.00689 ];
ABS_KEFF                  (idx, [1:   2]) = [  3.87565E-01 0.00628 ];
ABS_KINF                  (idx, [1:   2]) = [  1.89113E+00 0.00162 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% Forward-weighted delayed neutron parameters:

FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  1.73758E-02 0.03399  7.17489E-04 0.18969  2.35200E-03 0.06077  2.80326E-03 0.10362  8.12680E-03 0.05116  2.52121E-03 0.11491  8.55088E-04 0.17934 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  7.92462E-01 0.09233  7.91069E-03 0.14129  3.18241E-02 0.0E+00  1.09375E-01 2.7E-09  3.16990E-01 0.0E+00  1.26372E+00 0.04963  6.04546E+00 0.12157 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  1.03210E-02 0.07259  3.31689E-04 0.28728  1.51378E-03 0.21016  1.58096E-03 0.20461  4.74562E-03 0.08191  1.53341E-03 0.17066  6.15519E-04 0.22380 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  9.15950E-01 0.15156  1.24906E-02 0.0E+00  3.18241E-02 0.0E+00  1.09375E-01 2.7E-09  3.16990E-01 0.0E+00  1.35398E+00 4.7E-09  8.63638E+00 0.0E+00 ];

% Adjoint weighted time constants using Nauchi's method:

ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  2.04669E-05 0.02959  2.04499E-05 0.03026  1.90221E-05 0.30994 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  7.94484E-06 0.02891  7.94008E-06 0.02990  7.25017E-06 0.29965 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  1.05509E-02 0.07256  4.17413E-04 0.39281  1.18702E-03 0.16089  1.53045E-03 0.22225  5.59593E-03 0.09775  1.32662E-03 0.21967  4.93435E-04 0.32453 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  8.16848E-01 0.20111  1.24906E-02 1.0E-08  3.18241E-02 0.0E+00  1.09375E-01 6.7E-09  3.16990E-01 2.7E-09  1.35398E+00 3.9E-09  8.63638E+00 0.0E+00 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  1.87150E-05 0.07425  1.88170E-05 0.07498  3.27818E-06 0.68690 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  7.27676E-06 0.07303  7.31795E-06 0.07388  1.20370E-06 0.66795 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  8.93233E-03 0.26316  0.00000E+00 0.0E+00  1.68785E-03 0.54380  1.65309E-03 0.68553  4.72540E-03 0.42353  6.56197E-04 1.00000  2.09791E-04 1.00000 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  4.06071E-01 0.34933  0.00000E+00 0.0E+00  3.18241E-02 0.0E+00  1.09375E-01 9.1E-09  3.16990E-01 0.0E+00  1.35398E+00 0.0E+00  8.63638E+00 0.0E+00 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  9.16342E-03 0.24916  0.00000E+00 0.0E+00  1.69854E-03 0.57041  1.94331E-03 0.66782  4.60095E-03 0.38894  6.67342E-04 1.00000  2.53281E-04 1.00000 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  4.01006E-01 0.35971  0.00000E+00 0.0E+00  3.18241E-02 1.2E-08  1.09375E-01 0.0E+00  3.16990E-01 0.0E+00  1.35398E+00 0.0E+00  8.63638E+00 0.0E+00 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -5.60107E+02 0.37078 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  2.02278E-05 0.01003 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  7.85462E-06 0.00915 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  1.11870E-02 0.05946 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -5.55259E+02 0.06144 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  3.77902E-07 0.00973 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  1.56953E-06 0.00542  1.57015E-06 0.00551  1.45701E-06 0.05599 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  3.85100E-05 0.00836  3.85304E-05 0.00867  3.58298E-05 0.09802 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  1.52000E-01 0.00750  1.52949E-01 0.00751  1.00010E-01 0.07669 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.10269E+01 0.08901 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  1.00182E+01 0.00641  1.35439E+01 0.01204 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '0' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  1.00000E-11  5.00000E-09  1.00000E-08  1.50000E-08  2.00000E-08  2.50000E-08  3.00000E-08  3.50000E-08  4.20000E-08  5.00000E-08  5.80000E-08  6.70000E-08  8.00000E-08  1.00000E-07  1.40000E-07  1.80000E-07  2.20000E-07  2.50000E-07  2.80000E-07  3.00000E-07  3.20000E-07  3.50000E-07  4.00000E-07  5.00000E-07  6.25000E-07  7.80000E-07  8.50000E-07  9.10000E-07  9.50000E-07  9.72000E-07  9.96000E-07  1.02000E-06  1.04500E-06  1.07100E-06  1.09700E-06  1.12300E-06  1.15000E-06  1.30000E-06  1.50000E-06  1.85500E-06  2.10000E-06  2.60000E-06  3.30000E-06  4.00000E-06  9.87700E-06  1.59680E-05  2.77000E-05  4.80520E-05  7.55014E-05  1.48728E-04  3.67262E-04  9.06898E-04  1.42510E-03  2.23945E-03  3.51910E-03  5.50000E-03  9.11800E-03  1.50300E-02  2.47800E-02  4.08500E-02  6.74300E-02  1.11000E-01  1.83000E-01  3.02500E-01  5.00000E-01  8.21000E-01  1.35300E+00  2.23100E+00  3.67900E+00  6.06550E+00  2.00000E+01 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  1.84093E+04 0.0E+00  7.58022E+04 0.0E+00  1.48585E+05 0.0E+00  1.64916E+05 0.0E+00  1.45318E+05 0.0E+00  1.16838E+05 0.0E+00  8.08483E+04 0.0E+00  5.71023E+04 0.0E+00  3.99437E+04 0.0E+00  2.92448E+04 0.0E+00  2.40129E+04 0.0E+00  1.95627E+04 0.0E+00  1.78722E+04 0.0E+00  1.61079E+04 0.0E+00  1.47977E+04 0.0E+00  1.22700E+04 0.0E+00  1.22395E+04 0.0E+00  1.12517E+04 0.0E+00  1.07345E+04 0.0E+00  1.97348E+04 0.0E+00  1.79167E+04 0.0E+00  1.26147E+04 0.0E+00  8.01663E+03 0.0E+00  9.12318E+03 0.0E+00  8.60523E+03 0.0E+00  7.20367E+03 0.0E+00  1.32054E+04 0.0E+00  2.66889E+03 0.0E+00  3.40442E+03 0.0E+00  3.05588E+03 0.0E+00  1.65922E+03 0.0E+00  2.84452E+03 0.0E+00  1.88106E+03 0.0E+00  1.62186E+03 0.0E+00  2.92401E+02 0.0E+00  2.96122E+02 0.0E+00  2.92359E+02 0.0E+00  2.92932E+02 0.0E+00  3.35514E+02 0.0E+00  3.29761E+02 0.0E+00  3.44354E+02 0.0E+00  2.92394E+02 0.0E+00  5.73509E+02 0.0E+00  9.35188E+02 0.0E+00  1.13555E+03 0.0E+00  2.94238E+03 0.0E+00  3.00455E+03 0.0E+00  3.42293E+03 0.0E+00  2.51264E+03 0.0E+00  1.91852E+03 0.0E+00  1.51093E+03 0.0E+00  1.91676E+03 0.0E+00  3.77349E+03 0.0E+00  5.44781E+03 0.0E+00  1.08287E+04 0.0E+00  1.76873E+04 0.0E+00  2.74646E+04 0.0E+00  1.82467E+04 0.0E+00  1.34965E+04 0.0E+00  9.90911E+03 0.0E+00  8.97454E+03 0.0E+00  8.98130E+03 0.0E+00  7.58311E+03 0.0E+00  5.12572E+03 0.0E+00  4.85451E+03 0.0E+00  4.33160E+03 0.0E+00  3.64526E+03 0.0E+00  2.91352E+03 0.0E+00  1.95197E+03 0.0E+00  7.30304E+02 0.0E+00 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.89513E+00 0.0E+00 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  7.18938E+17 0.0E+00  1.07595E+17 0.0E+00 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  4.01348E-01 0.0E+00  1.99825E+00 0.0E+00 ];
INF_CAPT                  (idx, [1:   4]) = [  1.82421E-03 0.0E+00  1.26274E-02 0.0E+00 ];
INF_ABS                   (idx, [1:   4]) = [  1.01292E-02 0.0E+00  2.89663E-02 0.0E+00 ];
INF_FISS                  (idx, [1:   4]) = [  8.30501E-03 0.0E+00  1.63389E-02 0.0E+00 ];
INF_NSF                   (idx, [1:   4]) = [  2.13573E-02 0.0E+00  3.98130E-02 0.0E+00 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.57162E+00 0.0E+00  2.43670E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.02270E+02 0.0E+00  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  2.43469E-08 0.0E+00  2.70082E-06 0.0E+00 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  3.91251E-01 0.0E+00  1.96912E+00 0.0E+00 ];
INF_SCATT1                (idx, [1:   4]) = [  2.17088E-01 0.0E+00  5.39426E-01 0.0E+00 ];
INF_SCATT2                (idx, [1:   4]) = [  8.62517E-02 0.0E+00  1.26050E-01 0.0E+00 ];
INF_SCATT3                (idx, [1:   4]) = [  7.19947E-03 0.0E+00  3.74255E-02 0.0E+00 ];
INF_SCATT4                (idx, [1:   4]) = [ -8.86991E-03 0.0E+00 -1.26719E-02 0.0E+00 ];
INF_SCATT5                (idx, [1:   4]) = [  8.75931E-04 0.0E+00  9.89512E-03 0.0E+00 ];
INF_SCATT6                (idx, [1:   4]) = [  4.96134E-03 0.0E+00 -2.48431E-02 0.0E+00 ];
INF_SCATT7                (idx, [1:   4]) = [  6.34136E-04 0.0E+00  9.65934E-04 0.0E+00 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  3.91294E-01 0.0E+00  1.96912E+00 0.0E+00 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.17085E-01 0.0E+00  5.39426E-01 0.0E+00 ];
INF_SCATTP2               (idx, [1:   4]) = [  8.62530E-02 0.0E+00  1.26050E-01 0.0E+00 ];
INF_SCATTP3               (idx, [1:   4]) = [  7.19858E-03 0.0E+00  3.74255E-02 0.0E+00 ];
INF_SCATTP4               (idx, [1:   4]) = [ -8.86930E-03 0.0E+00 -1.26719E-02 0.0E+00 ];
INF_SCATTP5               (idx, [1:   4]) = [  8.73774E-04 0.0E+00  9.89512E-03 0.0E+00 ];
INF_SCATTP6               (idx, [1:   4]) = [  4.95638E-03 0.0E+00 -2.48431E-02 0.0E+00 ];
INF_SCATTP7               (idx, [1:   4]) = [  6.36634E-04 0.0E+00  9.65934E-04 0.0E+00 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  1.16421E-01 0.0E+00  1.25905E+00 0.0E+00 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  2.86318E+00 0.0E+00  2.64750E-01 0.0E+00 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  1.00863E-02 0.0E+00  2.89663E-02 0.0E+00 ];
INF_REMXS                 (idx, [1:   4]) = [  2.08030E-02 0.0E+00  2.98233E-02 0.0E+00 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  3.80545E-01 0.0E+00  1.07059E-02 0.0E+00  6.93764E-04 0.0E+00  1.96843E+00 0.0E+00 ];
INF_S1                    (idx, [1:   8]) = [  2.13888E-01 0.0E+00  3.19992E-03 0.0E+00  3.86239E-04 0.0E+00  5.39039E-01 0.0E+00 ];
INF_S2                    (idx, [1:   8]) = [  8.72213E-02 0.0E+00 -9.69579E-04 0.0E+00  2.13576E-04 0.0E+00  1.25836E-01 0.0E+00 ];
INF_S3                    (idx, [1:   8]) = [  8.33162E-03 0.0E+00 -1.13214E-03 0.0E+00  4.22012E-05 0.0E+00  3.73833E-02 0.0E+00 ];
INF_S4                    (idx, [1:   8]) = [ -8.50508E-03 0.0E+00 -3.64828E-04 0.0E+00 -5.34314E-06 0.0E+00 -1.26666E-02 0.0E+00 ];
INF_S5                    (idx, [1:   8]) = [  8.48680E-04 0.0E+00  2.72517E-05 0.0E+00 -4.91095E-05 0.0E+00  9.94423E-03 0.0E+00 ];
INF_S6                    (idx, [1:   8]) = [  5.05062E-03 0.0E+00 -8.92799E-05 0.0E+00 -3.73302E-05 0.0E+00 -2.48058E-02 0.0E+00 ];
INF_S7                    (idx, [1:   8]) = [  7.58040E-04 0.0E+00 -1.23904E-04 0.0E+00 -4.27697E-05 0.0E+00  1.00870E-03 0.0E+00 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.80588E-01 0.0E+00  1.07059E-02 0.0E+00  6.93764E-04 0.0E+00  1.96843E+00 0.0E+00 ];
INF_SP1                   (idx, [1:   8]) = [  2.13885E-01 0.0E+00  3.19992E-03 0.0E+00  3.86239E-04 0.0E+00  5.39039E-01 0.0E+00 ];
INF_SP2                   (idx, [1:   8]) = [  8.72226E-02 0.0E+00 -9.69579E-04 0.0E+00  2.13576E-04 0.0E+00  1.25836E-01 0.0E+00 ];
INF_SP3                   (idx, [1:   8]) = [  8.33073E-03 0.0E+00 -1.13214E-03 0.0E+00  4.22012E-05 0.0E+00  3.73833E-02 0.0E+00 ];
INF_SP4                   (idx, [1:   8]) = [ -8.50447E-03 0.0E+00 -3.64828E-04 0.0E+00 -5.34314E-06 0.0E+00 -1.26666E-02 0.0E+00 ];
INF_SP5                   (idx, [1:   8]) = [  8.46523E-04 0.0E+00  2.72517E-05 0.0E+00 -4.91095E-05 0.0E+00  9.94423E-03 0.0E+00 ];
INF_SP6                   (idx, [1:   8]) = [  5.04566E-03 0.0E+00 -8.92799E-05 0.0E+00 -3.73302E-05 0.0E+00 -2.48058E-02 0.0E+00 ];
INF_SP7                   (idx, [1:   8]) = [  7.60537E-04 0.0E+00 -1.23904E-04 0.0E+00 -4.27697E-05 0.0E+00  1.00870E-03 0.0E+00 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  1.03210E-02 0.07259  3.31689E-04 0.28728  1.51378E-03 0.21016  1.58096E-03 0.20461  4.74562E-03 0.08191  1.53341E-03 0.17066  6.15519E-04 0.22380 ];
LAMBDA                    (idx, [1:  14]) = [  9.15950E-01 0.15156  1.24906E-02 0.0E+00  3.18241E-02 0.0E+00  1.09375E-01 2.7E-09  3.16990E-01 0.0E+00  1.35398E+00 4.7E-09  8.63638E+00 0.0E+00 ];

