
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
TITLE                     (idx, [1: 29])  = 'BWR Temperature feedback test' ;
INPUT_FILE_NAME           (idx, [1:  5])  = 'input' ;
WORKING_DIRECTORY         (idx, [1: 56])  = '/Users/kerblm/Desktop/Serpent_2.1.26/Coupling/2d_coolifc' ;
HOSTNAME                  (idx, [1: 15])  = 'inl604169.local' ;
CPU_TYPE                  (idx, [1:  7])  = 'Unknown' ;
START_DATE                (idx, [1: 24])  = 'Tue May 10 12:42:06 2016' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue May 10 12:43:09 2016' ;

% Run parameters:

POP                       (idx, 1)        = 3000 ;
CYCLES                    (idx, 1)        = 100 ;
SKIP                      (idx, 1)        = 10 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1462905726 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 0 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 0 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 1 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 0 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 1 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
SHARE_BUF_ARRAY           (idx, 1)        = 1 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 60])  = '/Users/kerblm/Desktop/c757mnyws00/xsdata/sss_jeff311u.xsdata' ;
DECAY_DATA_FILE_PATH      (idx, [1:  3])  = 'N/A' ;
SFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
NFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  1.88692E-01 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  2.46554E-02 0.00181  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.75345E-01 4.6E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  6.89656E-01 0.00020  0.00000E+00 0.0E+00 ];
IFC_COL_EFF               (idx, [1:   4]) = [  5.05473E-01 0.00013  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  3.51457E-01 0.00046  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  5.39535E+00 0.00116  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  4.30473E+01 0.00138  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  4.28797E+01 0.00139  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.92590E+01 0.00162  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  4.47683E+01 0.00114  0.00000E+00 0.0E+00 ];

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 100 ;
SOURCE_POPULATION         (idx, 1)        = 300326 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  3.00326E+03 0.00397 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  3.00326E+03 0.00397 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.04609E+00 ;
RUNNING_TIME              (idx, 1)        =  1.05492E+00 ;
INIT_TIME                 (idx, [1:  2])  = [  5.28167E-02  5.28167E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  1.83333E-03  1.83333E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  1.00025E+00  1.00025E+00  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  1.05487E+00  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 0.99164 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  9.92181E-01 0.00138 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  9.21605E-01 ;

% Memory usage:

ALLOC_MEMSIZE             (idx, 1)        = 329.41;
MEMSIZE                   (idx, 1)        = 288.62;
XS_MEMSIZE                (idx, 1)        = 39.28;
MAT_MEMSIZE               (idx, 1)        = 1.30;
RES_MEMSIZE               (idx, 1)        = 228.96;
MISC_MEMSIZE              (idx, 1)        = 19.08;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 40.79;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 2 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 167033 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 3 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 12 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 12 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 263 ;
TOT_TRANSMU_REA           (idx, 1)        = 0 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 0 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 0 ;
TMS_MODE                  (idx, 1)        = 0 ;
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

NORM_COEF                 (idx, [1:   4]) = [  1.59514E+11 0.00255  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  9.66585E-01 0.00474 ];
U235_FISS                 (idx, [1:   4]) = [  1.66678E+14 0.00096  9.01900E-01 0.00094 ];
U238_FISS                 (idx, [1:   4]) = [  1.81294E+13 0.00865  9.81001E-02 0.00866 ];
U235_CAPT                 (idx, [1:   4]) = [  3.75752E+13 0.00772  1.75387E-01 0.00600 ];
U238_CAPT                 (idx, [1:   4]) = [  1.61109E+14 0.00455  7.52194E-01 0.00146 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  6.00000E+03 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  1.10603E+00 5.8E-09 ];
TOT_GENRATE               (idx, [1:   2]) = [  4.57193E+14 0.00014 ];
TOT_FISSRATE              (idx, [1:   2]) = [  1.84808E+14 1.6E-05 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  2.14196E+14 0.00443 ];
TOT_ABSRATE               (idx, [1:   2]) = [  3.99004E+14 0.00238 ];
TOT_SRCRATE               (idx, [1:   2]) = [  4.78541E+14 0.00255 ];
TOT_FLUX                  (idx, [1:   2]) = [  2.61073E+16 0.00251 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  8.03586E+13 0.00631 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  4.79362E+14 0.00256 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  1.27653E+16 0.00242 ];
INI_FMASS                 (idx, 1)        =  5.42481E-03 ;
TOT_FMASS                 (idx, 1)        =  5.42481E-03 ;

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.47389E+00 0.00015 ];
FISSE                     (idx, [1:   2]) = [  2.02638E+02 1.6E-05 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  9.56515E-01 0.00259  9.48353E-01 0.00256  7.65011E-03 0.02592 ];
IMP_KEFF                  (idx, [1:   2]) = [  9.56003E-01 0.00255 ];
COL_KEFF                  (idx, [1:   2]) = [  9.56003E-01 0.00255 ];
ABS_KEFF                  (idx, [1:   2]) = [  9.56003E-01 0.00255 ];
ABS_KINF                  (idx, [1:   2]) = [  1.14884E+00 0.00237 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% Forward-weighted delayed neutron parameters:

FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  8.22388E-03 0.02208  2.37446E-04 0.12064  1.17983E-03 0.05573  7.10953E-04 0.06716  1.62466E-03 0.04549  2.47081E-03 0.03900  9.18969E-04 0.06262  8.16564E-04 0.07407  2.64650E-04 0.13220 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  4.69387E-01 0.03749  6.23335E-03 0.10050  2.71600E-02 0.02052  3.78467E-02 0.03533  1.31712E-01 0.01010  2.92467E-01 0.0E+00  6.06504E-01 0.03161  1.45496E+00 0.03533  1.77730E+00 0.10050 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Adjoint weighted time constants using Nauchi's method:

ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  3.42648E-05 0.00463  3.42262E-05 0.00463  3.90602E-05 0.04581 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  3.27569E-05 0.00411  3.27196E-05 0.00410  3.73722E-05 0.04607 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  7.90212E-03 0.02702  2.38998E-04 0.18165  1.20893E-03 0.08505  6.24706E-04 0.11715  1.54764E-03 0.07631  2.40280E-03 0.06071  9.36372E-04 0.10228  6.99852E-04 0.13156  2.42819E-04 0.21433 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  4.49579E-01 0.06063  1.24667E-02 3.8E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 4.6E-09  2.92467E-01 0.0E+00  6.66488E-01 4.2E-09  1.63478E+00 0.0E+00  3.55460E+00 4.7E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  3.34714E-05 0.02669  3.34404E-05 0.02669  2.80172E-05 0.10549 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  3.19816E-05 0.02667  3.19523E-05 0.02668  2.67113E-05 0.10490 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  9.20961E-03 0.10152  4.58129E-04 0.55017  1.38961E-03 0.24248  3.84206E-04 0.50024  2.35989E-03 0.21152  2.52834E-03 0.18547  8.56924E-04 0.31145  8.70346E-04 0.35630  3.62166E-04 0.48351 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  4.31788E-01 0.17222  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 3.9E-09  1.63478E+00 0.0E+00  3.55460E+00 0.0E+00 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  3.45157E-05 0.00314 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  3.29937E-05 0.00192 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  8.32410E-03 0.01578 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -2.41561E+02 0.01642 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  3.73549E-07 0.00236 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  4.53758E-06 0.00203  4.53769E-06 0.00202  4.57167E-06 0.02287 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  3.88988E-05 0.00259  3.89004E-05 0.00259  3.85741E-05 0.02892 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  4.59050E-01 0.00203  4.58748E-01 0.00206  5.23290E-01 0.03716 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.19471E+01 0.05472 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  2.66313E+01 0.00128  3.38011E+01 0.00152 ];

