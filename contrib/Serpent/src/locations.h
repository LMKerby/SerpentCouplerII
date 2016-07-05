#ifndef LOCATIONS_H
#define LOCATIONS_H

#ifdef __cplusplus
extern "C" {
#endif
/*****************************************************************************/

/***** Main data array (direct pointers) *************************************/

/* Data sizes */

#define DATA_ASCII_DATA_SIZE              0
#define DATA_ALLOC_MAIN_SIZE              1
#define DATA_REAL_MAIN_SIZE               3
#define DATA_ALLOC_ACE_SIZE               4
#define DATA_REAL_ACE_SIZE                5
#define DATA_ALLOC_PRIVA_SIZE             6
#define DATA_REAL_PRIVA_SIZE              7
#define DATA_ALLOC_RES1_SIZE              8
#define DATA_REAL_RES1_SIZE               9
#define DATA_ALLOC_RES2_SIZE             10
#define DATA_REAL_RES2_SIZE              11
#define DATA_ALLOC_BUF_SIZE              12
#define DATA_REAL_BUF_SIZE               13

#define DATA_TOTAL_BYTES                 14
#define DATA_REAL_BYTES                  15
#define DATA_BYTE_COUNT                  16

#define DATA_TOT_XS_BYTES                17
#define DATA_TOT_MAT_BYTES               18
#define DATA_TOT_RES_BYTES               19
#define DATA_TOT_MISC_BYTES              20

/* CPU name and memory size */

#define DATA_CPU_MEM                     21
#define DATA_PTR_CPU_NAME                22
#define DATA_CPU_MHZ                     23

/* Date */

#define DATA_PTR_DATE                    24

/* Title */

#define DATA_PTR_TITLE                   25

/* Host name, working directory and compile date */

#define DATA_PTR_HOSTNAME                26
#define DATA_PTR_WORKDIR                 27
#define DATA_PTR_COMPILE_DATE            28

/* CPU time */

#define DATA_CPU_T0                      30
#define DATA_CPU_TIME                    35

/* List pointers */

#define DATA_PTR_S0                      50
#define DATA_PTR_T0                      51
#define DATA_PTR_M0                      52
#define DATA_PTR_NUC0                    53
#define DATA_PTR_ACE0                    54
#define DATA_PTR_DECAY_ACE0              55
#define DATA_PTR_TRANSMU_ACE0            56
#define DATA_PTR_C0                      57
#define DATA_PTR_NST0                    58
#define DATA_PTR_TR0                     59
#define DATA_PTR_GPL0                    60
#define DATA_PTR_MPL0                    61
#define DATA_PTR_PB0                     62
#define DATA_PTR_L0                      63
#define DATA_PTR_LVL0                    64
#define DATA_PTR_SCORE0                  65
#define DATA_PTR_U0                      66
#define DATA_PTR_SRC0                    67
#define DATA_PTR_DET0                    68
#define DATA_PTR_ENE0                    69
#define DATA_PTR_GCU0                    70
#define DATA_PTR_IFC0                    71
#define DATA_PTR_DIV0                    72
#define DATA_PTR_MVOL0                   73
#define DATA_PTR_SYM0                    74
#define DATA_PTR_REP0                    75
#define DATA_PTR_TME0                    76
#define DATA_PTR_FUN0                    77
#define DATA_PTR_RIA0                    80
#define DATA_PTR_ADF0                    82
#define DATA_PTR_PPW0                    83
#define DATA_PTR_MORA0                   84
#define DATA_PTR_ICM0                    85
#define DATA_PTR_CSM0                    86
#define DATA_PTR_UMSH0                   87
#define DATA_PTR_STL0                    88
#define DATA_PTR_FIN0                    89
#define DATA_PTR_BRA0                    90
#define DATA_PTR_COEF0                   91
#define DATA_PTR_ALB0                    92
#define DATA_PTR_MFLOW0                  93
#define DATA_PTR_WWD0                    94
#define DATA_PTR_FINBIN0                 95

/* File names and paths */

#define DATA_PTR_ACEDATA_FNAME_LIST     132
#define DATA_PTR_INPUT_FNAME            133
#define DATA_PTR_DECDATA_FNAME_LIST     134
#define DATA_PTR_NFYDATA_FNAME_LIST     135
#define DATA_PTR_SFYDATA_FNAME_LIST     136
#define DATA_PTR_BRADATA_FNAME_LIST     137
#define DATA_PTR_XSTEST_FNAME           138

/* Stuff for burnup calculation */

#define DATA_PTR_COEF_BU_PT             139
#define DATA_MORE_COEF_CALC             140
#define DATA_BU_SPECTRUM_COLLAPSE       141
#define DATA_BU_URES_EMIN               144
#define DATA_BU_URES_EMAX               145
#define DATA_BU_ACT_MIN_Z               146
#define DATA_BU_ACT_MAX_Z               147
#define DATA_MAX_DIV_SEP_LVL            148

/* Nuclide counters */

#define DATA_N_TOT_NUCLIDES             150
#define DATA_N_TRANSPORT_NUCLIDES       151
#define DATA_N_DOSIMETRY_NUCLIDES       152
#define DATA_N_DECAY_NUCLIDES           153
#define DATA_N_PHOTON_NUCLIDES          154

/* List of outer boundaries */

#define DATA_PTR_OUTER_BOUNDS           156

/* Geometry level data */

#define DATA_PTR_ROOT_UNIVERSE          158
#define DATA_GEOM_LEVELS                159
#define DATA_PTR_COLLISION_UNI          160

/* Global zone index */

#define DATA_PTR_ZONE_IDX               161

/* Geometry dimensions */

#define DATA_GEOM_MINX                  162
#define DATA_GEOM_MAXX                  163
#define DATA_GEOM_MINY                  164
#define DATA_GEOM_MAXY                  165
#define DATA_GEOM_MINZ                  166
#define DATA_GEOM_MAXZ                  167

#define DATA_GEOM_DIM                   168

/* ACE data pointers */

#define DATA_PTR_ACE_NFY_DATA           170
#define DATA_PTR_ACE_SFY_DATA           171

/* Fission product lists. (These are used to find the  */
/* ace data for FP nuclides in CombineFissionYields()) */

#define DATA_TOT_FP_NUCLIDES            180
#define DATA_PTR_FP_LIB_ID_LIST         181
#define DATA_PTR_FP_ZAI_LIST            182
#define DATA_PTR_AC_ZAI_LIST            183

/* Pointer to lost nuclide data */

#define DATA_PTR_NUCLIDE_LOST           185

/* Reaction counters */

#define DATA_N_TRANSPORT_REA            190
#define DATA_N_SPECIAL_REA              191
#define DATA_N_DECAY_REA                192
#define DATA_N_TRANSMUTATION_REA        193
#define DATA_N_TRANSPORT_BRANCH         194
#define DATA_N_DECAY_BRANCH             195
#define DATA_N_DEAD_PATH                196

/* Delayed neutron precursor groups */

#define DATA_PRECURSOR_GROUPS           200

/* Counters */

#define DATA_N_GEOM_PLOTS               210
#define DATA_N_PBED                     211
#define DATA_N_MATERIALS                212
#define DATA_N_BURN_MATERIALS           213
#define DATA_N_TOT_CELLS                214
#define DATA_N_UNION_CELLS              215

/* Global and unionized energy grid data */

#define DATA_ERG_TOL                    220
#define DATA_ERG_INITIAL_PTS            221
#define DATA_ERG_IMPORTANT_PTS          222
#define DATA_ERG_PTR_UNIONIZED_NGRID    223
#define DATA_ERG_PTR_UNIONIZED_PGRID    224

/* Minimum and maximum energy allowed in transport calculation */

#define DATA_NEUTRON_EMIN               240
#define DATA_NEUTRON_EMAX               241

#define DATA_PHOTON_EMIN                242
#define DATA_PHOTON_EMAX                243

#define DATA_NEUTRON_XS_EMIN            244
#define DATA_NEUTRON_XS_EMAX            245

#define DATA_PHOTON_XS_EMIN             246
#define DATA_PHOTON_XS_EMAX             247

/* Minimum macroscopic cross section */

#define DATA_MIN_NMACROXS               248
#define DATA_MIN_PMACROXS               249

/* Boundary condition and albedo */

#define DATA_STOP_AT_BOUNDARY           250
#define DATA_GEOM_BC0                   251
#define DATA_GEOM_BC1                   252
#define DATA_GEOM_BC2                   253
#define DATA_GEOM_BC3                   254
#define DATA_GEOM_ALBEDO1               255
#define DATA_GEOM_ALBEDO2               256
#define DATA_GEOM_ALBEDO3               257

/* Cut-offs */

#define DATA_DEP_TTA_CUTOFF             260
#define DATA_DEP_HALF_LIFE_CUTOFF       261
#define DATA_DEP_FP_YIELD_CUTOFF        262
#define DATA_MIN_TOTXS                  263
#define DATA_URES_DILU_CUT              264
#define DATA_TIME_CUT_TMIN              265
#define DATA_TIME_CUT_TMAX              266
#define DATA_GEN_CUT                    267
#define DATA_MAX_PROMPT_CHAIN_LENGTH    268

/* Equilibrium Xe-135 calculation */

#define DATA_XE135_DC                   270
#define DATA_I135_DC                    271
#define DATA_PM149_DC                   272

#define DATA_XENON_EQUILIBRIUM_MODE     273
#define DATA_SAMARIUM_EQUILIBRIUM_MODE  274
#define DATA_PTR_XENON_MAT_LIST         275
#define DATA_PTR_SAMARIUM_MAT_LIST      276

/* Xenon entropy */

#define DATA_XENON_ENTROPY              279

/* Warning messages */

#define DATA_WARN_NFY_SUM               280
#define DATA_WARN_SFY_SUM               281
#define DATA_WARN_ERROR_AWR             282
#define DATA_WARN_ERROR_BRANCH          283

/* k-eff iteration */

#define DATA_ITER_MODE                  284
#define DATA_ITER_KEFF                  285
#define DATA_ITER_NCYC                  286
#define DATA_ITER_VAL                   287
#define DATA_ITER_FIX                   288

#define DATA_ITER_ALB_F1                289
#define DATA_ITER_ALB_F2                290
#define DATA_ITER_ALB_F3                291

/* Dummy variable (used by GetText(), etc.) */

#define DATA_DUMMY                      300

/* Geometry plotter */

#define DATA_STOP_AFTER_PLOT            310
#define DATA_PLOTTER_MODE               311
#define DATA_QUICK_PLOT_MODE            312

/* Monte Carlo volume calculation */

#define DATA_VOLUME_MC_NMAX             320
#define DATA_VOLUME_MC_TMAX             321
#define DATA_VOLUME_MC_EMAX             322

/* URES variables */

#define DATA_URES_AVAIL                 340
#define DATA_USE_URES                   341
#define DATA_URES_PTR_USE_LIST          342
#define DATA_URES_USED                  343
#define DATA_URES_EMIN                  344
#define DATA_URES_EMAX                  345

/* DBRC and TMS (majorant generation) */

#define DATA_QPARAM_DBRC                360
#define DATA_QPARAM_TMS                 361
#define DATA_USE_DBRC                   362
#define DATA_PTR_DBRC                   363
#define DATA_DBRC_EMIN                  364
#define DATA_DBRC_EMAX                  365
#define DATA_PTR_DBRC_COUNT             366
#define DATA_PTR_DBRC_EXCEED_COUNT      367

/* Normalization */

#define DATA_PTR_NORM                   380
#define DATA_NORM_U235_FISSE            381
#define DATA_NORM_INCLUDE_DH            382
#define DATA_NORM_BURN                  383
#define DATA_NORM_PTR_FISSH             384
#define DATA_INI_FMASS                  385
#define DATA_INI_BURN_FMASS             386
#define DATA_TOT_FMASS                  387
#define DATA_TOT_BURN_FMASS             388


/* Isomeric branching ratio data */

#define DATA_PTR_BRA_LIST               400

/* Majorants */

#define DATA_PTR_MAJORANT               401
#define DATA_PTR_PHOTON_MAJORANT        402
#define DATA_MAJORANT_PTR_EXTRA_XS      403

/* Optimization and memory/data options */

#define DATA_OPTI_MODE                  420
#define DATA_OPTI_UNIONIZE_GRID         421
#define DATA_OPTI_RECONSTRUCT_MICROXS   422
#define DATA_OPTI_RECONSTRUCT_MACROXS   423
#define DATA_OPTI_INCLUDE_SPECIALS      424
#define DATA_OPTI_MODE0_INCLUDE_TOTAL   425
#define DATA_OPTI_IMPLICIT_RR           426
#define DATA_OPTI_GC_CALC               427
#define DATA_OPTI_MG_MODE               428
#define DATA_OPTI_SHARED_BUF            429
#define DATA_OPTI_SHARED_RES2           430
#define DATA_OPTI_OMP_REPRODUCIBILITY   431
#define DATA_OPTI_REPLAY                432
#define DATA_OPTI_ENTROPY_CALC          433
#define DATA_OPTI_PRINT_HIS             434
#define DATA_OPTI_MPI_REPRODUCIBILITY   435
#define DATA_OPTI_POISON_CALC           436
#define DATA_OPTI_MPI_BATCH_SIZE        437

#define DATA_POISON_XS_VOL_RAT          438

#define DATA_SOURCE_PT_ANIM             439
#define DATA_SOURCE_PT_ANIM_F           440
#define DATA_SOURCE_PT_ANIM_PALETTE     441

/* Delta-tracking */

#define DATA_OPT_USE_DT                 450
#define DATA_DT_NTHRESH                 451
#define DATA_DT_PTHRESH                 452
#define DATA_DT_PTR_BLOCK_LIST          453
#define DATA_DT_PTR_FORCE_LIST          454
#define DATA_DT_ENFORCE_NEXT_TRACK      455

/* Ignore void cells */

#define DATA_IGNORE_VOID_CELLS          460
#define DATA_PTR_VOID_CELL              461

/* Implicit Monte Carlo (TODO: ota toi OPT pois nimest�) */

#define DATA_OPT_IMPL_CAPT              470
#define DATA_OPT_IMPL_FISS              471
#define DATA_OPT_IMPL_FISS_NUBAR        472
#define DATA_OPT_IMPL_NXN               473
#define DATA_OPT_ROULETTE_W0            474
#define DATA_OPT_ROULETTE_P0            475
#define DATA_USE_WEIGHT_WINDOWS         476
#define DATA_WWD_LOWER_BOUND            477
#define DATA_WWD_UPPER_BOUND            478

/* Cross section plotter */

#define DATA_XSPLOT_NE                  480
#define DATA_XSPLOT_EMIN                481
#define DATA_XSPLOT_EMAX                482

/* Run parameters */

#define DATA_NEUTRON_MAX_TRACK_LOOP_ERR 496
#define DATA_PHOTON_MAX_TRACK_LOOP_ERR  497
#define DATA_NEUTRON_MAX_TRACK_LOOP     498
#define DATA_PHOTON_MAX_TRACK_LOOP      499
#define DATA_PRINT_INTERVAL             500
#define DATA_CRIT_POP                   501
#define DATA_CRIT_CYCLES                502
#define DATA_CRIT_SKIP                  503
#define DATA_SRC_POP                    504
#define DATA_SRC_BATCHES                505
#define DATA_SIMUL_BATCH_SIZE           506
#define DATA_SIMULATION_MODE            507
#define DATA_CYCLE_IDX                  508
#define DATA_SIMULATION_COMPLETED       509
#define DATA_BATCH_INTERVAL             510
#define DATA_BATCH_COUNT                511
#define DATA_DYN_PTR_TIME_BINS          512
#define DATA_DYN_POP_MIN                513
#define DATA_DYN_POP_MAX                514
#define DATA_DYN_TMIN                   515
#define DATA_DYN_TMAX                   516
#define DATA_DYN_NB                     517
#define DATA_DYN_TB                     518
#define DATA_DYN_WGT0                   519
#define DATA_TERMINATE_ON_DIE           520

#define DATA_OMP_MAX_THREADS            521
#define DATA_PTR_OMP_HISTORY_COUNT      522

/* History index for debugging */

#define DATA_PTR_PRIVA_HIS_IDX          523

/* Running modes */

#define DATA_NEUTRON_TRANSPORT_MODE     550
#define DATA_PHOTON_TRANSPORT_MODE      551
#define DATA_BURNUP_CALCULATION_MODE    552
#define DATA_VOLUME_CALCULATION_MODE    553
#define DATA_PARTICLE_DISPERSER_MODE    554
#define DATA_PARTICLE_REDEPLETE_MODE    555
#define DATA_USE_DECAY_SRC              556
#define DATA_PHOTON_PRODUCTION          557

/* Storage space for neutrons and gammas */

#define DATA_PART_PTR_NSTACK            560
#define DATA_PART_PTR_GSTACK            561
#define DATA_PART_PTR_QUE               562
#define DATA_PART_PTR_SOURCE            563
#define DATA_PART_PTR_BANK              564
#define DATA_PART_PTR_TRK_BANK          565
#define DATA_PART_ALLOC_N               566
#define DATA_PART_ALLOC_G               567
#define DATA_PART_NBUF_FACTOR           568
#define DATA_PART_GBUF_FACTOR           569
#define DATA_PART_MIN_NSTACK            570
#define DATA_PART_MIN_GSTACK            571

/* Recorded events */

#define DATA_EVENT_RECORD_FLAGS         572
#define DATA_PTR_EVENT_BANK             573
#define DATA_EVENT_BANK_SZ              574
#define DATA_EVENT_MAX_GEN              575

/* Size of history list */

#define DATA_HIST_LIST_SIZE             578

/* Run-time variables */

#define DATA_CYCLE_PROMPT_WGT           580
#define DATA_CYCLE_DELAYED_WGT          581
#define DATA_CYCLE_KEFF                 582
#define DATA_NHIST_TOT                  583
#define DATA_N_POP_EIG                  584
#define DATA_POP_EIG                    585
#define DATA_WIELANDT_MODE              586
#define DATA_WIELANDT_KEFF              587
#define DATA_WIELANDT_P                 588
#define DATA_WIELANDT_KP                589

/* N�� on ehk� sama asia ? */

#define DATA_NHIST_CYCLE                590
#define DATA_CYCLE_BATCH_SIZE           591

/* Estimated running times */

#define DATA_ESTIM_CYCLE_TIME           592
#define DATA_ESTIM_TOT_TIME             593

/* Radioactivity data */

#define DATA_TOT_ING_TOX                601
#define DATA_TOT_INH_TOX                602
#define DATA_TOT_ACTIVITY               603
#define DATA_TOT_SFRATE                 604
#define DATA_TOT_DECAY_HEAT             605
#define DATA_BURN_SFRATE                606
#define DATA_BURN_DECAY_HEAT            607
#define DATA_ACT_ACTIVITY               608
#define DATA_ACT_DECAY_HEAT             609
#define DATA_FP_ACTIVITY                610
#define DATA_FP_DECAY_HEAT              611
#define DATA_SR90_ACTIVITY              612
#define DATA_TE132_ACTIVITY             613
#define DATA_I131_ACTIVITY              614
#define DATA_I132_ACTIVITY              615
#define DATA_CS134_ACTIVITY             616
#define DATA_CS137_ACTIVITY             617
#define DATA_TOT_PHOTON_SRC_RATE        618
#define DATA_TOT_PHOTON_SRC_VOL         619
#define DATA_PHOTON_SRC_MAX_I           620

/* Group constant generation */

#define DATA_ERG_FG_NG                  650
#define DATA_ERG_FG_PTR_GRID            651
#define DATA_ERG_FG_PTR_PREDEF          652
#define DATA_GCU_PTR_UNI                653
#define DATA_HOMOFLUX_SOLVER            654
#define DATA_ADF_TRAPZ_PT               655

/* Core power distribution */

#define DATA_CORE_PDE_DEPTH             660
#define DATA_CORE_PDE_NZ                661
#define DATA_CORE_PDE_ZMIN              662
#define DATA_CORE_PDE_ZMAX              663
#define DATA_CORE_PDE_PTR_CORE          664
#define DATA_CORE_PDE_PTR_ASS           665
#define DATA_CORE_PDE_PTR_RES0          666
#define DATA_CORE_PDE_PTR_RES1          667
#define DATA_CORE_PDE_PTR_RES2          668
#define DATA_CORE_PDE_N0                669
#define DATA_CORE_PDE_N1                670
#define DATA_CORE_PDE_N2                671

/* Uniform fission source */

#define DATA_UFS_MODE                   700
#define DATA_UFS_PTR_FACTORS            701
#define DATA_UFS_PTR_SRC_MESH           702
#define DATA_UFS_ORDER                  703
#define DATA_UFS_MIN                    704
#define DATA_UFS_MAX                    705
#define DATA_UFS_PTR_LAT                706
#define DATA_UFS_NX                     707
#define DATA_UFS_XMIN                   708
#define DATA_UFS_XMAX                   709
#define DATA_UFS_NY                     710
#define DATA_UFS_YMIN                   711
#define DATA_UFS_YMAX                   712
#define DATA_UFS_NZ                     713
#define DATA_UFS_ZMIN                   714
#define DATA_UFS_ZMAX                   715

/* Fission source entropy */

#define DATA_ENTROPY_NX                 751
#define DATA_ENTROPY_NY                 752
#define DATA_ENTROPY_NZ                 753
#define DATA_ENTROPY_XMIN               754
#define DATA_ENTROPY_XMAX               755
#define DATA_ENTROPY_YMIN               756
#define DATA_ENTROPY_YMAX               757
#define DATA_ENTROPY_ZMIN               758
#define DATA_ENTROPY_ZMAX               759
#define DATA_ENTROPY_PTR_SPT_STAT       760
#define DATA_ENTROPY_PTR_SWG_STAT       761

/* Burnup calculation stuff */

#define DATA_BURN_CALC_INI_MASS         800
#define DATA_BURN_DECAY_CALC            801
#define DATA_BURN_STEP_PC               802
#define DATA_BURN_BUMODE                803
#define DATA_BURN_CRAM_K                804
#define DATA_BURN_STEP                  805
#define DATA_BURN_TIME_INTERVAL         806
#define DATA_BURN_BURNUP_INTERVAL       807
#define DATA_BURN_PTR_DEP               808
#define DATA_BURN_CUM_BURNTIME          809
#define DATA_BURN_CUM_BURNUP            810
#define DATA_BURN_CUM_REAL_BURNUP       811
#define DATA_BURN_PTR_INVENTORY         812
#define DATA_BURN_INVENTORY_NUCLIDES    813
#define DATA_BURN_PRINT_DEPMTX          814
#define DATA_BURN_ENECUT                815
#define DATA_BURN_TOT_STEPS             816
#define DATA_BURN_PRED_STEP             817
#define DATA_BURN_CORR_STEP             818
#define DATA_BURN_PRINT_COMP            819
#define DATA_BURN_PRINT_COMP_LIM        820
#define DATA_BURN_STEP_TYPE             821
#define DATA_BURN_MAT_OUTPUT            822
#define DATA_BURN_PRINT_INTERMEDIATE    823
#define DATA_BURN_MATERIALS_FLAG        824
#define DATA_BURN_PREV_KEFF             825
#define DATA_BURN_PREV_DKEFF            826
#define DATA_BURN_ACTI_PTR_DEP          827
#define DATA_BURN_ACTI_PREV_F           828

/* Counters needed for burnup calculation */

#define DATA_TOT_NUCLIDES               850
#define DATA_PRED_TRANSPORT_TIME        851
#define DATA_CORR_TRANSPORT_TIME        852
#define DATA_COEF_TRANSPORT_TIME        853

/* Parameters for top inventory calculation */

#define DATA_BURN_INV_TOP_MASS          900
#define DATA_BURN_INV_TOP_ACTIVITY      901
#define DATA_BURN_INV_TOP_SF            902
#define DATA_BURN_INV_TOP_GSRC          903
#define DATA_BURN_INV_TOP_DECAY_HEAT    904
#define DATA_BURN_INV_TOP_ING_TOX       905
#define DATA_BURN_INV_TOP_INH_TOX       906

/* Critical spectrum calculation */

#define DATA_B1_CALC                    910
#define DATA_B1_BURNUP_CORR             911
#define DATA_B1_ERR_LIMIT               912
#define DATA_B1_REPEATED                913
#define DATA_B1_CONVERGED               914

/* Micro-group structure */

#define DATA_MICRO_CALC_BATCH_COUNT     920
#define DATA_MICRO_CALC_BATCH_SIZE      921
#define DATA_MICRO_CALC_NORM            922
#define DATA_MICRO_PTR_EGRID            923
#define DATA_MICRO_PTR_IDX_MAP          924

/* lengths of the previous step and preceding predictor (AIs) */

#define DATA_BURN_PS1_LENGTH           1000
#define DATA_BURN_PRED_LENGTH          1001

/* Old total powers, used in SetDepStepSize only. (AIs) */

#define DATA_BURN_POW_PS1              1010
#define DATA_BURN_POW_BOS              1011
#define DATA_BURN_POW_EOS              1012

/* weights in the coefficients of 2. order poly fit for burnup calculation */
/* see depletionpolyfit.c for more (AIs)*/

#define DATA_BURN_FIT_C2W1             1020
#define DATA_BURN_FIT_C2W2             1021
#define DATA_BURN_FIT_C2W3             1022
#define DATA_BURN_FIT_C1W1             1023
#define DATA_BURN_FIT_C1W2             1024
#define DATA_BURN_FIT_C1W3             1025
#define DATA_BURN_FIT_C0W1             1026
#define DATA_BURN_FIT_C0W2             1027
#define DATA_BURN_FIT_C0W3             1028
#define DATA_BURN_FIT_TYPE             1029

#define DATA_BURN_PRED_TYPE            1030
#define DATA_BURN_PRED_NSS             1031
#define DATA_BURN_CORR_TYPE            1032
#define DATA_BURN_CORR_NSS             1033

/* Restart file */

#define DATA_WRITE_RESTART_FILE        1040
#define DATA_RESTART_WRITE_PTR_FNAME   1041
#define DATA_READ_RESTART_FILE         1042
#define DATA_RESTART_READ_PTR_FNAME    1043
#define DATA_RESTART_READ_POINT       1044
#define DATA_RESTART_READ_IDX          1045

/* Coefficient calculation */

#define DATA_COEF_CALC_IDX             1050
#define DATA_TOT_COEF_CALC             1051
#define DATA_COEF_CALC_BU_IDX          1052
#define DATA_TOT_COEF_BU               1053
#define DATA_COEF_CALC_TOT_RUNS        1054
#define DATA_COEF_CALC_RUN_IDX         1055
#define DATA_COEF_CALC_PTR_PARAM_LIST  1056
#define DATA_COEF_CALC_SPECIAL_MODE    1057
#define DATA_COEF_CALC_INCLUDE_ERRORS  1058

/* Energy grid for coarse multi-group cross sections */

#define DATA_COARSE_MG_NE              1060
#define DATA_COARSE_MG_PTR_GRID        1061

/* Delayed nubar option */

#define DATA_USE_DELNU                 1070

/* Doppler-broadening mode and temperature feedback */

#define DATA_USE_DOPPLER_PREPROCESSOR  1071
#define DATA_TMS_MODE                  1072
#define DATA_USE_DENSITY_FACTOR        1073

/* Tracking collison counter */

#define DATA_PTR_COLLISION_COUNT       1081

/* Buffer and RES2 reduced flag */

#define DATA_BUF_REDUCED               1082
#define DATA_RES2_REDUCED              1083

/* Perform statistical tests on group constants */

#define DATA_GC_STAT_TESTS             1084
#define DATA_RUN_STAT_TESTS            1085

/* Include scattering production in removal xs */

#define DATA_GC_REMXS_MULT             1088

/* Calculation of analog reaction rates */

#define DATA_ANA_RR_NCALC              1089
#define DATA_ANA_RR_PCALC              1090

/* Heat and gamma production cross sections */

#define DATA_INCLUDE_HEAT_PROD_XS      1091
#define DATA_INCLUDE_PHOT_PROD_XS      1092

/* Source points for track plotter */

#define DATA_TRACK_PLOTTER_HIS         1093

/* Track plotter time interval options */

#define DATA_TRACK_PLOT_HIS_LENGTH     1100
#define DATA_TRACK_PLOT_FRAMES         1101
#define DATA_TRACK_PLOT_TMIN           1102
#define DATA_TRACK_PLOT_TMAX           1103
#define DATA_TRACK_PLOT_NHIS           1104
#define DATA_TRACK_PLOT_ANIM           1105

/* Pointers to pre-allocated work arrays */

#define DATA_PTR_WORK_GRID1            1110
#define DATA_PTR_WORK_GRID2            1111
#define DATA_PTR_WORK_GRID3            1112

#define DATA_PTR_WORK_PRIVA_GRID1      1113
#define DATA_PTR_WORK_PRIVA_GRID2      1114

/* Reaction sampling */

#define DATA_NPHYS_SAMPLE_FISS         1120
#define DATA_NPHYS_SAMPLE_CAPT         1121
#define DATA_NPHYS_SAMPLE_SCATT        1122

/* Memory operations */

#define DATA_ALLOW_MEM_OP              1123
#define DATA_PRIVA_MEM_READY           1124

/* Other stuff */

#define DATA_PRINT_PREV_COMPLETE       1130
#define DATA_PTR_CRIT_SRC_DET          1131
#define DATA_PTR_RIA_SRC               1132
#define DATA_ALPHA_EIG                 1133
#define DATA_NORM_COEF_N               1134
#define DATA_NORM_COEF_G               1135
#define DATA_SORT_COUNT                1136
#define DATA_NORM_PTR_RAD_SRC_MAT      1137

/* Minimum xs for CFE */

#define DATA_CFE_N_MIN_L               1140
#define DATA_CFE_N_MIN_T               1141
#define DATA_CFE_G_MIN_L               1142
#define DATA_CFE_G_MIN_T               1143

/* Cache-optimized xs block */

#define DATA_PTR_CACHE_OPTI_XS         1160
#define DATA_CACHE_OPTI_EMAX           1161
#define DATA_CACHE_OPTI_NE             1162
#define DATA_CACHE_OPTI_NREA           1163

/* Fission matrix */

#define DATA_PTR_FMTX                  1170
#define DATA_FMTX_TYPE                 1171

/* Corrector iteration */

#define DATA_BURN_SIE                  1180
#define DATA_BURN_CI_TYPE              1181 /* see below (AIs) */
#define DATA_BURN_CI_MAXI              1182 /* max # of iterations (AIs)*/
#define DATA_BURN_CI_I                 1183 /* curent iter. (AIs) */
#define DATA_BURN_CI_LAST              1184 /* YES/NO flags last iteration */

#define DATA_BURN_CI_NBATCH            1185
#define DATA_BURN_CI_CYCLES            1186
#define DATA_BURN_CI_SKIP              1187
#define DATA_BURN_STEP_TOT             1188
#define DATA_BURN_CI_FLAG              1189

/* Number of progenies for beta-eff and prompt lifetime calculation */

#define DATA_IFP_OPT_PRINT_ALL         1190
#define DATA_IFP_CHAIN_LENGTH          1191
#define DATA_PERT_VAR_A                1192
#define DATA_PERT_VAR_C                1193
#define DATA_PERT_N_BATCH              1194

/* Group constant generation at multiple levels */

#define DATA_MULTI_LEVEL_GCU           1197

/* Interface current method */

#define DATA_ICM_CALC                  1201
#define DATA_ICM_PTR_OUTFILE           1202
#define DATA_ICM_NSEG                  1203
#define DATA_ICM_NSUB                  1204
#define DATA_ICM_NMU0                  1205
#define DATA_ICM_NMU1                  1206
#define DATA_ICM_NMU2                  1207
#define DATA_ICM_PTR_SUB               1208
#define DATA_ICM_PTR_MU0               1209
#define DATA_ICM_PTR_MU1               1210
#define DATA_ICM_PTR_MU2               1211
#define DATA_ICM_NG0                   1212
#define DATA_ICM_NG1                   1213
#define DATA_ICM_PTR_ENE0              1214
#define DATA_ICM_PTR_ENE1              1215

/* For external coupling */

#define DATA_WAITING                   1230
#define DATA_PPID                      1231
#define DATA_RUN_CC                    1232

/* Material composition files */

#define DATA_PTR_COMP_FILE             1240

/* Counters for line number calculation */

#define DATA_LINE_NUM_N0               1251
#define DATA_LINE_NUM_NL0              1252

/* Temporary array for storing items in the generation of adaptive meshes */

#define DATA_ADA_MESH_BUF              1260

/* STL geometry stuff */

#define DATA_STL_TEMP_ARRAY_SIZE       1270
#define DATA_STL_TEST_N_PTS            1271
#define DATA_STL_TEST_N_DIR            1272
#define DATA_STL_GEOM_TEST_MODE        1273
#define DATA_STL_FACET_EXD             1274
#define DATA_STL_ENFORCE_DT            1275

/* Global solution relaxation stuff */

#define DATA_SOL_REL_NTOT              1276
#define DATA_SOL_REL_N1                1277
#define DATA_SOL_REL_NCUR              1278
#define DATA_SOL_REL_CYCLES            1279
#define DATA_SOL_REL_ALPHA             1280
#define DATA_SOL_REL_ITER              1281
#define DATA_SOL_REL_FACT              1282
#define DATA_SOL_REL_MAX_ITER          1283
#define DATA_SOL_REL_MAX_POP           1284

#define DATA_BURN_CI_TOLER             1285

/* Photon transport data */

#define DATA_PHOTON_EKN                1290
#define DATA_PHOTON_USE_DOPPLER        1291
#define DATA_PHOTON_COMP_EANG          1292
#define DATA_PHOTON_BREM_TOT           1293
#define DATA_PHOTON_BREM_E             1294
#define DATA_PHOTON_USE_TTB            1295

#define DATA_PHOTON_COH_FNAME          1296
#define DATA_PHOTON_INCOH_FNAME        1297
#define DATA_PHOTON_RELAX_FNAME        1298
#define DATA_PHOTON_PESS_FNAME         1299
#define DATA_PHOTON_PETOT_FNAME        1300
#define DATA_PHOTON_CP_FNAME           1301
#define DATA_PHOTON_ELSP_FNAME         1302
#define DATA_PHOTON_ELBR_FNAME         1303
#define DATA_PHOTON_DATA_DIR           1304
#define DATA_PHOTON_TTBPM              1305
#define DATA_PHOTON_TTBEC              1306

/* Flag for iteration */

#define DATA_ITERATE                   1310

/* Fission source passing between iterations / depsteps */

#define DATA_USE_FSP                   1311

/* Number of inactive cycles to use with FSP */

#define DATA_FSP_CRIT_SKIP             1312

/* Number of additional materials in openFOAM interface */

#define DATA_OF_N_EXTRAMAT             1313

/* Communications filenames (inwards and outwards) */

#define DATA_PTR_COM_IN                1314
#define DATA_PTR_COM_OUT               1315

/* Signalling mode in coupled calculation */
/* none/POSIX/file */

#define DATA_CC_SIG_MODE               1316

/* Particle counts for different batches in dynamic mode */

#define DATA_PTR_DYN_PARTCOUNT         1317
#define DATA_PART_PTR_BOI_STORE        1318
#define DATA_PART_PTR_EOI_STORE        1319

/* Fraction of total memory to use at max */

#define DATA_CPU_MEM_FRAC              1320

/* FINIX input file name lists */

#define DATA_PTR_FINROD_FNAME          1321
#define DATA_PTR_FINOPTI_FNAME         1322
#define DATA_PTR_FINSCEN_FNAME         1323
#define DATA_PTR_FININIT_FNAME         1324

/* Pointers for precursor tracking and delayed neutrons */

#define DATA_PART_PTR_PSTACK           1325
#define DATA_PART_PTR_PSOURCE          1326
#define DATA_PART_ALLOC_P              1327
#define DATA_PART_MIN_PSTACK           1328
#define DATA_PART_PBUF_FACTOR          1329

#define DATA_PTR_PREC_DET              1330
#define DATA_PRECURSOR_TRANSPORT_MODE  1331

#define DATA_PREC_SRC_FACT             1332
#define DATA_PREC_STORE_TRESH          1333

/* Last global statistical variable */

#define DATA_LAST_GLOBAL_STAT          1334

/* Last value in data block */

#define DATA_LAST_VALUE                1340

/*****************************************************************************/

/***** Statistical variables *************************************************/

#define RES_TOT_NEUTRON_LEAKRATE        (DATA_LAST_VALUE + 1)
#define RES_TOT_NEUTRON_LOSSRATE        (DATA_LAST_VALUE + 2)
#define RES_TOT_NEUTRON_SRCRATE         (DATA_LAST_VALUE + 3)
#define RES_TOT_NEUTRON_CUTRATE         (DATA_LAST_VALUE + 4)
#define RES_TOT_NEUTRON_RR              (DATA_LAST_VALUE + 5)
#define RES_TOT_NEUTRON_FLUX            (DATA_LAST_VALUE + 6)
#define RES_TOT_NEUTRON_POWER           (DATA_LAST_VALUE + 7)
#define RES_TOT_PHOTON_LEAKRATE         (DATA_LAST_VALUE + 8)
#define RES_TOT_PHOTON_LOSSRATE         (DATA_LAST_VALUE + 9)
#define RES_TOT_PHOTON_SRCRATE          (DATA_LAST_VALUE + 10)
#define RES_TOT_PHOTON_CUTRATE          (DATA_LAST_VALUE + 11)
#define RES_TOT_PHOTON_RR               (DATA_LAST_VALUE + 12)
#define RES_TOT_PHOTON_FLUX             (DATA_LAST_VALUE + 13)
#define RES_TOT_PHOTON_HEATRATE         (DATA_LAST_VALUE + 14)
#define RES_TOT_FISSRATE                (DATA_LAST_VALUE + 15)
#define RES_TOT_NSF                     (DATA_LAST_VALUE + 16)
#define RES_TOT_NUBAR                   (DATA_LAST_VALUE + 17)
#define RES_TOT_FISSE                   (DATA_LAST_VALUE + 18)
#define RES_TOT_CAPTRATE                (DATA_LAST_VALUE + 19)
#define RES_TOT_INLPRODRATE             (DATA_LAST_VALUE + 20)
#define RES_TOT_ELARATE                 (DATA_LAST_VALUE + 21)
#define RES_TOT_ABSRATE                 (DATA_LAST_VALUE + 22)
#define RES_TOT_POWDENS                 (DATA_LAST_VALUE + 23)
#define RES_IMP_KEFF                    (DATA_LAST_VALUE + 24)
#define RES_IMP_KINF                    (DATA_LAST_VALUE + 25)
#define RES_ANA_KEFF                    (DATA_LAST_VALUE + 26)
#define RES_COL_KEFF                    (DATA_LAST_VALUE + 27)
#define RES_MEAN_POP_SIZE               (DATA_LAST_VALUE + 28)
#define RES_MEAN_POP_WGT                (DATA_LAST_VALUE + 29)
#define RES_SRC_MULT                    (DATA_LAST_VALUE + 30)
#define RES_TOT_GENRATE                 (DATA_LAST_VALUE + 31)
#define RES_TOT_RECIPVEL                (DATA_LAST_VALUE + 32)
#define RES_WIELANDT_K                  (DATA_LAST_VALUE + 33)
#define RES_WIELANDT_P                  (DATA_LAST_VALUE + 34)
#define RES_ANA_FISS_FRAC               (DATA_LAST_VALUE + 40)
#define RES_ANA_CAPT_FRAC               (DATA_LAST_VALUE + 41)
#define RES_ANA_CONV_RATIO              (DATA_LAST_VALUE + 42)
#define RES_CYCLE_RUNTIME               (DATA_LAST_VALUE + 43)
#define RES_CPU_USAGE                   (DATA_LAST_VALUE + 44)
#define RES_INI_SRC_WGT                 (DATA_LAST_VALUE + 57)
#define RES_NEW_SRC_WGT                 (DATA_LAST_VALUE + 58)
#define RES_EXT_K                       (DATA_LAST_VALUE + 59)
#define RES_NORM_COEF                   (DATA_LAST_VALUE + 60)
#define RES_ITER_VAL                    (DATA_LAST_VALUE + 68)
#define RES_ALB_NEUTRON_LEAKRATE        (DATA_LAST_VALUE + 69)
#define RES_GEOM_ALBEDO                 (DATA_LAST_VALUE + 70)

/* Forward-weighted time constants */

#define RES_FWD_ANA_BETA_ZERO           (DATA_LAST_VALUE + 73)
#define RES_FWD_ANA_LAMBDA              (DATA_LAST_VALUE + 74)

/* Adjoint-weighted time constants */

#define RES_ADJ_MEULEKAMP_BETA_EFF      (DATA_LAST_VALUE + 77)
#define RES_ADJ_MEULEKAMP_LAMBDA        (DATA_LAST_VALUE + 78)
#define RES_ADJ_NAUCHI_GEN_TIME         (DATA_LAST_VALUE + 79)
#define RES_ADJ_NAUCHI_LIFETIME         (DATA_LAST_VALUE + 80)
#define RES_ADJ_NAUCHI_BETA_EFF         (DATA_LAST_VALUE + 81)
#define RES_ADJ_NAUCHI_LAMBDA           (DATA_LAST_VALUE + 82)
#define RES_ADJ_IFP_GEN_TIME            (DATA_LAST_VALUE + 83)
#define RES_ADJ_IFP_LIFETIME            (DATA_LAST_VALUE + 84)
#define RES_ADJ_IFP_IMP_BETA_EFF        (DATA_LAST_VALUE + 85)
#define RES_ADJ_IFP_IMP_LAMBDA          (DATA_LAST_VALUE + 86)
#define RES_ADJ_IFP_ANA_BETA_EFF        (DATA_LAST_VALUE + 87)
#define RES_ADJ_IFP_ANA_LAMBDA          (DATA_LAST_VALUE + 88)
#define RES_ADJ_IFP_ROSSI_ALPHA         (DATA_LAST_VALUE + 89)
#define RES_ADJ_PERT_GEN_TIME           (DATA_LAST_VALUE + 90)
#define RES_ADJ_PERT_LIFETIME           (DATA_LAST_VALUE + 91)
#define RES_ADJ_PERT_BETA_EFF           (DATA_LAST_VALUE + 92)
#define RES_ADJ_PERT_ROSSI_ALPHA        (DATA_LAST_VALUE + 93)

/* Misc. analog time constants */

#define RES_ANA_MEAN_NCOL               (DATA_LAST_VALUE + 94)
#define RES_ANA_PHOTON_LIFETIME         (DATA_LAST_VALUE + 95)
#define RES_ANA_DELAYED_EMTIME          (DATA_LAST_VALUE + 96)
#define RES_ANA_SLOW_TIME               (DATA_LAST_VALUE + 97)
#define RES_ANA_THERM_TIME              (DATA_LAST_VALUE + 98)
#define RES_ANA_THERM_FRAC              (DATA_LAST_VALUE + 99)

#define RES_ST_TRACK_FRAC               (DATA_LAST_VALUE + 100)
#define RES_DT_TRACK_FRAC               (DATA_LAST_VALUE + 101)
#define RES_DT_TRACK_EFF                (DATA_LAST_VALUE + 102)
#define RES_IFC_COL_EFF                 (DATA_LAST_VALUE + 103)
#define RES_TOT_COL_EFF                 (DATA_LAST_VALUE + 104)
#define RES_REA_SAMPLING_EFF            (DATA_LAST_VALUE + 105)
#define RES_REA_SAMPLING_FAIL           (DATA_LAST_VALUE + 106)
#define RES_AVG_TRACKS                  (DATA_LAST_VALUE + 107)
#define RES_AVG_SURF_CROSS              (DATA_LAST_VALUE + 108)
#define RES_AVG_REAL_COL                (DATA_LAST_VALUE + 109)
#define RES_AVG_VIRT_COL                (DATA_LAST_VALUE + 110)
#define RES_SRC_SAMPLING_EFF            (DATA_LAST_VALUE + 111)
#define RES_SRC_MEAN_WGT                (DATA_LAST_VALUE + 112)
#define RES_AVG_TRACK_LOOPS             (DATA_LAST_VALUE + 113)

#define RES_MEAN_NGEN                   (DATA_LAST_VALUE + 114)

#define RES_DYN_PERIOD                  (DATA_LAST_VALUE + 115)
#define RES_DYN_POP                     (DATA_LAST_VALUE + 116)
#define RES_PROMPT_GEN_POP              (DATA_LAST_VALUE + 117)
#define RES_PROMPT_GEN_CUMU             (DATA_LAST_VALUE + 118)
#define RES_PROMPT_GEN_TIMES            (DATA_LAST_VALUE + 119)
#define RES_PROMPT_CHAIN_LENGTH         (DATA_LAST_VALUE + 120)

/* Myrkyt */

#define RES_I135_EQUIL_CONC             (DATA_LAST_VALUE + 121)
#define RES_XE135_EQUIL_CONC            (DATA_LAST_VALUE + 122)
#define RES_XE135_ABSRATE               (DATA_LAST_VALUE + 123)

#define RES_PM149_EQUIL_CONC            (DATA_LAST_VALUE + 124)
#define RES_SM149_EQUIL_CONC            (DATA_LAST_VALUE + 125)
#define RES_SM149_ABSRATE               (DATA_LAST_VALUE + 126)

#define RES_TMS_SAMPLING_EFF            (DATA_LAST_VALUE + 127)
#define RES_TMS_FAIL_STAT               (DATA_LAST_VALUE + 128)
#define RES_MIN_MACROXS                 (DATA_LAST_VALUE + 129)

#define RES_STL_RAY_TEST                (DATA_LAST_VALUE + 130)

#define RES_PHOTOELE_CAPT_RATE          (DATA_LAST_VALUE + 132)
#define RES_PAIRPROD_CAPT_RATE          (DATA_LAST_VALUE + 133)

#define RES_STAT_VARIABLES              134

/* Size of fixed data block and minimum acceptable pointer (NOTE: tossa */
/* pit�� nyt olla + 1 koska noi viittaa samaan blokkiin. Sin�ns� turhan */
/* monimutkaisesti tehty. */

#define DATA_FIXED_BLOCK_SIZE       (DATA_LAST_VALUE + RES_STAT_VARIABLES + 1)
#define VALID_PTR                   (DATA_FIXED_BLOCK_SIZE - 1)

/*****************************************************************************/

/***** List data *************************************************************/

/* This is data stored for every item */

#define LIST_DATA_SIZE   4

#define LIST_PTR_NEXT    0
#define LIST_PTR_PREV    1
#define LIST_PTR_COMMON  2
#define LIST_PTR_DIRECT  3

/* Simplified one-way LIFO list (nolla-alkio pit�� olla toi next) */

#define LIFO_LIST_DATA_SIZE   1

#define LIFO_LIST_PTR_NEXT    0

/* This is common data in a separate structure */

#define LIST_COMMON_DATA_SIZE      5

#define LIST_COMMON_ITEM_SIZE      0
#define LIST_COMMON_PTR_ROOT       1
#define LIST_COMMON_N_ITEMS        2
#define LIST_COMMON_PTR_FIRST      3
#define LIST_COMMON_PTR_LAST       4

/* Listable value pair block */

#define VALUE_PAIR_BLOCK_SIZE    (LIST_DATA_SIZE + 2)

#define VALUE_PAIR_VAL1          (LIST_DATA_SIZE + 0)
#define VALUE_PAIR_VAL2          (LIST_DATA_SIZE + 1)

/*****************************************************************************/

/***** Common variables for input parameters *********************************/

/* PARAM_N_COMMON:iin ei lis�t� LIST_DATA_SIZE:a */

#define PARAM_N_COMMON   3

#define PARAM_PTR_NAME   (LIST_DATA_SIZE + 0)
#define PARAM_PTR_FNAME  (LIST_DATA_SIZE + 1)
#define PARAM_LINE       (LIST_DATA_SIZE + 2)

/*****************************************************************************/

/***** Material block ********************************************************/

/* TODO: N�it� nimi� pit�� seriously j�rkev�itt�� !!! */

#define MATERIAL_BLOCK_SIZE            (LIST_DATA_SIZE + PARAM_N_COMMON + 142)

#define MATERIAL_OPTIONS               (LIST_DATA_SIZE + PARAM_N_COMMON + 0)
#define MATERIAL_PTR_NAME              (LIST_DATA_SIZE + PARAM_N_COMMON + 1)
#define MATERIAL_PTR_COMP              (LIST_DATA_SIZE + PARAM_N_COMMON + 2)
#define MATERIAL_PTR_MIX               (LIST_DATA_SIZE + PARAM_N_COMMON + 3)
#define MATERIAL_RGB                   (LIST_DATA_SIZE + PARAM_N_COMMON + 5)
#define MATERIAL_VOLUME                (LIST_DATA_SIZE + PARAM_N_COMMON + 6)
#define MATERIAL_VOLUME_GIVEN          (LIST_DATA_SIZE + PARAM_N_COMMON + 7)
#define MATERIAL_MASS                  (LIST_DATA_SIZE + PARAM_N_COMMON + 8)
#define MATERIAL_MASS_GIVEN            (LIST_DATA_SIZE + PARAM_N_COMMON + 9)
#define MATERIAL_INI_FISS_MDENS        (LIST_DATA_SIZE + PARAM_N_COMMON + 10)
#define MATERIAL_INI_FMASS             (LIST_DATA_SIZE + PARAM_N_COMMON + 11)
#define MATERIAL_PTR_MC_VOLUME         (LIST_DATA_SIZE + PARAM_N_COMMON + 12)
#define MATERIAL_PTR_MC_DENSITY        (LIST_DATA_SIZE + PARAM_N_COMMON + 13)
#define MATERIAL_BURN_RINGS            (LIST_DATA_SIZE + PARAM_N_COMMON + 14)
#define MATERIAL_COLOUR_IDX            (LIST_DATA_SIZE + PARAM_N_COMMON + 15)
#define MATERIAL_ADENS                 (LIST_DATA_SIZE + PARAM_N_COMMON + 16)
#define MATERIAL_MDENS                 (LIST_DATA_SIZE + PARAM_N_COMMON + 17)
#define MATERIAL_PTR_TOT_REA_LIST      (LIST_DATA_SIZE + PARAM_N_COMMON + 18)
#define MATERIAL_PTR_ELA_REA_LIST      (LIST_DATA_SIZE + PARAM_N_COMMON + 19)
#define MATERIAL_PTR_ABS_REA_LIST      (LIST_DATA_SIZE + PARAM_N_COMMON + 20)
#define MATERIAL_PTR_FISS_REA_LIST     (LIST_DATA_SIZE + PARAM_N_COMMON + 21)
#define MATERIAL_PTR_HEATT_REA_LIST    (LIST_DATA_SIZE + PARAM_N_COMMON + 22)
#define MATERIAL_PTR_PHOTP_REA_LIST    (LIST_DATA_SIZE + PARAM_N_COMMON + 23)
#define MATERIAL_PTR_INLP_REA_LIST     (LIST_DATA_SIZE + PARAM_N_COMMON + 24)
#define MATERIAL_PTR_PHOT_TOT_LIST     (LIST_DATA_SIZE + PARAM_N_COMMON + 25)
#define MATERIAL_PTR_PHOT_HEAT_LIST    (LIST_DATA_SIZE + PARAM_N_COMMON + 26)
#define MATERIAL_PTR_TOT_URES_LIST     (LIST_DATA_SIZE + PARAM_N_COMMON + 27)
#define MATERIAL_PTR_ABS_URES_LIST     (LIST_DATA_SIZE + PARAM_N_COMMON + 28)
#define MATERIAL_PTR_ELA_URES_LIST     (LIST_DATA_SIZE + PARAM_N_COMMON + 29)
#define MATERIAL_PTR_FISS_URES_LIST    (LIST_DATA_SIZE + PARAM_N_COMMON + 30)
#define MATERIAL_PTR_HEAT_URES_LIST    (LIST_DATA_SIZE + PARAM_N_COMMON + 31)
#define MATERIAL_PTR_TMP_MAJORANT_LIST (LIST_DATA_SIZE + PARAM_N_COMMON + 32)
#define MATERIAL_PTR_TOTXS             (LIST_DATA_SIZE + PARAM_N_COMMON + 33)
#define MATERIAL_PTR_ELAXS             (LIST_DATA_SIZE + PARAM_N_COMMON + 34)
#define MATERIAL_PTR_ABSXS             (LIST_DATA_SIZE + PARAM_N_COMMON + 35)
#define MATERIAL_PTR_FISSXS            (LIST_DATA_SIZE + PARAM_N_COMMON + 36)
#define MATERIAL_PTR_INLPXS            (LIST_DATA_SIZE + PARAM_N_COMMON + 37)
#define MATERIAL_PTR_FISSE             (LIST_DATA_SIZE + PARAM_N_COMMON + 38)
#define MATERIAL_PTR_NSF               (LIST_DATA_SIZE + PARAM_N_COMMON + 39)
#define MATERIAL_PTR_HEATTXS           (LIST_DATA_SIZE + PARAM_N_COMMON + 40)
#define MATERIAL_PTR_PHOTPXS           (LIST_DATA_SIZE + PARAM_N_COMMON + 41)
#define MATERIAL_PTR_TOTPHOTXS         (LIST_DATA_SIZE + PARAM_N_COMMON + 42)
#define MATERIAL_PTR_HEATPHOTXS        (LIST_DATA_SIZE + PARAM_N_COMMON + 43)
#define MATERIAL_PTR_TMP_MAJORANTXS    (LIST_DATA_SIZE + PARAM_N_COMMON + 44)
#define MATERIAL_MEM_SIZE              (LIST_DATA_SIZE + PARAM_N_COMMON + 45)
#define MATERIAL_TOT_DIV_MEM_SIZE      (LIST_DATA_SIZE + PARAM_N_COMMON + 46)
#define MATERIAL_PTR_SAB               (LIST_DATA_SIZE + PARAM_N_COMMON + 47)
#define MATERIAL_PTR_DEP_TRA_LIST      (LIST_DATA_SIZE + PARAM_N_COMMON + 48)
#define MATERIAL_PTR_DEP_FISS_LIST     (LIST_DATA_SIZE + PARAM_N_COMMON + 49)
#define MATERIAL_ACTIVITY              (LIST_DATA_SIZE + PARAM_N_COMMON + 50)
#define MATERIAL_PHOTON_SRC_RATE       (LIST_DATA_SIZE + PARAM_N_COMMON + 51)
#define MATERIAL_SFRATE                (LIST_DATA_SIZE + PARAM_N_COMMON + 52)
#define MATERIAL_DECAY_HEAT            (LIST_DATA_SIZE + PARAM_N_COMMON + 53)
#define MATERIAL_PTR_FLUX_SPEC         (LIST_DATA_SIZE + PARAM_N_COMMON + 54)
#define MATERIAL_PTR_FLUX_SPEC_SUM     (LIST_DATA_SIZE + PARAM_N_COMMON + 55)
#define MATERIAL_PTR_BURN_FLUX         (LIST_DATA_SIZE + PARAM_N_COMMON + 56)
#define MATERIAL_BURN_FLUX_PS1         (LIST_DATA_SIZE + PARAM_N_COMMON + 57)
#define MATERIAL_BURN_FLUX_BOS         (LIST_DATA_SIZE + PARAM_N_COMMON + 58)
#define MATERIAL_BURN_FLUX_EOS         (LIST_DATA_SIZE + PARAM_N_COMMON + 59)
#define MATERIAL_BURN_FLUX_SSA         (LIST_DATA_SIZE + PARAM_N_COMMON + 60)
#define MATERIAL_BURN_POW_PS1          (LIST_DATA_SIZE + PARAM_N_COMMON + 61)
#define MATERIAL_BURN_POW_BOS          (LIST_DATA_SIZE + PARAM_N_COMMON + 62)
#define MATERIAL_BURN_POW_EOS          (LIST_DATA_SIZE + PARAM_N_COMMON + 63)
#define MATERIAL_BURNUP                (LIST_DATA_SIZE + PARAM_N_COMMON + 64)
#define MATERIAL_OMP_ID                (LIST_DATA_SIZE + PARAM_N_COMMON + 65)
#define MATERIAL_MPI_ID                (LIST_DATA_SIZE + PARAM_N_COMMON + 66)
#define MATERIAL_PTR_DATA_BLOCK        (LIST_DATA_SIZE + PARAM_N_COMMON + 67)
#define MATERIAL_DATA_BLOCK_SIZE       (LIST_DATA_SIZE + PARAM_N_COMMON + 68)
#define MATERIAL_PTR_IFC               (LIST_DATA_SIZE + PARAM_N_COMMON + 69)
#define MATERIAL_PROC_IDX              (LIST_DATA_SIZE + PARAM_N_COMMON + 70)
#define MATERIAL_BURN_IDX              (LIST_DATA_SIZE + PARAM_N_COMMON + 71)
#define MATERIAL_PTR_GCU               (LIST_DATA_SIZE + PARAM_N_COMMON + 72)
#define MATERIAL_PTR_DIV               (LIST_DATA_SIZE + PARAM_N_COMMON + 73)
#define MATERIAL_DIV_TYPE              (LIST_DATA_SIZE + PARAM_N_COMMON + 74)
#define MATERIAL_DIV_PTR_PARENT        (LIST_DATA_SIZE + PARAM_N_COMMON + 75)
#define MATERIAL_DIV_N_TOT_ZONES       (LIST_DATA_SIZE + PARAM_N_COMMON + 76)
#define MATERIAL_DIV_N_SUB_ZONES       (LIST_DATA_SIZE + PARAM_N_COMMON + 77)
#define MATERIAL_DIV_N_ZONES           (LIST_DATA_SIZE + PARAM_N_COMMON + 78)
#define MATERIAL_DIV_ZONE_IDX          (LIST_DATA_SIZE + PARAM_N_COMMON + 79)
#define MATERIAL_BURN_PRINT_OUTPUT     (LIST_DATA_SIZE + PARAM_N_COMMON + 80)
#define MATERIAL_DT_MODE               (LIST_DATA_SIZE + PARAM_N_COMMON + 81)
#define MATERIAL_URES_EMIN             (LIST_DATA_SIZE + PARAM_N_COMMON + 82)
#define MATERIAL_URES_EMAX             (LIST_DATA_SIZE + PARAM_N_COMMON + 83)
#define MATERIAL_VOL_COUNT             (LIST_DATA_SIZE + PARAM_N_COMMON + 84)
#define MATERIAL_BURN_SORT_FLAG        (LIST_DATA_SIZE + PARAM_N_COMMON + 85)
#define MATERIAL_FMTX_IDX              (LIST_DATA_SIZE + PARAM_N_COMMON + 87)
#define MATERIAL_PTR_I135_ISO          (LIST_DATA_SIZE + PARAM_N_COMMON + 88)
#define MATERIAL_PTR_XE135_ISO         (LIST_DATA_SIZE + PARAM_N_COMMON + 89)
#define MATERIAL_PTR_PM147_ISO         (LIST_DATA_SIZE + PARAM_N_COMMON + 90)
#define MATERIAL_PTR_PM148_ISO         (LIST_DATA_SIZE + PARAM_N_COMMON + 91)
#define MATERIAL_PTR_PM148M_ISO        (LIST_DATA_SIZE + PARAM_N_COMMON + 92)
#define MATERIAL_PTR_PM149_ISO         (LIST_DATA_SIZE + PARAM_N_COMMON + 93)
#define MATERIAL_PTR_SM149_ISO         (LIST_DATA_SIZE + PARAM_N_COMMON + 94)
#define MATERIAL_PTR_I135_PROD_RATE    (LIST_DATA_SIZE + PARAM_N_COMMON + 95)
#define MATERIAL_PTR_XE135_PROD_RATE   (LIST_DATA_SIZE + PARAM_N_COMMON + 96)
#define MATERIAL_PTR_PM149_PROD_RATE   (LIST_DATA_SIZE + PARAM_N_COMMON + 97)
#define MATERIAL_PTR_SM149_PROD_RATE   (LIST_DATA_SIZE + PARAM_N_COMMON + 98)
#define MATERIAL_PTR_I135_ABS_RATE     (LIST_DATA_SIZE + PARAM_N_COMMON + 99)
#define MATERIAL_PTR_XE135_ABS_RATE    (LIST_DATA_SIZE + PARAM_N_COMMON + 100)
#define MATERIAL_PTR_PM149_ABS_RATE    (LIST_DATA_SIZE + PARAM_N_COMMON + 101)
#define MATERIAL_PTR_SM149_ABS_RATE    (LIST_DATA_SIZE + PARAM_N_COMMON + 102)
#define MATERIAL_XENON_EQUIL_CALC      (LIST_DATA_SIZE + PARAM_N_COMMON + 103)
#define MATERIAL_SAMARIUM_EQUIL_CALC   (LIST_DATA_SIZE + PARAM_N_COMMON + 104)
#define MATERIAL_PTR_I135_CONC         (LIST_DATA_SIZE + PARAM_N_COMMON + 105)
#define MATERIAL_PTR_XE135_CONC        (LIST_DATA_SIZE + PARAM_N_COMMON + 106)
#define MATERIAL_PTR_PM149_CONC        (LIST_DATA_SIZE + PARAM_N_COMMON + 107)
#define MATERIAL_PTR_SM149_CONC        (LIST_DATA_SIZE + PARAM_N_COMMON + 108)
#define MATERIAL_DEFAULT_PTR_LIB_ID    (LIST_DATA_SIZE + PARAM_N_COMMON + 109)
#define MATERIAL_DEFAULT_TMP           (LIST_DATA_SIZE + PARAM_N_COMMON + 110)
#define MATERIAL_DOPPLER_TEMP          (LIST_DATA_SIZE + PARAM_N_COMMON + 111)
#define MATERIAL_TMS_MODE              (LIST_DATA_SIZE + PARAM_N_COMMON + 112)
#define MATERIAL_TMS_TMIN              (LIST_DATA_SIZE + PARAM_N_COMMON + 113)
#define MATERIAL_TMS_TMAX              (LIST_DATA_SIZE + PARAM_N_COMMON + 114)
#define MATERIAL_USE_IFC               (LIST_DATA_SIZE + PARAM_N_COMMON + 115)
#define MATERIAL_PTR_NORM              (LIST_DATA_SIZE + PARAM_N_COMMON + 116)
#define MATERIAL_COEF_TEMP             (LIST_DATA_SIZE + PARAM_N_COMMON + 117)
#define MATERIAL_COEF_SAB              (LIST_DATA_SIZE + PARAM_N_COMMON + 118)
#define MATERIAL_BURN_FLUX_REL         (LIST_DATA_SIZE + PARAM_N_COMMON + 119)
#define MATERIAL_BURN_FLUX_AVE         (LIST_DATA_SIZE + PARAM_N_COMMON + 120)
#define MATERIAL_CI_BOS_ABSXS          (LIST_DATA_SIZE + PARAM_N_COMMON + 121)
#define MATERIAL_CI_EOS_ABSXS          (LIST_DATA_SIZE + PARAM_N_COMMON + 122)
#define MATERIAL_CI_AVE_ABSXS          (LIST_DATA_SIZE + PARAM_N_COMMON + 123)
#define MATERIAL_CI_BOS_ABSXS2         (LIST_DATA_SIZE + PARAM_N_COMMON + 124)
#define MATERIAL_CI_EOS_ABSXS2         (LIST_DATA_SIZE + PARAM_N_COMMON + 125)
#define MATERIAL_CI_AVE_ABSXS2         (LIST_DATA_SIZE + PARAM_N_COMMON + 126)
#define MATERIAL_CI_IDE                (LIST_DATA_SIZE + PARAM_N_COMMON + 127)
#define MATERIAL_PTR_DETBIN            (LIST_DATA_SIZE + PARAM_N_COMMON + 128)
#define MATERIAL_PTR_INFLOW            (LIST_DATA_SIZE + PARAM_N_COMMON + 129)
#define MATERIAL_PTR_OUTFLOW           (LIST_DATA_SIZE + PARAM_N_COMMON + 130)
#define MATERIAL_FLOW_IDX              (LIST_DATA_SIZE + PARAM_N_COMMON + 131)
#define MATERIAL_FLOW_N                (LIST_DATA_SIZE + PARAM_N_COMMON + 132)
#define MATERIAL_FLOW_PTR_FIRST        (LIST_DATA_SIZE + PARAM_N_COMMON + 133)
#define MATERIAL_PTR_ORIG_NUC_COMP     (LIST_DATA_SIZE + PARAM_N_COMMON + 134)
#define MATERIAL_PTR_DECAY_SRC         (LIST_DATA_SIZE + PARAM_N_COMMON + 135)
#define MATERIAL_PHOTON_ATT_NE         (LIST_DATA_SIZE + PARAM_N_COMMON + 136)
#define MATERIAL_PTR_PHOTON_ATT_E      (LIST_DATA_SIZE + PARAM_N_COMMON + 137)
#define MATERIAL_PTR_PHOTON_ATT_F      (LIST_DATA_SIZE + PARAM_N_COMMON + 138)
#define MATERIAL_SAMPLED_PHOTON_SRC    (LIST_DATA_SIZE + PARAM_N_COMMON + 139)
#define MATERIAL_MAX_ADENS             (LIST_DATA_SIZE + PARAM_N_COMMON + 140)
#define MATERIAL_PTR_TTB               (LIST_DATA_SIZE + PARAM_N_COMMON + 141)

/*****************************************************************************/

/***** Nuclide composition in material ***************************************/

#define COMPOSITION_BLOCK_SIZE   (LIST_DATA_SIZE + 4)

#define COMPOSITION_PTR_NUCLIDE  (LIST_DATA_SIZE + 0)
#define COMPOSITION_ADENS        (LIST_DATA_SIZE + 1)
#define COMPOSITION_ADENS_BOS    (LIST_DATA_SIZE + 2)
#define COMPOSITION_ADENS_AVE    (LIST_DATA_SIZE + 3)

/*****************************************************************************/

/***** Mixture composition ***************************************************/

/* TODO: muuta tän nimi ISO:ksi? */

#define MIXTURE_BLOCK_SIZE  (LIST_DATA_SIZE + 3)

#define MIXTURE_PTR_MAT     (LIST_DATA_SIZE + 0)
#define MIXTURE_VFRAC       (LIST_DATA_SIZE + 1)
#define MIXTURE_MFRAC       (LIST_DATA_SIZE + 2)

/*****************************************************************************/

/***** XS data array *********************************************************/

#define NUCLIDE_BLOCK_SIZE             (LIST_DATA_SIZE + 94)

#define NUCLIDE_PTR_NAME               (LIST_DATA_SIZE +  0)
#define NUCLIDE_TYPE                   (LIST_DATA_SIZE +  1)
#define NUCLIDE_OPTIONS                (LIST_DATA_SIZE +  2)
#define NUCLIDE_TYPE_FLAGS             (LIST_DATA_SIZE +  3)
#define NUCLIDE_ZAI                    (LIST_DATA_SIZE +  4)
#define NUCLIDE_ZA                     (LIST_DATA_SIZE +  5)
#define NUCLIDE_Z                      (LIST_DATA_SIZE +  6)
#define NUCLIDE_A                      (LIST_DATA_SIZE +  7)
#define NUCLIDE_I                      (LIST_DATA_SIZE +  8)
#define NUCLIDE_AW                     (LIST_DATA_SIZE +  9)
#define NUCLIDE_AWR                    (LIST_DATA_SIZE + 10)
#define NUCLIDE_PTR_LIB_ID             (LIST_DATA_SIZE + 11)
#define NUCLIDE_PTR_ACE                (LIST_DATA_SIZE + 12)
#define NUCLIDE_PTR_DECAY_ACE          (LIST_DATA_SIZE + 13)
#define NUCLIDE_PTR_PHOTON_ACE         (LIST_DATA_SIZE + 14)
#define NUCLIDE_PTR_REA                (LIST_DATA_SIZE + 15)
#define NUCLIDE_N_TRANSPORT_REA        (LIST_DATA_SIZE + 16)
#define NUCLIDE_N_SPECIAL_REA          (LIST_DATA_SIZE + 17)
#define NUCLIDE_N_DECAY_REA            (LIST_DATA_SIZE + 18)
#define NUCLIDE_N_TRANSPORT_BRANCH     (LIST_DATA_SIZE + 19)
#define NUCLIDE_N_DECAY_BRANCH         (LIST_DATA_SIZE + 20)
#define NUCLIDE_N_TRANSMUTATION_REA    (LIST_DATA_SIZE + 21)
#define NUCLIDE_N_DEAD_PATH            (LIST_DATA_SIZE + 22)
#define NUCLIDE_N_TRANSMUTATION_PATH   (LIST_DATA_SIZE + 23)
#define NUCLIDE_LAMBDA                 (LIST_DATA_SIZE + 24)
#define NUCLIDE_DECAY_E                (LIST_DATA_SIZE + 25)
#define NUCLIDE_SFBR                   (LIST_DATA_SIZE + 26)
#define NUCLIDE_PTR_NFY_DATA           (LIST_DATA_SIZE + 27)
#define NUCLIDE_PTR_SFY_DATA           (LIST_DATA_SIZE + 28)
#define NUCLIDE_NFY_NE                 (LIST_DATA_SIZE + 29)
#define NUCLIDE_PATH_LEVEL             (LIST_DATA_SIZE + 30)
#define NUCLIDE_PTR_EGRID              (LIST_DATA_SIZE + 31)
#define NUCLIDE_EGRID_NE               (LIST_DATA_SIZE + 32)
#define NUCLIDE_EMIN                   (LIST_DATA_SIZE + 33)
#define NUCLIDE_EMAX                   (LIST_DATA_SIZE + 34)
#define NUCLIDE_MAX_TOTXS              (LIST_DATA_SIZE + 35)
#define NUCLIDE_PTR_TOTXS              (LIST_DATA_SIZE + 36)
#define NUCLIDE_PTR_HEATPRODXS         (LIST_DATA_SIZE + 37)
#define NUCLIDE_PTR_PHOTPRODXS         (LIST_DATA_SIZE + 38)
#define NUCLIDE_PTR_SUM_ABSXS          (LIST_DATA_SIZE + 39)
#define NUCLIDE_PTR_ELAXS              (LIST_DATA_SIZE + 40)
#define NUCLIDE_PTR_FISSXS             (LIST_DATA_SIZE + 41)
#define NUCLIDE_PTR_NGAMMAXS           (LIST_DATA_SIZE + 42)
#define NUCLIDE_PTR_NGAMMAXS_ISO       (LIST_DATA_SIZE + 43)
#define NUCLIDE_PTR_PHOTON_INCOHEXS    (LIST_DATA_SIZE + 44)
#define NUCLIDE_PTR_PHOTON_COHEXS      (LIST_DATA_SIZE + 45)
#define NUCLIDE_PTR_PHOTON_PHOTOELXS   (LIST_DATA_SIZE + 46)
#define NUCLIDE_PTR_PHOTON_PAIRPRODXS  (LIST_DATA_SIZE + 47)
#define NUCLIDE_PTR_PHOTON_HEATPRODXS  (LIST_DATA_SIZE + 48)
#define NUCLIDE_PTR_SAMPLE_REA_LIST    (LIST_DATA_SIZE + 49)
#define NUCLIDE_URES_EMIN              (LIST_DATA_SIZE + 50)
#define NUCLIDE_URES_EMAX              (LIST_DATA_SIZE + 51)
#define NUCLIDE_PTR_URES_RND           (LIST_DATA_SIZE + 52)
#define NUCLIDE_PTR_MATRIX_IDX         (LIST_DATA_SIZE + 53)
#define NUCLIDE_MEMSIZE                (LIST_DATA_SIZE + 54)
#define NUCLIDE_PTR_PHOTON_LINE_SPEC   (LIST_DATA_SIZE + 55)
#define NUCLIDE_IDX                    (LIST_DATA_SIZE + 56)
#define NUCLIDE_INVENTORY_IDX          (LIST_DATA_SIZE + 57)
#define NUCLIDE_TMP_IDX                (LIST_DATA_SIZE + 58)
#define NUCLIDE_PTR_TOTFISS_REA        (LIST_DATA_SIZE + 59)
#define NUCLIDE_SPEC_ING_TOX           (LIST_DATA_SIZE + 60)
#define NUCLIDE_SPEC_INH_TOX           (LIST_DATA_SIZE + 61)
#define NUCLIDE_SPEC_PHOTON_I          (LIST_DATA_SIZE + 62)
#define NUCLIDE_ACE_PREC_GROUPS        (LIST_DATA_SIZE + 63)
#define NUCLIDE_PREV_COL_Z2            (LIST_DATA_SIZE + 64)
#define NUCLIDE_PREV_COL_COS           (LIST_DATA_SIZE + 65)
#define NUCLIDE_PREV_COL_ER            (LIST_DATA_SIZE + 66)
#define NUCLIDE_MIN_AFRAC              (LIST_DATA_SIZE + 67)
#define NUCLIDE_MAX_AFRAC              (LIST_DATA_SIZE + 68)
#define NUCLIDE_URES_SAMPLING          (LIST_DATA_SIZE + 69)
#define NUCLIDE_PTR_PHOTON_DATA        (LIST_DATA_SIZE + 70)
#define NUCLIDE_TMP_ADENS              (LIST_DATA_SIZE + 71)
#define NUCLIDE_TOP_MASS               (LIST_DATA_SIZE + 72)
#define NUCLIDE_TOP_ACTIVITY           (LIST_DATA_SIZE + 73)
#define NUCLIDE_TOP_SF                 (LIST_DATA_SIZE + 74)
#define NUCLIDE_TOP_GSRC               (LIST_DATA_SIZE + 75)
#define NUCLIDE_TOP_DECAY_HEAT         (LIST_DATA_SIZE + 76)
#define NUCLIDE_TOP_ING_TOX            (LIST_DATA_SIZE + 77)
#define NUCLIDE_TOP_INH_TOX            (LIST_DATA_SIZE + 78)
#define NUCLIDE_TEMP                   (LIST_DATA_SIZE + 79)
#define NUCLIDE_XS_TEMP                (LIST_DATA_SIZE + 80)
#define NUCLIDE_MAJORANT_TEMP          (LIST_DATA_SIZE + 81)
#define NUCLIDE_TMS_MIN_TEMP           (LIST_DATA_SIZE + 82)
#define NUCLIDE_TMS_MAX_TEMP           (LIST_DATA_SIZE + 83)
#define NUCLIDE_DBRC_MAX_TEMP          (LIST_DATA_SIZE + 84)
#define NUCLIDE_PTR_SAB_NUC            (LIST_DATA_SIZE + 85)
#define NUCLIDE_BRA_TYPE               (LIST_DATA_SIZE + 86)
#define NUCLIDE_PTR_SAB                (LIST_DATA_SIZE + 87)
#define NUCLIDE_SAB_EMAX               (LIST_DATA_SIZE + 88)
#define NUCLIDE_SAB_EMAXLOW            (LIST_DATA_SIZE + 89)
#define NUCLIDE_SAB_PTR_FREE           (LIST_DATA_SIZE + 90)
#define NUCLIDE_PTR_RELAX              (LIST_DATA_SIZE + 91)
#define NUCLIDE_PTR_PHOTON_PROD        (LIST_DATA_SIZE + 92)
#define NUCLIDE_FISSE                  (LIST_DATA_SIZE + 93)

/*****************************************************************************/

/****** FP identifier ********************************************************/

#define FP_IDENT_BLOCK_SIZE  (LIST_DATA_SIZE + 3)

#define FP_IDENT_PTR_ID      (LIST_DATA_SIZE + 0)
#define FP_IDENT_TEMP        (LIST_DATA_SIZE + 1)
#define FP_IDENT_TMS         (LIST_DATA_SIZE + 2)

/*****************************************************************************/

/***** Photon spectrum *******************************************************/

#define PHOTON_LINE_SPEC_BLOCK_SIZE  (LIST_DATA_SIZE + 2)

#define PHOTON_LINE_SPEC_E           (LIST_DATA_SIZE + 0)
#define PHOTON_LINE_SPEC_RI          (LIST_DATA_SIZE + 1)

/*****************************************************************************/

/***** Statistical variable **************************************************/

/* T�h�n ei voi laittaa LIST_DATA_SIZE:a */

#define STAT_BLOCK_SIZE              3

#define STAT_N                       0
#define STAT_X                       1
#define STAT_X2                      2

#define BUF_BLOCK_SIZE               3

#define BUF_VAL                      0
#define BUF_WGT                      1
#define BUF_N                        2

/*****************************************************************************/

/***** ACE data array ********************************************************/

/* T�ss� ei tarvita list dataa koska pointterit viittaa ACE taulukkoon */

#define ACE_BLOCK_SIZE                23

#define ACE_PTR_ALIAS                  0
#define ACE_PTR_NAME                   1
#define ACE_TYPE                       2
#define ACE_AW                         3
#define ACE_AWR                        4
#define ACE_ZAI                        5
#define ACE_ZA                         6
#define ACE_I                          7
#define ACE_TEMP                       8
#define ACE_PTR_LIB_ID                 9
#define ACE_PTR_NXS                   10
#define ACE_PTR_JXS                   11
#define ACE_PTR_XSS                   12
#define ACE_PTR_FILE                  13
#define ACE_LAMBDA                    14
#define ACE_DECAY_E                   15
#define ACE_PTR_DECAY_LIST            16
#define ACE_SFBR                      17
#define ACE_DECAY_NDK                 18
#define ACE_BOUND_ZA                  19
#define ACE_PTR_RAD_SPEC              20
#define ACE_DELNU_PREC                21
#define ACE_PTR_NEXT                  22

/*****************************************************************************/

/***** Energy grid structure *************************************************/

/* NOTE: T�m� on bin��ripuurakenne joka ei k�yt� linkitetyn listan */
/*       pointtereita tai rutiineja. */

#define ENERGY_GRID_BLOCK_SIZE    16

#define ENERGY_GRID_NE             0
#define ENERGY_GRID_I0             1
#define ENERGY_GRID_EMIN           2
#define ENERGY_GRID_EMAX           3
#define ENERGY_GRID_EMID           4
#define ENERGY_GRID_PTR_LOW        5
#define ENERGY_GRID_PTR_HIGH       6
#define ENERGY_GRID_PTR_DATA       7
#define ENERGY_GRID_NB             8
#define ENERGY_GRID_PTR_BINS       9
#define ENERGY_GRID_LOG_EMIN      10
#define ENERGY_GRID_LOG_EMAX      11
#define ENERGY_GRID_TYPE          12
#define ENERGY_GRID_PTR_PREV_VAL  13
#define ENERGY_GRID_INTERP_MODE   14
#define ENERGY_GRID_ALLOC_NE      15

/*****************************************************************************/

/***** Reaction data array ***************************************************/

#define REACTION_BLOCK_SIZE             (LIST_DATA_SIZE + 68)

#define REACTION_TYPE                   (LIST_DATA_SIZE +  0)
#define REACTION_NR                     (LIST_DATA_SIZE +  1)
#define REACTION_MT                     (LIST_DATA_SIZE +  2)
#define REACTION_RTYP2                  (LIST_DATA_SIZE +  3)
#define REACTION_RTYP3                  (LIST_DATA_SIZE +  4)
#define REACTION_RTYP4                  (LIST_DATA_SIZE +  5)
#define REACTION_RTYP5                  (LIST_DATA_SIZE +  6)
#define REACTION_BR                     (LIST_DATA_SIZE +  7)
#define REACTION_RFS                    (LIST_DATA_SIZE +  8)
#define REACTION_AWR                    (LIST_DATA_SIZE +  9)
#define REACTION_Q                      (LIST_DATA_SIZE + 10)
#define REACTION_TY                     (LIST_DATA_SIZE + 11)
#define REACTION_WGT_F                  (LIST_DATA_SIZE + 12)
#define REACTION_EMIN                   (LIST_DATA_SIZE + 13)
#define REACTION_EMAX                   (LIST_DATA_SIZE + 14)
#define REACTION_PTR_EGRID              (LIST_DATA_SIZE + 15)
#define REACTION_PTR_XS                 (LIST_DATA_SIZE + 16)
#define REACTION_XS_I0                  (LIST_DATA_SIZE + 17)
#define REACTION_XS_NE                  (LIST_DATA_SIZE + 18)
#define REACTION_TGT_ZAI                (LIST_DATA_SIZE + 19)
#define REACTION_PTR_TGT                (LIST_DATA_SIZE + 20)
#define REACTION_PTR_BRANCH_PARENT      (LIST_DATA_SIZE + 21)
#define REACTION_BRANCH_MT              (LIST_DATA_SIZE + 22)
#define REACTION_PTR_FISSY              (LIST_DATA_SIZE + 23)
#define REACTION_PTR_PARTIAL_LIST       (LIST_DATA_SIZE + 24)
#define REACTION_PTR_NUCLIDE            (LIST_DATA_SIZE + 26)
#define REACTION_PTR_MAT                (LIST_DATA_SIZE + 27)
#define REACTION_PTR_PREV_XS            (LIST_DATA_SIZE + 28)
#define REACTION_PTR_PREV_URES_XS0      (LIST_DATA_SIZE + 29)
#define REACTION_PTR_MAJORANT_XS        (LIST_DATA_SIZE + 30)
#define REACTION_PTR_URES               (LIST_DATA_SIZE + 31)
#define REACTION_PTR_ANG                (LIST_DATA_SIZE + 32)
#define REACTION_PTR_PHOTON_DIST        (LIST_DATA_SIZE + 33)
#define REACTION_URES_EMIN              (LIST_DATA_SIZE + 34)
#define REACTION_URES_EMAX              (LIST_DATA_SIZE + 35)
#define REACTION_PTR_ANA_RATE           (LIST_DATA_SIZE + 36)
#define REACTION_ITP                    (LIST_DATA_SIZE + 37)
#define REACTION_SAB_EMAX               (LIST_DATA_SIZE + 38)
#define REACTION_SAB_FRAC               (LIST_DATA_SIZE + 39)
#define REACTION_PTR_TNUBAR             (LIST_DATA_SIZE + 40)
#define REACTION_PTR_DNUBAR             (LIST_DATA_SIZE + 41)
#define REACTION_PTR_PREC_LIST          (LIST_DATA_SIZE + 42)
#define REACTION_PTR_ERG                (LIST_DATA_SIZE + 43)
#define REACTION_PTR_TRANSMUXS          (LIST_DATA_SIZE + 44)
#define REACTION_FISSY_IE0              (LIST_DATA_SIZE + 45)
#define REACTION_FISSY_IE1              (LIST_DATA_SIZE + 46)
#define REACTION_FISSY_IE2              (LIST_DATA_SIZE + 47)
#define REACTION_I135_YIELD             (LIST_DATA_SIZE + 48)
#define REACTION_XE135_YIELD            (LIST_DATA_SIZE + 49)
#define REACTION_PM147_YIELD            (LIST_DATA_SIZE + 50)
#define REACTION_PM148_YIELD            (LIST_DATA_SIZE + 51)
#define REACTION_PM148M_YIELD           (LIST_DATA_SIZE + 52)
#define REACTION_PM149_YIELD            (LIST_DATA_SIZE + 53)
#define REACTION_SM149_YIELD            (LIST_DATA_SIZE + 54)
#define REACTION_MODE                   (LIST_DATA_SIZE + 55)
#define REACTION_PTR_MGXS               (LIST_DATA_SIZE + 56)
#define REACTION_PTR_URES_MAX           (LIST_DATA_SIZE + 57)
#define REACTION_URES_MAX_N0            (LIST_DATA_SIZE + 58)
#define REACTION_URES_MAX_NP            (LIST_DATA_SIZE + 59)
#define REACTION_SAB_MIN_EM_E           (LIST_DATA_SIZE + 60)
#define REACTION_SAB_MAX_EM_E           (LIST_DATA_SIZE + 61)
#define REACTION_CACHE_OPTI_IDX         (LIST_DATA_SIZE + 62)
#define REACTION_PTR_TMP_MAJORANT       (LIST_DATA_SIZE + 63)
#define REACTION_PTR_0K_DATA            (LIST_DATA_SIZE + 64)
#define REACTION_PTR_BRA_STATE          (LIST_DATA_SIZE + 65)
#define REACTION_PTR_BRA_ISOMER         (LIST_DATA_SIZE + 66)
#define REACTION_PTR_MULT               (LIST_DATA_SIZE + 67)

/*****************************************************************************/

/***** Reaction list array ***************************************************/

#define SAMPLE_LIST_BLOCK_SIZE  (LIST_DATA_SIZE + 3)

#define SAMPLE_LIST_PTR_COMP    (LIST_DATA_SIZE + 0)
#define SAMPLE_LIST_PTR_REA     (LIST_DATA_SIZE + 1)
#define SAMPLE_LIST_PTR_COUNT   (LIST_DATA_SIZE + 2)

#define RLS_BLOCK_SIZE          (LIST_DATA_SIZE + 4)

#define RLS_PTR_MAT             (LIST_DATA_SIZE + 0)
#define RLS_REA_MODE            (LIST_DATA_SIZE + 1)
#define RLS_PTR_REA0            (LIST_DATA_SIZE + 2)
#define RLS_PTR_NEXT            (LIST_DATA_SIZE + 3)

#define RLS_DATA_BLOCK_SIZE     (LIST_DATA_SIZE + 8)

#define RLS_DATA_PTR_NUCLIDE    (LIST_DATA_SIZE + 0)
#define RLS_DATA_PTR_REA        (LIST_DATA_SIZE + 1)
#define RLS_DATA_COMP_IDX       (LIST_DATA_SIZE + 2)
#define RLS_DATA_EMIN           (LIST_DATA_SIZE + 3)
#define RLS_DATA_EMAX           (LIST_DATA_SIZE + 4)
#define RLS_DATA_PTR_COUNT      (LIST_DATA_SIZE + 5)
#define RLS_DATA_MAX_ADENS      (LIST_DATA_SIZE + 6)
#define RLS_DATA_CUT            (LIST_DATA_SIZE + 7)

/*****************************************************************************/

/***** Nubar data ************************************************************/

#define NUBAR_BLOCK_SIZE              (LIST_DATA_SIZE + 5)

#define NUBAR_DATA_TYPE               (LIST_DATA_SIZE + 0)
#define NUBAR_PTR_POLY_DATA           (LIST_DATA_SIZE + 1)
#define NUBAR_PTR_EGRID               (LIST_DATA_SIZE + 2)
#define NUBAR_PTR_PTS                 (LIST_DATA_SIZE + 3)
#define NUBAR_PTR_PREV_VAL            (LIST_DATA_SIZE + 4)

/*****************************************************************************/

/***** Precursor data ********************************************************/

#define PREC_BLOCK_SIZE               (LIST_DATA_SIZE + 5)

#define PREC_IDX                      (LIST_DATA_SIZE + 0)
#define PREC_LAMBDA                   (LIST_DATA_SIZE + 1)
#define PREC_PTR_EGRID                (LIST_DATA_SIZE + 2)
#define PREC_PTR_PTS                  (LIST_DATA_SIZE + 3)
#define PREC_PTR_ERG                  (LIST_DATA_SIZE + 4)

/*****************************************************************************/

/***** Precursor detector data ***********************************************/

/* There is probably only one of these (so no list needed) */

#define PRECDET_BLOCK_SIZE              (LIST_DATA_SIZE + 27)

#define PRECDET_PTR_STAT                (LIST_DATA_SIZE +  0)
#define PRECDET_PTR_DIM                 (LIST_DATA_SIZE +  1)
#define PRECDET_PTR_MESH                (LIST_DATA_SIZE +  2)
#define PRECDET_PTR_PREC_ARRAY          (LIST_DATA_SIZE +  3)
#define PRECDET_PTR_LAM_ARRAY           (LIST_DATA_SIZE +  4)
#define PRECDET_PTR_REA_ARRAY           (LIST_DATA_SIZE +  5)
#define PRECDET_NT                      (LIST_DATA_SIZE +  6)
#define PRECDET_NG                      (LIST_DATA_SIZE +  7)
#define PRECDET_W_EMIT                  (LIST_DATA_SIZE +  8)
#define PRECDET_MAX_WGT                 (LIST_DATA_SIZE +  9)
#define PRECDET_W_LIVE                  (LIST_DATA_SIZE + 10)
#define PRECDET_W_AVE                   (LIST_DATA_SIZE + 11)
#define PRECDET_N_EMIT                  (LIST_DATA_SIZE + 12)
#define PRECDET_N_LIVE                  (LIST_DATA_SIZE + 13)
#define PRECDET_PTR_LIVE_DET            (LIST_DATA_SIZE + 14)
#define PRECDET_PTR_FILE_DET            (LIST_DATA_SIZE + 15)
#define PRECDET_PTR_OUT_FNAME           (LIST_DATA_SIZE + 16)
#define PRECDET_PTR_IN_FNAME            (LIST_DATA_SIZE + 17)
#define PRECDET_N0                      (LIST_DATA_SIZE + 18)
#define PRECDET_N1                      (LIST_DATA_SIZE + 19)
#define PRECDET_N2                      (LIST_DATA_SIZE + 20)
#define PRECDET_PTR_PREC_DET            (LIST_DATA_SIZE + 21)
#define PRECDET_PTR_PREC_SRC            (LIST_DATA_SIZE + 22)
#define PRECDET_AVE_EMIT                (LIST_DATA_SIZE + 23)
#define PRECDET_SAVE_FRAC_LIVE          (LIST_DATA_SIZE + 24)
#define PRECDET_SAVE_FRAC_PREC          (LIST_DATA_SIZE + 25)
#define PRECDET_PTR_MESH_LIST           (LIST_DATA_SIZE + 26)

/*****************************************************************************/

/***** Photon production data ************************************************/

#define PHOTON_PROD_BLOCK_SIZE         (LIST_DATA_SIZE + 5)

#define PHOTON_PROD_MT                 (LIST_DATA_SIZE + 0)
#define PHOTON_PROD_PTR_EGRID          (LIST_DATA_SIZE + 1)
#define PHOTON_PROD_PTR_PRODXS         (LIST_DATA_SIZE + 2)
#define PHOTON_PROD_PTR_ANG            (LIST_DATA_SIZE + 3)
#define PHOTON_PROD_PTR_ERG            (LIST_DATA_SIZE + 4)

/*****************************************************************************/

/***** RIA simulation ********************************************************/

#define RIA_BLOCK_SIZE                (LIST_DATA_SIZE + PARAM_N_COMMON + 2)

#define RIA_PTR_NAME                  (LIST_DATA_SIZE + PARAM_N_COMMON + 0)
#define RIA_PTR_TME                   (LIST_DATA_SIZE + PARAM_N_COMMON + 1)

/*****************************************************************************/

/***** Energy distribution data **********************************************/

#define ERG_BLOCK_SIZE                (LIST_DATA_SIZE + 8)

#define ERG_PTR_NUCLIDE               (LIST_DATA_SIZE + 0)
#define ERG_PTR_EGRID                 (LIST_DATA_SIZE + 1)
#define ERG_PTR_PROB                  (LIST_DATA_SIZE + 2)
#define ERG_LAW                       (LIST_DATA_SIZE + 3)
#define ERG_INTERP                    (LIST_DATA_SIZE + 4)
#define ERG_PTR_DATA                  (LIST_DATA_SIZE + 5)
#define ERG_NR                        (LIST_DATA_SIZE + 6)
#define ERG_PTR_INTERP                (LIST_DATA_SIZE + 7)

/*****************************************************************************/

/***** Photon distribution data **********************************************/

/* TODO: pit�is poistaa turhat */

#define PHOTON_DIST_BLOCK_SIZE        (LIST_DATA_SIZE + 32)

#define PHOTON_DIST_RAYL_N            (LIST_DATA_SIZE + 0)
#define PHOTON_DIST_RAYL_X2           (LIST_DATA_SIZE + 1)
#define PHOTON_DIST_RAYL_FF2INT       (LIST_DATA_SIZE + 2)
#define PHOTON_DIST_RAYL_FF2          (LIST_DATA_SIZE + 3)
#define PHOTON_DIST_MINC              (LIST_DATA_SIZE + 4)
#define PHOTON_DIST_PTR_VIC           (LIST_DATA_SIZE + 5)
#define PHOTON_DIST_PTR_INC_FF        (LIST_DATA_SIZE + 6)
#define PHOTON_DIST_RAYL_LX2          (LIST_DATA_SIZE + 7)
#define PHOTON_DIST_RAYL_C            (LIST_DATA_SIZE + 8)

#define PHOTON_DIST_N_UOCCUP          (LIST_DATA_SIZE + 9)
#define PHOTON_DIST_PTR_UOCCUP        (LIST_DATA_SIZE + 10)
#define PHOTON_DIST_N_INC_CP          (LIST_DATA_SIZE + 11)
#define PHOTON_DIST_PTR_INC_CPPZ      (LIST_DATA_SIZE + 12)
#define PHOTON_DIST_PTR_INC_CP        (LIST_DATA_SIZE + 13)
#define PHOTON_DIST_PTR_INC_CPINT     (LIST_DATA_SIZE + 14)
#define PHOTON_DIST_INC_PZMAX         (LIST_DATA_SIZE + 15)
#define PHOTON_DIST_PTR_UI            (LIST_DATA_SIZE + 16)
#define PHOTON_DIST_INC_PZMINIDX      (LIST_DATA_SIZE + 17)
#define PHOTON_DIST_PTR_INC_ELNCDF    (LIST_DATA_SIZE + 18)
#define PHOTON_DIST_PTR_INC_EXTA      (LIST_DATA_SIZE + 19)

#define PHOTON_DIST_PTR_PP_MDXSFC     (LIST_DATA_SIZE + 20)
#define PHOTON_DIST_PP_FC             (LIST_DATA_SIZE + 21)
#define PHOTON_DIST_PTR_PP_F0         (LIST_DATA_SIZE + 22)
#define PHOTON_DIST_PP_G0             (LIST_DATA_SIZE + 23)

#define PHOTON_DIST_N_PE_SS           (LIST_DATA_SIZE + 24)
#define PHOTON_DIST_PTR_PE_NXS        (LIST_DATA_SIZE + 25)
#define PHOTON_DIST_PTR_PE_EB         (LIST_DATA_SIZE + 26)
#define PHOTON_DIST_PTR_PE_SSE        (LIST_DATA_SIZE + 27)
#define PHOTON_DIST_PTR_PE_SSXS       (LIST_DATA_SIZE + 28)
#define PHOTON_DIST_N_PE_TOT          (LIST_DATA_SIZE + 29)
#define PHOTON_DIST_PTR_PE_ETOT       (LIST_DATA_SIZE + 30)
#define PHOTON_DIST_PTR_PE_XSTOT      (LIST_DATA_SIZE + 31)

/*****************************************************************************/

/***** Thick-target bremsstrahlung data **************************************/

#define TTB_BLOCK_SIZE                (LIST_DATA_SIZE + 9)

#define TTB_NE                        (LIST_DATA_SIZE + 0)
#define TTB_E                         (LIST_DATA_SIZE + 1)
#define TTB_LE                        (LIST_DATA_SIZE + 2)
#define TTB_BRECDF                    (LIST_DATA_SIZE + 3)
#define TTB_BRPCDF                    (LIST_DATA_SIZE + 4)
#define TTB_BREPDF                    (LIST_DATA_SIZE + 5)
#define TTB_BRPPDF                    (LIST_DATA_SIZE + 6)
#define TTB_LYE                       (LIST_DATA_SIZE + 7)
#define TTB_LYP                       (LIST_DATA_SIZE + 8)

/*****************************************************************************/

/***** Atomic relaxation data ************************************************/

#define RELAX_BLOCK_SIZE              (LIST_DATA_SIZE + 12)

#define RELAX_NSS                     (LIST_DATA_SIZE + 0)
#define RELAX_NTR                     (LIST_DATA_SIZE + 1)
#define RELAX_EBI                     (LIST_DATA_SIZE + 2)
#define RELAX_ELN                     (LIST_DATA_SIZE + 3)
#define RELAX_ETR                     (LIST_DATA_SIZE + 4)
#define RELAX_FTR                     (LIST_DATA_SIZE + 5)
#define RELAX_SUBI                    (LIST_DATA_SIZE + 6)
#define RELAX_SUBJ                    (LIST_DATA_SIZE + 7)
#define RELAX_SUBK                    (LIST_DATA_SIZE + 8)
#define RELAX_D2IMAP                  (LIST_DATA_SIZE + 9)
#define RELAX_MAXVAC                  (LIST_DATA_SIZE + 10)
#define RELAX_ELNCDF                  (LIST_DATA_SIZE + 11)

/*****************************************************************************/

/***** Isomeric branching ratio data *****************************************/

#define BRA_LIST_BLOCK_SIZE           (LIST_DATA_SIZE + 4)

#define BRA_LIST_ZAI                  (LIST_DATA_SIZE + 0)
#define BRA_LIST_MT                   (LIST_DATA_SIZE + 1)
#define BRA_LIST_NS                   (LIST_DATA_SIZE + 2)
#define BRA_LIST_PTR_STATES           (LIST_DATA_SIZE + 3)

#define BRA_STATE_BLOCK_SIZE          (LIST_DATA_SIZE + 5)

#define BRA_STATE_LFS                 (LIST_DATA_SIZE + 0)
#define BRA_STATE_INTT                (LIST_DATA_SIZE + 1)
#define BRA_STATE_NP                  (LIST_DATA_SIZE + 2)
#define BRA_STATE_PTR_ERG             (LIST_DATA_SIZE + 3)
#define BRA_STATE_PTR_FRAC            (LIST_DATA_SIZE + 4)

/*****************************************************************************/

/***** Depletion transmutation list ******************************************/

#define DEP_TRA_BLOCK_SIZE            (LIST_DATA_SIZE + 10)

#define DEP_TRA_PTR_REA               (LIST_DATA_SIZE + 0)
#define DEP_TRA_E0                    (LIST_DATA_SIZE + 1)
#define DEP_TRA_PTR_RESU              (LIST_DATA_SIZE + 2)
#define DEP_TRA_PS1                   (LIST_DATA_SIZE + 3)
#define DEP_TRA_BOS                   (LIST_DATA_SIZE + 4)
#define DEP_TRA_EOS                   (LIST_DATA_SIZE + 5)
#define DEP_TRA_REL                   (LIST_DATA_SIZE + 6)
#define DEP_TRA_AV0                   (LIST_DATA_SIZE + 7)
#define DEP_TRA_AV1                   (LIST_DATA_SIZE + 8)
#define DEP_TRA_AVE                   (LIST_DATA_SIZE + 9)

/*****************************************************************************/

/***** Fission yield data array **********************************************/

/* Data block */

#define FISSION_YIELD_BLOCK_SIZE       (LIST_DATA_SIZE + 6)

#define FISSION_YIELD_PARENT_ZAI       (LIST_DATA_SIZE + 0)
#define FISSION_YIELD_E                (LIST_DATA_SIZE + 1)
#define FISSION_YIELD_NFP              (LIST_DATA_SIZE + 2)
#define FISSION_YIELD_ORIG_NFP         (LIST_DATA_SIZE + 3)
#define FISSION_YIELD_PTR_DISTR        (LIST_DATA_SIZE + 4)
#define FISSION_YIELD_PTR_NEXT         (LIST_DATA_SIZE + 5)

/* Yield entry */

#define FY_BLOCK_SIZE                  (LIST_DATA_SIZE + 4)

#define FY_TGT_ZAI                     (LIST_DATA_SIZE + 0)
#define FY_PTR_TGT                     (LIST_DATA_SIZE + 1)
#define FY_INDEPENDENT_FRAC            (LIST_DATA_SIZE + 2)
#define FY_CUMULATIVE_FRAC             (LIST_DATA_SIZE + 3)

/*****************************************************************************/

/***** Decay array ***********************************************************/

#define DECAY_BLOCK_SIZE               (LIST_DATA_SIZE + 8)

#define DECAY_RTYP1                    (LIST_DATA_SIZE + 0)
#define DECAY_RTYP2                    (LIST_DATA_SIZE + 1)
#define DECAY_RTYP3                    (LIST_DATA_SIZE + 2)
#define DECAY_RTYP4                    (LIST_DATA_SIZE + 3)
#define DECAY_RTYP5                    (LIST_DATA_SIZE + 4)
#define DECAY_RFS                      (LIST_DATA_SIZE + 5)
#define DECAY_Q                        (LIST_DATA_SIZE + 6)
#define DECAY_BR                       (LIST_DATA_SIZE + 7)

/*****************************************************************************/

/***** Radiation spectrum block size *****************************************/

#define RAD_SPEC_BLOCK_SIZE            (LIST_DATA_SIZE + 10)

#define RAD_SPEC_TYPE                  (LIST_DATA_SIZE + 0)
#define RAD_SPEC_AVG_E                 (LIST_DATA_SIZE + 1)
#define RAD_SPEC_DISC_NORM             (LIST_DATA_SIZE + 2)
#define RAD_SPEC_DISC_NE               (LIST_DATA_SIZE + 3)
#define RAD_SPEC_PTR_DISC_E            (LIST_DATA_SIZE + 4)
#define RAD_SPEC_PTR_DISC_RI           (LIST_DATA_SIZE + 5)
#define RAD_SPEC_CONT_NORM             (LIST_DATA_SIZE + 6)
#define RAD_SPEC_CONT_NE               (LIST_DATA_SIZE + 7)
#define RAD_SPEC_PTR_CONT_E            (LIST_DATA_SIZE + 8)
#define RAD_SPEC_PTR_CONT_RI           (LIST_DATA_SIZE + 9)

/*****************************************************************************/

/***** Thermal scattering library data ***************************************/

#define THERM_BLOCK_SIZE              (LIST_DATA_SIZE + PARAM_N_COMMON +  10)

#define THERM_OPTIONS                 (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define THERM_PTR_ALIAS               (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define THERM_ZA                      (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define THERM_T                       (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define THERM_PTR_SAB                 (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define THERM_PTR_THERM               (LIST_DATA_SIZE + PARAM_N_COMMON +  5)
#define THERM_PTR_COMP                (LIST_DATA_SIZE + PARAM_N_COMMON +  6)
#define THERM_INTERP_MODE             (LIST_DATA_SIZE + PARAM_N_COMMON +  7)
#define THERM_OTF_MIN_TEMP            (LIST_DATA_SIZE + PARAM_N_COMMON +  8)
#define THERM_OTF_MAX_TEMP            (LIST_DATA_SIZE + PARAM_N_COMMON +  9)

/*****************************************************************************/

/***** S(a,b) nuclide data  **************************************************/

#define SAB_BLOCK_SIZE               (LIST_DATA_SIZE + 6)

#define SAB_PTR_NAME                 (LIST_DATA_SIZE + 0)
#define SAB_T                        (LIST_DATA_SIZE + 1)
#define SAB_PTR_ISO                  (LIST_DATA_SIZE + 2)
#define SAB_FRAC                     (LIST_DATA_SIZE + 3)
#define SAB_PTR_PREV_FRAC            (LIST_DATA_SIZE + 4)
#define SAB_PTR_PREV_SAB1            (LIST_DATA_SIZE + 5)

/*****************************************************************************/

/***** Unresolved resonance probability table data ***************************/

#define URES_BLOCK_SIZE               (LIST_DATA_SIZE + 10)

#define URES_PTR_EGRID                (LIST_DATA_SIZE +  0)
#define URES_NP                       (LIST_DATA_SIZE +  1)
#define URES_INT                      (LIST_DATA_SIZE +  2)
#define URES_PTR_PROB                 (LIST_DATA_SIZE +  3)
#define URES_IFF                      (LIST_DATA_SIZE +  4)
#define URES_PTR_FACT                 (LIST_DATA_SIZE +  5)
#define URES_PTR_MAXF                 (LIST_DATA_SIZE +  6)
#define URES_PTR_RND                  (LIST_DATA_SIZE +  7)
#define URES_PTR_RND_CHK              (LIST_DATA_SIZE +  8)
#define URES_PTR_PREV_FACT            (LIST_DATA_SIZE +  9)

/*****************************************************************************/

/***** Angular distribution **************************************************/

#define ANG_BLOCK_SIZE                (LIST_DATA_SIZE + 5)

#define ANG_PTR_EGRID                 (LIST_DATA_SIZE + 0)
#define ANG_TYPE                      (LIST_DATA_SIZE + 1)
#define ANG_PTR_D0                    (LIST_DATA_SIZE + 2)
#define ANG_BINS                      (LIST_DATA_SIZE + 3)
#define ANG_INTT                      (LIST_DATA_SIZE + 4)

/*****************************************************************************/

/***** Surface ***************************************************************/

/* Surface types */

#define SURFACE_TYPES   39

#define SURF_CYL         1
#define SURF_PX          2
#define SURF_PY          3
#define SURF_PZ          4
#define SURF_INF         5
#define SURF_SQC         6
#define SURF_HEXYC       7
#define SURF_HEXXC       8
#define SURF_SPH         9
#define SURF_CROSS      10
#define SURF_PAD        11
#define SURF_CUBE       12
#define SURF_CONE       13
#define SURF_SVC        14
#define SURF_CUBOID     15
#define SURF_HEXYPRISM  16
#define SURF_HEXXPRISM  17
#define SURF_DODE       18
#define SURF_OCTA       19
#define SURF_ASTRA      20
#define SURF_PLANE      21
#define SURF_QUADRATIC  22
#define SURF_CYLX       23
#define SURF_CYLY       24
#define SURF_CYLZ       25
#define SURF_GCROSS     26
#define SURF_PPD        27
#define SURF_RECT       28
#define SURF_CKX        29
#define SURF_CKY        30
#define SURF_CKZ        31
#define SURF_X          32
#define SURF_Y          33
#define SURF_Z          34
#define SURF_TORX       35
#define SURF_TORY       36
#define SURF_TORZ       37
#define SURF_CYLV       38
#define SURF_USER       39

/* Data block */

#define SURFACE_BLOCK_SIZE  (LIST_DATA_SIZE + PARAM_N_COMMON + 6)

#define SURFACE_PTR_NAME    (LIST_DATA_SIZE + PARAM_N_COMMON + 0)
#define SURFACE_OPTIONS     (LIST_DATA_SIZE + PARAM_N_COMMON + 1)
#define SURFACE_TYPE        (LIST_DATA_SIZE + PARAM_N_COMMON + 2)
#define SURFACE_PTR_PARAMS  (LIST_DATA_SIZE + PARAM_N_COMMON + 3)
#define SURFACE_N_PARAMS    (LIST_DATA_SIZE + PARAM_N_COMMON + 4)
#define SURFACE_PTR_TRANS   (LIST_DATA_SIZE + PARAM_N_COMMON + 5)

/* List of outer boundaries */

#define BOUNDS_BLOCK_SIZE   (LIST_DATA_SIZE + 1)

#define BOUNDS_PTR_SURF     (LIST_DATA_SIZE + 0)

/*****************************************************************************/

/***** Cell ******************************************************************/

/* Cell types */

#define CELL_TYPE_MAT      1
#define CELL_TYPE_VOID     2
#define CELL_TYPE_OUTSIDE  3
#define CELL_TYPE_FILL     4

/* Data block (t�n koko voi olla merkitt�v� tekij� unstructured */
/* mesh -tyyppisiss� geometrioissa. */

#define CELL_BLOCK_SIZE          (LIST_DATA_SIZE + PARAM_N_COMMON + 19)

#define CELL_PTR_NAME            (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define CELL_TYPE                (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define CELL_OPTIONS             (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define CELL_PTR_UNI             (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define CELL_PTR_MAT             (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define CELL_PTR_FILL            (LIST_DATA_SIZE + PARAM_N_COMMON +  5)
#define CELL_PTR_SURF_LIST       (LIST_DATA_SIZE + PARAM_N_COMMON +  6)
#define CELL_PTR_SURF_COMP       (LIST_DATA_SIZE + PARAM_N_COMMON +  7)
#define CELL_PTR_SURF_INSC       (LIST_DATA_SIZE + PARAM_N_COMMON +  8)
#define CELL_COL_COUNT           (LIST_DATA_SIZE + PARAM_N_COMMON +  9)
#define CELL_VOLUME              (LIST_DATA_SIZE + PARAM_N_COMMON + 10)
#define CELL_VOL_COUNT           (LIST_DATA_SIZE + PARAM_N_COMMON + 11)
#define CELL_PTR_REG_MAT         (LIST_DATA_SIZE + PARAM_N_COMMON + 12)
#define CELL_UMSH_IDX            (LIST_DATA_SIZE + PARAM_N_COMMON + 13)
#define CELL_UMSH_DET_BIN        (LIST_DATA_SIZE + PARAM_N_COMMON + 14)
#define CELL_PTR_BC_SURF         (LIST_DATA_SIZE + PARAM_N_COMMON + 15)
#define CELL_PTR_DETBIN          (LIST_DATA_SIZE + PARAM_N_COMMON + 16)
#define CELL_PTR_PREV_TET        (LIST_DATA_SIZE + PARAM_N_COMMON + 17)
#define CELL_PTR_TRANS           (LIST_DATA_SIZE + PARAM_N_COMMON + 18)

/* Universe cell list */

#define CELL_LIST_BLOCK_SIZE  (LIST_DATA_SIZE + 3)

#define CELL_LIST_PTR_CELL    (LIST_DATA_SIZE + 0)
#define CELL_LIST_PTR_COUNT   (LIST_DATA_SIZE + 1)
#define CELL_LIST_REG_IDX     (LIST_DATA_SIZE + 2)

/* List of intersections */

#define CELL_INSC_BLOCK_SIZE           (LIST_DATA_SIZE + 4)

#define CELL_INSC_PTR_SURF             (LIST_DATA_SIZE + 0)
#define CELL_INSC_SIDE                 (LIST_DATA_SIZE + 1)
#define CELL_INSC_PTR_OUT_COUNT        (LIST_DATA_SIZE + 2)
#define CELL_INSC_PTR_NEXT_TET_CELL    (LIST_DATA_SIZE + 3)

/*****************************************************************************/

/***** Super-imposed cell mesh ***********************************************/

#define CELL_MESH_BLOCK_SIZE      (LIST_DATA_SIZE + 12)

#define CELL_MESH_PTR_UNI         (LIST_DATA_SIZE +  0)
#define CELL_MESH_PTR_MESH        (LIST_DATA_SIZE +  1)
#define CELL_MESH_TYPE            (LIST_DATA_SIZE +  2)
#define CELL_MESH_N0              (LIST_DATA_SIZE +  3)
#define CELL_MESH_N1              (LIST_DATA_SIZE +  4)
#define CELL_MESH_N2              (LIST_DATA_SIZE +  5)
#define CELL_MESH_MIN0            (LIST_DATA_SIZE +  6)
#define CELL_MESH_MAX0            (LIST_DATA_SIZE +  7)
#define CELL_MESH_MIN1            (LIST_DATA_SIZE +  8)
#define CELL_MESH_MAX1            (LIST_DATA_SIZE +  9)
#define CELL_MESH_MIN2            (LIST_DATA_SIZE + 10)
#define CELL_MESH_MAX2            (LIST_DATA_SIZE + 11)

/*****************************************************************************/

/***** Nest ******************************************************************/

#define NEST_BLOCK_SIZE         (LIST_DATA_SIZE + PARAM_N_COMMON + 8)

#define NEST_PTR_NAME           (LIST_DATA_SIZE + PARAM_N_COMMON + 0)
#define NEST_OPTIONS            (LIST_DATA_SIZE + PARAM_N_COMMON + 1)
#define NEST_PTR_UNI            (LIST_DATA_SIZE + PARAM_N_COMMON + 2)
#define NEST_PTR_REGIONS        (LIST_DATA_SIZE + PARAM_N_COMMON + 3)
#define NEST_TYPE               (LIST_DATA_SIZE + PARAM_N_COMMON + 4)
#define NEST_PTR_COL_REG        (LIST_DATA_SIZE + PARAM_N_COMMON + 5)
#define NEST_COUNT              (LIST_DATA_SIZE + PARAM_N_COMMON + 6)
#define NEST_PTR_BC_SURF        (LIST_DATA_SIZE + PARAM_N_COMMON + 7)

#define NEST_REG_BLOCK_SIZE     (LIST_DATA_SIZE + 6)

#define NEST_REG_PTR_FILL       (LIST_DATA_SIZE + 0)
#define NEST_REG_TMP_PTR        (LIST_DATA_SIZE + 1)
#define NEST_REG_PTR_SURF_IN    (LIST_DATA_SIZE + 2)
#define NEST_REG_PTR_SURF_OUT   (LIST_DATA_SIZE + 3)
#define NEST_REG_PTR_CELL       (LIST_DATA_SIZE + 4)
#define NEST_REG_IDX            (LIST_DATA_SIZE + 5)

/*****************************************************************************/

/***** Coordinate transformation *********************************************/

#define TRANS_BLOCK_SIZE  (LIST_DATA_SIZE + PARAM_N_COMMON + 38)

#define TRANS_PTR_NAME    (LIST_DATA_SIZE + PARAM_N_COMMON + 0)
#define TRANS_TYPE        (LIST_DATA_SIZE + PARAM_N_COMMON + 1)
#define TRANS_PTR_UNI     (LIST_DATA_SIZE + PARAM_N_COMMON + 2)
#define TRANS_PTR_SURF    (LIST_DATA_SIZE + PARAM_N_COMMON + 3)
#define TRANS_PTR_CELL    (LIST_DATA_SIZE + PARAM_N_COMMON + 4)
#define TRANS_OPTIONS     (LIST_DATA_SIZE + PARAM_N_COMMON + 5)
#define TRANS_ROT         (LIST_DATA_SIZE + PARAM_N_COMMON + 6)
#define TRANS_X0          (LIST_DATA_SIZE + PARAM_N_COMMON + 7)
#define TRANS_Y0          (LIST_DATA_SIZE + PARAM_N_COMMON + 8)
#define TRANS_Z0          (LIST_DATA_SIZE + PARAM_N_COMMON + 9)
#define TRANS_RX1         (LIST_DATA_SIZE + PARAM_N_COMMON + 10)
#define TRANS_RX2         (LIST_DATA_SIZE + PARAM_N_COMMON + 11)
#define TRANS_RX3         (LIST_DATA_SIZE + PARAM_N_COMMON + 12)
#define TRANS_RX4         (LIST_DATA_SIZE + PARAM_N_COMMON + 13)
#define TRANS_RX5         (LIST_DATA_SIZE + PARAM_N_COMMON + 14)
#define TRANS_RX6         (LIST_DATA_SIZE + PARAM_N_COMMON + 15)
#define TRANS_RX7         (LIST_DATA_SIZE + PARAM_N_COMMON + 16)
#define TRANS_RX8         (LIST_DATA_SIZE + PARAM_N_COMMON + 17)
#define TRANS_RX9         (LIST_DATA_SIZE + PARAM_N_COMMON + 18)
#define TRANS_RY1         (LIST_DATA_SIZE + PARAM_N_COMMON + 19)
#define TRANS_RY2         (LIST_DATA_SIZE + PARAM_N_COMMON + 20)
#define TRANS_RY3         (LIST_DATA_SIZE + PARAM_N_COMMON + 21)
#define TRANS_RY4         (LIST_DATA_SIZE + PARAM_N_COMMON + 22)
#define TRANS_RY5         (LIST_DATA_SIZE + PARAM_N_COMMON + 23)
#define TRANS_RY6         (LIST_DATA_SIZE + PARAM_N_COMMON + 24)
#define TRANS_RY7         (LIST_DATA_SIZE + PARAM_N_COMMON + 25)
#define TRANS_RY8         (LIST_DATA_SIZE + PARAM_N_COMMON + 26)
#define TRANS_RY9         (LIST_DATA_SIZE + PARAM_N_COMMON + 27)
#define TRANS_RZ1         (LIST_DATA_SIZE + PARAM_N_COMMON + 28)
#define TRANS_RZ2         (LIST_DATA_SIZE + PARAM_N_COMMON + 29)
#define TRANS_RZ3         (LIST_DATA_SIZE + PARAM_N_COMMON + 30)
#define TRANS_RZ4         (LIST_DATA_SIZE + PARAM_N_COMMON + 31)
#define TRANS_RZ5         (LIST_DATA_SIZE + PARAM_N_COMMON + 32)
#define TRANS_RZ6         (LIST_DATA_SIZE + PARAM_N_COMMON + 33)
#define TRANS_RZ7         (LIST_DATA_SIZE + PARAM_N_COMMON + 34)
#define TRANS_RZ8         (LIST_DATA_SIZE + PARAM_N_COMMON + 35)
#define TRANS_RZ9         (LIST_DATA_SIZE + PARAM_N_COMMON + 36)
#define TRANS_PTR_LVL     (LIST_DATA_SIZE + PARAM_N_COMMON + 37)

/*****************************************************************************/

/***** Geometry plotter ******************************************************/

#define GPL_BLOCK_SIZE  (LIST_DATA_SIZE + PARAM_N_COMMON + 15)

#define GPL_IDX         (LIST_DATA_SIZE + PARAM_N_COMMON + 0)
#define GPL_GEOM        (LIST_DATA_SIZE + PARAM_N_COMMON + 1)
#define GPL_PIX_X       (LIST_DATA_SIZE + PARAM_N_COMMON + 2)
#define GPL_PIX_Y       (LIST_DATA_SIZE + PARAM_N_COMMON + 3)
#define GPL_POS         (LIST_DATA_SIZE + PARAM_N_COMMON + 4)
#define GPL_XMIN        (LIST_DATA_SIZE + PARAM_N_COMMON + 5)
#define GPL_XMAX        (LIST_DATA_SIZE + PARAM_N_COMMON + 6)
#define GPL_YMIN        (LIST_DATA_SIZE + PARAM_N_COMMON + 7)
#define GPL_YMAX        (LIST_DATA_SIZE + PARAM_N_COMMON + 8)
#define GPL_ZMIN        (LIST_DATA_SIZE + PARAM_N_COMMON + 9)
#define GPL_ZMAX        (LIST_DATA_SIZE + PARAM_N_COMMON + 10)
#define GPL_PLOT_BOUND  (LIST_DATA_SIZE + PARAM_N_COMMON + 11)
#define GPL_IMP_SCALE   (LIST_DATA_SIZE + PARAM_N_COMMON + 12)
#define GPL_IMP_MIN     (LIST_DATA_SIZE + PARAM_N_COMMON + 13)
#define GPL_IMP_MAX     (LIST_DATA_SIZE + PARAM_N_COMMON + 14)

/*****************************************************************************/

/***** Pebble bed geometry ***************************************************/

#define PBED_BLOCK_SIZE               (LIST_DATA_SIZE + PARAM_N_COMMON + 12)

#define PBED_OPTIONS                  (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define PBED_PTR_FNAME                (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define PBED_N_PEBBLES                (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define PBED_PTR_NAME                 (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define PBED_PTR_UNI                  (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define PBED_PTR_PEBBLES              (LIST_DATA_SIZE + PARAM_N_COMMON +  5)
#define PBED_PTR_SEARCH_MESH          (LIST_DATA_SIZE + PARAM_N_COMMON +  6)
#define PBED_PTR_BG_UNIV              (LIST_DATA_SIZE + PARAM_N_COMMON +  7)
#define PBED_CALC_RESULTS             (LIST_DATA_SIZE + PARAM_N_COMMON +  8)
#define PBED_PTR_COL_PEBBLE           (LIST_DATA_SIZE + PARAM_N_COMMON +  9)
#define PBED_PTR_PEBBLE_TYPES         (LIST_DATA_SIZE + PARAM_N_COMMON + 10)
#define PBED_PTR_POW                  (LIST_DATA_SIZE + PARAM_N_COMMON + 11)

#define PEBBLE_BLOCK_SIZE              (LIST_DATA_SIZE + 6)

#define PEBBLE_PTR_UNIV                (LIST_DATA_SIZE + 0)
#define PEBBLE_X0                      (LIST_DATA_SIZE + 1)
#define PEBBLE_Y0                      (LIST_DATA_SIZE + 2)
#define PEBBLE_Z0                      (LIST_DATA_SIZE + 3)
#define PEBBLE_RAD                     (LIST_DATA_SIZE + 4)
#define PEBBLE_IDX                     (LIST_DATA_SIZE + 5)

#define PEBTYPE_BLOCK_SIZE             (LIST_DATA_SIZE + 2)

#define PEBTYPE_PTR_UNIV               (LIST_DATA_SIZE + 0)
#define PEBTYPE_COUNT                  (LIST_DATA_SIZE + 1)

/*****************************************************************************/

/***** Unstructured mesh geometry ********************************************/

#define UMSH_BLOCK_SIZE               (LIST_DATA_SIZE + PARAM_N_COMMON + 13)

#define UMSH_PTR_NAME                 (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define UMSH_PTR_UNI                  (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define UMSH_PTR_BG_UNIV              (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define UMSH_N_ORIG_CELLS             (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define UMSH_N_CELLS                  (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define UMSH_OPTIONS                  (LIST_DATA_SIZE + PARAM_N_COMMON +  5)
#define UMSH_PTR_POINTS_FNAME         (LIST_DATA_SIZE + PARAM_N_COMMON +  6)
#define UMSH_PTR_FACES_FNAME          (LIST_DATA_SIZE + PARAM_N_COMMON +  7)
#define UMSH_PTR_OWNER_FNAME          (LIST_DATA_SIZE + PARAM_N_COMMON +  8)
#define UMSH_PTR_NEIGHBOUR_FNAME      (LIST_DATA_SIZE + PARAM_N_COMMON +  9)
#define UMSH_PTR_MATERIALS_FNAME      (LIST_DATA_SIZE + PARAM_N_COMMON + 10)
#define UMSH_PTR_IFC                  (LIST_DATA_SIZE + PARAM_N_COMMON + 11)
#define UMSH_PTR_FNAME                (LIST_DATA_SIZE + PARAM_N_COMMON + 12)

/*****************************************************************************/

/***** STL based geometry ****************************************************/

#define STL_BLOCK_SIZE            (LIST_DATA_SIZE + PARAM_N_COMMON + 20)

#define STL_PTR_NAME              (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define STL_PTR_UNI               (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define STL_OPTIONS               (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define STL_PTR_BG_UNIV           (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define STL_PTR_FILES             (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define STL_PTR_BODIES            (LIST_DATA_SIZE + PARAM_N_COMMON +  5)
#define STL_PTR_SOLIDS            (LIST_DATA_SIZE + PARAM_N_COMMON +  6)
#define STL_SEARCH_MESH_ADA_SPLIT (LIST_DATA_SIZE + PARAM_N_COMMON +  7)
#define STL_SEARCH_MESH_PTR_SZ    (LIST_DATA_SIZE + PARAM_N_COMMON +  8)
#define STL_PTR_FACET_MESH        (LIST_DATA_SIZE + PARAM_N_COMMON +  9)
#define STL_PTR_SOLID_MESH        (LIST_DATA_SIZE + PARAM_N_COMMON + 10)
#define STL_XMIN                  (LIST_DATA_SIZE + PARAM_N_COMMON + 11)
#define STL_XMAX                  (LIST_DATA_SIZE + PARAM_N_COMMON + 12)
#define STL_YMIN                  (LIST_DATA_SIZE + PARAM_N_COMMON + 13)
#define STL_YMAX                  (LIST_DATA_SIZE + PARAM_N_COMMON + 14)
#define STL_ZMIN                  (LIST_DATA_SIZE + PARAM_N_COMMON + 15)
#define STL_ZMAX                  (LIST_DATA_SIZE + PARAM_N_COMMON + 16)
#define STL_MERGE_RAD             (LIST_DATA_SIZE + PARAM_N_COMMON + 17)
#define STL_SEARCH_MESH_CELLS     (LIST_DATA_SIZE + PARAM_N_COMMON + 18)
#define STL_SEARCH_MODE           (LIST_DATA_SIZE + PARAM_N_COMMON + 19)

#define STL_FILE_BLOCK_SIZE       (LIST_DATA_SIZE + 6)

#define STL_FILE_PTR_FNAME        (LIST_DATA_SIZE + 0)
#define STL_FILE_PTR_SNAME        (LIST_DATA_SIZE + 1)
#define STL_FILE_SCALING          (LIST_DATA_SIZE + 2)
#define STL_FILE_X0               (LIST_DATA_SIZE + 3)
#define STL_FILE_Y0               (LIST_DATA_SIZE + 4)
#define STL_FILE_Z0               (LIST_DATA_SIZE + 5)

#define STL_BODY_BLOCK_SIZE       (LIST_DATA_SIZE + 8)

#define STL_BODY_PTR_BNAME        (LIST_DATA_SIZE + 0)
#define STL_BODY_PTR_CNAME        (LIST_DATA_SIZE + 1)
#define STL_BODY_PTR_MNAME        (LIST_DATA_SIZE + 2)
#define STL_BODY_PTR_FILL         (LIST_DATA_SIZE + 3)
#define STL_BODY_PTR_CELL         (LIST_DATA_SIZE + 4)
#define STL_BODY_N_PARTS          (LIST_DATA_SIZE + 5)
#define STL_BODY_N_POINTS         (LIST_DATA_SIZE + 6)
#define STL_BODY_N_FACETS         (LIST_DATA_SIZE + 7)

#define STL_SOLID_BLOCK_SIZE      (LIST_DATA_SIZE + 14)

#define STL_SOLID_PTR_STL_NAME    (LIST_DATA_SIZE +  0)
#define STL_SOLID_PTR_FNAME       (LIST_DATA_SIZE +  1)
#define STL_SOLID_PTR_CELL        (LIST_DATA_SIZE +  2)
#define STL_SOLID_N_POINTS        (LIST_DATA_SIZE +  3)
#define STL_SOLID_N_FACETS        (LIST_DATA_SIZE +  4)
#define STL_SOLID_PTR_FACETS      (LIST_DATA_SIZE +  5)
#define STL_SOLID_XMIN            (LIST_DATA_SIZE +  6)
#define STL_SOLID_XMAX            (LIST_DATA_SIZE +  7)
#define STL_SOLID_YMIN            (LIST_DATA_SIZE +  8)
#define STL_SOLID_YMAX            (LIST_DATA_SIZE +  9)
#define STL_SOLID_ZMIN            (LIST_DATA_SIZE + 10)
#define STL_SOLID_ZMAX            (LIST_DATA_SIZE + 11)
#define STL_SOLID_PTR_BODY        (LIST_DATA_SIZE + 12)
#define STL_SOLID_REG_IDX         (LIST_DATA_SIZE + 13)

#define STL_POINT_BLOCK_SIZE       (LIST_DATA_SIZE + 3)

#define STL_POINT_X                (LIST_DATA_SIZE + 0)
#define STL_POINT_Y                (LIST_DATA_SIZE + 1)
#define STL_POINT_Z                (LIST_DATA_SIZE + 2)

#define STL_FACET_BLOCK_SIZE      13

#define STL_FACET_NORM_U           0
#define STL_FACET_NORM_V           1
#define STL_FACET_NORM_W           2
#define STL_FACET_PTR_PT1          3
#define STL_FACET_PTR_PT2          4
#define STL_FACET_PTR_PT3          5
#define STL_FACET_PTR_SOLID        6
#define STL_FACET_XMIN             7
#define STL_FACET_XMAX             8
#define STL_FACET_YMIN             9
#define STL_FACET_YMAX            10
#define STL_FACET_ZMIN            11
#define STL_FACET_ZMAX            12

/* STL Search mesh content (pointers to facet and solid lists) */

#define STL_SEARCH_MESH_CONTENT_BLOCK_SIZE  2

#define STL_SEARCH_MESH_CONTENT_PTR_FACETS  0
#define STL_SEARCH_MESH_CONTENT_PTR_SOLIDS  1

/*****************************************************************************/

/***** Search mesh ***********************************************************/

/* List for search mesh */

#define SEARCH_MESH_CELL_BLOCK_SIZE    (LIST_DATA_SIZE + 2)

#define SEARCH_MESH_CELL_CONTENT       (LIST_DATA_SIZE + 0)
#define SEARCH_MESH_PTR_CELL_COUNT     (LIST_DATA_SIZE + 1)

/*****************************************************************************/

/***** Multi-physics interface ***********************************************/

#define IFC_BLOCK_SIZE                (LIST_DATA_SIZE + PARAM_N_COMMON + 78)

#define IFC_IDX                       (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define IFC_DIM                       (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define IFC_PTR_INPUT_FNAME           (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define IFC_PTR_OUTPUT_FNAME          (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define IFC_TYPE                      (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define IFC_EXCL_RAD                  (LIST_DATA_SIZE + PARAM_N_COMMON +  5)
#define IFC_EXP                       (LIST_DATA_SIZE + PARAM_N_COMMON +  6)
#define IFC_PTR_MAT                   (LIST_DATA_SIZE + PARAM_N_COMMON +  7)
#define IFC_MIN_DENSITY               (LIST_DATA_SIZE + PARAM_N_COMMON +  8)
#define IFC_MAX_DENSITY               (LIST_DATA_SIZE + PARAM_N_COMMON +  9)
#define IFC_MIN_TEMP                  (LIST_DATA_SIZE + PARAM_N_COMMON + 10)
#define IFC_MAX_TEMP                  (LIST_DATA_SIZE + PARAM_N_COMMON + 11)
#define IFC_SEARCH_MESH_NX            (LIST_DATA_SIZE + PARAM_N_COMMON + 12)
#define IFC_SEARCH_MESH_NY            (LIST_DATA_SIZE + PARAM_N_COMMON + 13)
#define IFC_SEARCH_MESH_NZ            (LIST_DATA_SIZE + PARAM_N_COMMON + 14)
#define IFC_SEARCH_MESH_ADA_SPLIT     (LIST_DATA_SIZE + PARAM_N_COMMON + 15)
#define IFC_SEARCH_MESH_ADA_PTR_SZ    (LIST_DATA_SIZE + PARAM_N_COMMON + 16)
#define IFC_PTR_SEARCH_MESH           (LIST_DATA_SIZE + PARAM_N_COMMON + 17)
#define IFC_NP                        (LIST_DATA_SIZE + PARAM_N_COMMON + 18)
#define IFC_PTR_OUTPUT                (LIST_DATA_SIZE + PARAM_N_COMMON + 19)
#define IFC_NZ                        (LIST_DATA_SIZE + PARAM_N_COMMON + 20)
#define IFC_ZMIN                      (LIST_DATA_SIZE + PARAM_N_COMMON + 21)
#define IFC_ZMAX                      (LIST_DATA_SIZE + PARAM_N_COMMON + 22)
#define IFC_NR                        (LIST_DATA_SIZE + PARAM_N_COMMON + 23)
#define IFC_PTR_OUT                   (LIST_DATA_SIZE + PARAM_N_COMMON + 24)
#define IFC_PTR_POINTS                (LIST_DATA_SIZE + PARAM_N_COMMON + 25)
#define IFC_FUNC_NP                   (LIST_DATA_SIZE + PARAM_N_COMMON + 26)
#define IFC_FUNC_PTR_PARAM            (LIST_DATA_SIZE + PARAM_N_COMMON + 27)
#define IFC_CALC_OUTPUT               (LIST_DATA_SIZE + PARAM_N_COMMON + 28)
#define IFC_PTR_TET_MSH               (LIST_DATA_SIZE + PARAM_N_COMMON + 29)
#define IFC_PTR_TET_MSH_PRNTS         (LIST_DATA_SIZE + PARAM_N_COMMON + 30)
#define IFC_MESH_XMIN                 (LIST_DATA_SIZE + PARAM_N_COMMON + 31)
#define IFC_MESH_XMAX                 (LIST_DATA_SIZE + PARAM_N_COMMON + 32)
#define IFC_MESH_YMIN                 (LIST_DATA_SIZE + PARAM_N_COMMON + 33)
#define IFC_MESH_YMAX                 (LIST_DATA_SIZE + PARAM_N_COMMON + 34)
#define IFC_MESH_ZMIN                 (LIST_DATA_SIZE + PARAM_N_COMMON + 35)
#define IFC_MESH_ZMAX                 (LIST_DATA_SIZE + PARAM_N_COMMON + 36)
#define IFC_PTR_SCORE                 (LIST_DATA_SIZE + PARAM_N_COMMON + 37)
#define IFC_PTR_STAT                  (LIST_DATA_SIZE + PARAM_N_COMMON + 38)
#define IFC_PTR_STAT_REL              (LIST_DATA_SIZE + PARAM_N_COMMON + 39)
#define IFC_PTR_STAT_GRD              (LIST_DATA_SIZE + PARAM_N_COMMON + 40)
#define IFC_PTR_STAT_VOL              (LIST_DATA_SIZE + PARAM_N_COMMON + 41)
#define IFC_STAT_NREG                 (LIST_DATA_SIZE + PARAM_N_COMMON + 42)
#define IFC_PTR_FUEP                  (LIST_DATA_SIZE + PARAM_N_COMMON + 43)
#define IFC_PTR_PREV_CELL             (LIST_DATA_SIZE + PARAM_N_COMMON + 44)
#define IFC_PTR_PREV_COL_CELL         (LIST_DATA_SIZE + PARAM_N_COMMON + 45)
#define IFC_OUT_PTR_LIM               (LIST_DATA_SIZE + PARAM_N_COMMON + 46)
#define IFC_OUT_PTR_Z                 (LIST_DATA_SIZE + PARAM_N_COMMON + 47)
#define IFC_OUT_PTR_PHI               (LIST_DATA_SIZE + PARAM_N_COMMON + 48)
#define IFC_OUT_PTR_R2                (LIST_DATA_SIZE + PARAM_N_COMMON + 49)
#define IFC_PTR_OF_PFILE              (LIST_DATA_SIZE + PARAM_N_COMMON + 50)
#define IFC_PTR_OF_FFILE              (LIST_DATA_SIZE + PARAM_N_COMMON + 51)
#define IFC_PTR_OF_OFILE              (LIST_DATA_SIZE + PARAM_N_COMMON + 52)
#define IFC_PTR_OF_NFILE              (LIST_DATA_SIZE + PARAM_N_COMMON + 53)
#define IFC_PTR_OF_RFILE              (LIST_DATA_SIZE + PARAM_N_COMMON + 54)
#define IFC_PTR_OF_TFILE              (LIST_DATA_SIZE + PARAM_N_COMMON + 55)
#define IFC_PTR_OF_BATCHES            (LIST_DATA_SIZE + PARAM_N_COMMON + 56)
#define IFC_MEM_SIZE                  (LIST_DATA_SIZE + PARAM_N_COMMON + 57)
#define IFC_PTR_OF_MFILE              (LIST_DATA_SIZE + PARAM_N_COMMON + 58)
#define IFC_MSH_N_CELLS               (LIST_DATA_SIZE + PARAM_N_COMMON + 59)
#define IFC_MSH_N_CHILD_CELLS         (LIST_DATA_SIZE + PARAM_N_COMMON + 60)
#define IFC_PTR_MNAMES                (LIST_DATA_SIZE + PARAM_N_COMMON + 61)
#define IFC_PTR_SURF_LIST             (LIST_DATA_SIZE + PARAM_N_COMMON + 62)
#define IFC_PTR_OWNR_LIST             (LIST_DATA_SIZE + PARAM_N_COMMON + 63)
#define IFC_PTR_NBR_LIST              (LIST_DATA_SIZE + PARAM_N_COMMON + 64)
#define IFC_PTR_GCELL_LIST            (LIST_DATA_SIZE + PARAM_N_COMMON + 65)
#define IFC_PTR_POINT_LIST            (LIST_DATA_SIZE + PARAM_N_COMMON + 66)
#define IFC_PTR_SURF_LIST_PRNTS       (LIST_DATA_SIZE + PARAM_N_COMMON + 67)
#define IFC_PTR_OWNR_LIST_PRNTS       (LIST_DATA_SIZE + PARAM_N_COMMON + 68)
#define IFC_PTR_NBR_LIST_PRNTS        (LIST_DATA_SIZE + PARAM_N_COMMON + 69)
#define IFC_PTR_POINT_LIST_PRNTS      (LIST_DATA_SIZE + PARAM_N_COMMON + 70)
#define IFC_PTR_PRNT_CELL_CP_LIST     (LIST_DATA_SIZE + PARAM_N_COMMON + 71)
#define IFC_PTR_PRNT_FACE_CP_LIST     (LIST_DATA_SIZE + PARAM_N_COMMON + 72)
#define IFC_PTR_CELL_CP_LIST          (LIST_DATA_SIZE + PARAM_N_COMMON + 73)
#define IFC_NC_PRNTS                  (LIST_DATA_SIZE + PARAM_N_COMMON + 74)
#define IFC_NF_PRNTS                  (LIST_DATA_SIZE + PARAM_N_COMMON + 75)
#define IFC_NC                        (LIST_DATA_SIZE + PARAM_N_COMMON + 76)
#define IFC_NF                        (LIST_DATA_SIZE + PARAM_N_COMMON + 77)

/* Points */

#define IFC_PT_LIST_BLOCK_SIZE        (LIST_DATA_SIZE + 5)

#define IFC_PT_X                      (LIST_DATA_SIZE + 0)
#define IFC_PT_Y                      (LIST_DATA_SIZE + 1)
#define IFC_PT_Z                      (LIST_DATA_SIZE + 2)
#define IFC_PT_DF                     (LIST_DATA_SIZE + 3)
#define IFC_PT_TMP                    (LIST_DATA_SIZE + 4)

/* Type 2: Regular mesh */

#define IFC_REG_MSH_LIST_BLOCK_SIZE    (LIST_DATA_SIZE + 2)

#define IFC_REG_MSH_DF                 (LIST_DATA_SIZE + 0)
#define IFC_REG_MSH_TMP                (LIST_DATA_SIZE + 1)

/* Type 4: Unstructured tetrahedral mesh based distribution */

#define IFC_TET_MSH_LIST_BLOCK_SIZE    (LIST_DATA_SIZE + 15)

#define IFC_TET_MSH_IDX                (LIST_DATA_SIZE +  0)
#define IFC_TET_MSH_PTR_CELL           (LIST_DATA_SIZE +  1)
#define IFC_TET_MSH_PTR_PARENT         (LIST_DATA_SIZE +  2)
#define IFC_TET_MSH_DF                 (LIST_DATA_SIZE +  3)
#define IFC_TET_MSH_TMP                (LIST_DATA_SIZE +  4)
#define IFC_TET_MSH_XMIN               (LIST_DATA_SIZE +  5)
#define IFC_TET_MSH_XMAX               (LIST_DATA_SIZE +  6)
#define IFC_TET_MSH_YMIN               (LIST_DATA_SIZE +  7)
#define IFC_TET_MSH_YMAX               (LIST_DATA_SIZE +  8)
#define IFC_TET_MSH_ZMIN               (LIST_DATA_SIZE +  9)
#define IFC_TET_MSH_ZMAX               (LIST_DATA_SIZE + 10)
#define IFC_TET_MSH_STAT_IDX           (LIST_DATA_SIZE + 11)
#define IFC_TET_MSH_NF                 (LIST_DATA_SIZE + 12)
#define IFC_TET_MSH_PTR_FACES          (LIST_DATA_SIZE + 13)
#define IFC_TET_MSH_PTR_SIDES          (LIST_DATA_SIZE + 14)

/* OpenFOAM batches */

#define IFC_OF_BATCH_BLOCK_SIZE        (LIST_DATA_SIZE + 1)

#define IFC_OF_BATCH_PTR_NAME          (LIST_DATA_SIZE + 0)

/* Type 5: Interface for fuel performance codes */

#define IFC_FUEP_LIST_BLOCK_SIZE       (LIST_DATA_SIZE + 25)

#define IFC_FUEP_LIST_BLOCK_SIZE       (LIST_DATA_SIZE + 25)
#define IFC_FUEP_PTR_UNI               (LIST_DATA_SIZE +  0)
#define IFC_FUEP_OUT_PTR_LIM           (LIST_DATA_SIZE +  1)
#define IFC_FUEP_OUT_PTR_Z             (LIST_DATA_SIZE +  2)
#define IFC_FUEP_OUT_PTR_R2            (LIST_DATA_SIZE +  3)
#define IFC_FUEP_OUT_PTR_PHI           (LIST_DATA_SIZE +  4)
#define IFC_FUEP_PTR_POWER             (LIST_DATA_SIZE +  5)
#define IFC_FUEP_OUT_PTR_FLIM          (LIST_DATA_SIZE +  6)
#define IFC_FUEP_OUT_PTR_FZ            (LIST_DATA_SIZE +  7)
#define IFC_FUEP_OUT_PTR_FR2           (LIST_DATA_SIZE +  8)
#define IFC_FUEP_OUT_PTR_FPHI          (LIST_DATA_SIZE +  9)
#define IFC_FUEP_PTR_FLUX              (LIST_DATA_SIZE + 10)
#define IFC_FUEP_TYPE                  (LIST_DATA_SIZE + 11)
#define IFC_FUEP_PTR_UNI_LIST          (LIST_DATA_SIZE + 12)
#define IFC_FUEP_N_UNI                 (LIST_DATA_SIZE + 13)
#define IFC_FUEP_OUT_PTR_TB            (LIST_DATA_SIZE + 14)
#define IFC_FUEP_OUT_PTR_FTB           (LIST_DATA_SIZE + 15)
#define IFC_FUEP_N_T                   (LIST_DATA_SIZE + 16)
#define IFC_FUEP_TMIN                  (LIST_DATA_SIZE + 17)
#define IFC_FUEP_TMAX                  (LIST_DATA_SIZE + 18)
#define IFC_FUEP_PTR_T                 (LIST_DATA_SIZE + 19)
#define IFC_FUEP_PTR_POWER_REL         (LIST_DATA_SIZE + 20)
#define IFC_FUEP_PTR_POWER_GRD         (LIST_DATA_SIZE + 21)
#define IFC_FUEP_PTR_POWER_PREV        (LIST_DATA_SIZE + 22)
#define IFC_FUEP_LIM_PTR_T             (LIST_DATA_SIZE + 23)
#define IFC_FUEP_PTR_FINIX             (LIST_DATA_SIZE + 24)

#define IFC_FUEP_T_BLOCK_SIZE          (LIST_DATA_SIZE + 9)

#define IFC_FUEP_T_TMIN                (LIST_DATA_SIZE + 0)
#define IFC_FUEP_T_TMAX                (LIST_DATA_SIZE + 1)
#define IFC_FUEP_T_N_AX                (LIST_DATA_SIZE + 2)
#define IFC_FUEP_T_PTR_AX              (LIST_DATA_SIZE + 3)
#define IFC_FUEP_FINIX_PTR_ROD         (LIST_DATA_SIZE + 4)
#define IFC_FUEP_FINIX_PTR_SCENARIO    (LIST_DATA_SIZE + 5)
#define IFC_FUEP_FINIX_PTR_BC          (LIST_DATA_SIZE + 6)
#define IFC_FUEP_FINIX_PTR_RESULTS     (LIST_DATA_SIZE + 7)
#define IFC_FUEP_FINIX_PTR_OPTIONS     (LIST_DATA_SIZE + 8)

#define IFC_FUEP_AX_BLOCK_SIZE         (LIST_DATA_SIZE +  5)

#define IFC_FUEP_AX_ZMIN               (LIST_DATA_SIZE +  1)
#define IFC_FUEP_AX_ZMAX               (LIST_DATA_SIZE +  2)
#define IFC_FUEP_AX_N_ANG	       (LIST_DATA_SIZE +  3)
#define IFC_FUEP_AX_PTR_ANG            (LIST_DATA_SIZE +  4)

#define IFC_FUEP_ANG_BLOCK_SIZE         (LIST_DATA_SIZE +  10)

#define IFC_FUEP_ANG_AMIN               (LIST_DATA_SIZE +  0)
#define IFC_FUEP_ANG_AMAX               (LIST_DATA_SIZE +  1)
#define IFC_FUEP_ANG_CMIN               (LIST_DATA_SIZE +  2)
#define IFC_FUEP_ANG_CMAX               (LIST_DATA_SIZE +  3)
#define IFC_FUEP_ANG_N_RAD 	        (LIST_DATA_SIZE +  4)
#define IFC_FUEP_ANG_PTR_DF             (LIST_DATA_SIZE +  5)
#define IFC_FUEP_ANG_PTR_TEMP0          (LIST_DATA_SIZE +  6)
#define IFC_FUEP_ANG_PTR_TEMP1          (LIST_DATA_SIZE +  7)
#define IFC_FUEP_ANG_PTR_COLD_R2        (LIST_DATA_SIZE +  8)
#define IFC_FUEP_ANG_PTR_HOT_R2         (LIST_DATA_SIZE +  9)

/* Scoring regions (this is needed to keep list sorted based on geometry */
/* region index, while the output list is sorted based on coordinates).  */

#define IFC_SCORE_LIST_BLOCK_SIZE     (LIST_DATA_SIZE + 3)

#define IFC_SCORE_REG_IDX             (LIST_DATA_SIZE + 0)
#define IFC_SCORE_STAT_IDX            (LIST_DATA_SIZE + 1)
#define IFC_SCORE_PTR_OUT             (LIST_DATA_SIZE + 2)

/* Output for power distribution (other types than tet mesh) */

#define IFC_OUT_LIST_BLOCK_SIZE       (LIST_DATA_SIZE + 5)

#define IFC_OUT_X0                    (LIST_DATA_SIZE + 0)
#define IFC_OUT_Y0                    (LIST_DATA_SIZE + 1)
#define IFC_OUT_R                     (LIST_DATA_SIZE + 2)
#define IFC_OUT_PTR_IFC               (LIST_DATA_SIZE + 3)
#define IFC_OUT_PTR_SCORE             (LIST_DATA_SIZE + 4)

/* Pointers for FUEP parameters (NOTE: N�it� ei voi varate NewItem():ll�, */
/* eik� niihin voi k�ytt�� NextItem() etc. funktioita koska tossa ei ole  */
/* (LIST_DATA_SIZE:a mukana. */

#define FUEP_LIM_SIZE                14

#define FUEP_NZ                       0
#define FUEP_ZMIN                     1
#define FUEP_ZMAX                     2
#define FUEP_NA                       3
#define FUEP_AMIN                     4
#define FUEP_AMAX                     5
#define FUEP_NR                       6
#define FUEP_RMIN                     7
#define FUEP_RMAX                     8
#define FUEP_EMIN                     9
#define FUEP_EMAX                    10
#define FUEP_NT                      11
#define FUEP_TMIN                    12
#define FUEP_TMAX                    13

/*****************************************************************************/

/***** Lattice ***************************************************************/

#define LATTICE_TYPES    13

#define LAT_TYPE_S       1
#define LAT_TYPE_HX      2
#define LAT_TYPE_HY      3
#define LAT_TYPE_CLU     4
#define LAT_TYPE_RND     5
#define LAT_TYPE_INFS    6
#define LAT_TYPE_INFHY   7
#define LAT_TYPE_INFHX   8
#define LAT_TYPE_ZSTACK  9
#define LAT_TYPE_CUBOID  11
#define LAT_TYPE_XPRISM  12
#define LAT_TYPE_YPRISM  13

#define LAT_BLOCK_SIZE                (LIST_DATA_SIZE + PARAM_N_COMMON + 19)

#define LAT_PTR_NAME                  (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define LAT_PTR_UNI                   (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define LAT_OPTIONS                   (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define LAT_TYPE                      (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define LAT_ORIG_X0                   (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define LAT_ORIG_Y0                   (LIST_DATA_SIZE + PARAM_N_COMMON +  5)
#define LAT_ORIG_Z0                   (LIST_DATA_SIZE + PARAM_N_COMMON +  6)
#define LAT_PITCH                     (LIST_DATA_SIZE + PARAM_N_COMMON +  7)
#define LAT_PITCHX                    (LIST_DATA_SIZE + PARAM_N_COMMON +  8)
#define LAT_PITCHY                    (LIST_DATA_SIZE + PARAM_N_COMMON +  9)
#define LAT_PITCHZ                    (LIST_DATA_SIZE + PARAM_N_COMMON + 10)
#define LAT_NX                        (LIST_DATA_SIZE + PARAM_N_COMMON + 11)
#define LAT_NY                        (LIST_DATA_SIZE + PARAM_N_COMMON + 12)
#define LAT_NZ                        (LIST_DATA_SIZE + PARAM_N_COMMON + 13)
#define LAT_NTOT                      (LIST_DATA_SIZE + PARAM_N_COMMON + 14)
#define LAT_N_RINGS                   (LIST_DATA_SIZE + PARAM_N_COMMON + 15)
#define LAT_PTR_FILL                  (LIST_DATA_SIZE + PARAM_N_COMMON + 16)
#define LAT_PTR_Z                     (LIST_DATA_SIZE + PARAM_N_COMMON + 17)
#define LAT_COL_CELL_IDX              (LIST_DATA_SIZE + PARAM_N_COMMON + 18)

#define RING_BLOCK_SIZE                (LIST_DATA_SIZE + 5)

#define RING_N_SEC                     (LIST_DATA_SIZE + 0)
#define RING_RAD                       (LIST_DATA_SIZE + 1)
#define RING_RLIM                      (LIST_DATA_SIZE + 2)
#define RING_TILT                      (LIST_DATA_SIZE + 3)
#define RING_PTR_FILL                  (LIST_DATA_SIZE + 4)

/*****************************************************************************/

/***** Depletion history *****************************************************/

/* Depletion step types */

#define DEP_STEP_BU_STEP   1
#define DEP_STEP_BU_TOT    2
#define DEP_STEP_DAY_STEP  3
#define DEP_STEP_DAY_TOT   4
#define DEP_STEP_DEC_STEP  5
#define DEP_STEP_DEC_TOT   6
#define DEP_STEP_ACT_STEP  7
#define DEP_STEP_ACT_TOT   8

/* predictor and corrector types available for burnup calculations (AIs) */

#define PRED_TYPE_CONSTANT   10
#define PRED_TYPE_LINEAR     11
#define CORR_TYPE_NONE       19
#define CORR_TYPE_CONSTANT   20
#define CORR_TYPE_LINEAR     21
#define CORR_TYPE_QUADRATIC  22

/* NOTE: predictor =/= the last corrector iteration, keep as NO on pred. steps*/

/* the types */

#define CI_TYPE_NONE                   -1
#define CI_TYPE_INNER                   0
#define CI_TYPE_OUTER                   1

/* inner iteration means that pEOS flux is iterated as in Dufek's method,
   i.e., using backwards constant extrapolation, after which the final round
   uses the selected corrector scheme. Outer iteration means that each
   iteration uses the selected corrector scheme. Default is OUTER with MAXI=1
   which results in the old (no iteration) behavior. (AIs) */


#define DEP_HIS_BLOCK_SIZE            (LIST_DATA_SIZE + PARAM_N_COMMON + 10)

#define DEP_HIS_STEP_TYPE             (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define DEP_HIS_N_STEPS               (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define DEP_HIS_PTR_STEPS             (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define DEP_HIS_PTR_NORM              (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define DEP_HIS_PRED_TYPE             (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define DEP_HIS_PRED_NSS              (LIST_DATA_SIZE + PARAM_N_COMMON +  5)
#define DEP_HIS_CORR_TYPE             (LIST_DATA_SIZE + PARAM_N_COMMON +  6)
#define DEP_HIS_CORR_NSS              (LIST_DATA_SIZE + PARAM_N_COMMON +  7)
#define DEP_HIS_PTR_REPROC            (LIST_DATA_SIZE + PARAM_N_COMMON +  8)
#define DEP_HIS_PTR_BRANCH            (LIST_DATA_SIZE + PARAM_N_COMMON +  9)

/*****************************************************************************/

/***** Source array **********************************************************/

#define SRC_BLOCK_SIZE                (LIST_DATA_SIZE + PARAM_N_COMMON + 37)

#define SRC_PTR_NAME                  (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define SRC_TYPE                      (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define SRC_WGT                       (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define SRC_E                         (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define SRC_PTR_XSDATA                (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define SRC_PTR_REA                   (LIST_DATA_SIZE + PARAM_N_COMMON +  5)
#define SRC_PTR_MAT                   (LIST_DATA_SIZE + PARAM_N_COMMON +  6)
#define SRC_PTR_CELL                  (LIST_DATA_SIZE + PARAM_N_COMMON +  7)
#define SRC_PTR_UNIV                  (LIST_DATA_SIZE + PARAM_N_COMMON +  8)
#define SRC_PTR_EBINS                 (LIST_DATA_SIZE + PARAM_N_COMMON +  9)
#define SRC_PTR_SURF                  (LIST_DATA_SIZE + PARAM_N_COMMON + 10)
#define SRC_SURF_SIDE                 (LIST_DATA_SIZE + PARAM_N_COMMON + 11)
#define SRC_X0                        (LIST_DATA_SIZE + PARAM_N_COMMON + 12)
#define SRC_Y0                        (LIST_DATA_SIZE + PARAM_N_COMMON + 13)
#define SRC_Z0                        (LIST_DATA_SIZE + PARAM_N_COMMON + 14)
#define SRC_U0                        (LIST_DATA_SIZE + PARAM_N_COMMON + 15)
#define SRC_V0                        (LIST_DATA_SIZE + PARAM_N_COMMON + 16)
#define SRC_W0                        (LIST_DATA_SIZE + PARAM_N_COMMON + 17)
#define SRC_XMIN                      (LIST_DATA_SIZE + PARAM_N_COMMON + 18)
#define SRC_XMAX                      (LIST_DATA_SIZE + PARAM_N_COMMON + 19)
#define SRC_YMIN                      (LIST_DATA_SIZE + PARAM_N_COMMON + 20)
#define SRC_YMAX                      (LIST_DATA_SIZE + PARAM_N_COMMON + 21)
#define SRC_ZMIN                      (LIST_DATA_SIZE + PARAM_N_COMMON + 22)
#define SRC_ZMAX                      (LIST_DATA_SIZE + PARAM_N_COMMON + 23)
#define SRC_TMIN                      (LIST_DATA_SIZE + PARAM_N_COMMON + 24)
#define SRC_TMAX                      (LIST_DATA_SIZE + PARAM_N_COMMON + 25)
#define SRC_PTR_RAD_SRC_MAT           (LIST_DATA_SIZE + PARAM_N_COMMON + 26)
#define SRC_RAD_SRC_MODE              (LIST_DATA_SIZE + PARAM_N_COMMON + 27)
#define SRC_READ_FILE_TYPE            (LIST_DATA_SIZE + PARAM_N_COMMON + 28)
#define SRC_READ_PTR_FILE             (LIST_DATA_SIZE + PARAM_N_COMMON + 29)
#define SRC_READ_PTR_BUF              (LIST_DATA_SIZE + PARAM_N_COMMON + 30)
#define SRC_READ_BUF_SZ               (LIST_DATA_SIZE + PARAM_N_COMMON + 31)
#define SRC_READ_BUF_IDX              (LIST_DATA_SIZE + PARAM_N_COMMON + 32)
#define SRC_READ_FILE_POS             (LIST_DATA_SIZE + PARAM_N_COMMON + 33)
#define SRC_PTR_USR                   (LIST_DATA_SIZE + PARAM_N_COMMON + 34)
#define SRC_PTR_PLASMA_SRC            (LIST_DATA_SIZE + PARAM_N_COMMON + 35)
#define SRC_READ_BINARY               (LIST_DATA_SIZE + PARAM_N_COMMON + 36)

/* Source energy bin */

#define SRC_EBIN_BLOCK_SIZE           (LIST_DATA_SIZE + 2)

#define SRC_EBIN_EMAX                 (LIST_DATA_SIZE + 0)
#define SRC_EBIN_WGT                  (LIST_DATA_SIZE + 1)

/* Source buffer entry (toi ei oo linkitetty lista) */

#define SRC_BUF_BLOCK_SIZE            9

#define SRC_BUF_X                     0
#define SRC_BUF_Y                     1
#define SRC_BUF_Z                     2
#define SRC_BUF_U                     3
#define SRC_BUF_V                     4
#define SRC_BUF_W                     5
#define SRC_BUF_E                     6
#define SRC_BUF_WGT                   7
#define SRC_BUF_T                     8

/* User-defined routine */

#define SRC_USR_BLOCK_SIZE            (LIST_DATA_SIZE + 2)
#define SRC_USR_NP                    (LIST_DATA_SIZE + 0)
#define SRC_USR_PTR_PARAMS            (LIST_DATA_SIZE + 1)

/* Decay source */

#define SRC_DECCAY_BLOCK_SIZE          (LIST_DATA_SIZE + 5)

#define SRC_DECCAY_PTR_NUCLIDE         (LIST_DATA_SIZE + 0)
#define SRC_DECCAY_PTR_SPEC            (LIST_DATA_SIZE + 1)
#define SRC_DECCAY_I                   (LIST_DATA_SIZE + 2)
#define SRC_DECCAY_CUM_P               (LIST_DATA_SIZE + 3)
#define SRC_DECCAY_WGT                 (LIST_DATA_SIZE + 4)

/* Fusion plasma source (only thermal included) */

#define SRC_PLASMA_BLOCK_SIZE          (LIST_DATA_SIZE + 13)

#define SRC_PLASMA_NTOT                (LIST_DATA_SIZE + 0)
#define SRC_PLASMA_X0                  (LIST_DATA_SIZE + 1)
#define SRC_PLASMA_Y0                  (LIST_DATA_SIZE + 2)
#define SRC_PLASMA_CO                  (LIST_DATA_SIZE + 3)
#define SRC_PLASMA_XYB                 (LIST_DATA_SIZE + 4)
#define SRC_PLASMA_PTR_XB              (LIST_DATA_SIZE + 5)
#define SRC_PLASMA_PTR_YB              (LIST_DATA_SIZE + 6)
#define SRC_PLASMA_N_ISO               (LIST_DATA_SIZE + 7)
#define SRC_PLASMA_PTR_PHI             (LIST_DATA_SIZE + 8)
#define SRC_PLASMA_PTR_P_ISO           (LIST_DATA_SIZE + 9)
#define SRC_PLASMA_REACTIONS           (LIST_DATA_SIZE + 10)
#define SRC_PLASMA_PTR_REACT           (LIST_DATA_SIZE + 11)
#define SRC_PLASMA_DATA_SIZE           (LIST_DATA_SIZE + 12)

#define SRC_PLASMA_REA_BLOCK_SIZE      (LIST_DATA_SIZE + 6)

#define SRC_PLASMA_REA_PROB            (LIST_DATA_SIZE + 0)
#define SRC_PLASMA_REA_T               (LIST_DATA_SIZE + 1)
#define SRC_PLASMA_REA_N_E             (LIST_DATA_SIZE + 2)
#define SRC_PLASMA_REA_PTR_E           (LIST_DATA_SIZE + 3)
#define SRC_PLASMA_REA_PTR_P_E         (LIST_DATA_SIZE + 4)
#define SRC_PLASMA_REA_PTR_GEOM        (LIST_DATA_SIZE + 5)

#define SRC_PLASMA_REA_GEOM_BLOCK_SIZE (LIST_DATA_SIZE + 14)

/* RHO-THETA coordinate system */

#define SRC_PLASMA_REA_N_RHO           (LIST_DATA_SIZE + 0)
#define SRC_PLASMA_REA_PTR_RHO         (LIST_DATA_SIZE + 1)
#define SRC_PLASMA_REA_PTR_P_RHO       (LIST_DATA_SIZE + 2)
#define SRC_PLASMA_REA_XYB             (LIST_DATA_SIZE + 3)
#define SRC_PLASMA_REA_PTR_XB          (LIST_DATA_SIZE + 4)
#define SRC_PLASMA_REA_PTR_YB          (LIST_DATA_SIZE + 5)
#define SRC_PLASMA_REA_PTR_P_XY        (LIST_DATA_SIZE + 6)

/* RZ coordinate system */

#define SRC_PLASMA_REA_N_RZ            (LIST_DATA_SIZE + 7)
#define SRC_PLASMA_REA_PTR_R           (LIST_DATA_SIZE + 8)
#define SRC_PLASMA_REA_PTR_Z           (LIST_DATA_SIZE + 9)
#define SRC_PLASMA_REA_PTR_P_RZ        (LIST_DATA_SIZE + 10)
#define SRC_PLASMA_REA_N_PHI           (LIST_DATA_SIZE + 11)
#define SRC_PLASMA_REA_PTR_PHI         (LIST_DATA_SIZE + 12)
#define SRC_PLASMA_REA_PTR_P_PHI       (LIST_DATA_SIZE + 13)

/*****************************************************************************/

/***** Detector array ********************************************************/

#define DET_BLOCK_SIZE                (LIST_DATA_SIZE + PARAM_N_COMMON + 37)

#define DET_PTR_NAME                  (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define DET_TYPE                      (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define DET_PTR_EGRID                 (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define DET_PTR_RBINS                 (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define DET_PTR_UBINS                 (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define DET_PTR_LBINS                 (LIST_DATA_SIZE + PARAM_N_COMMON +  5)
#define DET_PTR_MBINS                 (LIST_DATA_SIZE + PARAM_N_COMMON +  6)
#define DET_PTR_CBINS                 (LIST_DATA_SIZE + PARAM_N_COMMON +  7)
#define DET_PTR_SBINS                 (LIST_DATA_SIZE + PARAM_N_COMMON +  8)
#define DET_PTR_TME                   (LIST_DATA_SIZE + PARAM_N_COMMON + 10)
#define DET_N_EBINS                   (LIST_DATA_SIZE + PARAM_N_COMMON + 11)
#define DET_N_UBINS                   (LIST_DATA_SIZE + PARAM_N_COMMON + 12)
#define DET_N_CBINS                   (LIST_DATA_SIZE + PARAM_N_COMMON + 13)
#define DET_N_MBINS                   (LIST_DATA_SIZE + PARAM_N_COMMON + 14)
#define DET_N_LBINS                   (LIST_DATA_SIZE + PARAM_N_COMMON + 15)
#define DET_N_RBINS                   (LIST_DATA_SIZE + PARAM_N_COMMON + 16)
#define DET_N_TBINS                   (LIST_DATA_SIZE + PARAM_N_COMMON + 18)
#define DET_N_TOT_BINS                (LIST_DATA_SIZE + PARAM_N_COMMON + 19)
#define DET_PTR_STAT                  (LIST_DATA_SIZE + PARAM_N_COMMON + 20)
#define DET_VOL                       (LIST_DATA_SIZE + PARAM_N_COMMON + 21)
#define DET_PTR_MUL                   (LIST_DATA_SIZE + PARAM_N_COMMON + 22)
#define DET_PTR_ADJOINT               (LIST_DATA_SIZE + PARAM_N_COMMON + 23)
#define DET_PTR_MESH                  (LIST_DATA_SIZE + PARAM_N_COMMON + 24)
#define DET_PARTICLE                  (LIST_DATA_SIZE + PARAM_N_COMMON + 25)
#define DET_WRITE_PTR_FILE            (LIST_DATA_SIZE + PARAM_N_COMMON + 26)
#define DET_WRITE_PTR_BUF             (LIST_DATA_SIZE + PARAM_N_COMMON + 27)
#define DET_WRITE_BUF_SZ              (LIST_DATA_SIZE + PARAM_N_COMMON + 28)
#define DET_WRITE_BUF_IDX             (LIST_DATA_SIZE + PARAM_N_COMMON + 29)
#define DET_WRITE_PROB                (LIST_DATA_SIZE + PARAM_N_COMMON + 30)
#define DET_WRITE_HIS                 (LIST_DATA_SIZE + PARAM_N_COMMON + 31)
#define DET_SCORE_FISS_REG_ONLY       (LIST_DATA_SIZE + PARAM_N_COMMON + 32)
#define DET_DIRVEC_U                  (LIST_DATA_SIZE + PARAM_N_COMMON + 33)
#define DET_DIRVEC_V                  (LIST_DATA_SIZE + PARAM_N_COMMON + 34)
#define DET_DIRVEC_W                  (LIST_DATA_SIZE + PARAM_N_COMMON + 35)
#define DET_WRITE_BINARY              (LIST_DATA_SIZE + PARAM_N_COMMON + 36)

/* Detector reaction bin */

#define DET_RBIN_BLOCK_SIZE           (LIST_DATA_SIZE + 10)

#define DET_RBIN_PTR_MAT              (LIST_DATA_SIZE + 0)
#define DET_RBIN_MT                   (LIST_DATA_SIZE + 1)
#define DET_RBIN_PTR_REA              (LIST_DATA_SIZE + 2)
#define DET_RBIN_VOID_MODE            (LIST_DATA_SIZE + 3)
#define DET_RBIN_PTR_FUN              (LIST_DATA_SIZE + 4)
#define DET_RBIN_ATTN_NP              (LIST_DATA_SIZE + 5)
#define DET_RBIN_ATTN_PTR_E           (LIST_DATA_SIZE + 6)
#define DET_RBIN_ATTN_PTR_F           (LIST_DATA_SIZE + 7)
#define DET_RBIN_CONVERT              (LIST_DATA_SIZE + 8)
#define DET_RBIN_PTR_PULSE_DATA       (LIST_DATA_SIZE + 9)

/* Pulse data for pulse-height detector (PRIVA-array, no LIST_DATA_SIZE) */

#define DET_PULSE_BLOCK_SIZE          3

#define DET_PULSE_PHOTON_IDX          0
#define DET_PULSE_EDEP                1
#define DET_PULSE_WGT                 2

/* Detector universe bin */

#define DET_UBIN_BLOCK_SIZE           (LIST_DATA_SIZE + 1)

#define DET_UBIN_PTR_UNI              (LIST_DATA_SIZE + 0)

/* Detector lattice bin */

#define DET_LBIN_BLOCK_SIZE           (LIST_DATA_SIZE + 1)

#define DET_LBIN_PTR_LAT              (LIST_DATA_SIZE + 0)

/* Detector material bin */

#define DET_MBIN_BLOCK_SIZE           (LIST_DATA_SIZE + 1)

#define DET_MBIN_PTR_MAT              (LIST_DATA_SIZE + 0)

/* Detector cell bin */

#define DET_CBIN_BLOCK_SIZE           (LIST_DATA_SIZE + 5)

#define DET_CBIN_PTR_CELL             (LIST_DATA_SIZE + 0)
#define DET_CBIN_SUPER_CELL           (LIST_DATA_SIZE + 1)
#define DET_CBIN_UMSH_PTR_UMSH        (LIST_DATA_SIZE + 2)
#define DET_CBIN_UMSH_PTR_CELLS       (LIST_DATA_SIZE + 3)
#define DET_CBIN_UMSH_PTR_BINS        (LIST_DATA_SIZE + 4)

/* Bins for super-imposed detector */

#define DET_SBIN_BLOCK_SIZE           (LIST_DATA_SIZE + 3)

#define DET_SBIN_TYPE                 (LIST_DATA_SIZE + 0)
#define DET_SBIN_PTR_SURF             (LIST_DATA_SIZE + 1)
#define DET_SBIN_SURF_NORM            (LIST_DATA_SIZE + 2)

/* Bins list for material detectors */

#define DETBIN_BLOCK_SIZE             (LIST_DATA_SIZE + 2)

#define DETBIN_PTR_DET                (LIST_DATA_SIZE + 0)
#define DETBIN_BIN                    (LIST_DATA_SIZE + 1)

/*****************************************************************************/

/***** User-defined energy grid **********************************************/

#define ENE_BLOCK_SIZE                (LIST_DATA_SIZE + PARAM_N_COMMON +  7)

#define ENE_PTR_NAME                  (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define ENE_TYPE                      (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define ENE_NB                        (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define ENE_EMIN                      (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define ENE_EMAX                      (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define ENE_PTR_GRID                  (LIST_DATA_SIZE + PARAM_N_COMMON +  5)
#define ENE_PTR_PREDEF                (LIST_DATA_SIZE + PARAM_N_COMMON +  6)

/*****************************************************************************/

/***** User-defined time binning *********************************************/

#define TME_BLOCK_SIZE                (LIST_DATA_SIZE + PARAM_N_COMMON +  6)

#define TME_PTR_NAME                  (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define TME_TYPE                      (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define TME_NB                        (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define TME_TMIN                      (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define TME_TMAX                      (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define TME_PTR_BINS                  (LIST_DATA_SIZE + PARAM_N_COMMON +  5)

/*****************************************************************************/

/***** User-defined response function ****************************************/

#define FUN_BLOCK_SIZE                (LIST_DATA_SIZE + PARAM_N_COMMON +  6)

#define FUN_PTR_NAME                  (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define FUN_TYPE                      (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define FUN_INT                       (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define FUN_NE                        (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define FUN_PTR_E                     (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define FUN_PTR_F                     (LIST_DATA_SIZE + PARAM_N_COMMON +  5)

/*****************************************************************************/

/***** Weight window array ***************************************************/

#define WWD_BLOCK_SIZE                (LIST_DATA_SIZE + PARAM_N_COMMON + 10)

#define WWD_NORM_FACT                 (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define WWD_NORM_X                    (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define WWD_NORM_Y                    (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define WWD_NORM_Z                    (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define WWD_MESH_MIN                  (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define WWD_MESH_MAX                  (LIST_DATA_SIZE + PARAM_N_COMMON +  5)
#define WWD_PTR_MESH                  (LIST_DATA_SIZE + PARAM_N_COMMON +  6)
#define WWD_PTR_MESH_DATA             (LIST_DATA_SIZE + PARAM_N_COMMON +  7)
#define WWD_POW                       (LIST_DATA_SIZE + PARAM_N_COMMON +  8)
#define WWD_PTR_ERG                   (LIST_DATA_SIZE + PARAM_N_COMMON +  9)

#define WWD_MESH_BLOCK_SIZE           (LIST_DATA_SIZE + 7)

#define WWD_MESH_IMP                  (LIST_DATA_SIZE + 0)
#define WWD_MESH_PTR_CURR             (LIST_DATA_SIZE + 1)
#define WWD_MESH_RES_OUT_CURR         (LIST_DATA_SIZE + 2)
#define WWD_MESH_RES_IN_CURR          (LIST_DATA_SIZE + 3)
#define WWD_MESH_RES_SRC_RATE         (LIST_DATA_SIZE + 4)

#define WWD_ERG_BLOCK_SIZE            (LIST_DATA_SIZE + 3)

#define WWD_ERG_EMIN                  (LIST_DATA_SIZE + 0)
#define WWD_ERG_EMAX                  (LIST_DATA_SIZE + 1)
#define WWD_ERG_IMP                   (LIST_DATA_SIZE + 2)

#define WWD_MESH_CURR_BLOCK_SIZE      (LIST_DATA_SIZE + 1)

#define WWD_MESH_CURR_PTR_NEIGHBOUR   (LIST_DATA_SIZE + 0)

/*****************************************************************************/

/***** Geometry levels *******************************************************/

/* Data stored in common array */

#define LVL_BLOCK_SIZE        (LIST_DATA_SIZE + 5)

#define LVL_NUMBER            (LIST_DATA_SIZE + 0)
#define LVL_PTR_PRIVATE_DATA  (LIST_DATA_SIZE + 1)
#define LVL_MAX_REGIONS       (LIST_DATA_SIZE + 2)
#define LVL_CUM_MAX_REGIONS   (LIST_DATA_SIZE + 3)
#define LVL_ZONE_IDX_MULT     (LIST_DATA_SIZE + 4)

/* Data stored in private array */

#define LVL_PRIV_BLOCK_SIZE     26

#define LVL_PRIV_TYPE            0
#define LVL_PRIV_PTR_NEST_REG    1
#define LVL_PRIV_PTR_LAT         2
#define LVL_PRIV_PTR_CELL        3
#define LVL_PRIV_PTR_PBED        4
#define LVL_PRIV_PTR_PEBBLE      5
#define LVL_PRIV_PTR_UMSH        6
#define LVL_PRIV_PTR_STL         7
#define LVL_PRIV_X               8
#define LVL_PRIV_Y               9
#define LVL_PRIV_Z              10
#define LVL_PRIV_U              11
#define LVL_PRIV_V              12
#define LVL_PRIV_W              13
#define LVL_PRIV_LAST           14
#define LVL_PRIV_PTR_MAT        15
#define LVL_PRIV_LAT_SURF_TYPE  16
#define LVL_PRIV_LAT_SURF_NP    17
#define LVL_PRIV_LAT_SURF_C0    18
#define LVL_PRIV_LAT_SURF_C1    19
#define LVL_PRIV_LAT_SURF_C2    20
#define LVL_PRIV_LAT_SURF_C3    21
#define LVL_PRIV_LAT_SURF_C4    22
#define LVL_PRIV_LAT_SURF_C5    23
#define LVL_PRIV_PTR_UNIV       24
#define LVL_PRIV_ZONE_IDX       25

/*****************************************************************************/

/***** Score block ***********************************************************/

#define SCORE_BLOCK_SIZE  (LIST_DATA_SIZE + 7)

#define SCORE_PTR_NAME    (LIST_DATA_SIZE + 0)
#define SCORE_DIM         (LIST_DATA_SIZE + 1)
#define SCORE_PTR_NMAX    (LIST_DATA_SIZE + 2)
#define SCORE_PTR_DATA    (LIST_DATA_SIZE + 3)
#define SCORE_PTR_BUF     (LIST_DATA_SIZE + 4)
#define SCORE_STAT_SIZE   (LIST_DATA_SIZE + 5)
#define SCORE_PTR_HIS     (LIST_DATA_SIZE + 6)

/*****************************************************************************/

/***** Universe **************************************************************/

/* Universe types */

#define UNIVERSE_TYPE_CELL     1
#define UNIVERSE_TYPE_NEST     2
#define UNIVERSE_TYPE_LATTICE  3
#define UNIVERSE_TYPE_PBED     4
#define UNIVERSE_TYPE_SUPER    5
#define UNIVERSE_TYPE_UMSH     6
#define UNIVERSE_TYPE_STL      7

/* Data block */

#define UNIVERSE_BLOCK_SIZE      (LIST_DATA_SIZE + 32)

#define UNIVERSE_PTR_NAME        (LIST_DATA_SIZE +  0)
#define UNIVERSE_OPTIONS         (LIST_DATA_SIZE +  1)
#define UNIVERSE_TYPE            (LIST_DATA_SIZE +  2)
#define UNIVERSE_PTR_CELL_LIST   (LIST_DATA_SIZE +  3)
#define UNIVERSE_PTR_NEST        (LIST_DATA_SIZE +  4)
#define UNIVERSE_PTR_LAT         (LIST_DATA_SIZE +  5)
#define UNIVERSE_PTR_PBED        (LIST_DATA_SIZE +  6)
#define UNIVERSE_PTR_UMSH        (LIST_DATA_SIZE +  7)
#define UNIVERSE_PTR_STL         (LIST_DATA_SIZE +  8)
#define UNIVERSE_PTR_TRANS       (LIST_DATA_SIZE +  9)
#define UNIVERSE_PTR_SYM         (LIST_DATA_SIZE + 10)
#define UNIVERSE_MINX            (LIST_DATA_SIZE + 11)
#define UNIVERSE_MAXX            (LIST_DATA_SIZE + 12)
#define UNIVERSE_MINY            (LIST_DATA_SIZE + 13)
#define UNIVERSE_MAXY            (LIST_DATA_SIZE + 14)
#define UNIVERSE_MINZ            (LIST_DATA_SIZE + 15)
#define UNIVERSE_MAXZ            (LIST_DATA_SIZE + 16)
#define UNIVERSE_DIM             (LIST_DATA_SIZE + 17)
#define UNIVERSE_COL_COUNT       (LIST_DATA_SIZE + 18)
#define UNIVERSE_PTR_GCU         (LIST_DATA_SIZE + 19)
#define UNIVERSE_LEVEL           (LIST_DATA_SIZE + 20)
#define UNIVERSE_PTR_PRIVA_X     (LIST_DATA_SIZE + 21)
#define UNIVERSE_PTR_PRIVA_Y     (LIST_DATA_SIZE + 22)
#define UNIVERSE_PTR_PRIVA_Z     (LIST_DATA_SIZE + 23)
#define UNIVERSE_PTR_PRIVA_T     (LIST_DATA_SIZE + 24)
#define UNIVERSE_PTR_PREV_REG    (LIST_DATA_SIZE + 25)
#define UNIVERSE_FMTX_IDX        (LIST_DATA_SIZE + 26)
#define UNIVERSE_PTR_IFC_FUEP    (LIST_DATA_SIZE + 27)
#define UNIVERSE_GCU_IDX         (LIST_DATA_SIZE + 28)
#define UNIVERSE_WARN_MULTI_LVL  (LIST_DATA_SIZE + 29)
#define UNIVERSE_PTR_CELL_MESH   (LIST_DATA_SIZE + 30)
#define UNIVERSE_PTR_NEXT_CELL   (LIST_DATA_SIZE + 31)

/*****************************************************************************/

/***** Particle (neutron / photon) *******************************************/

/* Particle types (indeksit pit�� menn� noin ett� sorttaus tyypin mukaan */
/* menee oikein) */

#define PARTICLE_TYPE_DUMMY      0
#define PARTICLE_TYPE_GAMMA      1
#define PARTICLE_TYPE_NEUTRON    2
#define PARTICLE_TYPE_PRECURSOR  3

/* Data block */

#define PARTICLE_BLOCK_SIZE       (LIST_DATA_SIZE + 38)

#define PARTICLE_HISTORY_IDX      (LIST_DATA_SIZE +  0)
#define PARTICLE_RNG_IDX          (LIST_DATA_SIZE +  1)
#define PARTICLE_TYPE             (LIST_DATA_SIZE +  2)
#define PARTICLE_X                (LIST_DATA_SIZE +  3)
#define PARTICLE_Y                (LIST_DATA_SIZE +  4)
#define PARTICLE_Z                (LIST_DATA_SIZE +  5)
#define PARTICLE_U                (LIST_DATA_SIZE +  6)
#define PARTICLE_V                (LIST_DATA_SIZE +  7)
#define PARTICLE_W                (LIST_DATA_SIZE +  8)
#define PARTICLE_E                (LIST_DATA_SIZE +  9)
#define PARTICLE_WGT              (LIST_DATA_SIZE + 10)
#define PARTICLE_T0               (LIST_DATA_SIZE + 11)
#define PARTICLE_T                (LIST_DATA_SIZE + 12)
#define PARTICLE_TD               (LIST_DATA_SIZE + 13)
#define PARTICLE_TT               (LIST_DATA_SIZE + 14)
#define PARTICLE_FMTX_IDX         (LIST_DATA_SIZE + 15)
#define PARTICLE_PREV_IMP         (LIST_DATA_SIZE + 16)
#define PARTICLE_PREV_IDX         (LIST_DATA_SIZE + 17)
#define PARTICLE_PTR_GCU          (LIST_DATA_SIZE + 18)

#ifdef OLD_HIST

#define PARTICLE_PTR_HIST         (LIST_DATA_SIZE + 19)

#endif

#define PARTICLE_PTR_MAT          (LIST_DATA_SIZE + 20)
#define PARTICLE_MPI_ID           (LIST_DATA_SIZE + 21)
#define PARTICLE_GEN_IDX          (LIST_DATA_SIZE + 22)
#define PARTICLE_DN_GROUP         (LIST_DATA_SIZE + 23)
#define PARTICLE_DN_LAMBDA        (LIST_DATA_SIZE + 24)

#ifdef OLD_IFP

#define PARTICLE_PTR_FISS_PROG    (LIST_DATA_SIZE + 25)

#endif

#define PARTICLE_ICM_PTR_ICM      (LIST_DATA_SIZE + 26)
#define PARTICLE_ICM_IDX          (LIST_DATA_SIZE + 27)
#define PARTICLE_ICM_MUA          (LIST_DATA_SIZE + 28)
#define PARTICLE_ICM_MUS          (LIST_DATA_SIZE + 29)
#define PARTICLE_ICM_G            (LIST_DATA_SIZE + 30)
#define PARTICLE_ICM_WGT          (LIST_DATA_SIZE + 31)
#define PARTICLE_ALB_PTR_GCU      (LIST_DATA_SIZE + 32)
#define PARTICLE_ALB_SURF_IDX     (LIST_DATA_SIZE + 33)
#define PARTICLE_ALB_G            (LIST_DATA_SIZE + 34)
#define PARTICLE_PTR_EVENTS       (LIST_DATA_SIZE + 35)
#define PARTICLE_COL_IDX          (LIST_DATA_SIZE + 36)
#define PARTICLE_MULTIPLICITY     (LIST_DATA_SIZE + 37)

/* History data */

#define HIST_BLOCK_SIZE           (LIST_DATA_SIZE + 14)

#define HIST_X                    (LIST_DATA_SIZE +  0)
#define HIST_Y                    (LIST_DATA_SIZE +  1)
#define HIST_Z                    (LIST_DATA_SIZE +  2)
#define HIST_U                    (LIST_DATA_SIZE +  3)
#define HIST_V                    (LIST_DATA_SIZE +  4)
#define HIST_W                    (LIST_DATA_SIZE +  5)
#define HIST_E                    (LIST_DATA_SIZE +  6)
#define HIST_T                    (LIST_DATA_SIZE +  7)
#define HIST_WGT                  (LIST_DATA_SIZE +  8)
#define HIST_FLX                  (LIST_DATA_SIZE +  9)
#define HIST_PTR_REA              (LIST_DATA_SIZE + 10)
#define HIST_PTR_MAT              (LIST_DATA_SIZE + 11)
#define HIST_TRK                  (LIST_DATA_SIZE + 12)
#define HIST_IDX                  (LIST_DATA_SIZE + 13)

/* Fission progenies */

#ifdef OLD_IFP

#define FISS_PROG_BLOCK_SIZE      (LIST_DATA_SIZE + 3)

#define FISS_PROG_DN_GROUP        (LIST_DATA_SIZE + 0)
#define FISS_PROG_LIFETIME        (LIST_DATA_SIZE + 1)
#define FISS_PROG_LAMBDA          (LIST_DATA_SIZE + 2)

#endif

/*****************************************************************************/

/***** Event array ***********************************************************/

#define EVENT_BLOCK_SIZE          (LIFO_LIST_DATA_SIZE + 10)

/* Common to all */

#define EVENT_HIS_COUNT           (LIFO_LIST_DATA_SIZE + 0)
#define EVENT_TYPE                (LIFO_LIST_DATA_SIZE + 1)

/* Fissions */

#define EVENT_DN_GROUP            (LIFO_LIST_DATA_SIZE + 2)
#define EVENT_LIFETIME            (LIFO_LIST_DATA_SIZE + 3)
#define EVENT_LAMBDA              (LIFO_LIST_DATA_SIZE + 4)

/* Collision points */

#define EVENT_X                   (LIFO_LIST_DATA_SIZE + 2)
#define EVENT_Y                   (LIFO_LIST_DATA_SIZE + 3)
#define EVENT_Z                   (LIFO_LIST_DATA_SIZE + 4)

/* Energy and time */

#define EVENT_E                   (LIFO_LIST_DATA_SIZE + 5)
#define EVENT_T                   (LIFO_LIST_DATA_SIZE + 6)

/* Weight and flux */

#define EVENT_WGT                 (LIFO_LIST_DATA_SIZE + 7)
#define EVENT_FLX                 (LIFO_LIST_DATA_SIZE + 8)

/* Material */

#define EVENT_PTR_MAT             (LIFO_LIST_DATA_SIZE + 9)

/*****************************************************************************/

/***** Mesh plot *************************************************************/

#define MPL_BLOCK_SIZE          (LIST_DATA_SIZE + PARAM_N_COMMON + 25)

#define MPL_NDIST               (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define MPL_DIV                 (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define MPL_AX                  (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define MPL_NX                  (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define MPL_NY                  (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define MPL_SYM                 (LIST_DATA_SIZE + PARAM_N_COMMON +  5)
#define MPL_XMIN                (LIST_DATA_SIZE + PARAM_N_COMMON +  6)
#define MPL_XMAX                (LIST_DATA_SIZE + PARAM_N_COMMON +  7)
#define MPL_YMIN                (LIST_DATA_SIZE + PARAM_N_COMMON +  8)
#define MPL_YMAX                (LIST_DATA_SIZE + PARAM_N_COMMON +  9)
#define MPL_ZMIN                (LIST_DATA_SIZE + PARAM_N_COMMON + 10)
#define MPL_ZMAX                (LIST_DATA_SIZE + PARAM_N_COMMON + 11)
#define MPL_MIN1                (LIST_DATA_SIZE + PARAM_N_COMMON + 12)
#define MPL_MAX1                (LIST_DATA_SIZE + PARAM_N_COMMON + 13)
#define MPL_MIN2                (LIST_DATA_SIZE + PARAM_N_COMMON + 14)
#define MPL_MAX2                (LIST_DATA_SIZE + PARAM_N_COMMON + 15)
#define MPL_PTR_VAL1            (LIST_DATA_SIZE + PARAM_N_COMMON + 16)
#define MPL_PTR_VAL2            (LIST_DATA_SIZE + PARAM_N_COMMON + 17)
#define MPL_PTR_DIV1            (LIST_DATA_SIZE + PARAM_N_COMMON + 18)
#define MPL_PTR_DIV2            (LIST_DATA_SIZE + PARAM_N_COMMON + 19)
#define MPL_PTR_FNAME           (LIST_DATA_SIZE + PARAM_N_COMMON + 20)
#define MPL_TYPE                (LIST_DATA_SIZE + PARAM_N_COMMON + 21)
#define MPL_COLMAP              (LIST_DATA_SIZE + PARAM_N_COMMON + 22)
#define MPL_PTR_DET             (LIST_DATA_SIZE + PARAM_N_COMMON + 23)
#define MPL_COLOR_SCALE         (LIST_DATA_SIZE + PARAM_N_COMMON + 24)

/*****************************************************************************/

/* FINIX block */

#define FINIX_BLOCK_SIZE            (LIST_DATA_SIZE + PARAM_N_COMMON  +  41)

#define FINIX_IDX                   (LIST_DATA_SIZE + PARAM_N_COMMON  +  0)
#define FINIX_PTR_UNI               (LIST_DATA_SIZE + PARAM_N_COMMON  +  1)
#define FINIX_NZ_POW                (LIST_DATA_SIZE + PARAM_N_COMMON  +  2)
#define FINIX_NR_POW                (LIST_DATA_SIZE + PARAM_N_COMMON  +  3)
#define FINIX_PTR_AX                (LIST_DATA_SIZE + PARAM_N_COMMON  +  4)
#define FINIX_PTR_POWER             (LIST_DATA_SIZE + PARAM_N_COMMON  +  5)
#define FINIX_NZ                    (LIST_DATA_SIZE + PARAM_N_COMMON  +  6)
#define FINIX_NR_FUEL               (LIST_DATA_SIZE + PARAM_N_COMMON  +  7)
#define FINIX_NR_CLAD               (LIST_DATA_SIZE + PARAM_N_COMMON  +  8)
#define FINIX_ZMIN                  (LIST_DATA_SIZE + PARAM_N_COMMON  +  9)
#define FINIX_ZMAX                  (LIST_DATA_SIZE + PARAM_N_COMMON  +  10)
#define FINIX_RMIN                  (LIST_DATA_SIZE + PARAM_N_COMMON  +  11)
#define FINIX_RMAX                  (LIST_DATA_SIZE + PARAM_N_COMMON  +  12)
#define FINIX_RODTYPE               (LIST_DATA_SIZE + PARAM_N_COMMON  +  13)
#define FINIX_PTR_ROD               (LIST_DATA_SIZE + PARAM_N_COMMON  +  14)
#define FINIX_PTR_BC                (LIST_DATA_SIZE + PARAM_N_COMMON  +  15)
#define FINIX_PTR_SCENARIO          (LIST_DATA_SIZE + PARAM_N_COMMON  +  16)
#define FINIX_PTR_RESULTS           (LIST_DATA_SIZE + PARAM_N_COMMON  +  17)
#define FINIX_PTR_PDEN              (LIST_DATA_SIZE + PARAM_N_COMMON  +  18)
#define FINIX_PTR_BU                (LIST_DATA_SIZE + PARAM_N_COMMON  +  19)
#define FINIX_PTR_PARAMS            (LIST_DATA_SIZE + PARAM_N_COMMON  +  20)
#define FINIX_PTR_BCOND             (LIST_DATA_SIZE + PARAM_N_COMMON  +  21)
#define FINIX_PTR_SRESULTS          (LIST_DATA_SIZE + PARAM_N_COMMON  +  22)
#define FINIX_PTR_VRESULTS          (LIST_DATA_SIZE + PARAM_N_COMMON  +  23)
#define FINIX_PTR_OPTIONS           (LIST_DATA_SIZE + PARAM_N_COMMON  +  24)
#define FINIX_BCTYPE                (LIST_DATA_SIZE + PARAM_N_COMMON  +  25)
#define FINIX_PTR_UNI_NAME          (LIST_DATA_SIZE + PARAM_N_COMMON  +  26)
#define FINIX_IFC_NREG              (LIST_DATA_SIZE + PARAM_N_COMMON  +  27)
#define FINIX_PTR_NEST              (LIST_DATA_SIZE + PARAM_N_COMMON  +  28)
#define FINIX_PTR_LHR               (LIST_DATA_SIZE + PARAM_N_COMMON  +  29)
#define FINIX_CMIN                  (LIST_DATA_SIZE + PARAM_N_COMMON  +  30)
#define FINIX_CMAX                  (LIST_DATA_SIZE + PARAM_N_COMMON  +  31)
#define FINIX_N_PIN                 (LIST_DATA_SIZE + PARAM_N_COMMON  +  32)
#define FINIX_PTR_FUEP              (LIST_DATA_SIZE + PARAM_N_COMMON  +  33)
#define FINIX_PTR_IFC               (LIST_DATA_SIZE + PARAM_N_COMMON  +  34)
#define FINIX_PTR_IFC_FNAME         (LIST_DATA_SIZE + PARAM_N_COMMON  +  35)
#define FINIX_PTR_RODNAME           (LIST_DATA_SIZE + PARAM_N_COMMON  +  36)
#define FINIX_PTR_OPTINAME          (LIST_DATA_SIZE + PARAM_N_COMMON  +  37)
#define FINIX_PTR_SCENNAME          (LIST_DATA_SIZE + PARAM_N_COMMON  +  38)
#define FINIX_PTR_POWMSH            (LIST_DATA_SIZE + PARAM_N_COMMON  +  39)
#define FINIX_N_RODS                (LIST_DATA_SIZE + PARAM_N_COMMON  +  40)

/* List of FINIX axial segments */

#define FINIX_AX_BLOCK_SIZE         (LIST_DATA_SIZE +  7)

#define FINIX_AX_IDX                (LIST_DATA_SIZE +  0)
#define FINIX_AX_ZMIN               (LIST_DATA_SIZE +  1)
#define FINIX_AX_ZMAX               (LIST_DATA_SIZE +  2)
#define FINIX_AX_CLADT              (LIST_DATA_SIZE +  3)
#define FINIX_AX_COOLT              (LIST_DATA_SIZE +  4)
#define FINIX_AX_HTCOE              (LIST_DATA_SIZE +  5)
#define FINIX_AX_HFLUX              (LIST_DATA_SIZE +  6)

/*****************************************************************************/

/***** Normalization *********************************************************/

/* NOTE: ton rakenteen yli loopataan setnormalization.c:ss� */

#define NORM_BLOCK_SIZE           (LIST_DATA_SIZE + 16)

#define NORM_POWER                (LIST_DATA_SIZE + 0)
#define NORM_POWDENS              (LIST_DATA_SIZE + 1)
#define NORM_GENRATE              (LIST_DATA_SIZE + 2)
#define NORM_FISSRATE             (LIST_DATA_SIZE + 3)
#define NORM_ABSRATE              (LIST_DATA_SIZE + 4)
#define NORM_LOSSRATE             (LIST_DATA_SIZE + 5)
#define NORM_FLUX                 (LIST_DATA_SIZE + 6)
#define NORM_SRCRATE              (LIST_DATA_SIZE + 7)
#define NORM_SFRATE               (LIST_DATA_SIZE + 8)
#define NORM_PTR_MAT              (LIST_DATA_SIZE + 9)
#define NORM_PTR_FISSRATE         (LIST_DATA_SIZE + 10)
#define NORM_PTR_NSF              (LIST_DATA_SIZE + 11)
#define NORM_PTR_FISSE            (LIST_DATA_SIZE + 12)
#define NORM_PTR_NEUTRON_FLUX     (LIST_DATA_SIZE + 13)
#define NORM_PTR_PHOTON_FLUX      (LIST_DATA_SIZE + 14)
#define NORM_PTR_PHOTON_HEATRATE  (LIST_DATA_SIZE + 15)

/*****************************************************************************/

/***** Mesh ******************************************************************/

#define MESH_BLOCK_SIZE         (LIST_DATA_SIZE + 22)

#define MESH_TYPE               (LIST_DATA_SIZE +  0)
#define MESH_CONTENT            (LIST_DATA_SIZE +  1)
#define MESH_CONTENT_DATA_TYPE  (LIST_DATA_SIZE +  2)
#define MESH_N0                 (LIST_DATA_SIZE +  3)
#define MESH_N1                 (LIST_DATA_SIZE +  4)
#define MESH_N2                 (LIST_DATA_SIZE +  5)
#define MESH_MIN0               (LIST_DATA_SIZE +  6)
#define MESH_MAX0               (LIST_DATA_SIZE +  7)
#define MESH_MIN1               (LIST_DATA_SIZE +  8)
#define MESH_MAX1               (LIST_DATA_SIZE +  9)
#define MESH_MIN2               (LIST_DATA_SIZE + 10)
#define MESH_MAX2               (LIST_DATA_SIZE + 11)
#define MESH_PTR_RES2           (LIST_DATA_SIZE + 12)
#define MESH_PTR_DATA           (LIST_DATA_SIZE + 13)
#define MESH_PTR_PTR            (LIST_DATA_SIZE + 14)
#define MESH_ORTHO_PTR_XLIM     (LIST_DATA_SIZE + 15)
#define MESH_ORTHO_PTR_YLIM     (LIST_DATA_SIZE + 16)
#define MESH_ORTHO_PTR_ZLIM     (LIST_DATA_SIZE + 17)
#define MESH_ADA_SPLIT          (LIST_DATA_SIZE + 18)
#define MESH_ADA_PTR_SZ         (LIST_DATA_SIZE + 19)
#define MESH_LOCAL_COORDS       (LIST_DATA_SIZE + 20)
#define MESH_PTR_NAME           (LIST_DATA_SIZE + 21)

/*****************************************************************************/

/***** Core power distribution ***********************************************/

#define CPD_BLOCK_SIZE          (LIST_DATA_SIZE + 2)

#define CPD_PTR_LAT             (LIST_DATA_SIZE + 0)
#define CPD_COL_COUNT           (LIST_DATA_SIZE + 1)

/*****************************************************************************/

/***** Fission matrix ********************************************************/

#define FMTX_BLOCK_SIZE  (LIST_DATA_SIZE + 7)

#define FMTX_PTR_MTX     (LIST_DATA_SIZE + 0)
#define FMTX_PTR_SRC     (LIST_DATA_SIZE + 1)
#define FMTX_SIZE        (LIST_DATA_SIZE + 2)
#define FMTX_PTR_MAT     (LIST_DATA_SIZE + 3)
#define FMTX_PTR_UNI     (LIST_DATA_SIZE + 4)
#define FMTX_LVL         (LIST_DATA_SIZE + 5)
#define FMTX_PTR_MESH    (LIST_DATA_SIZE + 6)

/*****************************************************************************/

/***** Material divider for burnup calculation *******************************/

#define DIV_BLOCK_SIZE           (LIST_DATA_SIZE + PARAM_N_COMMON + 21)

#define DIV_PTR_MAT              (LIST_DATA_SIZE + PARAM_N_COMMON + 0)
#define DIV_NX                   (LIST_DATA_SIZE + PARAM_N_COMMON + 1)
#define DIV_XMIN                 (LIST_DATA_SIZE + PARAM_N_COMMON + 2)
#define DIV_XMAX                 (LIST_DATA_SIZE + PARAM_N_COMMON + 3)
#define DIV_NY                   (LIST_DATA_SIZE + PARAM_N_COMMON + 4)
#define DIV_YMIN                 (LIST_DATA_SIZE + PARAM_N_COMMON + 5)
#define DIV_YMAX                 (LIST_DATA_SIZE + PARAM_N_COMMON + 6)
#define DIV_NZ                   (LIST_DATA_SIZE + PARAM_N_COMMON + 7)
#define DIV_ZMIN                 (LIST_DATA_SIZE + PARAM_N_COMMON + 8)
#define DIV_ZMAX                 (LIST_DATA_SIZE + PARAM_N_COMMON + 9)
#define DIV_NRAD                 (LIST_DATA_SIZE + PARAM_N_COMMON + 10)
#define DIV_RMIN                 (LIST_DATA_SIZE + PARAM_N_COMMON + 11)
#define DIV_RMAX                 (LIST_DATA_SIZE + PARAM_N_COMMON + 12)
#define DIV_NSEG                 (LIST_DATA_SIZE + PARAM_N_COMMON + 13)
#define DIV_SEG0                 (LIST_DATA_SIZE + PARAM_N_COMMON + 14)
#define DIV_SEP                  (LIST_DATA_SIZE + PARAM_N_COMMON + 15)
#define DIV_SEP_LVL              (LIST_DATA_SIZE + PARAM_N_COMMON + 16)
#define DIV_PTR_MAT_LIST         (LIST_DATA_SIZE + PARAM_N_COMMON + 17)
#define DIV_OUTPUT_FLAG          (LIST_DATA_SIZE + PARAM_N_COMMON + 18)
#define DIV_LVL_MAX              (LIST_DATA_SIZE + PARAM_N_COMMON + 19)
#define DIV_LIMS_CHECK           (LIST_DATA_SIZE + PARAM_N_COMMON + 20)

/* Divisor material lists */

#define DIV_MAT_LIST_BLOCK_SIZE  (LIST_DATA_SIZE + 3)

#define DIV_MAT_LIST_ZONE_IDX    (LIST_DATA_SIZE + 0)
#define DIV_MAT_LIST_PTR_UNIV    (LIST_DATA_SIZE + 1)
#define DIV_MAT_LIST_PTR_REG     (LIST_DATA_SIZE + 2)

/*****************************************************************************/

/***** Universe symmetry *****************************************************/

#define SYMMETRY_BLOCK_SIZE     (LIST_DATA_SIZE + PARAM_N_COMMON + 8)

#define SYMMETRY_PTR_UNI        (LIST_DATA_SIZE + PARAM_N_COMMON + 0)
#define SYMMETRY_AXIS           (LIST_DATA_SIZE + PARAM_N_COMMON + 1)
#define SYMMETRY_BC             (LIST_DATA_SIZE + PARAM_N_COMMON + 2)
#define SYMMETRY_THETA0         (LIST_DATA_SIZE + PARAM_N_COMMON + 3)
#define SYMMETRY_ROT            (LIST_DATA_SIZE + PARAM_N_COMMON + 4)
#define SYMMETRY_X0             (LIST_DATA_SIZE + PARAM_N_COMMON + 5)
#define SYMMETRY_Y0             (LIST_DATA_SIZE + PARAM_N_COMMON + 6)
#define SYMMETRY_SYM            (LIST_DATA_SIZE + PARAM_N_COMMON + 7)

/*****************************************************************************/

/***** Reprocessing **********************************************************/

#define REPROC_BLOCK_SIZE       (LIST_DATA_SIZE + PARAM_N_COMMON + 6)

#define REPROC_PTR_NAME         (LIST_DATA_SIZE + PARAM_N_COMMON + 0)
#define REPROC_OPTIONS          (LIST_DATA_SIZE + PARAM_N_COMMON + 1)
#define REPROC_PTR_SWAP_LIST    (LIST_DATA_SIZE + PARAM_N_COMMON + 2)
#define REPROC_PTR_RPL_LIST     (LIST_DATA_SIZE + PARAM_N_COMMON + 3)
#define REPROC_PTR_REM_LIST     (LIST_DATA_SIZE + PARAM_N_COMMON + 4)
#define REPROC_PTR_CON_LIST     (LIST_DATA_SIZE + PARAM_N_COMMON + 5)

/* Swap list */

#define REPROC_SWAP_BLOCK_SIZE  (LIST_DATA_SIZE + 2)

#define REPROC_SWAP_PTR_UNI1    (LIST_DATA_SIZE + 0)
#define REPROC_SWAP_PTR_UNI2    (LIST_DATA_SIZE + 1)

/* Replace list */

#define REPROC_RPL_BLOCK_SIZE  (LIST_DATA_SIZE + 2)

#define REPROC_RPL_PTR_UNI1    (LIST_DATA_SIZE + 0)
#define REPROC_RPL_PTR_UNI2    (LIST_DATA_SIZE + 1)

/* Removal list */

#define REPROC_REM_BLOCK_SIZE  (LIST_DATA_SIZE + 2)

#define REPROC_REM_PTR_MAT1    (LIST_DATA_SIZE + 0)
#define REPROC_REM_PTR_MAT2    (LIST_DATA_SIZE + 1)

/* Continuous reprocessing */

#define REPROC_CON_BLOCK_SIZE  (LIST_DATA_SIZE + 4)

#define REPROC_CON_PTR_MAT1    (LIST_DATA_SIZE + 0)
#define REPROC_CON_PTR_MAT2    (LIST_DATA_SIZE + 1)
#define REPROC_CON_PTR_MFLOW   (LIST_DATA_SIZE + 2)
#define REPROC_CON_MODE        (LIST_DATA_SIZE + 3)

/*****************************************************************************/

/***** Nuclide mass flow list ************************************************/

#define MFLOW_BLOCK_SIZE        (LIST_DATA_SIZE + PARAM_N_COMMON + 3)

#define MFLOW_PTR_NAME          (LIST_DATA_SIZE + PARAM_N_COMMON + 0)
#define MFLOW_OPTIONS           (LIST_DATA_SIZE + PARAM_N_COMMON + 1)
#define MFLOW_PTR_DATA          (LIST_DATA_SIZE + PARAM_N_COMMON + 2)

#define MFLOW_LIST_BLOCK_SIZE   (LIST_DATA_SIZE + 6)

#define MFLOW_LIST_ZAI          (LIST_DATA_SIZE + 0)
#define MFLOW_LIST_RATE         (LIST_DATA_SIZE + 1)
#define MFLOW_LIST_TOT          (LIST_DATA_SIZE + 2)
#define MFLOW_LIST_PTR_ORIG     (LIST_DATA_SIZE + 3)
#define MFLOW_LIST_PTR_ISO0     (LIST_DATA_SIZE + 4)
#define MFLOW_LIST_PTR_ISO1     (LIST_DATA_SIZE + 5)

/*****************************************************************************/

/***** Nuclide inventory list ************************************************/

#define INVENTORY_BLOCK_SIZE    (LIST_DATA_SIZE + 3)

#define INVENTORY_PTR_ENTRY     (LIST_DATA_SIZE + 0)
#define INVENTORY_PTR_NAME      (LIST_DATA_SIZE + 1)
#define INVENTORY_ZAI           (LIST_DATA_SIZE + 2)

/*****************************************************************************/

/***** Assembly discontinuity factors ****************************************/

#define ADF_BLOCK_SIZE          (LIST_DATA_SIZE + 9)

#define ADF_PTR_GCU             (LIST_DATA_SIZE + 0)
#define ADF_PTR_SURF            (LIST_DATA_SIZE + 1)
#define ADF_NSURF               (LIST_DATA_SIZE + 2)
#define ADF_NCORN               (LIST_DATA_SIZE + 3)
#define ADF_PTR_SURF_AREA       (LIST_DATA_SIZE + 4)
#define ADF_PTR_MID_AREA        (LIST_DATA_SIZE + 5)
#define ADF_PTR_CORN_AREA       (LIST_DATA_SIZE + 6)
#define ADF_VOL                 (LIST_DATA_SIZE + 7)
#define ADF_SYM                 (LIST_DATA_SIZE + 8)

/*****************************************************************************/

/***** Interface current method **********************************************/

#define ICM_BLOCK_SIZE          (LIST_DATA_SIZE + PARAM_N_COMMON + 27)

#define ICM_PTR_ID              (LIST_DATA_SIZE + PARAM_N_COMMON +  0)
#define ICM_PTR_SURF            (LIST_DATA_SIZE + PARAM_N_COMMON +  1)
#define ICM_PTR_LAT             (LIST_DATA_SIZE + PARAM_N_COMMON +  2)
#define ICM_NP                  (LIST_DATA_SIZE + PARAM_N_COMMON +  3)
#define ICM_BREAK_PTR_COUNT     (LIST_DATA_SIZE + PARAM_N_COMMON +  4)
#define ICM_X0                  (LIST_DATA_SIZE + PARAM_N_COMMON +  5)
#define ICM_Y0                  (LIST_DATA_SIZE + PARAM_N_COMMON +  6)
#define ICM_Z0                  (LIST_DATA_SIZE + PARAM_N_COMMON +  7)

#define ICM_RES_CURR0           (LIST_DATA_SIZE + PARAM_N_COMMON +  8)
#define ICM_RES_CC1             (LIST_DATA_SIZE + PARAM_N_COMMON +  9)
#define ICM_RES_CC2             (LIST_DATA_SIZE + PARAM_N_COMMON + 10)
#define ICM_RES_AFLX1           (LIST_DATA_SIZE + PARAM_N_COMMON + 11)
#define ICM_RES_AFLX2           (LIST_DATA_SIZE + PARAM_N_COMMON + 12)
#define ICM_RES_ASRC1           (LIST_DATA_SIZE + PARAM_N_COMMON + 13)
#define ICM_RES_ASRC2           (LIST_DATA_SIZE + PARAM_N_COMMON + 14)
#define ICM_RES_AFISS1          (LIST_DATA_SIZE + PARAM_N_COMMON + 15)
#define ICM_RES_AFISS2          (LIST_DATA_SIZE + PARAM_N_COMMON + 16)
#define ICM_RES_AABS1           (LIST_DATA_SIZE + PARAM_N_COMMON + 17)
#define ICM_RES_AABS2           (LIST_DATA_SIZE + PARAM_N_COMMON + 18)
#define ICM_RES_APOW1           (LIST_DATA_SIZE + PARAM_N_COMMON + 19)
#define ICM_RES_APOW2           (LIST_DATA_SIZE + PARAM_N_COMMON + 20)
#define ICM_RES_LEAK1           (LIST_DATA_SIZE + PARAM_N_COMMON + 21)
#define ICM_RES_LEAK2           (LIST_DATA_SIZE + PARAM_N_COMMON + 22)
#define ICM_RES_PPOW1           (LIST_DATA_SIZE + PARAM_N_COMMON + 23)
#define ICM_RES_PPOW2           (LIST_DATA_SIZE + PARAM_N_COMMON + 24)
#define ICM_RES_PFLX1           (LIST_DATA_SIZE + PARAM_N_COMMON + 25)
#define ICM_RES_PFLX2           (LIST_DATA_SIZE + PARAM_N_COMMON + 26)

/*****************************************************************************/

/***** Pin power distributions ***********************************************/

#define PPW_BLOCK_SIZE          (LIST_DATA_SIZE +  4)

#define PPW_PTR_GCU             (LIST_DATA_SIZE +  0)
#define PPW_PTR_LAT             (LIST_DATA_SIZE +  1)
#define PPW_LAT_TYPE            (LIST_DATA_SIZE +  2)
#define PPW_NP                  (LIST_DATA_SIZE +  3)

/*****************************************************************************/

/***** Albedo (group constant) ***********************************************/

#define ALB_BLOCK_SIZE          (LIST_DATA_SIZE +  4)

#define ALB_PTR_GCU             (LIST_DATA_SIZE +  0)
#define ALB_PTR_SURF            (LIST_DATA_SIZE +  1)
#define ALB_DIR                 (LIST_DATA_SIZE +  2)
#define ALB_NSURF               (LIST_DATA_SIZE +  3)

/*****************************************************************************/

/***** Material volumes list *************************************************/

#define MVOL_BLOCK_SIZE         (LIST_DATA_SIZE + 3)

#define MVOL_PTR_MAT            (LIST_DATA_SIZE + 0)
#define MVOL_REG_IDX            (LIST_DATA_SIZE + 1)
#define MVOL_VOL                (LIST_DATA_SIZE + 2)

/*****************************************************************************/

/***** Additional cross sections to majorant *********************************/

#define MAJORANT_EXTRA_BLOCK_SIZE     (LIST_DATA_SIZE + 4)

#define MAJORANT_EXTRA_TYPE           (LIST_DATA_SIZE + 0)
#define MAJORANT_EXTRA_PTR_NUC        (LIST_DATA_SIZE + 1)
#define MAJORANT_EXTRA_PTR_MAT        (LIST_DATA_SIZE + 2)
#define MAJORANT_EXTRA_FRAC           (LIST_DATA_SIZE + 3)

/*****************************************************************************/

/***** Mora output ***********************************************************/

#define MORA_BLOCK_SIZE         (LIST_DATA_SIZE + 16)

#define MORA_PTR_FNAME          (LIST_DATA_SIZE +  0)
#define MORA_PTR_UNIV           (LIST_DATA_SIZE +  1)
#define MORA_PTR_EG             (LIST_DATA_SIZE +  2)
#define MORA_N_EG               (LIST_DATA_SIZE +  3)
#define MORA_N_COS              (LIST_DATA_SIZE +  4)
#define MORA_PTR_TOT            (LIST_DATA_SIZE +  5)
#define MORA_PTR_CAPT           (LIST_DATA_SIZE +  6)
#define MORA_PTR_FISS           (LIST_DATA_SIZE +  7)
#define MORA_PTR_SCATTP         (LIST_DATA_SIZE +  8)
#define MORA_PTR_SCATTW         (LIST_DATA_SIZE +  9)
#define MORA_PTR_PNU            (LIST_DATA_SIZE + 10)
#define MORA_PTR_DNU            (LIST_DATA_SIZE + 11)
#define MORA_PTR_KAPPA          (LIST_DATA_SIZE + 12)
#define MORA_PTR_CHIP           (LIST_DATA_SIZE + 13)
#define MORA_PTR_CHID           (LIST_DATA_SIZE + 14)
#define MORA_PTR_FLX            (LIST_DATA_SIZE + 15)

/*****************************************************************************/

/***** Group constant generation *********************************************/

#define GCU_BLOCK_SIZE                 (LIST_DATA_SIZE + 365)

#define GCU_PTR_UNIV                   (LIST_DATA_SIZE +  0)
#define GCU_PTR_ADF                    (LIST_DATA_SIZE +  1)
#define GCU_B1_CONV                    (LIST_DATA_SIZE +  2)

/* T�� varmaan my�hemmin pois */

#define GCU_PTR_MORA                   (LIST_DATA_SIZE +  3)

/* Discontinuity factors, etc. */

#define GCU_RES_FG_DF_HET_SURF_FLUX    (LIST_DATA_SIZE +  4)
#define GCU_RES_FG_DF_HET_CORN_FLUX    (LIST_DATA_SIZE +  5)
#define GCU_RES_FG_DF_HET_VOL_FLUX     (LIST_DATA_SIZE +  6)
#define GCU_RES_FG_DF_SURF_DF          (LIST_DATA_SIZE +  7)
#define GCU_RES_FG_DF_CORN_DF          (LIST_DATA_SIZE +  8)
#define GCU_RES_FG_DF_SURF_IN_CURR     (LIST_DATA_SIZE +  9)
#define GCU_RES_FG_DF_SURF_OUT_CURR    (LIST_DATA_SIZE + 10)
#define GCU_RES_FG_DF_SURF_NET_CURR    (LIST_DATA_SIZE + 11)
#define GCU_RES_FG_DF_MID_IN_CURR      (LIST_DATA_SIZE + 12)
#define GCU_RES_FG_DF_MID_OUT_CURR     (LIST_DATA_SIZE + 13)
#define GCU_RES_FG_DF_MID_NET_CURR     (LIST_DATA_SIZE + 14)
#define GCU_RES_FG_DF_CORN_IN_CURR     (LIST_DATA_SIZE + 15)
#define GCU_RES_FG_DF_CORN_OUT_CURR    (LIST_DATA_SIZE + 16)
#define GCU_RES_FG_DF_CORN_NET_CURR    (LIST_DATA_SIZE + 17)
#define GCU_RES_FG_DF_HOM_SURF_FLUX    (LIST_DATA_SIZE + 18)
#define GCU_RES_FG_DF_HOM_CORN_FLUX    (LIST_DATA_SIZE + 19)
#define GCU_RES_FG_DF_HOM_VOL_FLUX     (LIST_DATA_SIZE + 20)

/* Diffusion coefficient for E.D */

#define GCU_RES_DIFFCOEF_ED            (LIST_DATA_SIZE + 21)

/* Pin-power distributions */

#define GCU_PTR_PPW                    (LIST_DATA_SIZE + 22)
#define GCU_RES_FG_PPW_POW             (LIST_DATA_SIZE + 23)
#define GCU_RES_FG_PPW_HOM_FLUX        (LIST_DATA_SIZE + 24)
#define GCU_RES_FG_PPW_FF              (LIST_DATA_SIZE + 25)
#define GCU_RES_FG_PPW_XYZ             (LIST_DATA_SIZE + 26)

/* Albedos */

#define GCU_PTR_ALB                    (LIST_DATA_SIZE + 27)
#define GCU_RES_FG_ALB_IN_CURR         (LIST_DATA_SIZE + 28)
#define GCU_RES_FG_ALB_OUT_CURR        (LIST_DATA_SIZE + 29)
#define GCU_RES_FG_TOT_ALB             (LIST_DATA_SIZE + 30)
#define GCU_RES_FG_PART_ALB            (LIST_DATA_SIZE + 31)

/* Micro-group data */

#define GCU_MICRO_FLX                  (LIST_DATA_SIZE +  32)
#define GCU_MICRO_FISS_FLX             (LIST_DATA_SIZE +  33)
#define GCU_MICRO_TOT                  (LIST_DATA_SIZE +  34)
#define GCU_MICRO_ABS                  (LIST_DATA_SIZE +  35)
#define GCU_MICRO_FISS                 (LIST_DATA_SIZE +  36)
#define GCU_MICRO_NSF                  (LIST_DATA_SIZE +  37)
#define GCU_MICRO_FISSE                (LIST_DATA_SIZE +  38)
#define GCU_MICRO_INV_V                (LIST_DATA_SIZE +  39)
#define GCU_MICRO_CHIT                 (LIST_DATA_SIZE +  40)
#define GCU_MICRO_CHIP                 (LIST_DATA_SIZE +  41)
#define GCU_MICRO_CHID                 (LIST_DATA_SIZE +  42)
#define GCU_MICRO_SCATT0               (LIST_DATA_SIZE +  43)
#define GCU_MICRO_SCATT1               (LIST_DATA_SIZE +  44)
#define GCU_MICRO_SCATT2               (LIST_DATA_SIZE +  45)
#define GCU_MICRO_SCATT3               (LIST_DATA_SIZE +  46)
#define GCU_MICRO_SCATT4               (LIST_DATA_SIZE +  47)
#define GCU_MICRO_SCATT5               (LIST_DATA_SIZE +  48)
#define GCU_MICRO_SCATT6               (LIST_DATA_SIZE +  49)
#define GCU_MICRO_SCATT7               (LIST_DATA_SIZE +  50)
#define GCU_MICRO_SCATTP0              (LIST_DATA_SIZE +  51)
#define GCU_MICRO_SCATTP1              (LIST_DATA_SIZE +  52)
#define GCU_MICRO_SCATTP2              (LIST_DATA_SIZE +  53)
#define GCU_MICRO_SCATTP3              (LIST_DATA_SIZE +  54)
#define GCU_MICRO_SCATTP4              (LIST_DATA_SIZE +  55)
#define GCU_MICRO_SCATTP5              (LIST_DATA_SIZE +  56)
#define GCU_MICRO_SCATTP6              (LIST_DATA_SIZE +  57)
#define GCU_MICRO_SCATTP7              (LIST_DATA_SIZE +  58)
#define GCU_MICRO_I135_YIELD           (LIST_DATA_SIZE +  59)
#define GCU_MICRO_XE135_YIELD          (LIST_DATA_SIZE +  60)
#define GCU_MICRO_PM147_YIELD          (LIST_DATA_SIZE +  61)
#define GCU_MICRO_PM148_YIELD          (LIST_DATA_SIZE +  62)
#define GCU_MICRO_PM148M_YIELD         (LIST_DATA_SIZE +  63)
#define GCU_MICRO_PM149_YIELD          (LIST_DATA_SIZE +  64)
#define GCU_MICRO_SM149_YIELD          (LIST_DATA_SIZE +  65)
#define GCU_MICRO_I135_ABS             (LIST_DATA_SIZE +  66)
#define GCU_MICRO_XE135_ABS            (LIST_DATA_SIZE +  67)
#define GCU_MICRO_PM147_ABS            (LIST_DATA_SIZE +  68)
#define GCU_MICRO_PM148_ABS            (LIST_DATA_SIZE +  69)
#define GCU_MICRO_PM148M_ABS           (LIST_DATA_SIZE +  70)
#define GCU_MICRO_PM149_ABS            (LIST_DATA_SIZE +  71)
#define GCU_MICRO_SM149_ABS            (LIST_DATA_SIZE +  72)
#define GCU_MICRO_XE135_MACRO_ABS      (LIST_DATA_SIZE +  73)
#define GCU_MICRO_SM149_MACRO_ABS      (LIST_DATA_SIZE +  74)
#define GCU_MICRO_ADF_SURF_FLUX        (LIST_DATA_SIZE +  75)
#define GCU_MICRO_ADF_CORN_FLUX        (LIST_DATA_SIZE +  76)
#define GCU_MICRO_ADF_CELL_FLUX        (LIST_DATA_SIZE +  77)
#define GCU_MICRO_ADF_SURF_IN_CURR     (LIST_DATA_SIZE +  78)
#define GCU_MICRO_ADF_SURF_OUT_CURR    (LIST_DATA_SIZE +  79)
#define GCU_MICRO_ADF_MID_IN_CURR      (LIST_DATA_SIZE +  80)
#define GCU_MICRO_ADF_MID_OUT_CURR     (LIST_DATA_SIZE +  81)
#define GCU_MICRO_ADF_CORN_IN_CURR     (LIST_DATA_SIZE +  82)
#define GCU_MICRO_ADF_CORN_OUT_CURR    (LIST_DATA_SIZE +  83)
#define GCU_MICRO_B1_FLX               (LIST_DATA_SIZE +  84)
#define GCU_MICRO_B1_DIFFCOEF          (LIST_DATA_SIZE +  85)
#define GCU_MICRO_DIFFCOEF_ED          (LIST_DATA_SIZE +  86)
#define GCU_MICRO_PPW_POW              (LIST_DATA_SIZE +  87)
#define GCU_MICRO_PPW_XYZ              (LIST_DATA_SIZE +  88)
#define GCU_MICRO_ALB_IN_CURR          (LIST_DATA_SIZE +  89)
#define GCU_MICRO_ALB_OUT_CURR         (LIST_DATA_SIZE +  90)
#define GCU_INF_FLX                    (LIST_DATA_SIZE +  91)
#define GCU_INF_MICRO_FLX              (LIST_DATA_SIZE +  92)
#define GCU_INF_FISS_FLX               (LIST_DATA_SIZE +  93)
#define GCU_INF_KINF                   (LIST_DATA_SIZE +  94)
#define GCU_INF_TOT                    (LIST_DATA_SIZE +  97)
#define GCU_INF_CAPT                   (LIST_DATA_SIZE +  98)
#define GCU_INF_FISS                   (LIST_DATA_SIZE +  99)
#define GCU_INF_ABS                    (LIST_DATA_SIZE + 100)
#define GCU_INF_NSF                    (LIST_DATA_SIZE + 101)
#define GCU_INF_KAPPA                  (LIST_DATA_SIZE + 102)
#define GCU_INF_INVV                   (LIST_DATA_SIZE + 103)
#define GCU_INF_NUBAR                  (LIST_DATA_SIZE + 104)
#define GCU_INF_RABSXS                 (LIST_DATA_SIZE + 105)
#define GCU_INF_REMXS                  (LIST_DATA_SIZE + 106)
#define GCU_INF_CHIT                   (LIST_DATA_SIZE + 107)
#define GCU_INF_CHIP                   (LIST_DATA_SIZE + 108)
#define GCU_INF_CHID                   (LIST_DATA_SIZE + 109)
#define GCU_INF_S0                     (LIST_DATA_SIZE + 110)
#define GCU_INF_S1                     (LIST_DATA_SIZE + 111)
#define GCU_INF_S2                     (LIST_DATA_SIZE + 112)
#define GCU_INF_S3                     (LIST_DATA_SIZE + 113)
#define GCU_INF_S4                     (LIST_DATA_SIZE + 114)
#define GCU_INF_S5                     (LIST_DATA_SIZE + 115)
#define GCU_INF_S6                     (LIST_DATA_SIZE + 116)
#define GCU_INF_S7                     (LIST_DATA_SIZE + 117)
#define GCU_INF_SP0                    (LIST_DATA_SIZE + 118)
#define GCU_INF_SP1                    (LIST_DATA_SIZE + 119)
#define GCU_INF_SP2                    (LIST_DATA_SIZE + 120)
#define GCU_INF_SP3                    (LIST_DATA_SIZE + 121)
#define GCU_INF_SP4                    (LIST_DATA_SIZE + 122)
#define GCU_INF_SP5                    (LIST_DATA_SIZE + 123)
#define GCU_INF_SP6                    (LIST_DATA_SIZE + 124)
#define GCU_INF_SP7                    (LIST_DATA_SIZE + 125)
#define GCU_INF_SCATT0                 (LIST_DATA_SIZE + 126)
#define GCU_INF_SCATT1                 (LIST_DATA_SIZE + 127)
#define GCU_INF_SCATT2                 (LIST_DATA_SIZE + 128)
#define GCU_INF_SCATT3                 (LIST_DATA_SIZE + 129)
#define GCU_INF_SCATT4                 (LIST_DATA_SIZE + 130)
#define GCU_INF_SCATT5                 (LIST_DATA_SIZE + 131)
#define GCU_INF_SCATT6                 (LIST_DATA_SIZE + 132)
#define GCU_INF_SCATT7                 (LIST_DATA_SIZE + 133)
#define GCU_INF_SCATTP0                (LIST_DATA_SIZE + 134)
#define GCU_INF_SCATTP1                (LIST_DATA_SIZE + 135)
#define GCU_INF_SCATTP2                (LIST_DATA_SIZE + 136)
#define GCU_INF_SCATTP3                (LIST_DATA_SIZE + 137)
#define GCU_INF_SCATTP4                (LIST_DATA_SIZE + 138)
#define GCU_INF_SCATTP5                (LIST_DATA_SIZE + 139)
#define GCU_INF_SCATTP6                (LIST_DATA_SIZE + 140)
#define GCU_INF_SCATTP7                (LIST_DATA_SIZE + 141)
#define GCU_INF_I135_YIELD             (LIST_DATA_SIZE + 142)
#define GCU_INF_XE135_YIELD            (LIST_DATA_SIZE + 143)
#define GCU_INF_PM147_YIELD            (LIST_DATA_SIZE + 144)
#define GCU_INF_PM148_YIELD            (LIST_DATA_SIZE + 145)
#define GCU_INF_PM148M_YIELD           (LIST_DATA_SIZE + 146)
#define GCU_INF_PM149_YIELD            (LIST_DATA_SIZE + 147)
#define GCU_INF_SM149_YIELD            (LIST_DATA_SIZE + 148)
#define GCU_INF_I135_ABS               (LIST_DATA_SIZE + 149)
#define GCU_INF_XE135_ABS              (LIST_DATA_SIZE + 150)
#define GCU_INF_PM147_ABS              (LIST_DATA_SIZE + 151)
#define GCU_INF_PM148_ABS              (LIST_DATA_SIZE + 152)
#define GCU_INF_PM148M_ABS             (LIST_DATA_SIZE + 153)
#define GCU_INF_PM149_ABS              (LIST_DATA_SIZE + 154)
#define GCU_INF_SM149_ABS              (LIST_DATA_SIZE + 155)
#define GCU_INF_XE135_MACRO_ABS        (LIST_DATA_SIZE + 157)
#define GCU_INF_SM149_MACRO_ABS        (LIST_DATA_SIZE + 159)
#define GCU_INF_TRANSPXS               (LIST_DATA_SIZE + 160)
#define GCU_INF_DIFFCOEF               (LIST_DATA_SIZE + 161)
#define GCU_B1_KEFF                    (LIST_DATA_SIZE + 162)
#define GCU_B1_KINF                    (LIST_DATA_SIZE + 163)
#define GCU_B1_B2                      (LIST_DATA_SIZE + 164)
#define GCU_B1_ERR                     (LIST_DATA_SIZE + 165)
#define GCU_B1_FISS_FLX                (LIST_DATA_SIZE + 166)
#define GCU_B1_MICRO_FLX               (LIST_DATA_SIZE + 167)
#define GCU_B1_FLX                     (LIST_DATA_SIZE + 168)
#define GCU_B1_TOT                     (LIST_DATA_SIZE + 169)
#define GCU_B1_CAPT                    (LIST_DATA_SIZE + 170)
#define GCU_B1_FISS                    (LIST_DATA_SIZE + 171)
#define GCU_B1_ABS                     (LIST_DATA_SIZE + 172)
#define GCU_B1_NSF                     (LIST_DATA_SIZE + 173)
#define GCU_B1_KAPPA                   (LIST_DATA_SIZE + 174)
#define GCU_B1_INVV                    (LIST_DATA_SIZE + 175)
#define GCU_B1_NUBAR                   (LIST_DATA_SIZE + 176)
#define GCU_B1_RABSXS                  (LIST_DATA_SIZE + 177)
#define GCU_B1_REMXS                   (LIST_DATA_SIZE + 178)
#define GCU_B1_CHIT                    (LIST_DATA_SIZE + 179)
#define GCU_B1_CHIP                    (LIST_DATA_SIZE + 180)
#define GCU_B1_CHID                    (LIST_DATA_SIZE + 183)
#define GCU_B1_S0                      (LIST_DATA_SIZE + 184)
#define GCU_B1_S1                      (LIST_DATA_SIZE + 185)
#define GCU_B1_S2                      (LIST_DATA_SIZE + 186)
#define GCU_B1_S3                      (LIST_DATA_SIZE + 187)
#define GCU_B1_S4                      (LIST_DATA_SIZE + 188)
#define GCU_B1_S5                      (LIST_DATA_SIZE + 189)
#define GCU_B1_S6                      (LIST_DATA_SIZE + 190)
#define GCU_B1_S7                      (LIST_DATA_SIZE + 191)
#define GCU_B1_SP0                     (LIST_DATA_SIZE + 192)
#define GCU_B1_SP1                     (LIST_DATA_SIZE + 193)
#define GCU_B1_SP2                     (LIST_DATA_SIZE + 194)
#define GCU_B1_SP3                     (LIST_DATA_SIZE + 195)
#define GCU_B1_SP4                     (LIST_DATA_SIZE + 196)
#define GCU_B1_SP5                     (LIST_DATA_SIZE + 197)
#define GCU_B1_SP6                     (LIST_DATA_SIZE + 198)
#define GCU_B1_SP7                     (LIST_DATA_SIZE + 199)
#define GCU_B1_SCATT0                  (LIST_DATA_SIZE + 200)
#define GCU_B1_SCATT1                  (LIST_DATA_SIZE + 201)
#define GCU_B1_SCATT2                  (LIST_DATA_SIZE + 202)
#define GCU_B1_SCATT3                  (LIST_DATA_SIZE + 203)
#define GCU_B1_SCATT4                  (LIST_DATA_SIZE + 204)
#define GCU_B1_SCATT5                  (LIST_DATA_SIZE + 205)
#define GCU_B1_SCATT6                  (LIST_DATA_SIZE + 206)
#define GCU_B1_SCATT7                  (LIST_DATA_SIZE + 207)
#define GCU_B1_SCATTP0                 (LIST_DATA_SIZE + 208)
#define GCU_B1_SCATTP1                 (LIST_DATA_SIZE + 209)
#define GCU_B1_SCATTP2                 (LIST_DATA_SIZE + 210)
#define GCU_B1_SCATTP3                 (LIST_DATA_SIZE + 211)
#define GCU_B1_SCATTP4                 (LIST_DATA_SIZE + 212)
#define GCU_B1_SCATTP5                 (LIST_DATA_SIZE + 213)
#define GCU_B1_SCATTP6                 (LIST_DATA_SIZE + 214)
#define GCU_B1_SCATTP7                 (LIST_DATA_SIZE + 215)
#define GCU_B1_I135_YIELD              (LIST_DATA_SIZE + 216)
#define GCU_B1_XE135_YIELD             (LIST_DATA_SIZE + 217)
#define GCU_B1_PM147_YIELD             (LIST_DATA_SIZE + 218)
#define GCU_B1_PM148_YIELD             (LIST_DATA_SIZE + 219)
#define GCU_B1_PM148M_YIELD            (LIST_DATA_SIZE + 220)
#define GCU_B1_PM149_YIELD             (LIST_DATA_SIZE + 221)
#define GCU_B1_SM149_YIELD             (LIST_DATA_SIZE + 222)
#define GCU_B1_I135_ABS                (LIST_DATA_SIZE + 223)
#define GCU_B1_XE135_ABS               (LIST_DATA_SIZE + 224)
#define GCU_B1_PM147_ABS               (LIST_DATA_SIZE + 225)
#define GCU_B1_PM148_ABS               (LIST_DATA_SIZE + 226)
#define GCU_B1_PM148M_ABS              (LIST_DATA_SIZE + 227)
#define GCU_B1_PM149_ABS               (LIST_DATA_SIZE + 228)
#define GCU_B1_SM149_ABS               (LIST_DATA_SIZE + 229)
#define GCU_B1_XE135_MACRO_ABS         (LIST_DATA_SIZE + 230)
#define GCU_B1_SM149_MACRO_ABS         (LIST_DATA_SIZE + 231)
#define GCU_B1_TRANSPXS                (LIST_DATA_SIZE + 232)
#define GCU_B1_DIFFCOEF                (LIST_DATA_SIZE + 233)

#define GCU_MEULEKAMP_TOT_FISS         (LIST_DATA_SIZE + 234)
#define GCU_MEULEKAMP_BETA_EFF         (LIST_DATA_SIZE + 235)
#define GCU_MEULEKAMP_LAMBDA           (LIST_DATA_SIZE + 236)

#define GCU_PTR_FIRST_STAT             (LIST_DATA_SIZE + 237)
#define GCU_PTR_LAST_STAT              (LIST_DATA_SIZE + 238)

/*****************************************************************************/

/***** Depletion branches ****************************************************/

/* Vai olisko "deviation", tms. parempi nimi? */

#define DEP_BRA_BLOCK_SIZE              (LIST_DATA_SIZE + PARAM_N_COMMON + 8)

#define DEP_BRA_PTR_NAME                (LIST_DATA_SIZE + PARAM_N_COMMON + 0)
#define DEP_BRA_PTR_STP                 (LIST_DATA_SIZE + PARAM_N_COMMON + 1)
#define DEP_BRA_PTR_VAR                 (LIST_DATA_SIZE + PARAM_N_COMMON + 2)
#define DEP_BRA_PTR_REPLACE_MAT         (LIST_DATA_SIZE + PARAM_N_COMMON + 3)
#define DEP_BRA_PTR_REPLACE_UNI         (LIST_DATA_SIZE + PARAM_N_COMMON + 4)
#define DEP_BRA_PTR_TRANS               (LIST_DATA_SIZE + PARAM_N_COMMON + 5)
#define DEP_BRA_PTR_GCU                 (LIST_DATA_SIZE + PARAM_N_COMMON + 6)
#define DEP_BRA_NORM                    (LIST_DATA_SIZE + PARAM_N_COMMON + 7)

/* Change in material state */

#define DEP_BRA_STP_BLOCK_SIZE          (LIST_DATA_SIZE + 4)

#define DEP_BRA_STP_PTR_MAT             (LIST_DATA_SIZE + 0)
#define DEP_BRA_STP_DENSITY             (LIST_DATA_SIZE + 1)
#define DEP_BRA_STP_TEMP                (LIST_DATA_SIZE + 2)
#define DEP_BRA_STP_PTR_SAB             (LIST_DATA_SIZE + 3)

#define DEP_BRA_STP_SAB_BLOCK_SIZE      (LIST_DATA_SIZE + 3)

#define DEP_BRA_STP_SAB_PTR_THERM       (LIST_DATA_SIZE + 0)
#define DEP_BRA_STP_SAB_PTR_LIB1        (LIST_DATA_SIZE + 1)
#define DEP_BRA_STP_SAB_PTR_LIB2        (LIST_DATA_SIZE + 2)

/* Variable (to be passed into output) */

#define DEP_BRA_VAR_BLOCK_SIZE          (LIST_DATA_SIZE + 2)

#define DEP_BRA_VAR_PTR_NAME            (LIST_DATA_SIZE + 0)
#define DEP_BRA_VAR_PTR_VALUE           (LIST_DATA_SIZE + 1)

/* Replace material */

#define DEP_BRA_REPLACE_MAT_BLOCK_SIZE  (LIST_DATA_SIZE + 2)

#define DEP_BRA_REPLACE_MAT_PTR_MAT1    (LIST_DATA_SIZE + 0)
#define DEP_BRA_REPLACE_MAT_PTR_MAT2    (LIST_DATA_SIZE + 1)

/* Replace universe */

#define DEP_BRA_REPLACE_UNI_BLOCK_SIZE  (LIST_DATA_SIZE + 2)

#define DEP_BRA_REPLACE_UNI_PTR_UNI1    (LIST_DATA_SIZE + 0)
#define DEP_BRA_REPLACE_UNI_PTR_UNI2    (LIST_DATA_SIZE + 1)

/* Transformation */

#define DEP_BRA_TRANS_BLOCK_SIZE        (LIST_DATA_SIZE + 2)

#define DEP_BRA_TRANS_PTR_UNI           (LIST_DATA_SIZE + 0)
#define DEP_BRA_TRANS_PTR_TRANS         (LIST_DATA_SIZE + 1)

/*****************************************************************************/

/***** Coefficient calculations **********************************************/

#define COEF_BLOCK_SIZE                 (LIST_DATA_SIZE + PARAM_N_COMMON + 4)

#define COEF_N_BU                       (LIST_DATA_SIZE + PARAM_N_COMMON + 0)
#define COEF_PTR_BU_PTS                 (LIST_DATA_SIZE + PARAM_N_COMMON + 1)
#define COEF_N_TOT                      (LIST_DATA_SIZE + PARAM_N_COMMON + 2)
#define COEF_PTR_MTX                    (LIST_DATA_SIZE + PARAM_N_COMMON + 3)

#define COEF_MTX_BLOCK_SIZE             (LIST_DATA_SIZE + 4)

#define COEF_MTX_N_BRA                  (LIST_DATA_SIZE + 0)
#define COEF_MTX_N_CUMU                 (LIST_DATA_SIZE + 1)
#define COEF_MTX_PTR_BRA                (LIST_DATA_SIZE + 2)
#define COEF_MTX_PTR_VAR                (LIST_DATA_SIZE + 3)

/*****************************************************************************/

#ifdef __cplusplus
} /* closing curly bracket */
#endif

#endif
