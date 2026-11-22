MODULE FileNames_Module
   USE Numeric_Kinds_Module
   USE FileUnits_Module
   IMPLICIT NONE
   SAVE
   CHARACTER(LEN=70) filename
   INTEGER(KIND=IW_Kind)           ::  angeigv_unit_pid   !
   INTEGER(KIND=IW_Kind)           ::  angeigf_unit_pid   !
   INTEGER(KIND=IW_Kind)           ::  overlap_unit_pid   !
   INTEGER(KIND=IW_Kind)           ::  potmatsvd_unit_pid !
   INTEGER(KIND=IW_Kind)           ::  potmatdbs_unit_pid !
   CHARACTER(LEN=7), PARAMETER::    BinOut_Dir  = 'BinOut/' 
   CHARACTER(LEN=5), PARAMETER::      done_dir  = 'done/'                   !
   CHARACTER(LEN=8), PARAMETER::      work_dir  = 'working/'                !
   CHARACTER(LEN=6), PARAMETER::     trash_dir  = 'trash/'                  !
   CHARACTER(LEN=15), PARAMETER:: potmatdbs_dir = 'PotMatDBS.bin/'          !
   CHARACTER(LEN=15), PARAMETER:: potmatsvd_dir = 'PotMatSVD.bin/'          !
   CHARACTER(LEN=13), PARAMETER::   overlap_dir = 'Overlap.bin/'            !
   CHARACTER(LEN=13), PARAMETER::    angeig_dir = 'AngEigen.bin/'           !
   CHARACTER(LEN=16), PARAMETER::      prob_dir = 'Probabilities/'          !
   CHARACTER(LEN=10), PARAMETER::      opac_dir = 'Opacities/'              !
   CHARACTER(LEN=10), PARAMETER::     aph3d_dir = 'Aph3d/'                  !
   CHARACTER(LEN=24), PARAMETER::      Psi_name = 'graph_psi.txt'           !
   CHARACTER(LEN=24), PARAMETER::    Wrons_name = 'out_wrons.txt'           !
   CHARACTER(LEN=24), PARAMETER::       Qn_name = 'out_quantum_njl.txt'     !
   CHARACTER(LEN=24), PARAMETER::     Angf_name = 'out_aafns.txt'           !
   CHARACTER(LEN=24), PARAMETER::   Level0_name = 'wg_level0.txt'              !
   CHARACTER(LEN=24), PARAMETER::   Levels_name = 'wg_levels.txt'           !
   CHARACTER(LEN=24), PARAMETER::potmatdbs_name = 'PotMatDBS.bin'           !
   CHARACTER(LEN=24), PARAMETER::potmatsvd_name = 'PotMatSVD.bin'           !
   CHARACTER(LEN=24), PARAMETER::  overlap_name = 'Overlap.bin'             !
   CHARACTER(LEN=24), PARAMETER::  angeigf_file = 'AngEigf'                        !
   CHARACTER(LEN=24), PARAMETER::     rmat_file = 'R'                        !
   CHARACTER(LEN=24), PARAMETER::     prob_file = 'Prob'                        !
   CHARACTER(LEN=24), PARAMETER::     opac_file = 'Opac'                        !
   CHARACTER(LEN=24), PARAMETER::   dia2ad_name = 'Dia2Ad.bin'              !
   CHARACTER(LEN=24), PARAMETER::     Prob_name = 'wg_probs.txt'            !
   CHARACTER(LEN=24), PARAMETER::    Probf_name = 'wg_probf.txt'            !
   CHARACTER(LEN=24), PARAMETER::        BCname = 'out_inp_bc.txt'          !
   CHARACTER(LEN=24), PARAMETER::    HfOut_name = 'out_hf.txt'              !
   CHARACTER(LEN=24), PARAMETER::    HfPlt_name = 'plt_hf.txt'              !
   CHARACTER(LEN=24), PARAMETER::     Diad_name = 'out_adb_dia.txt'         !
   CHARACTER(LEN=24), PARAMETER::    treal_name = 'tmatrix_real.txt'        !
   CHARACTER(LEN=24), PARAMETER::    timag_name = 'tmatrix_imag.txt'        !
   CHARACTER(LEN=24), PARAMETER::     opac_name = 'opacities.txt'           !
   CHARACTER(LEN=24), PARAMETER::       nj_name = 'out_quantum_nj.txt'      !
   CHARACTER(LEN=24), PARAMETER::       cs_name = 'out_cross_sections.txt'  !
   CHARACTER(LEN=24), PARAMETER::    probe_name = 'probf_vs_e.txt'          !
   CHARACTER(LEN=24), PARAMETER::   probqb_name = 'probqb.txt'              !
   CHARACTER(LEN=24), PARAMETER::  probqbf_name = 'probqbf.txt'             !
   CHARACTER(LEN=24), PARAMETER::    angbf_name = 'ang_eigenf.bin'          !
   CHARACTER(LEN=24), PARAMETER::    nonun_name = 'unitarity.txt'           !
   CHARACTER(LEN=24), PARAMETER::  mainrec_name = 'main.rec'                !
   CHARACTER(LEN=24), PARAMETER::  chanrec_name = 'channels.rec'            !
   CHARACTER(LEN=24), PARAMETER::  enerrec_name = 'energies.rec'            !
   CHARACTER(LEN=24), PARAMETER:: energies_data = 'energies.data'           !
   CHARACTER(LEN=24), PARAMETER:: channels_data = 'channels.data'           !
   CHARACTER(LEN=24), PARAMETER::  idnames_data = 'idnames.data'            !
   CHARACTER(LEN=24), PARAMETER::   method_data = 'method.data'             !
   CHARACTER(LEN=24), PARAMETER::     main_data = 'main.data'               !
   CHARACTER(LEN=24), PARAMETER::     prop_data = 'propagation.data'        !
   CHARACTER(LEN=24), PARAMETER::   masses_data = 'masses.data'             !
   CHARACTER(LEN=24), PARAMETER::   arranges_data = 'arrangements.data'     !
   CHARACTER(LEN=24), PARAMETER:: asymptotic_data = 'asymptotic.data'       !
   CHARACTER(LEN=24), PARAMETER::    rnrec_name = 'rn.rec'                  !
   CHARACTER(LEN=24), PARAMETER::    tnrec_name = 'tn.rec'                  !
   CHARACTER(LEN=24), PARAMETER::vecoldrec_name = 'vecold.rec'              !
   CHARACTER(LEN=24), PARAMETER::    Debug_name = 'debug.data'              !
   CHARACTER(LEN=24), PARAMETER:: basisout_name = 'basis.txt'               !
   CHARACTER(LEN=24), PARAMETER::  propout_name = 'propagation.txt'         !
   CHARACTER(LEN=24), PARAMETER::  asymptout_name = 'asympt.txt'            !
   CHARACTER(LEN=24), PARAMETER::     timing_name = 'Timing.txt'            !
   CHARACTER(LEN=24), PARAMETER::       rmat_name = 'R'                     !
   CHARACTER(LEN=24), PARAMETER::  ReadPot_name = 'PotMat.bin'              !
   CHARACTER(LEN=24), PARAMETER::  angeigv_name = 'ang_eigenv.bin'          !
   CHARACTER(LEN=24), PARAMETER::  angeigf_name = 'ang_eigenf.bin'          !
   CHARACTER(LEN=24), PARAMETER::  ReadOvr_name = 'OverDer.bin'             !
   CHARACTER(LEN=24), PARAMETER:: OpacityJ_Name= 'OpacityJ.txt'             !
   CHARACTER(LEN=24), PARAMETER::    Prob_grace_name = 'wg_probs.grace'     !
   CHARACTER(LEN=24), PARAMETER:: psi_grace_name = 'graph_psi_grace'        !
   CHARACTER(LEN=24), PARAMETER:: potmatsvd_file = 'PotmatSVD.bin'          !
   CHARACTER(LEN=24), PARAMETER:: potmatdbs_file = 'PotmatDBS.bin'          !
   CHARACTER(LEN=21), PARAMETER:: URU_File      =BinOut_dir//'URU_Matrix.bin'
   CHARACTER(LEN=22), PARAMETER:: UREU_File     =BinOut_dir//'UREU_Matrix.bin'
   CHARACTER(LEN=19), PARAMETER:: K_File        =BinOut_dir//'K_Matrix.bin'
   CHARACTER(LEN=20), PARAMETER:: KE_File       =BinOut_dir//'KE_Matrix.bin'
END MODULE FileNames_Module
