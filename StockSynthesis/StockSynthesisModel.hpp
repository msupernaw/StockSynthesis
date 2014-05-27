/* 
 * File:   StockSynthesisModel.hpp
 * Author: matthewsupernaw
 *
 * Created on May 27, 2014, 1:32 PM
 */

#ifndef STOCKSYNTHESISMODEL_HPP
#define	STOCKSYNTHESISMODEL_HPP
#include "Model.hpp"
#include <string>
#include <vector>
#include <valarray>

namespace ss {

    template<class REAL_T>
    class StockSynthesisData : public ss::DataModule<REAL_T> {
        int z;
        int z1;
        int z2;
        int L1;
        int L2;
        int A2;
        int a1;
        int f;
        int g;
        int gg;
        int a;
        int b;
        int p;
        int p1;
        int p2;
        int i;
        int y;
        int yz;
        int s;
        int s2;
        int smid;
        int t;
        int j;
        int j1;
        int j2;
        int k;
        int s_off;
        int Fishon;
        int NP;
        int Ip;
        int firstseas;
        int t_base;
        int niter;
        int loop;
        int TG_t;
        int Fcast_catch_start;
        int ParCount;
        int N_SC;
        int N_DC;
        int N_CC;
        int N_FC;
        int icycle;
        int Ncycle;
        int No_Report;
        REAL_T mcmcFlag;
        REAL_T temp;
        REAL_T temp1;
        std::string datfilename;
        std::string ctlfilename;
        int readparfile;
        int rundetail;
        int reportdetail;
        int docheckup;
        int Do_ParmTrace;
        int Do_CumReport;
        int Do_all_priors;
        int SoftBound;
        int N_nudata;
        int Turn_off_phase;
        int burn_intvl;
        int thin_intvl;
        REAL_T jitter;
        int STD_Yr_min;
        int STD_Yr_max;
        int N_STD_Yr_RD;
        int N_STD_Yr;
        std::vector<int> STD_Yr_RD;
        int save_for_report;
        int save_gparm;
        int save_gparm_print;
        int N_warn;
        int mceval_counter;
        int mceval_header;
        int mcmc_counter;
        int done_run;
        std::vector<REAL_T> func_eval;
        std::vector<REAL_T> func_conv;
        REAL_T final_conv;
        int retro_yr;
        int fishery_on_off;
        int Smry_Age;
        int depletion_basis;
        REAL_T depletion_level;
        int SPR_reporting;
        int F_reporting;
        std::vector<int> F_reporting_ages;
        int F_report_basis;
        int finish_starter;
        REAL_T pi;
        REAL_T neglog19;
        REAL_T NilNumbers;
        REAL_T dummy_datum;
        int dummy_phase;
        int runnumber;
        int N_prof_var;
        int prof_var_cnt;
        int prof_junk;
        std::vector<REAL_T> prof_var;
        int styr;
        int endyr;
        int nseas;
        std::vector<REAL_T> seasdur;
        int TimeMax;
        int TimeMax_Fcast_std;
        int YrMax;
        int eq_yr;
        int bio_yr;
        REAL_T sumseas;
        std::vector<REAL_T> seasdur2;
        int spawn_seas;
        int Nfleet;
        int Nsurvey;
        int Ntypes;
        int pop;
        std::valarray<std::valarray<int> > pfleetname;
        std::string fleetnameread;
        std::vector<REAL_T> surveytime;
        std::vector<int> fleet_area;
        std::vector<REAL_T> catchunits;
        std::vector<REAL_T> catch_se_rd;
        std::valarray<std::valarray<REAL_T> > catch_se;
        int gender;
        int nages;
        std::vector<int> age_vector;
        std::vector<REAL_T> r_ages;
        std::vector<REAL_T> frac_ages;
        std::vector<int> years;
        std::vector<REAL_T> r_years;
        std::valarray<std::valarray<REAL_T> > yr_cr2;
        std::valarray<std::valarray<REAL_T> > yr_disc2;
        std::vector<REAL_T> yr_mnwt2;
        std::valarray<std::valarray<REAL_T> > have_data;
        std::vector<REAL_T> obs_equ_catch;
        int N_ReadCatch;
        std::valarray<std::valarray<REAL_T> > catch_bioT;
        std::valarray<std::valarray<REAL_T> > catch_ret_obs;
        std::vector<int> have_catch;
        std::vector<REAL_T> catch_seas_area;
        std::valarray<std::valarray<REAL_T> > totcatch_byarea;
        std::vector<REAL_T> totcat;
        int first_catch_yr;
        int nobs_cr_rd;
        int nobs_cr;
        std::valarray<std::valarray<int> > cr_units_rd;
        std::vector<int> cr_units;
        std::vector<int> cr_errtype;
        std::valarray<std::valarray<REAL_T> > indexdata;
        std::vector<int> nyr_cr;
        int in_superperiod;
        std::vector<int> N_suprper_cr;
        std::valarray<std::valarray<int> > yr_cr;
        std::valarray<std::valarray<int> > yr_cr_y;
        std::valarray<std::valarray<int> > yr_cr_s;
        std::valarray<std::valarray<int> > yr_cr_use;
        std::valarray<std::valarray<REAL_T> > obs_cr;
        std::valarray<std::valarray<REAL_T> > Ln_obs_cr;
        std::valarray<std::valarray<REAL_T> > se_cr_obs;
        std::valarray<std::valarray<REAL_T> > vul_bio;
        std::valarray<std::valarray<int> > suprper_cr1;
        std::valarray<std::valarray<int> > suprper_cr2;
        int Ndisc_fleets;
        int nobs_disc;
        std::valarray<std::valarray<int> > disc_units_rd;
        std::vector<int> disc_units;
        std::vector<int> disc_errtype;
        std::vector<REAL_T> disc_errtype_r;
        int nobs_disc_rd;
        std::valarray<std::valarray<REAL_T> > discdata;
        std::vector<int> nyr_disc;
        std::vector<int> N_suprper_disc;
        std::valarray<std::valarray<int> > yr_disc;
        std::valarray<std::valarray<int> > yr_disc_y;
        std::valarray<std::valarray<int> > yr_disc_s;
        std::valarray<std::valarray<int> > yr_disc_use;
        std::valarray<std::valarray<REAL_T> > obs_disc;
        std::valarray<std::valarray<REAL_T> > cv_disc;
        std::valarray<std::valarray<REAL_T> > sd_disc;
        std::valarray<std::valarray<int> > suprper_disc1;
        std::valarray<std::valarray<int> > suprper_disc2;
        int nobs_mnwt_rd;
        int nobs_mnwt;
        REAL_T DF_bodywt;
        std::valarray<std::valarray<REAL_T> > mnwtdata1;
        std::valarray<std::valarray<REAL_T> > mnwtdata;
        REAL_T binwidth2;
        REAL_T minLread;
        REAL_T maxLread;
        int nlen_bin2;
        int nlen_binP;
        REAL_T minL;
        REAL_T minL_m;
        REAL_T maxL;
        int nlength;
        int nlength1;
        int nlength2;
        REAL_T startbin;
        int LenBin_option;
        std::vector<REAL_T> PopBin_Read;
        std::vector<REAL_T> len_bins_rd;
        REAL_T min_tail;
        REAL_T min_comp;
        int CombGender_l;
        int nlen_bin;
        std::vector<REAL_T> len_bins_dat;
        std::vector<REAL_T> len_bins_dat2;
        std::vector<REAL_T> len_bins_dat_m;
        std::vector<REAL_T> len_bins;
        std::vector<REAL_T> len_bins2;
        std::vector<REAL_T> binwidth;
        std::vector<REAL_T> len_bins_m;
        std::vector<REAL_T> len_bins_m2;
        std::vector<REAL_T> len_bins_sq;
        std::vector<REAL_T> male_offset;
        std::valarray<std::valarray<REAL_T> > make_len_bin;
        int ibin;
        int ibinsave;
        int fini;
        REAL_T topbin;
        REAL_T botbin;
        int nobsl_rd;
        int Nobs_l_tot;
        std::valarray<std::valarray<REAL_T> > lendata;
        std::vector<int> Nobs_l;
        std::vector<int> N_suprper_l;
        std::valarray<std::valarray<int> > yr_l_t;
        std::valarray<std::valarray<int> > yr_l_y;
        std::vector<REAL_T> obs_l;
        std::valarray<std::valarray<REAL_T> > obs_l_all;
        std::valarray<std::valarray<REAL_T> > nsamp_l;
        std::valarray<std::valarray<REAL_T> > nsamp_l_read;
        std::valarray<std::valarray<int> > gen_l;
        std::valarray<std::valarray<int> > mkt_l;
        std::vector<REAL_T> header_l;
        std::vector<REAL_T> tails_l;
        std::vector<int> tails_w;
        std::valarray<std::valarray<int> > suprper_l1;
        std::valarray<std::valarray<int> > suprper_l2;
        int floop;
        int tloop;
        int n_abins;
        int n_abins1;
        int n_abins2;
        int Use_AgeKeyZero;
        int AgeKeyParm;
        std::vector<REAL_T> age_bins1;
        int N_ageerr;
        std::vector<REAL_T> age_err_rd;
        std::vector<REAL_T> age_bins;
        int nobsa_rd;
        int Nobs_a_tot;
        int Lbin_method;
        int CombGender_a;
        std::valarray<std::valarray<REAL_T> > agedata;
        std::vector<int> Nobs_a;
        std::vector<int> N_suprper_a;
        std::valarray<std::valarray<int> > yr_a_t;
        std::valarray<std::valarray<int> > yr_a_y;
        std::vector<REAL_T> obs_a;
        std::valarray<std::valarray<REAL_T> > nsamp_a;
        std::valarray<std::valarray<REAL_T> > nsamp_a_read;
        std::valarray<std::valarray<int> > ageerr_type_a;
        std::valarray<std::valarray<int> > gen_a;
        std::valarray<std::valarray<int> > mkt_a;
        std::vector<REAL_T> Lbin_filter;
        std::valarray<std::valarray<int> > use_Lbin_filter;
        std::valarray<std::valarray<int> > Lbin_lo;
        std::valarray<std::valarray<int> > Lbin_hi;
        std::vector<REAL_T> tails_a;
        std::vector<REAL_T> header_a;
        std::valarray<std::valarray<int> > suprper_a1;
        std::valarray<std::valarray<int> > suprper_a2;
        int nobs_ms_rd;
        int nobs_ms_tot;
        std::valarray<std::valarray<REAL_T> > sizeagedata;
        std::vector<int> Nobs_ms;
        std::vector<int> N_suprper_ms;
        std::valarray<std::valarray<int> > yr_ms_t;
        std::valarray<std::valarray<int> > yr_ms_y;
        std::vector<REAL_T> obs_ms;
        std::vector<REAL_T> obs_ms_n;
        std::vector<REAL_T> obs_ms_n_read;
        std::valarray<std::valarray<int> > ageerr_type_ms;
        std::valarray<std::valarray<int> > gen_ms;
        std::valarray<std::valarray<int> > mkt_ms;
        std::valarray<std::valarray<int> > use_ms;
        std::vector<REAL_T> header_ms;
        std::valarray<std::valarray<int> > suprper_ms1;
        std::valarray<std::valarray<int> > suprper_ms2;
        int N_envvar;
        int N_envdata;
        std::valarray<std::valarray<REAL_T> > env_data_RD;
        std::valarray<std::valarray<REAL_T> > env_temp;
        int SzFreqMethod;
        int iobs;
        int SzFreq_Nmeth;
        std::vector<REAL_T> SzFreq_HaveObs;
        std::valarray<std::valarray<int> > SzFreq_HaveObs2;
        std::vector<int> SzFreq_Nbins;
        std::vector<int> SzFreq_units;
        std::vector<int> SzFreq_scale;
        std::vector<REAL_T> SzFreq_mincomp;
        std::vector<int> SzFreq_nobs;
        std::vector<int> SzFreq_Nbins_seas_g;
        std::vector<int> SzFreq_Nbins3;
        int SzFreqMethod_seas;
        std::valarray<std::valarray<REAL_T> > SzFreq_bins1;
        std::valarray<std::valarray<REAL_T> > SzFreq_bins;
        std::valarray<std::valarray<REAL_T> > SzFreq_bins2;
        std::vector<int> SzFreq_Omit_Small;
        int SzFreq_totobs;
        int SzFreq_N_Like;
        std::vector<int> SzFreq_Setup;
        std::vector<int> SzFreq_Setup2;
        std::valarray<std::valarray<REAL_T> > SzFreq_obs1;
        std::valarray<std::valarray<int> > SzFreq_obs_hdr;
        std::vector<REAL_T> SzFreq_sampleN;
        std::vector<REAL_T> SzFreq_effN;
        std::vector<REAL_T> SzFreq_eachlike;
        std::valarray<std::valarray<REAL_T> > SzFreq_obs;
        std::valarray<std::valarray<int> > SzFreq_LikeComponent;
        std::vector<int> N_suprper_SzFreq;
        std::vector<REAL_T> SzFreq_like_base;
        std::valarray<std::valarray<int> > suprper_start_SzFreq;
        std::valarray<std::valarray<int> > suprper_end_SzFreq;
        int Do_TG;
        std::vector<REAL_T> TG_temp;
        int TG;
        int N_TG;
        int N_TG2;
        int TG_timestart;
        int N_TG_recap;
        int TG_mixperiod;
        int TG_maxperiods;
        std::vector<int> TG_endtime;
        std::valarray<std::valarray<REAL_T> > TG_release;
        std::valarray<std::valarray<REAL_T> > TG_recap_data;
        std::vector<REAL_T> TG_recap_obs;
        int Do_Morphcomp;
        std::vector<REAL_T> mc_temp;
        int Morphcomp_nobs;
        int Morphcomp_nmorph;
        REAL_T Morphcomp_mincomp;
        std::valarray<std::valarray<REAL_T> > Morphcomp_obs;
        std::vector<REAL_T> Morphcomp_havedata;
        int fid;
        int Do_Benchmark;
        int Do_MSY;
        int did_MSY;
        int show_MSY;
        REAL_T SPR_target;
        REAL_T BTGT_target;
        std::vector<int> Bmark_Yr;
        std::vector<int> Bmark_t;
        std::vector<int> Bmark_Yr_rd;
        int Bmark_RelF_Basis;
        int Do_Forecast;
        std::vector<REAL_T> Fcast_Input;
        int N_Fcast_Yrs;
        std::vector<int> Fcast_yr;
        int Fcast_Sel_yr1;
        int Fcast_Sel_yr2;
        int Fcast_RelF_yr1;
        int Fcast_RelF_yr2;
        int Fcast_RelF_Basis;
        REAL_T Fcast_Flevel;
        int Do_Rebuilder;
        int Rebuild_Ydecl;
        int Rebuild_Yinit;
        int HarvestPolicy;
        REAL_T H4010_top;
        REAL_T H4010_bot;
        REAL_T H4010_scale;
        int Do_Impl_Error;
        REAL_T Impl_Error_Std;
        std::vector<int> Fcast_Loop_Control;
        int N_Fcast_catches;
        int Fcast_InputCatch_Basis;
        int Fcast_Catch_Basis;
        int Fcast_Catch_Allocation_Groups;
        int Fcast_Do_Fleet_Cap;
        int Fcast_Do_Area_Cap;
        int Fcast_Cap_FirstYear;
        std::valarray<std::valarray<REAL_T> > Fcast_RelF;
        std::vector<REAL_T> Fcast_Input_rd;
        std::valarray<std::valarray<REAL_T> > Fcast_RelF_Input;
        std::vector<REAL_T> Fcast_MaxFleetCatch;
        std::vector<REAL_T> Max_Fcast_Catch;
        std::vector<int> Allocation_Fleet_Assignments;
        std::vector<REAL_T> Fcast_Catch_Allocation;
        std::vector<int> more_Fcast_input;
        std::valarray<std::valarray<REAL_T> > Fcast_InputCatch;
        std::valarray<std::valarray<REAL_T> > Fcast_InputCatch_1;
        REAL_T fif;
        std::vector<int> STD_Yr_Reverse;
        std::vector<int> STD_Yr_Reverse_Dep;
        std::vector<int> STD_Yr_Reverse_Ofish;
        std::vector<int> STD_Yr_Reverse_F;
        int N_STD_Yr_Dep;
        int N_STD_Yr_Ofish;
        int N_STD_Yr_F;
        int N_STD_Mgmt_Quant;
        int N_GP;
        int N_GP_sub;
        REAL_T sd_ratio;
        REAL_T sd_withinmorph;
        REAL_T sd_betweenmorph;
        std::vector<int> ishadow;
        std::vector<REAL_T> shadow;
        std::vector<REAL_T> submorphdist;
        int gmorph;
        int gp;
        int gp2;
        int birthseas;
        std::vector<int> sx;
        std::vector<int> GP;
        std::vector<int> GP3;
        std::vector<int> GP4;
        std::vector<int> Bseas;
        std::vector<int> GP2;
        std::vector<int> use_morph;
        std::valarray<std::valarray<int> > TG_use_morph;
        std::valarray<std::valarray<int> > recr_dist_pattern_2;
        std::vector<REAL_T> azero_seas;
        std::vector<REAL_T> azero_G;
        std::vector<REAL_T> curr_age1;
        int first_grow_age;
        std::vector<REAL_T> curr_age2;
        std::vector<REAL_T> lin_grow1;
        std::vector<REAL_T> lin_grow2;
        std::valarray<std::valarray<int> > first_grow_age2;
        std::vector<REAL_T> recr_dist_pattern;
        int recr_dist_inx;
        int recr_dist_read;
        std::vector<int> recr_dist_input;
        std::valarray<std::valarray<int> > recr_dist_pattern_1;
        std::valarray<std::valarray<int> > GP_finder;
        int do_migration;
        REAL_T migr_firstage;
        std::valarray<std::valarray<REAL_T> > migr_start;
        std::valarray<std::valarray<REAL_T> > move_def;
        std::vector<REAL_T> move_pattern;
        int do_migr2;
        std::vector<int> firstBseas;
        std::valarray<std::valarray<REAL_T> > move_def2;
        int k1;
        int k2;
        int k3;
        int N_Block_Designs;
        std::vector<int> Nblk;
        std::vector<int> Nblk2;
        std::valarray<std::valarray<int> > Block_Design;
        int N_MGparm;
        int N_natMparms;
        int N_growparms;
        int MGparm_per_def;
        int recr_dist_parms;
        REAL_T natM_amin;
        REAL_T natM_amax;
        REAL_T fracfemale;
        int natM_type;
        std::vector<REAL_T> tempvec4;
        std::vector<REAL_T> NatM_break;
        std::valarray<std::valarray<REAL_T> > Age_NatMort;
        int Grow_type;
        REAL_T AFIX;
        REAL_T AFIX2;
        REAL_T AFIX2_forCV;
        REAL_T AFIX_delta;
        REAL_T AFIX_plus;
        std::vector<REAL_T> tempvec5;
        std::valarray<std::valarray<REAL_T> > Len_At_Age_rd;
        REAL_T SD_add_to_LAA;
        int CV_depvar;
        int CV_depvar_a;
        int CV_depvar_b;
        int Maturity_Option;
        int WTage_rd;
        std::valarray<std::valarray<REAL_T> > Age_Maturity;
        int First_Mature_Age;
        int Fecund_Option;
        int Hermaphro_Option;
        int MGparm_Hermaphro;
        std::vector<int> Hermaphro_more;
        int Hermaphro_seas;
        int Hermaphro_maleSPB;
        int MGparm_def;
        int MG_adjust_method;
        std::valarray<std::valarray<int> > time_vary_MG;
        std::vector<int> MG_active;
        int do_once;
        int doit;
        std::vector<REAL_T> femfrac;
        int MGP_CGD;
        int CGD;
        std::valarray<std::valarray<REAL_T> > MGparm_1;
        std::vector<int> MGparm_offset;
        int N_MGparm_env;
        int customMGenvsetup;
        std::vector<int> MGparm_env;
        std::vector<int> MGparm_envuse;
        std::vector<int> MGparm_envtype;
        std::vector<int> mgp_type;
        std::valarray<std::valarray<REAL_T> > MGparm_env_1;
        int N_MGparm_blk;
        std::valarray<std::valarray<int> > Block_Defs_MG;
        int N_MGparm_trend;
        int N_MGparm_trend2;
        std::vector<int> MGparm_trend_point;
        int customblocksetup_MG;
        std::valarray<std::valarray<REAL_T> > MGparm_blk_1;
        std::valarray<std::valarray<REAL_T> > MGparm_trend_1;
        std::vector<int> MGparm_trend_rev;
        std::vector<int> MGparm_trend_rev_1;
        std::vector<int> MGparm_seas_effects;
        int MGparm_doseas;
        int N_MGparm_seas;
        std::valarray<std::valarray<REAL_T> > MGparm_seas_1;
        int N_MGparm2;
        std::vector<REAL_T> MGparm_LO;
        std::vector<REAL_T> MGparm_HI;
        std::vector<REAL_T> MGparm_RD;
        std::vector<REAL_T> MGparm_PR;
        std::vector<int> MGparm_PRtype;
        std::vector<REAL_T> MGparm_CV;
        std::vector<int> MGparm_PH;
        int N_MGparm_dev;
        std::vector<int> MGparm_dev_minyr;
        std::vector<int> MGparm_dev_maxyr;
        std::vector<REAL_T> MGparm_dev_stddev;
        std::vector<int> MGparm_dev_type;
        std::vector<int> MGparm_dev_select;
        int MGparm_dev_PH;
        int SR_fxn;
        int N_SRparm;
        int N_SRparm2;
        std::valarray<std::valarray<REAL_T> > SR_parm_1;
        int SR_env_link;
        int SR_env_target_RD;
        int SR_env_target;
        int SR_autocorr;
        std::vector<REAL_T> SRvec_LO;
        std::vector<REAL_T> SRvec_HI;
        std::vector<int> SRvec_PH;
        int do_recdev;
        int recdev_start;
        int recdev_end;
        int recdev_PH_rd;
        int recdev_PH;
        int recdev_adv;
        std::vector<REAL_T> recdev_options_rd;
        std::vector<REAL_T> recdev_options;
        int recdev_early_start_rd;
        int recdev_early_start;
        int recdev_early_end;
        int recdev_first;
        int recdev_early_PH;
        int Fcast_recr_PH;
        int Fcast_recr_PH2;
        REAL_T Fcast_recr_lambda;
        std::vector<REAL_T> recdev_adj;
        int recdev_cycle;
        int recdev_do_early;
        int recdev_read;
        REAL_T recdev_LO;
        REAL_T recdev_HI;
        std::vector<int> recdev_doit;
        std::valarray<std::valarray<REAL_T> > recdev_cycle_parm_RD;
        std::vector<REAL_T> recdev_cycle_LO;
        std::vector<REAL_T> recdev_cycle_HI;
        std::vector<int> recdev_cycle_PH;
        std::valarray<std::valarray<REAL_T> > recdev_input;
        REAL_T F_ballpark;
        int F_ballpark_yr;
        int F_Method;
        int F_Method_use;
        REAL_T max_harvest_rate;
        REAL_T Equ_F_joiner;
        std::vector<REAL_T> F_setup;
        int F_detail;
        int F_Tune;
        std::valarray<std::valarray<REAL_T> > F_setup2;
        std::valarray<std::valarray<REAL_T> > init_F_parm_1;
        std::vector<REAL_T> init_F_LO;
        std::vector<REAL_T> init_F_HI;
        std::vector<REAL_T> init_F_RD;
        std::vector<REAL_T> init_F_PR;
        std::vector<REAL_T> init_F_PRtype;
        std::vector<REAL_T> init_F_CV;
        std::vector<int> init_F_PH;
        std::valarray<std::valarray<REAL_T> > Q_setup;
        int Q_Npar2;
        int Q_Npar;
        int ask_detail;
        std::valarray<std::valarray<int> > Q_setup_parms;
        std::vector<REAL_T> Q_parm_LO;
        std::vector<REAL_T> Q_parm_HI;
        std::vector<int> Q_parm_PH;
        int Q_parm_detail;
        std::valarray<std::valarray<REAL_T> > Q_parm_1;
        std::valarray<std::valarray<REAL_T> > Q_parm_2;
        std::vector<int> seltype_Nparam;
        std::valarray<std::valarray<int> > seltype;
        int N_selparm;
        int N_selparm2;
        std::vector<int> N_selparmvec;
        std::vector<int> Maleselparm;
        std::vector<int> RetainParm;
        std::vector<int> dolen;
        int blkparm;
        int firstselparm;
        int Do_Retain;
        std::valarray<std::valarray<REAL_T> > selparm_1;
        std::valarray<std::valarray<int> > time_vary_sel;
        std::valarray<std::valarray<int> > time_vary_makefishsel;
        int makefishsel_yr;
        int ALK_yr;
        int N_selparm_env;
        int customenvsetup;
        std::vector<int> selparm_env;
        std::vector<int> selparm_envuse;
        std::vector<int> selparm_envtype;
        std::valarray<std::valarray<REAL_T> > selparm_env_1;
        int N_selparm_blk;
        std::valarray<std::valarray<int> > Block_Defs_Sel;
        int customblocksetup;
        int N_selparm_trend;
        int N_selparm_trend2;
        std::vector<int> selparm_trend_point;
        std::valarray<std::valarray<REAL_T> > selparm_blk_1;
        std::valarray<std::valarray<REAL_T> > selparm_trend_1;
        std::vector<int> selparm_trend_rev;
        std::vector<int> selparm_trend_rev_1;
        int N_selparm_dev;
        int N_selparm_dev_tot;
        std::vector<int> selparm_dev_minyr;
        std::vector<int> selparm_dev_maxyr;
        std::vector<REAL_T> selparm_dev_stddev;
        std::vector<int> selparm_dev_type;
        std::vector<int> selparm_dev_select;
        REAL_T selparm_dev_PH;
        int selparm_adjust_method;
        std::vector<REAL_T> selparm_LO;
        std::vector<REAL_T> selparm_HI;
        std::vector<REAL_T> selparm_RD;
        std::vector<REAL_T> selparm_PR;
        std::vector<REAL_T> selparm_PRtype;
        std::vector<REAL_T> selparm_CV;
        std::vector<int> selparm_PH;
        int TG_custom;
        std::valarray<std::valarray<REAL_T> > TG_parm1;
        std::valarray<std::valarray<REAL_T> > TG_parm2;
        std::vector<REAL_T> TG_parm_LO;
        std::vector<REAL_T> TG_parm_HI;
        std::vector<int> TG_parm_PH;
        int Do_Var_adjust;
        std::valarray<std::valarray<REAL_T> > var_adjust1;
        std::valarray<std::valarray<REAL_T> > var_adjust;
        REAL_T max_lambda_phase;
        REAL_T sd_offset;
        std::valarray<std::valarray<REAL_T> > surv_lambda;
        std::valarray<std::valarray<REAL_T> > disc_lambda;
        std::valarray<std::valarray<REAL_T> > mnwt_lambda;
        std::valarray<std::valarray<REAL_T> > length_lambda;
        std::valarray<std::valarray<REAL_T> > age_lambda;
        std::valarray<std::valarray<REAL_T> > sizeage_lambda;
        std::vector<REAL_T> init_equ_lambda;
        std::valarray<std::valarray<REAL_T> > catch_lambda;
        std::vector<REAL_T> recrdev_lambda;
        std::vector<REAL_T> parm_prior_lambda;
        std::vector<REAL_T> parm_dev_lambda;
        std::vector<REAL_T> CrashPen_lambda;
        std::vector<REAL_T> Morphcomp_lambda;
        std::valarray<std::valarray<REAL_T> > SzFreq_lambda;
        std::valarray<std::valarray<REAL_T> > TG_lambda1;
        std::valarray<std::valarray<REAL_T> > TG_lambda2;
        int N_lambda_changes;
        std::valarray<std::valarray<REAL_T> > Lambda_changes;
        int Do_More_Std;
        std::vector<int> More_Std_Input;
        int Do_Selex_Std;
        int Selex_Std_AL;
        int Selex_Std_Year;
        int Selex_Std_Cnt;
        int Do_Growth_Std;
        int Growth_Std_Cnt;
        int Do_NatAge_Std;
        int NatAge_Std_Year;
        int NatAge_Std_Cnt;
        int Extra_Std_N;
        std::vector<int> Selex_Std_Pick;
        std::vector<int> Growth_Std_Pick;
        std::vector<int> NatAge_Std_Pick;
        int fim;
        int N_WTage_rd;
        int N_WTage_maxage;
        int y2;
        std::vector<REAL_T> junkvec;
        std::valarray<std::valarray<REAL_T> > WTage_in;
        std::vector<REAL_T> junkvec2;
        std::vector<REAL_T> WTage_emp;
        int CoVar_Count;
        int active_count;
        int active_parms;
        std::vector<int> active_parm;
        int max_phase;
        int report_phase;
        std::valarray<std::valarray<REAL_T> > save_G_parm;
        std::valarray<std::valarray<REAL_T> > save_seas_parm;
        std::valarray<std::valarray<REAL_T> > CoVar;


    public:

        int GetA2() const {
            return A2;
        }

        void SetA2(int A2) {
            this->A2 = A2;
        }

        REAL_T GetAFIX() const {
            return AFIX;
        }

        void SetAFIX(REAL_T AFIX) {
            this->AFIX = AFIX;
        }

        REAL_T GetAFIX2() const {
            return AFIX2;
        }

        void SetAFIX2(REAL_T AFIX2) {
            this->AFIX2 = AFIX2;
        }

        REAL_T GetAFIX2_forCV() const {
            return AFIX2_forCV;
        }

        void SetAFIX2_forCV(REAL_T AFIX2_forCV) {
            this->AFIX2_forCV = AFIX2_forCV;
        }

        REAL_T GetAFIX_delta() const {
            return AFIX_delta;
        }

        void SetAFIX_delta(REAL_T AFIX_delta) {
            this->AFIX_delta = AFIX_delta;
        }

        REAL_T GetAFIX_plus() const {
            return AFIX_plus;
        }

        void SetAFIX_plus(REAL_T AFIX_plus) {
            this->AFIX_plus = AFIX_plus;
        }

        int GetALK_yr() const {
            return ALK_yr;
        }

        void SetALK_yr(int ALK_yr) {
            this->ALK_yr = ALK_yr;
        }

        int GetAgeKeyParm() const {
            return AgeKeyParm;
        }

        void SetAgeKeyParm(int AgeKeyParm) {
            this->AgeKeyParm = AgeKeyParm;
        }

        std::valarray<std::valarray<REAL_T> > GetAge_Maturity() const {
            return Age_Maturity;
        }

        void SetAge_Maturity(std::valarray<std::valarray<REAL_T> > Age_Maturity) {
            this->Age_Maturity = Age_Maturity;
        }

        std::valarray<std::valarray<REAL_T> > GetAge_NatMort() const {
            return Age_NatMort;
        }

        void SetAge_NatMort(std::valarray<std::valarray<REAL_T> > Age_NatMort) {
            this->Age_NatMort = Age_NatMort;
        }

        std::vector<int> GetAllocation_Fleet_Assignments() const {
            return Allocation_Fleet_Assignments;
        }

        void SetAllocation_Fleet_Assignments(std::vector<int> Allocation_Fleet_Assignments) {
            this->Allocation_Fleet_Assignments = Allocation_Fleet_Assignments;
        }

        REAL_T GetBTGT_target() const {
            return BTGT_target;
        }

        void SetBTGT_target(REAL_T BTGT_target) {
            this->BTGT_target = BTGT_target;
        }

        std::valarray<std::valarray<int> > GetBlock_Defs_MG() const {
            return Block_Defs_MG;
        }

        void SetBlock_Defs_MG(std::valarray<std::valarray<int> > Block_Defs_MG) {
            this->Block_Defs_MG = Block_Defs_MG;
        }

        std::valarray<std::valarray<int> > GetBlock_Defs_Sel() const {
            return Block_Defs_Sel;
        }

        void SetBlock_Defs_Sel(std::valarray<std::valarray<int> > Block_Defs_Sel) {
            this->Block_Defs_Sel = Block_Defs_Sel;
        }

        std::valarray<std::valarray<int> > GetBlock_Design() const {
            return Block_Design;
        }

        void SetBlock_Design(std::valarray<std::valarray<int> > Block_Design) {
            this->Block_Design = Block_Design;
        }

        int GetBmark_RelF_Basis() const {
            return Bmark_RelF_Basis;
        }

        void SetBmark_RelF_Basis(int Bmark_RelF_Basis) {
            this->Bmark_RelF_Basis = Bmark_RelF_Basis;
        }

        std::vector<int> GetBmark_Yr() const {
            return Bmark_Yr;
        }

        void SetBmark_Yr(std::vector<int> Bmark_Yr) {
            this->Bmark_Yr = Bmark_Yr;
        }

        std::vector<int> GetBmark_Yr_rd() const {
            return Bmark_Yr_rd;
        }

        void SetBmark_Yr_rd(std::vector<int> Bmark_Yr_rd) {
            this->Bmark_Yr_rd = Bmark_Yr_rd;
        }

        std::vector<int> GetBmark_t() const {
            return Bmark_t;
        }

        void SetBmark_t(std::vector<int> Bmark_t) {
            this->Bmark_t = Bmark_t;
        }

        std::vector<int> GetBseas() const {
            return Bseas;
        }

        void SetBseas(std::vector<int> Bseas) {
            this->Bseas = Bseas;
        }

        int GetCGD() const {
            return CGD;
        }

        void SetCGD(int CGD) {
            this->CGD = CGD;
        }

        int GetCV_depvar() const {
            return CV_depvar;
        }

        void SetCV_depvar(int CV_depvar) {
            this->CV_depvar = CV_depvar;
        }

        int GetCV_depvar_a() const {
            return CV_depvar_a;
        }

        void SetCV_depvar_a(int CV_depvar_a) {
            this->CV_depvar_a = CV_depvar_a;
        }

        int GetCV_depvar_b() const {
            return CV_depvar_b;
        }

        void SetCV_depvar_b(int CV_depvar_b) {
            this->CV_depvar_b = CV_depvar_b;
        }

        std::valarray<std::valarray<REAL_T> > GetCoVar() const {
            return CoVar;
        }

        void SetCoVar(std::valarray<std::valarray<REAL_T> > CoVar) {
            this->CoVar = CoVar;
        }

        int GetCoVar_Count() const {
            return CoVar_Count;
        }

        void SetCoVar_Count(int CoVar_Count) {
            this->CoVar_Count = CoVar_Count;
        }

        int GetCombGender_a() const {
            return CombGender_a;
        }

        void SetCombGender_a(int CombGender_a) {
            this->CombGender_a = CombGender_a;
        }

        int GetCombGender_l() const {
            return CombGender_l;
        }

        void SetCombGender_l(int CombGender_l) {
            this->CombGender_l = CombGender_l;
        }

        std::vector<REAL_T> GetCrashPen_lambda() const {
            return CrashPen_lambda;
        }

        void SetCrashPen_lambda(std::vector<REAL_T> CrashPen_lambda) {
            this->CrashPen_lambda = CrashPen_lambda;
        }

        REAL_T GetDF_bodywt() const {
            return DF_bodywt;
        }

        void SetDF_bodywt(REAL_T DF_bodywt) {
            this->DF_bodywt = DF_bodywt;
        }

        int GetDo_Benchmark() const {
            return Do_Benchmark;
        }

        void SetDo_Benchmark(int Do_Benchmark) {
            this->Do_Benchmark = Do_Benchmark;
        }

        int GetDo_CumReport() const {
            return Do_CumReport;
        }

        void SetDo_CumReport(int Do_CumReport) {
            this->Do_CumReport = Do_CumReport;
        }

        int GetDo_Forecast() const {
            return Do_Forecast;
        }

        void SetDo_Forecast(int Do_Forecast) {
            this->Do_Forecast = Do_Forecast;
        }

        int GetDo_Growth_Std() const {
            return Do_Growth_Std;
        }

        void SetDo_Growth_Std(int Do_Growth_Std) {
            this->Do_Growth_Std = Do_Growth_Std;
        }

        int GetDo_Impl_Error() const {
            return Do_Impl_Error;
        }

        void SetDo_Impl_Error(int Do_Impl_Error) {
            this->Do_Impl_Error = Do_Impl_Error;
        }

        int GetDo_MSY() const {
            return Do_MSY;
        }

        void SetDo_MSY(int Do_MSY) {
            this->Do_MSY = Do_MSY;
        }

        int GetDo_More_Std() const {
            return Do_More_Std;
        }

        void SetDo_More_Std(int Do_More_Std) {
            this->Do_More_Std = Do_More_Std;
        }

        int GetDo_Morphcomp() const {
            return Do_Morphcomp;
        }

        void SetDo_Morphcomp(int Do_Morphcomp) {
            this->Do_Morphcomp = Do_Morphcomp;
        }

        int GetDo_NatAge_Std() const {
            return Do_NatAge_Std;
        }

        void SetDo_NatAge_Std(int Do_NatAge_Std) {
            this->Do_NatAge_Std = Do_NatAge_Std;
        }

        int GetDo_ParmTrace() const {
            return Do_ParmTrace;
        }

        void SetDo_ParmTrace(int Do_ParmTrace) {
            this->Do_ParmTrace = Do_ParmTrace;
        }

        int GetDo_Rebuilder() const {
            return Do_Rebuilder;
        }

        void SetDo_Rebuilder(int Do_Rebuilder) {
            this->Do_Rebuilder = Do_Rebuilder;
        }

        int GetDo_Retain() const {
            return Do_Retain;
        }

        void SetDo_Retain(int Do_Retain) {
            this->Do_Retain = Do_Retain;
        }

        int GetDo_Selex_Std() const {
            return Do_Selex_Std;
        }

        void SetDo_Selex_Std(int Do_Selex_Std) {
            this->Do_Selex_Std = Do_Selex_Std;
        }

        int GetDo_TG() const {
            return Do_TG;
        }

        void SetDo_TG(int Do_TG) {
            this->Do_TG = Do_TG;
        }

        int GetDo_Var_adjust() const {
            return Do_Var_adjust;
        }

        void SetDo_Var_adjust(int Do_Var_adjust) {
            this->Do_Var_adjust = Do_Var_adjust;
        }

        int GetDo_all_priors() const {
            return Do_all_priors;
        }

        void SetDo_all_priors(int Do_all_priors) {
            this->Do_all_priors = Do_all_priors;
        }

        REAL_T GetEqu_F_joiner() const {
            return Equ_F_joiner;
        }

        void SetEqu_F_joiner(REAL_T Equ_F_joiner) {
            this->Equ_F_joiner = Equ_F_joiner;
        }

        int GetExtra_Std_N() const {
            return Extra_Std_N;
        }

        void SetExtra_Std_N(int Extra_Std_N) {
            this->Extra_Std_N = Extra_Std_N;
        }

        int GetF_Method() const {
            return F_Method;
        }

        void SetF_Method(int F_Method) {
            this->F_Method = F_Method;
        }

        int GetF_Method_use() const {
            return F_Method_use;
        }

        void SetF_Method_use(int F_Method_use) {
            this->F_Method_use = F_Method_use;
        }

        int GetF_Tune() const {
            return F_Tune;
        }

        void SetF_Tune(int F_Tune) {
            this->F_Tune = F_Tune;
        }

        REAL_T GetF_ballpark() const {
            return F_ballpark;
        }

        void SetF_ballpark(REAL_T F_ballpark) {
            this->F_ballpark = F_ballpark;
        }

        int GetF_ballpark_yr() const {
            return F_ballpark_yr;
        }

        void SetF_ballpark_yr(int F_ballpark_yr) {
            this->F_ballpark_yr = F_ballpark_yr;
        }

        int GetF_detail() const {
            return F_detail;
        }

        void SetF_detail(int F_detail) {
            this->F_detail = F_detail;
        }

        int GetF_report_basis() const {
            return F_report_basis;
        }

        void SetF_report_basis(int F_report_basis) {
            this->F_report_basis = F_report_basis;
        }

        int GetF_reporting() const {
            return F_reporting;
        }

        void SetF_reporting(int F_reporting) {
            this->F_reporting = F_reporting;
        }

        std::vector<int> GetF_reporting_ages() const {
            return F_reporting_ages;
        }

        void SetF_reporting_ages(std::vector<int> F_reporting_ages) {
            this->F_reporting_ages = F_reporting_ages;
        }

        std::vector<REAL_T> GetF_setup() const {
            return F_setup;
        }

        void SetF_setup(std::vector<REAL_T> F_setup) {
            this->F_setup = F_setup;
        }

        std::valarray<std::valarray<REAL_T> > GetF_setup2() const {
            return F_setup2;
        }

        void SetF_setup2(std::valarray<std::valarray<REAL_T> > F_setup2) {
            this->F_setup2 = F_setup2;
        }

        int GetFcast_Cap_FirstYear() const {
            return Fcast_Cap_FirstYear;
        }

        void SetFcast_Cap_FirstYear(int Fcast_Cap_FirstYear) {
            this->Fcast_Cap_FirstYear = Fcast_Cap_FirstYear;
        }

        std::vector<REAL_T> GetFcast_Catch_Allocation() const {
            return Fcast_Catch_Allocation;
        }

        void SetFcast_Catch_Allocation(std::vector<REAL_T> Fcast_Catch_Allocation) {
            this->Fcast_Catch_Allocation = Fcast_Catch_Allocation;
        }

        int GetFcast_Catch_Allocation_Groups() const {
            return Fcast_Catch_Allocation_Groups;
        }

        void SetFcast_Catch_Allocation_Groups(int Fcast_Catch_Allocation_Groups) {
            this->Fcast_Catch_Allocation_Groups = Fcast_Catch_Allocation_Groups;
        }

        int GetFcast_Catch_Basis() const {
            return Fcast_Catch_Basis;
        }

        void SetFcast_Catch_Basis(int Fcast_Catch_Basis) {
            this->Fcast_Catch_Basis = Fcast_Catch_Basis;
        }

        int GetFcast_Do_Area_Cap() const {
            return Fcast_Do_Area_Cap;
        }

        void SetFcast_Do_Area_Cap(int Fcast_Do_Area_Cap) {
            this->Fcast_Do_Area_Cap = Fcast_Do_Area_Cap;
        }

        int GetFcast_Do_Fleet_Cap() const {
            return Fcast_Do_Fleet_Cap;
        }

        void SetFcast_Do_Fleet_Cap(int Fcast_Do_Fleet_Cap) {
            this->Fcast_Do_Fleet_Cap = Fcast_Do_Fleet_Cap;
        }

        REAL_T GetFcast_Flevel() const {
            return Fcast_Flevel;
        }

        void SetFcast_Flevel(REAL_T Fcast_Flevel) {
            this->Fcast_Flevel = Fcast_Flevel;
        }

        std::vector<REAL_T> GetFcast_Input() const {
            return Fcast_Input;
        }

        void SetFcast_Input(std::vector<REAL_T> Fcast_Input) {
            this->Fcast_Input = Fcast_Input;
        }

        std::valarray<std::valarray<REAL_T> > GetFcast_InputCatch() const {
            return Fcast_InputCatch;
        }

        void SetFcast_InputCatch(std::valarray<std::valarray<REAL_T> > Fcast_InputCatch) {
            this->Fcast_InputCatch = Fcast_InputCatch;
        }

        std::valarray<std::valarray<REAL_T> > GetFcast_InputCatch_1() const {
            return Fcast_InputCatch_1;
        }

        void SetFcast_InputCatch_1(std::valarray<std::valarray<REAL_T> > Fcast_InputCatch_1) {
            this->Fcast_InputCatch_1 = Fcast_InputCatch_1;
        }

        int GetFcast_InputCatch_Basis() const {
            return Fcast_InputCatch_Basis;
        }

        void SetFcast_InputCatch_Basis(int Fcast_InputCatch_Basis) {
            this->Fcast_InputCatch_Basis = Fcast_InputCatch_Basis;
        }

        std::vector<REAL_T> GetFcast_Input_rd() const {
            return Fcast_Input_rd;
        }

        void SetFcast_Input_rd(std::vector<REAL_T> Fcast_Input_rd) {
            this->Fcast_Input_rd = Fcast_Input_rd;
        }

        std::vector<int> GetFcast_Loop_Control() const {
            return Fcast_Loop_Control;
        }

        void SetFcast_Loop_Control(std::vector<int> Fcast_Loop_Control) {
            this->Fcast_Loop_Control = Fcast_Loop_Control;
        }

        std::vector<REAL_T> GetFcast_MaxFleetCatch() const {
            return Fcast_MaxFleetCatch;
        }

        void SetFcast_MaxFleetCatch(std::vector<REAL_T> Fcast_MaxFleetCatch) {
            this->Fcast_MaxFleetCatch = Fcast_MaxFleetCatch;
        }

        std::valarray<std::valarray<REAL_T> > GetFcast_RelF() const {
            return Fcast_RelF;
        }

        void SetFcast_RelF(std::valarray<std::valarray<REAL_T> > Fcast_RelF) {
            this->Fcast_RelF = Fcast_RelF;
        }

        int GetFcast_RelF_Basis() const {
            return Fcast_RelF_Basis;
        }

        void SetFcast_RelF_Basis(int Fcast_RelF_Basis) {
            this->Fcast_RelF_Basis = Fcast_RelF_Basis;
        }

        std::valarray<std::valarray<REAL_T> > GetFcast_RelF_Input() const {
            return Fcast_RelF_Input;
        }

        void SetFcast_RelF_Input(std::valarray<std::valarray<REAL_T> > Fcast_RelF_Input) {
            this->Fcast_RelF_Input = Fcast_RelF_Input;
        }

        int GetFcast_RelF_yr1() const {
            return Fcast_RelF_yr1;
        }

        void SetFcast_RelF_yr1(int Fcast_RelF_yr1) {
            this->Fcast_RelF_yr1 = Fcast_RelF_yr1;
        }

        int GetFcast_RelF_yr2() const {
            return Fcast_RelF_yr2;
        }

        void SetFcast_RelF_yr2(int Fcast_RelF_yr2) {
            this->Fcast_RelF_yr2 = Fcast_RelF_yr2;
        }

        int GetFcast_Sel_yr1() const {
            return Fcast_Sel_yr1;
        }

        void SetFcast_Sel_yr1(int Fcast_Sel_yr1) {
            this->Fcast_Sel_yr1 = Fcast_Sel_yr1;
        }

        int GetFcast_Sel_yr2() const {
            return Fcast_Sel_yr2;
        }

        void SetFcast_Sel_yr2(int Fcast_Sel_yr2) {
            this->Fcast_Sel_yr2 = Fcast_Sel_yr2;
        }

        int GetFcast_catch_start() const {
            return Fcast_catch_start;
        }

        void SetFcast_catch_start(int Fcast_catch_start) {
            this->Fcast_catch_start = Fcast_catch_start;
        }

        int GetFcast_recr_PH() const {
            return Fcast_recr_PH;
        }

        void SetFcast_recr_PH(int Fcast_recr_PH) {
            this->Fcast_recr_PH = Fcast_recr_PH;
        }

        int GetFcast_recr_PH2() const {
            return Fcast_recr_PH2;
        }

        void SetFcast_recr_PH2(int Fcast_recr_PH2) {
            this->Fcast_recr_PH2 = Fcast_recr_PH2;
        }

        REAL_T GetFcast_recr_lambda() const {
            return Fcast_recr_lambda;
        }

        void SetFcast_recr_lambda(REAL_T Fcast_recr_lambda) {
            this->Fcast_recr_lambda = Fcast_recr_lambda;
        }

        std::vector<int> GetFcast_yr() const {
            return Fcast_yr;
        }

        void SetFcast_yr(std::vector<int> Fcast_yr) {
            this->Fcast_yr = Fcast_yr;
        }

        int GetFecund_Option() const {
            return Fecund_Option;
        }

        void SetFecund_Option(int Fecund_Option) {
            this->Fecund_Option = Fecund_Option;
        }

        int GetFirst_Mature_Age() const {
            return First_Mature_Age;
        }

        void SetFirst_Mature_Age(int First_Mature_Age) {
            this->First_Mature_Age = First_Mature_Age;
        }

        int GetFishon() const {
            return Fishon;
        }

        void SetFishon(int Fishon) {
            this->Fishon = Fishon;
        }

        std::vector<int> GetGP() const {
            return GP;
        }

        void SetGP(std::vector<int> GP) {
            this->GP = GP;
        }

        std::vector<int> GetGP2() const {
            return GP2;
        }

        void SetGP2(std::vector<int> GP2) {
            this->GP2 = GP2;
        }

        std::vector<int> GetGP3() const {
            return GP3;
        }

        void SetGP3(std::vector<int> GP3) {
            this->GP3 = GP3;
        }

        std::vector<int> GetGP4() const {
            return GP4;
        }

        void SetGP4(std::vector<int> GP4) {
            this->GP4 = GP4;
        }

        std::valarray<std::valarray<int> > GetGP_finder() const {
            return GP_finder;
        }

        void SetGP_finder(std::valarray<std::valarray<int> > GP_finder) {
            this->GP_finder = GP_finder;
        }

        int GetGrow_type() const {
            return Grow_type;
        }

        void SetGrow_type(int Grow_type) {
            this->Grow_type = Grow_type;
        }

        int GetGrowth_Std_Cnt() const {
            return Growth_Std_Cnt;
        }

        void SetGrowth_Std_Cnt(int Growth_Std_Cnt) {
            this->Growth_Std_Cnt = Growth_Std_Cnt;
        }

        std::vector<int> GetGrowth_Std_Pick() const {
            return Growth_Std_Pick;
        }

        void SetGrowth_Std_Pick(std::vector<int> Growth_Std_Pick) {
            this->Growth_Std_Pick = Growth_Std_Pick;
        }

        REAL_T GetH4010_bot() const {
            return H4010_bot;
        }

        void SetH4010_bot(REAL_T H4010_bot) {
            this->H4010_bot = H4010_bot;
        }

        REAL_T GetH4010_scale() const {
            return H4010_scale;
        }

        void SetH4010_scale(REAL_T H4010_scale) {
            this->H4010_scale = H4010_scale;
        }

        REAL_T GetH4010_top() const {
            return H4010_top;
        }

        void SetH4010_top(REAL_T H4010_top) {
            this->H4010_top = H4010_top;
        }

        int GetHarvestPolicy() const {
            return HarvestPolicy;
        }

        void SetHarvestPolicy(int HarvestPolicy) {
            this->HarvestPolicy = HarvestPolicy;
        }

        int GetHermaphro_Option() const {
            return Hermaphro_Option;
        }

        void SetHermaphro_Option(int Hermaphro_Option) {
            this->Hermaphro_Option = Hermaphro_Option;
        }

        int GetHermaphro_maleSPB() const {
            return Hermaphro_maleSPB;
        }

        void SetHermaphro_maleSPB(int Hermaphro_maleSPB) {
            this->Hermaphro_maleSPB = Hermaphro_maleSPB;
        }

        std::vector<int> GetHermaphro_more() const {
            return Hermaphro_more;
        }

        void SetHermaphro_more(std::vector<int> Hermaphro_more) {
            this->Hermaphro_more = Hermaphro_more;
        }

        int GetHermaphro_seas() const {
            return Hermaphro_seas;
        }

        void SetHermaphro_seas(int Hermaphro_seas) {
            this->Hermaphro_seas = Hermaphro_seas;
        }

        REAL_T GetImpl_Error_Std() const {
            return Impl_Error_Std;
        }

        void SetImpl_Error_Std(REAL_T Impl_Error_Std) {
            this->Impl_Error_Std = Impl_Error_Std;
        }

        int GetIp() const {
            return Ip;
        }

        void SetIp(int Ip) {
            this->Ip = Ip;
        }

        int GetL1() const {
            return L1;
        }

        void SetL1(int L1) {
            this->L1 = L1;
        }

        int GetL2() const {
            return L2;
        }

        void SetL2(int L2) {
            this->L2 = L2;
        }

        std::valarray<std::valarray<REAL_T> > GetLambda_changes() const {
            return Lambda_changes;
        }

        void SetLambda_changes(std::valarray<std::valarray<REAL_T> > Lambda_changes) {
            this->Lambda_changes = Lambda_changes;
        }

        std::vector<REAL_T> GetLbin_filter() const {
            return Lbin_filter;
        }

        void SetLbin_filter(std::vector<REAL_T> Lbin_filter) {
            this->Lbin_filter = Lbin_filter;
        }

        std::valarray<std::valarray<int> > GetLbin_hi() const {
            return Lbin_hi;
        }

        void SetLbin_hi(std::valarray<std::valarray<int> > Lbin_hi) {
            this->Lbin_hi = Lbin_hi;
        }

        std::valarray<std::valarray<int> > GetLbin_lo() const {
            return Lbin_lo;
        }

        void SetLbin_lo(std::valarray<std::valarray<int> > Lbin_lo) {
            this->Lbin_lo = Lbin_lo;
        }

        int GetLbin_method() const {
            return Lbin_method;
        }

        void SetLbin_method(int Lbin_method) {
            this->Lbin_method = Lbin_method;
        }

        int GetLenBin_option() const {
            return LenBin_option;
        }

        void SetLenBin_option(int LenBin_option) {
            this->LenBin_option = LenBin_option;
        }

        std::valarray<std::valarray<REAL_T> > GetLen_At_Age_rd() const {
            return Len_At_Age_rd;
        }

        void SetLen_At_Age_rd(std::valarray<std::valarray<REAL_T> > Len_At_Age_rd) {
            this->Len_At_Age_rd = Len_At_Age_rd;
        }

        std::valarray<std::valarray<REAL_T> > GetLn_obs_cr() const {
            return Ln_obs_cr;
        }

        void SetLn_obs_cr(std::valarray<std::valarray<REAL_T> > Ln_obs_cr) {
            this->Ln_obs_cr = Ln_obs_cr;
        }

        int GetMGP_CGD() const {
            return MGP_CGD;
        }

        void SetMGP_CGD(int MGP_CGD) {
            this->MGP_CGD = MGP_CGD;
        }

        std::vector<int> GetMG_active() const {
            return MG_active;
        }

        void SetMG_active(std::vector<int> MG_active) {
            this->MG_active = MG_active;
        }

        int GetMG_adjust_method() const {
            return MG_adjust_method;
        }

        void SetMG_adjust_method(int MG_adjust_method) {
            this->MG_adjust_method = MG_adjust_method;
        }

        std::valarray<std::valarray<REAL_T> > GetMGparm_1() const {
            return MGparm_1;
        }

        void SetMGparm_1(std::valarray<std::valarray<REAL_T> > MGparm_1) {
            this->MGparm_1 = MGparm_1;
        }

        std::vector<REAL_T> GetMGparm_CV() const {
            return MGparm_CV;
        }

        void SetMGparm_CV(std::vector<REAL_T> MGparm_CV) {
            this->MGparm_CV = MGparm_CV;
        }

        std::vector<REAL_T> GetMGparm_HI() const {
            return MGparm_HI;
        }

        void SetMGparm_HI(std::vector<REAL_T> MGparm_HI) {
            this->MGparm_HI = MGparm_HI;
        }

        int GetMGparm_Hermaphro() const {
            return MGparm_Hermaphro;
        }

        void SetMGparm_Hermaphro(int MGparm_Hermaphro) {
            this->MGparm_Hermaphro = MGparm_Hermaphro;
        }

        std::vector<REAL_T> GetMGparm_LO() const {
            return MGparm_LO;
        }

        void SetMGparm_LO(std::vector<REAL_T> MGparm_LO) {
            this->MGparm_LO = MGparm_LO;
        }

        std::vector<int> GetMGparm_PH() const {
            return MGparm_PH;
        }

        void SetMGparm_PH(std::vector<int> MGparm_PH) {
            this->MGparm_PH = MGparm_PH;
        }

        std::vector<REAL_T> GetMGparm_PR() const {
            return MGparm_PR;
        }

        void SetMGparm_PR(std::vector<REAL_T> MGparm_PR) {
            this->MGparm_PR = MGparm_PR;
        }

        std::vector<int> GetMGparm_PRtype() const {
            return MGparm_PRtype;
        }

        void SetMGparm_PRtype(std::vector<int> MGparm_PRtype) {
            this->MGparm_PRtype = MGparm_PRtype;
        }

        std::vector<REAL_T> GetMGparm_RD() const {
            return MGparm_RD;
        }

        void SetMGparm_RD(std::vector<REAL_T> MGparm_RD) {
            this->MGparm_RD = MGparm_RD;
        }

        std::valarray<std::valarray<REAL_T> > GetMGparm_blk_1() const {
            return MGparm_blk_1;
        }

        void SetMGparm_blk_1(std::valarray<std::valarray<REAL_T> > MGparm_blk_1) {
            this->MGparm_blk_1 = MGparm_blk_1;
        }

        int GetMGparm_def() const {
            return MGparm_def;
        }

        void SetMGparm_def(int MGparm_def) {
            this->MGparm_def = MGparm_def;
        }

        int GetMGparm_dev_PH() const {
            return MGparm_dev_PH;
        }

        void SetMGparm_dev_PH(int MGparm_dev_PH) {
            this->MGparm_dev_PH = MGparm_dev_PH;
        }

        std::vector<int> GetMGparm_dev_maxyr() const {
            return MGparm_dev_maxyr;
        }

        void SetMGparm_dev_maxyr(std::vector<int> MGparm_dev_maxyr) {
            this->MGparm_dev_maxyr = MGparm_dev_maxyr;
        }

        std::vector<int> GetMGparm_dev_minyr() const {
            return MGparm_dev_minyr;
        }

        void SetMGparm_dev_minyr(std::vector<int> MGparm_dev_minyr) {
            this->MGparm_dev_minyr = MGparm_dev_minyr;
        }

        std::vector<int> GetMGparm_dev_select() const {
            return MGparm_dev_select;
        }

        void SetMGparm_dev_select(std::vector<int> MGparm_dev_select) {
            this->MGparm_dev_select = MGparm_dev_select;
        }

        std::vector<REAL_T> GetMGparm_dev_stddev() const {
            return MGparm_dev_stddev;
        }

        void SetMGparm_dev_stddev(std::vector<REAL_T> MGparm_dev_stddev) {
            this->MGparm_dev_stddev = MGparm_dev_stddev;
        }

        std::vector<int> GetMGparm_dev_type() const {
            return MGparm_dev_type;
        }

        void SetMGparm_dev_type(std::vector<int> MGparm_dev_type) {
            this->MGparm_dev_type = MGparm_dev_type;
        }

        int GetMGparm_doseas() const {
            return MGparm_doseas;
        }

        void SetMGparm_doseas(int MGparm_doseas) {
            this->MGparm_doseas = MGparm_doseas;
        }

        std::vector<int> GetMGparm_env() const {
            return MGparm_env;
        }

        void SetMGparm_env(std::vector<int> MGparm_env) {
            this->MGparm_env = MGparm_env;
        }

        std::valarray<std::valarray<REAL_T> > GetMGparm_env_1() const {
            return MGparm_env_1;
        }

        void SetMGparm_env_1(std::valarray<std::valarray<REAL_T> > MGparm_env_1) {
            this->MGparm_env_1 = MGparm_env_1;
        }

        std::vector<int> GetMGparm_envtype() const {
            return MGparm_envtype;
        }

        void SetMGparm_envtype(std::vector<int> MGparm_envtype) {
            this->MGparm_envtype = MGparm_envtype;
        }

        std::vector<int> GetMGparm_envuse() const {
            return MGparm_envuse;
        }

        void SetMGparm_envuse(std::vector<int> MGparm_envuse) {
            this->MGparm_envuse = MGparm_envuse;
        }

        std::vector<int> GetMGparm_offset() const {
            return MGparm_offset;
        }

        void SetMGparm_offset(std::vector<int> MGparm_offset) {
            this->MGparm_offset = MGparm_offset;
        }

        int GetMGparm_per_def() const {
            return MGparm_per_def;
        }

        void SetMGparm_per_def(int MGparm_per_def) {
            this->MGparm_per_def = MGparm_per_def;
        }

        std::valarray<std::valarray<REAL_T> > GetMGparm_seas_1() const {
            return MGparm_seas_1;
        }

        void SetMGparm_seas_1(std::valarray<std::valarray<REAL_T> > MGparm_seas_1) {
            this->MGparm_seas_1 = MGparm_seas_1;
        }

        std::vector<int> GetMGparm_seas_effects() const {
            return MGparm_seas_effects;
        }

        void SetMGparm_seas_effects(std::vector<int> MGparm_seas_effects) {
            this->MGparm_seas_effects = MGparm_seas_effects;
        }

        std::valarray<std::valarray<REAL_T> > GetMGparm_trend_1() const {
            return MGparm_trend_1;
        }

        void SetMGparm_trend_1(std::valarray<std::valarray<REAL_T> > MGparm_trend_1) {
            this->MGparm_trend_1 = MGparm_trend_1;
        }

        std::vector<int> GetMGparm_trend_point() const {
            return MGparm_trend_point;
        }

        void SetMGparm_trend_point(std::vector<int> MGparm_trend_point) {
            this->MGparm_trend_point = MGparm_trend_point;
        }

        std::vector<int> GetMGparm_trend_rev() const {
            return MGparm_trend_rev;
        }

        void SetMGparm_trend_rev(std::vector<int> MGparm_trend_rev) {
            this->MGparm_trend_rev = MGparm_trend_rev;
        }

        std::vector<int> GetMGparm_trend_rev_1() const {
            return MGparm_trend_rev_1;
        }

        void SetMGparm_trend_rev_1(std::vector<int> MGparm_trend_rev_1) {
            this->MGparm_trend_rev_1 = MGparm_trend_rev_1;
        }

        std::vector<int> GetMaleselparm() const {
            return Maleselparm;
        }

        void SetMaleselparm(std::vector<int> Maleselparm) {
            this->Maleselparm = Maleselparm;
        }

        int GetMaturity_Option() const {
            return Maturity_Option;
        }

        void SetMaturity_Option(int Maturity_Option) {
            this->Maturity_Option = Maturity_Option;
        }

        std::vector<REAL_T> GetMax_Fcast_Catch() const {
            return Max_Fcast_Catch;
        }

        void SetMax_Fcast_Catch(std::vector<REAL_T> Max_Fcast_Catch) {
            this->Max_Fcast_Catch = Max_Fcast_Catch;
        }

        std::vector<int> GetMore_Std_Input() const {
            return More_Std_Input;
        }

        void SetMore_Std_Input(std::vector<int> More_Std_Input) {
            this->More_Std_Input = More_Std_Input;
        }

        std::vector<REAL_T> GetMorphcomp_havedata() const {
            return Morphcomp_havedata;
        }

        void SetMorphcomp_havedata(std::vector<REAL_T> Morphcomp_havedata) {
            this->Morphcomp_havedata = Morphcomp_havedata;
        }

        std::vector<REAL_T> GetMorphcomp_lambda() const {
            return Morphcomp_lambda;
        }

        void SetMorphcomp_lambda(std::vector<REAL_T> Morphcomp_lambda) {
            this->Morphcomp_lambda = Morphcomp_lambda;
        }

        REAL_T GetMorphcomp_mincomp() const {
            return Morphcomp_mincomp;
        }

        void SetMorphcomp_mincomp(REAL_T Morphcomp_mincomp) {
            this->Morphcomp_mincomp = Morphcomp_mincomp;
        }

        int GetMorphcomp_nmorph() const {
            return Morphcomp_nmorph;
        }

        void SetMorphcomp_nmorph(int Morphcomp_nmorph) {
            this->Morphcomp_nmorph = Morphcomp_nmorph;
        }

        int GetMorphcomp_nobs() const {
            return Morphcomp_nobs;
        }

        void SetMorphcomp_nobs(int Morphcomp_nobs) {
            this->Morphcomp_nobs = Morphcomp_nobs;
        }

        std::valarray<std::valarray<REAL_T> > GetMorphcomp_obs() const {
            return Morphcomp_obs;
        }

        void SetMorphcomp_obs(std::valarray<std::valarray<REAL_T> > Morphcomp_obs) {
            this->Morphcomp_obs = Morphcomp_obs;
        }

        int GetNP() const {
            return NP;
        }

        void SetNP(int NP) {
            this->NP = NP;
        }

        int GetN_Block_Designs() const {
            return N_Block_Designs;
        }

        void SetN_Block_Designs(int N_Block_Designs) {
            this->N_Block_Designs = N_Block_Designs;
        }

        int GetN_CC() const {
            return N_CC;
        }

        void SetN_CC(int N_CC) {
            this->N_CC = N_CC;
        }

        int GetN_DC() const {
            return N_DC;
        }

        void SetN_DC(int N_DC) {
            this->N_DC = N_DC;
        }

        int GetN_FC() const {
            return N_FC;
        }

        void SetN_FC(int N_FC) {
            this->N_FC = N_FC;
        }

        int GetN_Fcast_Yrs() const {
            return N_Fcast_Yrs;
        }

        void SetN_Fcast_Yrs(int N_Fcast_Yrs) {
            this->N_Fcast_Yrs = N_Fcast_Yrs;
        }

        int GetN_Fcast_catches() const {
            return N_Fcast_catches;
        }

        void SetN_Fcast_catches(int N_Fcast_catches) {
            this->N_Fcast_catches = N_Fcast_catches;
        }

        int GetN_GP() const {
            return N_GP;
        }

        void SetN_GP(int N_GP) {
            this->N_GP = N_GP;
        }

        int GetN_GP_sub() const {
            return N_GP_sub;
        }

        void SetN_GP_sub(int N_GP_sub) {
            this->N_GP_sub = N_GP_sub;
        }

        int GetN_MGparm() const {
            return N_MGparm;
        }

        void SetN_MGparm(int N_MGparm) {
            this->N_MGparm = N_MGparm;
        }

        int GetN_MGparm2() const {
            return N_MGparm2;
        }

        void SetN_MGparm2(int N_MGparm2) {
            this->N_MGparm2 = N_MGparm2;
        }

        int GetN_MGparm_blk() const {
            return N_MGparm_blk;
        }

        void SetN_MGparm_blk(int N_MGparm_blk) {
            this->N_MGparm_blk = N_MGparm_blk;
        }

        int GetN_MGparm_dev() const {
            return N_MGparm_dev;
        }

        void SetN_MGparm_dev(int N_MGparm_dev) {
            this->N_MGparm_dev = N_MGparm_dev;
        }

        int GetN_MGparm_env() const {
            return N_MGparm_env;
        }

        void SetN_MGparm_env(int N_MGparm_env) {
            this->N_MGparm_env = N_MGparm_env;
        }

        int GetN_MGparm_seas() const {
            return N_MGparm_seas;
        }

        void SetN_MGparm_seas(int N_MGparm_seas) {
            this->N_MGparm_seas = N_MGparm_seas;
        }

        int GetN_MGparm_trend() const {
            return N_MGparm_trend;
        }

        void SetN_MGparm_trend(int N_MGparm_trend) {
            this->N_MGparm_trend = N_MGparm_trend;
        }

        int GetN_MGparm_trend2() const {
            return N_MGparm_trend2;
        }

        void SetN_MGparm_trend2(int N_MGparm_trend2) {
            this->N_MGparm_trend2 = N_MGparm_trend2;
        }

        int GetN_ReadCatch() const {
            return N_ReadCatch;
        }

        void SetN_ReadCatch(int N_ReadCatch) {
            this->N_ReadCatch = N_ReadCatch;
        }

        int GetN_SC() const {
            return N_SC;
        }

        void SetN_SC(int N_SC) {
            this->N_SC = N_SC;
        }

        int GetN_SRparm() const {
            return N_SRparm;
        }

        void SetN_SRparm(int N_SRparm) {
            this->N_SRparm = N_SRparm;
        }

        int GetN_SRparm2() const {
            return N_SRparm2;
        }

        void SetN_SRparm2(int N_SRparm2) {
            this->N_SRparm2 = N_SRparm2;
        }

        int GetN_STD_Mgmt_Quant() const {
            return N_STD_Mgmt_Quant;
        }

        void SetN_STD_Mgmt_Quant(int N_STD_Mgmt_Quant) {
            this->N_STD_Mgmt_Quant = N_STD_Mgmt_Quant;
        }

        int GetN_STD_Yr() const {
            return N_STD_Yr;
        }

        void SetN_STD_Yr(int N_STD_Yr) {
            this->N_STD_Yr = N_STD_Yr;
        }

        int GetN_STD_Yr_Dep() const {
            return N_STD_Yr_Dep;
        }

        void SetN_STD_Yr_Dep(int N_STD_Yr_Dep) {
            this->N_STD_Yr_Dep = N_STD_Yr_Dep;
        }

        int GetN_STD_Yr_F() const {
            return N_STD_Yr_F;
        }

        void SetN_STD_Yr_F(int N_STD_Yr_F) {
            this->N_STD_Yr_F = N_STD_Yr_F;
        }

        int GetN_STD_Yr_Ofish() const {
            return N_STD_Yr_Ofish;
        }

        void SetN_STD_Yr_Ofish(int N_STD_Yr_Ofish) {
            this->N_STD_Yr_Ofish = N_STD_Yr_Ofish;
        }

        int GetN_STD_Yr_RD() const {
            return N_STD_Yr_RD;
        }

        void SetN_STD_Yr_RD(int N_STD_Yr_RD) {
            this->N_STD_Yr_RD = N_STD_Yr_RD;
        }

        int GetN_TG() const {
            return N_TG;
        }

        void SetN_TG(int N_TG) {
            this->N_TG = N_TG;
        }

        int GetN_TG2() const {
            return N_TG2;
        }

        void SetN_TG2(int N_TG2) {
            this->N_TG2 = N_TG2;
        }

        int GetN_TG_recap() const {
            return N_TG_recap;
        }

        void SetN_TG_recap(int N_TG_recap) {
            this->N_TG_recap = N_TG_recap;
        }

        int GetN_WTage_maxage() const {
            return N_WTage_maxage;
        }

        void SetN_WTage_maxage(int N_WTage_maxage) {
            this->N_WTage_maxage = N_WTage_maxage;
        }

        int GetN_WTage_rd() const {
            return N_WTage_rd;
        }

        void SetN_WTage_rd(int N_WTage_rd) {
            this->N_WTage_rd = N_WTage_rd;
        }

        int GetN_ageerr() const {
            return N_ageerr;
        }

        void SetN_ageerr(int N_ageerr) {
            this->N_ageerr = N_ageerr;
        }

        int GetN_envdata() const {
            return N_envdata;
        }

        void SetN_envdata(int N_envdata) {
            this->N_envdata = N_envdata;
        }

        int GetN_envvar() const {
            return N_envvar;
        }

        void SetN_envvar(int N_envvar) {
            this->N_envvar = N_envvar;
        }

        int GetN_growparms() const {
            return N_growparms;
        }

        void SetN_growparms(int N_growparms) {
            this->N_growparms = N_growparms;
        }

        int GetN_lambda_changes() const {
            return N_lambda_changes;
        }

        void SetN_lambda_changes(int N_lambda_changes) {
            this->N_lambda_changes = N_lambda_changes;
        }

        int GetN_natMparms() const {
            return N_natMparms;
        }

        void SetN_natMparms(int N_natMparms) {
            this->N_natMparms = N_natMparms;
        }

        int GetN_nudata() const {
            return N_nudata;
        }

        void SetN_nudata(int N_nudata) {
            this->N_nudata = N_nudata;
        }

        int GetN_prof_var() const {
            return N_prof_var;
        }

        void SetN_prof_var(int N_prof_var) {
            this->N_prof_var = N_prof_var;
        }

        int GetN_selparm() const {
            return N_selparm;
        }

        void SetN_selparm(int N_selparm) {
            this->N_selparm = N_selparm;
        }

        int GetN_selparm2() const {
            return N_selparm2;
        }

        void SetN_selparm2(int N_selparm2) {
            this->N_selparm2 = N_selparm2;
        }

        int GetN_selparm_blk() const {
            return N_selparm_blk;
        }

        void SetN_selparm_blk(int N_selparm_blk) {
            this->N_selparm_blk = N_selparm_blk;
        }

        int GetN_selparm_dev() const {
            return N_selparm_dev;
        }

        void SetN_selparm_dev(int N_selparm_dev) {
            this->N_selparm_dev = N_selparm_dev;
        }

        int GetN_selparm_dev_tot() const {
            return N_selparm_dev_tot;
        }

        void SetN_selparm_dev_tot(int N_selparm_dev_tot) {
            this->N_selparm_dev_tot = N_selparm_dev_tot;
        }

        int GetN_selparm_env() const {
            return N_selparm_env;
        }

        void SetN_selparm_env(int N_selparm_env) {
            this->N_selparm_env = N_selparm_env;
        }

        int GetN_selparm_trend() const {
            return N_selparm_trend;
        }

        void SetN_selparm_trend(int N_selparm_trend) {
            this->N_selparm_trend = N_selparm_trend;
        }

        int GetN_selparm_trend2() const {
            return N_selparm_trend2;
        }

        void SetN_selparm_trend2(int N_selparm_trend2) {
            this->N_selparm_trend2 = N_selparm_trend2;
        }

        std::vector<int> GetN_selparmvec() const {
            return N_selparmvec;
        }

        void SetN_selparmvec(std::vector<int> N_selparmvec) {
            this->N_selparmvec = N_selparmvec;
        }

        std::vector<int> GetN_suprper_SzFreq() const {
            return N_suprper_SzFreq;
        }

        void SetN_suprper_SzFreq(std::vector<int> N_suprper_SzFreq) {
            this->N_suprper_SzFreq = N_suprper_SzFreq;
        }

        std::vector<int> GetN_suprper_a() const {
            return N_suprper_a;
        }

        void SetN_suprper_a(std::vector<int> N_suprper_a) {
            this->N_suprper_a = N_suprper_a;
        }

        std::vector<int> GetN_suprper_cr() const {
            return N_suprper_cr;
        }

        void SetN_suprper_cr(std::vector<int> N_suprper_cr) {
            this->N_suprper_cr = N_suprper_cr;
        }

        std::vector<int> GetN_suprper_disc() const {
            return N_suprper_disc;
        }

        void SetN_suprper_disc(std::vector<int> N_suprper_disc) {
            this->N_suprper_disc = N_suprper_disc;
        }

        std::vector<int> GetN_suprper_l() const {
            return N_suprper_l;
        }

        void SetN_suprper_l(std::vector<int> N_suprper_l) {
            this->N_suprper_l = N_suprper_l;
        }

        std::vector<int> GetN_suprper_ms() const {
            return N_suprper_ms;
        }

        void SetN_suprper_ms(std::vector<int> N_suprper_ms) {
            this->N_suprper_ms = N_suprper_ms;
        }

        int GetN_warn() const {
            return N_warn;
        }

        void SetN_warn(int N_warn) {
            this->N_warn = N_warn;
        }

        int GetNatAge_Std_Cnt() const {
            return NatAge_Std_Cnt;
        }

        void SetNatAge_Std_Cnt(int NatAge_Std_Cnt) {
            this->NatAge_Std_Cnt = NatAge_Std_Cnt;
        }

        std::vector<int> GetNatAge_Std_Pick() const {
            return NatAge_Std_Pick;
        }

        void SetNatAge_Std_Pick(std::vector<int> NatAge_Std_Pick) {
            this->NatAge_Std_Pick = NatAge_Std_Pick;
        }

        int GetNatAge_Std_Year() const {
            return NatAge_Std_Year;
        }

        void SetNatAge_Std_Year(int NatAge_Std_Year) {
            this->NatAge_Std_Year = NatAge_Std_Year;
        }

        std::vector<REAL_T> GetNatM_break() const {
            return NatM_break;
        }

        void SetNatM_break(std::vector<REAL_T> NatM_break) {
            this->NatM_break = NatM_break;
        }

        std::vector<int> GetNblk() const {
            return Nblk;
        }

        void SetNblk(std::vector<int> Nblk) {
            this->Nblk = Nblk;
        }

        std::vector<int> GetNblk2() const {
            return Nblk2;
        }

        void SetNblk2(std::vector<int> Nblk2) {
            this->Nblk2 = Nblk2;
        }

        int GetNcycle() const {
            return Ncycle;
        }

        void SetNcycle(int Ncycle) {
            this->Ncycle = Ncycle;
        }

        int GetNdisc_fleets() const {
            return Ndisc_fleets;
        }

        void SetNdisc_fleets(int Ndisc_fleets) {
            this->Ndisc_fleets = Ndisc_fleets;
        }

        int GetNfleet() const {
            return Nfleet;
        }

        void SetNfleet(int Nfleet) {
            this->Nfleet = Nfleet;
        }

        REAL_T GetNilNumbers() const {
            return NilNumbers;
        }

        void SetNilNumbers(REAL_T NilNumbers) {
            this->NilNumbers = NilNumbers;
        }

        int GetNo_Report() const {
            return No_Report;
        }

        void SetNo_Report(int No_Report) {
            this->No_Report = No_Report;
        }

        std::vector<int> GetNobs_a() const {
            return Nobs_a;
        }

        void SetNobs_a(std::vector<int> Nobs_a) {
            this->Nobs_a = Nobs_a;
        }

        int GetNobs_a_tot() const {
            return Nobs_a_tot;
        }

        void SetNobs_a_tot(int Nobs_a_tot) {
            this->Nobs_a_tot = Nobs_a_tot;
        }

        std::vector<int> GetNobs_l() const {
            return Nobs_l;
        }

        void SetNobs_l(std::vector<int> Nobs_l) {
            this->Nobs_l = Nobs_l;
        }

        int GetNobs_l_tot() const {
            return Nobs_l_tot;
        }

        void SetNobs_l_tot(int Nobs_l_tot) {
            this->Nobs_l_tot = Nobs_l_tot;
        }

        std::vector<int> GetNobs_ms() const {
            return Nobs_ms;
        }

        void SetNobs_ms(std::vector<int> Nobs_ms) {
            this->Nobs_ms = Nobs_ms;
        }

        int GetNsurvey() const {
            return Nsurvey;
        }

        void SetNsurvey(int Nsurvey) {
            this->Nsurvey = Nsurvey;
        }

        int GetNtypes() const {
            return Ntypes;
        }

        void SetNtypes(int Ntypes) {
            this->Ntypes = Ntypes;
        }

        int GetParCount() const {
            return ParCount;
        }

        void SetParCount(int ParCount) {
            this->ParCount = ParCount;
        }

        std::vector<REAL_T> GetPopBin_Read() const {
            return PopBin_Read;
        }

        void SetPopBin_Read(std::vector<REAL_T> PopBin_Read) {
            this->PopBin_Read = PopBin_Read;
        }

        int GetQ_Npar() const {
            return Q_Npar;
        }

        void SetQ_Npar(int Q_Npar) {
            this->Q_Npar = Q_Npar;
        }

        int GetQ_Npar2() const {
            return Q_Npar2;
        }

        void SetQ_Npar2(int Q_Npar2) {
            this->Q_Npar2 = Q_Npar2;
        }

        std::valarray<std::valarray<REAL_T> > GetQ_parm_1() const {
            return Q_parm_1;
        }

        void SetQ_parm_1(std::valarray<std::valarray<REAL_T> > Q_parm_1) {
            this->Q_parm_1 = Q_parm_1;
        }

        std::valarray<std::valarray<REAL_T> > GetQ_parm_2() const {
            return Q_parm_2;
        }

        void SetQ_parm_2(std::valarray<std::valarray<REAL_T> > Q_parm_2) {
            this->Q_parm_2 = Q_parm_2;
        }

        std::vector<REAL_T> GetQ_parm_HI() const {
            return Q_parm_HI;
        }

        void SetQ_parm_HI(std::vector<REAL_T> Q_parm_HI) {
            this->Q_parm_HI = Q_parm_HI;
        }

        std::vector<REAL_T> GetQ_parm_LO() const {
            return Q_parm_LO;
        }

        void SetQ_parm_LO(std::vector<REAL_T> Q_parm_LO) {
            this->Q_parm_LO = Q_parm_LO;
        }

        std::vector<int> GetQ_parm_PH() const {
            return Q_parm_PH;
        }

        void SetQ_parm_PH(std::vector<int> Q_parm_PH) {
            this->Q_parm_PH = Q_parm_PH;
        }

        int GetQ_parm_detail() const {
            return Q_parm_detail;
        }

        void SetQ_parm_detail(int Q_parm_detail) {
            this->Q_parm_detail = Q_parm_detail;
        }

        std::valarray<std::valarray<REAL_T> > GetQ_setup() const {
            return Q_setup;
        }

        void SetQ_setup(std::valarray<std::valarray<REAL_T> > Q_setup) {
            this->Q_setup = Q_setup;
        }

        std::valarray<std::valarray<int> > GetQ_setup_parms() const {
            return Q_setup_parms;
        }

        void SetQ_setup_parms(std::valarray<std::valarray<int> > Q_setup_parms) {
            this->Q_setup_parms = Q_setup_parms;
        }

        int GetRebuild_Ydecl() const {
            return Rebuild_Ydecl;
        }

        void SetRebuild_Ydecl(int Rebuild_Ydecl) {
            this->Rebuild_Ydecl = Rebuild_Ydecl;
        }

        int GetRebuild_Yinit() const {
            return Rebuild_Yinit;
        }

        void SetRebuild_Yinit(int Rebuild_Yinit) {
            this->Rebuild_Yinit = Rebuild_Yinit;
        }

        std::vector<int> GetRetainParm() const {
            return RetainParm;
        }

        void SetRetainParm(std::vector<int> RetainParm) {
            this->RetainParm = RetainParm;
        }

        REAL_T GetSD_add_to_LAA() const {
            return SD_add_to_LAA;
        }

        void SetSD_add_to_LAA(REAL_T SD_add_to_LAA) {
            this->SD_add_to_LAA = SD_add_to_LAA;
        }

        int GetSPR_reporting() const {
            return SPR_reporting;
        }

        void SetSPR_reporting(int SPR_reporting) {
            this->SPR_reporting = SPR_reporting;
        }

        REAL_T GetSPR_target() const {
            return SPR_target;
        }

        void SetSPR_target(REAL_T SPR_target) {
            this->SPR_target = SPR_target;
        }

        int GetSR_autocorr() const {
            return SR_autocorr;
        }

        void SetSR_autocorr(int SR_autocorr) {
            this->SR_autocorr = SR_autocorr;
        }

        int GetSR_env_link() const {
            return SR_env_link;
        }

        void SetSR_env_link(int SR_env_link) {
            this->SR_env_link = SR_env_link;
        }

        int GetSR_env_target() const {
            return SR_env_target;
        }

        void SetSR_env_target(int SR_env_target) {
            this->SR_env_target = SR_env_target;
        }

        int GetSR_env_target_RD() const {
            return SR_env_target_RD;
        }

        void SetSR_env_target_RD(int SR_env_target_RD) {
            this->SR_env_target_RD = SR_env_target_RD;
        }

        int GetSR_fxn() const {
            return SR_fxn;
        }

        void SetSR_fxn(int SR_fxn) {
            this->SR_fxn = SR_fxn;
        }

        std::valarray<std::valarray<REAL_T> > GetSR_parm_1() const {
            return SR_parm_1;
        }

        void SetSR_parm_1(std::valarray<std::valarray<REAL_T> > SR_parm_1) {
            this->SR_parm_1 = SR_parm_1;
        }

        std::vector<REAL_T> GetSRvec_HI() const {
            return SRvec_HI;
        }

        void SetSRvec_HI(std::vector<REAL_T> SRvec_HI) {
            this->SRvec_HI = SRvec_HI;
        }

        std::vector<REAL_T> GetSRvec_LO() const {
            return SRvec_LO;
        }

        void SetSRvec_LO(std::vector<REAL_T> SRvec_LO) {
            this->SRvec_LO = SRvec_LO;
        }

        std::vector<int> GetSRvec_PH() const {
            return SRvec_PH;
        }

        void SetSRvec_PH(std::vector<int> SRvec_PH) {
            this->SRvec_PH = SRvec_PH;
        }

        std::vector<int> GetSTD_Yr_RD() const {
            return STD_Yr_RD;
        }

        void SetSTD_Yr_RD(std::vector<int> STD_Yr_RD) {
            this->STD_Yr_RD = STD_Yr_RD;
        }

        std::vector<int> GetSTD_Yr_Reverse() const {
            return STD_Yr_Reverse;
        }

        void SetSTD_Yr_Reverse(std::vector<int> STD_Yr_Reverse) {
            this->STD_Yr_Reverse = STD_Yr_Reverse;
        }

        std::vector<int> GetSTD_Yr_Reverse_Dep() const {
            return STD_Yr_Reverse_Dep;
        }

        void SetSTD_Yr_Reverse_Dep(std::vector<int> STD_Yr_Reverse_Dep) {
            this->STD_Yr_Reverse_Dep = STD_Yr_Reverse_Dep;
        }

        std::vector<int> GetSTD_Yr_Reverse_F() const {
            return STD_Yr_Reverse_F;
        }

        void SetSTD_Yr_Reverse_F(std::vector<int> STD_Yr_Reverse_F) {
            this->STD_Yr_Reverse_F = STD_Yr_Reverse_F;
        }

        std::vector<int> GetSTD_Yr_Reverse_Ofish() const {
            return STD_Yr_Reverse_Ofish;
        }

        void SetSTD_Yr_Reverse_Ofish(std::vector<int> STD_Yr_Reverse_Ofish) {
            this->STD_Yr_Reverse_Ofish = STD_Yr_Reverse_Ofish;
        }

        int GetSTD_Yr_max() const {
            return STD_Yr_max;
        }

        void SetSTD_Yr_max(int STD_Yr_max) {
            this->STD_Yr_max = STD_Yr_max;
        }

        int GetSTD_Yr_min() const {
            return STD_Yr_min;
        }

        void SetSTD_Yr_min(int STD_Yr_min) {
            this->STD_Yr_min = STD_Yr_min;
        }

        int GetSelex_Std_AL() const {
            return Selex_Std_AL;
        }

        void SetSelex_Std_AL(int Selex_Std_AL) {
            this->Selex_Std_AL = Selex_Std_AL;
        }

        int GetSelex_Std_Cnt() const {
            return Selex_Std_Cnt;
        }

        void SetSelex_Std_Cnt(int Selex_Std_Cnt) {
            this->Selex_Std_Cnt = Selex_Std_Cnt;
        }

        std::vector<int> GetSelex_Std_Pick() const {
            return Selex_Std_Pick;
        }

        void SetSelex_Std_Pick(std::vector<int> Selex_Std_Pick) {
            this->Selex_Std_Pick = Selex_Std_Pick;
        }

        int GetSelex_Std_Year() const {
            return Selex_Std_Year;
        }

        void SetSelex_Std_Year(int Selex_Std_Year) {
            this->Selex_Std_Year = Selex_Std_Year;
        }

        int GetSmry_Age() const {
            return Smry_Age;
        }

        void SetSmry_Age(int Smry_Age) {
            this->Smry_Age = Smry_Age;
        }

        int GetSoftBound() const {
            return SoftBound;
        }

        void SetSoftBound(int SoftBound) {
            this->SoftBound = SoftBound;
        }

        int GetSzFreqMethod() const {
            return SzFreqMethod;
        }

        void SetSzFreqMethod(int SzFreqMethod) {
            this->SzFreqMethod = SzFreqMethod;
        }

        int GetSzFreqMethod_seas() const {
            return SzFreqMethod_seas;
        }

        void SetSzFreqMethod_seas(int SzFreqMethod_seas) {
            this->SzFreqMethod_seas = SzFreqMethod_seas;
        }

        std::vector<REAL_T> GetSzFreq_HaveObs() const {
            return SzFreq_HaveObs;
        }

        void SetSzFreq_HaveObs(std::vector<REAL_T> SzFreq_HaveObs) {
            this->SzFreq_HaveObs = SzFreq_HaveObs;
        }

        std::valarray<std::valarray<int> > GetSzFreq_HaveObs2() const {
            return SzFreq_HaveObs2;
        }

        void SetSzFreq_HaveObs2(std::valarray<std::valarray<int> > SzFreq_HaveObs2) {
            this->SzFreq_HaveObs2 = SzFreq_HaveObs2;
        }

        std::valarray<std::valarray<int> > GetSzFreq_LikeComponent() const {
            return SzFreq_LikeComponent;
        }

        void SetSzFreq_LikeComponent(std::valarray<std::valarray<int> > SzFreq_LikeComponent) {
            this->SzFreq_LikeComponent = SzFreq_LikeComponent;
        }

        int GetSzFreq_N_Like() const {
            return SzFreq_N_Like;
        }

        void SetSzFreq_N_Like(int SzFreq_N_Like) {
            this->SzFreq_N_Like = SzFreq_N_Like;
        }

        std::vector<int> GetSzFreq_Nbins() const {
            return SzFreq_Nbins;
        }

        void SetSzFreq_Nbins(std::vector<int> SzFreq_Nbins) {
            this->SzFreq_Nbins = SzFreq_Nbins;
        }

        std::vector<int> GetSzFreq_Nbins3() const {
            return SzFreq_Nbins3;
        }

        void SetSzFreq_Nbins3(std::vector<int> SzFreq_Nbins3) {
            this->SzFreq_Nbins3 = SzFreq_Nbins3;
        }

        std::vector<int> GetSzFreq_Nbins_seas_g() const {
            return SzFreq_Nbins_seas_g;
        }

        void SetSzFreq_Nbins_seas_g(std::vector<int> SzFreq_Nbins_seas_g) {
            this->SzFreq_Nbins_seas_g = SzFreq_Nbins_seas_g;
        }

        int GetSzFreq_Nmeth() const {
            return SzFreq_Nmeth;
        }

        void SetSzFreq_Nmeth(int SzFreq_Nmeth) {
            this->SzFreq_Nmeth = SzFreq_Nmeth;
        }

        std::vector<int> GetSzFreq_Omit_Small() const {
            return SzFreq_Omit_Small;
        }

        void SetSzFreq_Omit_Small(std::vector<int> SzFreq_Omit_Small) {
            this->SzFreq_Omit_Small = SzFreq_Omit_Small;
        }

        std::vector<int> GetSzFreq_Setup() const {
            return SzFreq_Setup;
        }

        void SetSzFreq_Setup(std::vector<int> SzFreq_Setup) {
            this->SzFreq_Setup = SzFreq_Setup;
        }

        std::vector<int> GetSzFreq_Setup2() const {
            return SzFreq_Setup2;
        }

        void SetSzFreq_Setup2(std::vector<int> SzFreq_Setup2) {
            this->SzFreq_Setup2 = SzFreq_Setup2;
        }

        std::valarray<std::valarray<REAL_T> > GetSzFreq_bins() const {
            return SzFreq_bins;
        }

        void SetSzFreq_bins(std::valarray<std::valarray<REAL_T> > SzFreq_bins) {
            this->SzFreq_bins = SzFreq_bins;
        }

        std::valarray<std::valarray<REAL_T> > GetSzFreq_bins1() const {
            return SzFreq_bins1;
        }

        void SetSzFreq_bins1(std::valarray<std::valarray<REAL_T> > SzFreq_bins1) {
            this->SzFreq_bins1 = SzFreq_bins1;
        }

        std::valarray<std::valarray<REAL_T> > GetSzFreq_bins2() const {
            return SzFreq_bins2;
        }

        void SetSzFreq_bins2(std::valarray<std::valarray<REAL_T> > SzFreq_bins2) {
            this->SzFreq_bins2 = SzFreq_bins2;
        }

        std::vector<REAL_T> GetSzFreq_eachlike() const {
            return SzFreq_eachlike;
        }

        void SetSzFreq_eachlike(std::vector<REAL_T> SzFreq_eachlike) {
            this->SzFreq_eachlike = SzFreq_eachlike;
        }

        std::vector<REAL_T> GetSzFreq_effN() const {
            return SzFreq_effN;
        }

        void SetSzFreq_effN(std::vector<REAL_T> SzFreq_effN) {
            this->SzFreq_effN = SzFreq_effN;
        }

        std::valarray<std::valarray<REAL_T> > GetSzFreq_lambda() const {
            return SzFreq_lambda;
        }

        void SetSzFreq_lambda(std::valarray<std::valarray<REAL_T> > SzFreq_lambda) {
            this->SzFreq_lambda = SzFreq_lambda;
        }

        std::vector<REAL_T> GetSzFreq_like_base() const {
            return SzFreq_like_base;
        }

        void SetSzFreq_like_base(std::vector<REAL_T> SzFreq_like_base) {
            this->SzFreq_like_base = SzFreq_like_base;
        }

        std::vector<REAL_T> GetSzFreq_mincomp() const {
            return SzFreq_mincomp;
        }

        void SetSzFreq_mincomp(std::vector<REAL_T> SzFreq_mincomp) {
            this->SzFreq_mincomp = SzFreq_mincomp;
        }

        std::vector<int> GetSzFreq_nobs() const {
            return SzFreq_nobs;
        }

        void SetSzFreq_nobs(std::vector<int> SzFreq_nobs) {
            this->SzFreq_nobs = SzFreq_nobs;
        }

        std::valarray<std::valarray<REAL_T> > GetSzFreq_obs() const {
            return SzFreq_obs;
        }

        void SetSzFreq_obs(std::valarray<std::valarray<REAL_T> > SzFreq_obs) {
            this->SzFreq_obs = SzFreq_obs;
        }

        std::valarray<std::valarray<REAL_T> > GetSzFreq_obs1() const {
            return SzFreq_obs1;
        }

        void SetSzFreq_obs1(std::valarray<std::valarray<REAL_T> > SzFreq_obs1) {
            this->SzFreq_obs1 = SzFreq_obs1;
        }

        std::valarray<std::valarray<int> > GetSzFreq_obs_hdr() const {
            return SzFreq_obs_hdr;
        }

        void SetSzFreq_obs_hdr(std::valarray<std::valarray<int> > SzFreq_obs_hdr) {
            this->SzFreq_obs_hdr = SzFreq_obs_hdr;
        }

        std::vector<REAL_T> GetSzFreq_sampleN() const {
            return SzFreq_sampleN;
        }

        void SetSzFreq_sampleN(std::vector<REAL_T> SzFreq_sampleN) {
            this->SzFreq_sampleN = SzFreq_sampleN;
        }

        std::vector<int> GetSzFreq_scale() const {
            return SzFreq_scale;
        }

        void SetSzFreq_scale(std::vector<int> SzFreq_scale) {
            this->SzFreq_scale = SzFreq_scale;
        }

        int GetSzFreq_totobs() const {
            return SzFreq_totobs;
        }

        void SetSzFreq_totobs(int SzFreq_totobs) {
            this->SzFreq_totobs = SzFreq_totobs;
        }

        std::vector<int> GetSzFreq_units() const {
            return SzFreq_units;
        }

        void SetSzFreq_units(std::vector<int> SzFreq_units) {
            this->SzFreq_units = SzFreq_units;
        }

        int GetTG() const {
            return TG;
        }

        void SetTG(int TG) {
            this->TG = TG;
        }

        int GetTG_custom() const {
            return TG_custom;
        }

        void SetTG_custom(int TG_custom) {
            this->TG_custom = TG_custom;
        }

        std::vector<int> GetTG_endtime() const {
            return TG_endtime;
        }

        void SetTG_endtime(std::vector<int> TG_endtime) {
            this->TG_endtime = TG_endtime;
        }

        std::valarray<std::valarray<REAL_T> > GetTG_lambda1() const {
            return TG_lambda1;
        }

        void SetTG_lambda1(std::valarray<std::valarray<REAL_T> > TG_lambda1) {
            this->TG_lambda1 = TG_lambda1;
        }

        std::valarray<std::valarray<REAL_T> > GetTG_lambda2() const {
            return TG_lambda2;
        }

        void SetTG_lambda2(std::valarray<std::valarray<REAL_T> > TG_lambda2) {
            this->TG_lambda2 = TG_lambda2;
        }

        int GetTG_maxperiods() const {
            return TG_maxperiods;
        }

        void SetTG_maxperiods(int TG_maxperiods) {
            this->TG_maxperiods = TG_maxperiods;
        }

        int GetTG_mixperiod() const {
            return TG_mixperiod;
        }

        void SetTG_mixperiod(int TG_mixperiod) {
            this->TG_mixperiod = TG_mixperiod;
        }

        std::valarray<std::valarray<REAL_T> > GetTG_parm1() const {
            return TG_parm1;
        }

        void SetTG_parm1(std::valarray<std::valarray<REAL_T> > TG_parm1) {
            this->TG_parm1 = TG_parm1;
        }

        std::valarray<std::valarray<REAL_T> > GetTG_parm2() const {
            return TG_parm2;
        }

        void SetTG_parm2(std::valarray<std::valarray<REAL_T> > TG_parm2) {
            this->TG_parm2 = TG_parm2;
        }

        std::vector<REAL_T> GetTG_parm_HI() const {
            return TG_parm_HI;
        }

        void SetTG_parm_HI(std::vector<REAL_T> TG_parm_HI) {
            this->TG_parm_HI = TG_parm_HI;
        }

        std::vector<REAL_T> GetTG_parm_LO() const {
            return TG_parm_LO;
        }

        void SetTG_parm_LO(std::vector<REAL_T> TG_parm_LO) {
            this->TG_parm_LO = TG_parm_LO;
        }

        std::vector<int> GetTG_parm_PH() const {
            return TG_parm_PH;
        }

        void SetTG_parm_PH(std::vector<int> TG_parm_PH) {
            this->TG_parm_PH = TG_parm_PH;
        }

        std::valarray<std::valarray<REAL_T> > GetTG_recap_data() const {
            return TG_recap_data;
        }

        void SetTG_recap_data(std::valarray<std::valarray<REAL_T> > TG_recap_data) {
            this->TG_recap_data = TG_recap_data;
        }

        std::vector<REAL_T> GetTG_recap_obs() const {
            return TG_recap_obs;
        }

        void SetTG_recap_obs(std::vector<REAL_T> TG_recap_obs) {
            this->TG_recap_obs = TG_recap_obs;
        }

        std::valarray<std::valarray<REAL_T> > GetTG_release() const {
            return TG_release;
        }

        void SetTG_release(std::valarray<std::valarray<REAL_T> > TG_release) {
            this->TG_release = TG_release;
        }

        int GetTG_t() const {
            return TG_t;
        }

        void SetTG_t(int TG_t) {
            this->TG_t = TG_t;
        }

        std::vector<REAL_T> GetTG_temp() const {
            return TG_temp;
        }

        void SetTG_temp(std::vector<REAL_T> TG_temp) {
            this->TG_temp = TG_temp;
        }

        int GetTG_timestart() const {
            return TG_timestart;
        }

        void SetTG_timestart(int TG_timestart) {
            this->TG_timestart = TG_timestart;
        }

        std::valarray<std::valarray<int> > GetTG_use_morph() const {
            return TG_use_morph;
        }

        void SetTG_use_morph(std::valarray<std::valarray<int> > TG_use_morph) {
            this->TG_use_morph = TG_use_morph;
        }

        int GetTimeMax() const {
            return TimeMax;
        }

        void SetTimeMax(int TimeMax) {
            this->TimeMax = TimeMax;
        }

        int GetTimeMax_Fcast_std() const {
            return TimeMax_Fcast_std;
        }

        void SetTimeMax_Fcast_std(int TimeMax_Fcast_std) {
            this->TimeMax_Fcast_std = TimeMax_Fcast_std;
        }

        int GetTurn_off_phase() const {
            return Turn_off_phase;
        }

        void SetTurn_off_phase(int Turn_off_phase) {
            this->Turn_off_phase = Turn_off_phase;
        }

        int GetUse_AgeKeyZero() const {
            return Use_AgeKeyZero;
        }

        void SetUse_AgeKeyZero(int Use_AgeKeyZero) {
            this->Use_AgeKeyZero = Use_AgeKeyZero;
        }

        std::vector<REAL_T> GetWTage_emp() const {
            return WTage_emp;
        }

        void SetWTage_emp(std::vector<REAL_T> WTage_emp) {
            this->WTage_emp = WTage_emp;
        }

        std::valarray<std::valarray<REAL_T> > GetWTage_in() const {
            return WTage_in;
        }

        void SetWTage_in(std::valarray<std::valarray<REAL_T> > WTage_in) {
            this->WTage_in = WTage_in;
        }

        int GetWTage_rd() const {
            return WTage_rd;
        }

        void SetWTage_rd(int WTage_rd) {
            this->WTage_rd = WTage_rd;
        }

        int GetYrMax() const {
            return YrMax;
        }

        void SetYrMax(int YrMax) {
            this->YrMax = YrMax;
        }

        int GetA() const {
            return a;
        }

        void SetA(int a) {
            this->a = a;
        }

        int GetA1() const {
            return a1;
        }

        void SetA1(int a1) {
            this->a1 = a1;
        }

        int GetActive_count() const {
            return active_count;
        }

        void SetActive_count(int active_count) {
            this->active_count = active_count;
        }

        std::vector<int> GetActive_parm() const {
            return active_parm;
        }

        void SetActive_parm(std::vector<int> active_parm) {
            this->active_parm = active_parm;
        }

        int GetActive_parms() const {
            return active_parms;
        }

        void SetActive_parms(int active_parms) {
            this->active_parms = active_parms;
        }

        std::vector<REAL_T> GetAge_bins() const {
            return age_bins;
        }

        void SetAge_bins(std::vector<REAL_T> age_bins) {
            this->age_bins = age_bins;
        }

        std::vector<REAL_T> GetAge_bins1() const {
            return age_bins1;
        }

        void SetAge_bins1(std::vector<REAL_T> age_bins1) {
            this->age_bins1 = age_bins1;
        }

        std::vector<REAL_T> GetAge_err_rd() const {
            return age_err_rd;
        }

        void SetAge_err_rd(std::vector<REAL_T> age_err_rd) {
            this->age_err_rd = age_err_rd;
        }

        std::valarray<std::valarray<REAL_T> > GetAge_lambda() const {
            return age_lambda;
        }

        void SetAge_lambda(std::valarray<std::valarray<REAL_T> > age_lambda) {
            this->age_lambda = age_lambda;
        }

        std::vector<int> GetAge_vector() const {
            return age_vector;
        }

        void SetAge_vector(std::vector<int> age_vector) {
            this->age_vector = age_vector;
        }

        std::valarray<std::valarray<REAL_T> > GetAgedata() const {
            return agedata;
        }

        void SetAgedata(std::valarray<std::valarray<REAL_T> > agedata) {
            this->agedata = agedata;
        }

        std::valarray<std::valarray<int> > GetAgeerr_type_a() const {
            return ageerr_type_a;
        }

        void SetAgeerr_type_a(std::valarray<std::valarray<int> > ageerr_type_a) {
            this->ageerr_type_a = ageerr_type_a;
        }

        std::valarray<std::valarray<int> > GetAgeerr_type_ms() const {
            return ageerr_type_ms;
        }

        void SetAgeerr_type_ms(std::valarray<std::valarray<int> > ageerr_type_ms) {
            this->ageerr_type_ms = ageerr_type_ms;
        }

        int GetAsk_detail() const {
            return ask_detail;
        }

        void SetAsk_detail(int ask_detail) {
            this->ask_detail = ask_detail;
        }

        std::vector<REAL_T> GetAzero_G() const {
            return azero_G;
        }

        void SetAzero_G(std::vector<REAL_T> azero_G) {
            this->azero_G = azero_G;
        }

        std::vector<REAL_T> GetAzero_seas() const {
            return azero_seas;
        }

        void SetAzero_seas(std::vector<REAL_T> azero_seas) {
            this->azero_seas = azero_seas;
        }

        int GetB() const {
            return b;
        }

        void SetB(int b) {
            this->b = b;
        }

        std::vector<REAL_T> GetBinwidth() const {
            return binwidth;
        }

        void SetBinwidth(std::vector<REAL_T> binwidth) {
            this->binwidth = binwidth;
        }

        REAL_T GetBinwidth2() const {
            return binwidth2;
        }

        void SetBinwidth2(REAL_T binwidth2) {
            this->binwidth2 = binwidth2;
        }

        int GetBio_yr() const {
            return bio_yr;
        }

        void SetBio_yr(int bio_yr) {
            this->bio_yr = bio_yr;
        }

        int GetBirthseas() const {
            return birthseas;
        }

        void SetBirthseas(int birthseas) {
            this->birthseas = birthseas;
        }

        int GetBlkparm() const {
            return blkparm;
        }

        void SetBlkparm(int blkparm) {
            this->blkparm = blkparm;
        }

        REAL_T GetBotbin() const {
            return botbin;
        }

        void SetBotbin(REAL_T botbin) {
            this->botbin = botbin;
        }

        int GetBurn_intvl() const {
            return burn_intvl;
        }

        void SetBurn_intvl(int burn_intvl) {
            this->burn_intvl = burn_intvl;
        }

        std::valarray<std::valarray<REAL_T> > GetCatch_bioT() const {
            return catch_bioT;
        }

        void SetCatch_bioT(std::valarray<std::valarray<REAL_T> > catch_bioT) {
            this->catch_bioT = catch_bioT;
        }

        std::valarray<std::valarray<REAL_T> > GetCatch_lambda() const {
            return catch_lambda;
        }

        void SetCatch_lambda(std::valarray<std::valarray<REAL_T> > catch_lambda) {
            this->catch_lambda = catch_lambda;
        }

        std::valarray<std::valarray<REAL_T> > GetCatch_ret_obs() const {
            return catch_ret_obs;
        }

        void SetCatch_ret_obs(std::valarray<std::valarray<REAL_T> > catch_ret_obs) {
            this->catch_ret_obs = catch_ret_obs;
        }

        std::valarray<std::valarray<REAL_T> > GetCatch_se() const {
            return catch_se;
        }

        void SetCatch_se(std::valarray<std::valarray<REAL_T> > catch_se) {
            this->catch_se = catch_se;
        }

        std::vector<REAL_T> GetCatch_se_rd() const {
            return catch_se_rd;
        }

        void SetCatch_se_rd(std::vector<REAL_T> catch_se_rd) {
            this->catch_se_rd = catch_se_rd;
        }

        std::vector<REAL_T> GetCatch_seas_area() const {
            return catch_seas_area;
        }

        void SetCatch_seas_area(std::vector<REAL_T> catch_seas_area) {
            this->catch_seas_area = catch_seas_area;
        }

        std::vector<REAL_T> GetCatchunits() const {
            return catchunits;
        }

        void SetCatchunits(std::vector<REAL_T> catchunits) {
            this->catchunits = catchunits;
        }

        std::vector<int> GetCr_errtype() const {
            return cr_errtype;
        }

        void SetCr_errtype(std::vector<int> cr_errtype) {
            this->cr_errtype = cr_errtype;
        }

        std::vector<int> GetCr_units() const {
            return cr_units;
        }

        void SetCr_units(std::vector<int> cr_units) {
            this->cr_units = cr_units;
        }

        std::valarray<std::valarray<int> > GetCr_units_rd() const {
            return cr_units_rd;
        }

        void SetCr_units_rd(std::valarray<std::valarray<int> > cr_units_rd) {
            this->cr_units_rd = cr_units_rd;
        }

        std::string GetCtlfilename() const {
            return ctlfilename;
        }

        void SetCtlfilename(std::string ctlfilename) {
            this->ctlfilename = ctlfilename;
        }

        std::vector<REAL_T> GetCurr_age1() const {
            return curr_age1;
        }

        void SetCurr_age1(std::vector<REAL_T> curr_age1) {
            this->curr_age1 = curr_age1;
        }

        std::vector<REAL_T> GetCurr_age2() const {
            return curr_age2;
        }

        void SetCurr_age2(std::vector<REAL_T> curr_age2) {
            this->curr_age2 = curr_age2;
        }

        int GetCustomMGenvsetup() const {
            return customMGenvsetup;
        }

        void SetCustomMGenvsetup(int customMGenvsetup) {
            this->customMGenvsetup = customMGenvsetup;
        }

        int GetCustomblocksetup() const {
            return customblocksetup;
        }

        void SetCustomblocksetup(int customblocksetup) {
            this->customblocksetup = customblocksetup;
        }

        int GetCustomblocksetup_MG() const {
            return customblocksetup_MG;
        }

        void SetCustomblocksetup_MG(int customblocksetup_MG) {
            this->customblocksetup_MG = customblocksetup_MG;
        }

        int GetCustomenvsetup() const {
            return customenvsetup;
        }

        void SetCustomenvsetup(int customenvsetup) {
            this->customenvsetup = customenvsetup;
        }

        std::valarray<std::valarray<REAL_T> > GetCv_disc() const {
            return cv_disc;
        }

        void SetCv_disc(std::valarray<std::valarray<REAL_T> > cv_disc) {
            this->cv_disc = cv_disc;
        }

        std::string GetDatfilename() const {
            return datfilename;
        }

        void SetDatfilename(std::string datfilename) {
            this->datfilename = datfilename;
        }

        int GetDepletion_basis() const {
            return depletion_basis;
        }

        void SetDepletion_basis(int depletion_basis) {
            this->depletion_basis = depletion_basis;
        }

        REAL_T GetDepletion_level() const {
            return depletion_level;
        }

        void SetDepletion_level(REAL_T depletion_level) {
            this->depletion_level = depletion_level;
        }

        int GetDid_MSY() const {
            return did_MSY;
        }

        void SetDid_MSY(int did_MSY) {
            this->did_MSY = did_MSY;
        }

        std::vector<int> GetDisc_errtype() const {
            return disc_errtype;
        }

        void SetDisc_errtype(std::vector<int> disc_errtype) {
            this->disc_errtype = disc_errtype;
        }

        std::vector<REAL_T> GetDisc_errtype_r() const {
            return disc_errtype_r;
        }

        void SetDisc_errtype_r(std::vector<REAL_T> disc_errtype_r) {
            this->disc_errtype_r = disc_errtype_r;
        }

        std::valarray<std::valarray<REAL_T> > GetDisc_lambda() const {
            return disc_lambda;
        }

        void SetDisc_lambda(std::valarray<std::valarray<REAL_T> > disc_lambda) {
            this->disc_lambda = disc_lambda;
        }

        std::vector<int> GetDisc_units() const {
            return disc_units;
        }

        void SetDisc_units(std::vector<int> disc_units) {
            this->disc_units = disc_units;
        }

        std::valarray<std::valarray<int> > GetDisc_units_rd() const {
            return disc_units_rd;
        }

        void SetDisc_units_rd(std::valarray<std::valarray<int> > disc_units_rd) {
            this->disc_units_rd = disc_units_rd;
        }

        std::valarray<std::valarray<REAL_T> > GetDiscdata() const {
            return discdata;
        }

        void SetDiscdata(std::valarray<std::valarray<REAL_T> > discdata) {
            this->discdata = discdata;
        }

        int GetDo_migr2() const {
            return do_migr2;
        }

        void SetDo_migr2(int do_migr2) {
            this->do_migr2 = do_migr2;
        }

        int GetDo_migration() const {
            return do_migration;
        }

        void SetDo_migration(int do_migration) {
            this->do_migration = do_migration;
        }

        int GetDo_once() const {
            return do_once;
        }

        void SetDo_once(int do_once) {
            this->do_once = do_once;
        }

        int GetDo_recdev() const {
            return do_recdev;
        }

        void SetDo_recdev(int do_recdev) {
            this->do_recdev = do_recdev;
        }

        int GetDocheckup() const {
            return docheckup;
        }

        void SetDocheckup(int docheckup) {
            this->docheckup = docheckup;
        }

        int GetDoit() const {
            return doit;
        }

        void SetDoit(int doit) {
            this->doit = doit;
        }

        std::vector<int> GetDolen() const {
            return dolen;
        }

        void SetDolen(std::vector<int> dolen) {
            this->dolen = dolen;
        }

        int GetDone_run() const {
            return done_run;
        }

        void SetDone_run(int done_run) {
            this->done_run = done_run;
        }

        REAL_T GetDummy_datum() const {
            return dummy_datum;
        }

        void SetDummy_datum(REAL_T dummy_datum) {
            this->dummy_datum = dummy_datum;
        }

        int GetDummy_phase() const {
            return dummy_phase;
        }

        void SetDummy_phase(int dummy_phase) {
            this->dummy_phase = dummy_phase;
        }

        int GetEndyr() const {
            return endyr;
        }

        void SetEndyr(int endyr) {
            this->endyr = endyr;
        }

        std::valarray<std::valarray<REAL_T> > GetEnv_data_RD() const {
            return env_data_RD;
        }

        void SetEnv_data_RD(std::valarray<std::valarray<REAL_T> > env_data_RD) {
            this->env_data_RD = env_data_RD;
        }

        std::valarray<std::valarray<REAL_T> > GetEnv_temp() const {
            return env_temp;
        }

        void SetEnv_temp(std::valarray<std::valarray<REAL_T> > env_temp) {
            this->env_temp = env_temp;
        }

        int GetEq_yr() const {
            return eq_yr;
        }

        void SetEq_yr(int eq_yr) {
            this->eq_yr = eq_yr;
        }

        int GetF() const {
            return f;
        }

        void SetF(int f) {
            this->f = f;
        }

        std::vector<REAL_T> GetFemfrac() const {
            return femfrac;
        }

        void SetFemfrac(std::vector<REAL_T> femfrac) {
            this->femfrac = femfrac;
        }

        int GetFid() const {
            return fid;
        }

        void SetFid(int fid) {
            this->fid = fid;
        }

        REAL_T GetFif() const {
            return fif;
        }

        void SetFif(REAL_T fif) {
            this->fif = fif;
        }

        int GetFim() const {
            return fim;
        }

        void SetFim(int fim) {
            this->fim = fim;
        }

        REAL_T GetFinal_conv() const {
            return final_conv;
        }

        void SetFinal_conv(REAL_T final_conv) {
            this->final_conv = final_conv;
        }

        int GetFini() const {
            return fini;
        }

        void SetFini(int fini) {
            this->fini = fini;
        }

        int GetFinish_starter() const {
            return finish_starter;
        }

        void SetFinish_starter(int finish_starter) {
            this->finish_starter = finish_starter;
        }

        std::vector<int> GetFirstBseas() const {
            return firstBseas;
        }

        void SetFirstBseas(std::vector<int> firstBseas) {
            this->firstBseas = firstBseas;
        }

        int GetFirst_catch_yr() const {
            return first_catch_yr;
        }

        void SetFirst_catch_yr(int first_catch_yr) {
            this->first_catch_yr = first_catch_yr;
        }

        int GetFirst_grow_age() const {
            return first_grow_age;
        }

        void SetFirst_grow_age(int first_grow_age) {
            this->first_grow_age = first_grow_age;
        }

        std::valarray<std::valarray<int> > GetFirst_grow_age2() const {
            return first_grow_age2;
        }

        void SetFirst_grow_age2(std::valarray<std::valarray<int> > first_grow_age2) {
            this->first_grow_age2 = first_grow_age2;
        }

        int GetFirstseas() const {
            return firstseas;
        }

        void SetFirstseas(int firstseas) {
            this->firstseas = firstseas;
        }

        int GetFirstselparm() const {
            return firstselparm;
        }

        void SetFirstselparm(int firstselparm) {
            this->firstselparm = firstselparm;
        }

        int GetFishery_on_off() const {
            return fishery_on_off;
        }

        void SetFishery_on_off(int fishery_on_off) {
            this->fishery_on_off = fishery_on_off;
        }

        std::vector<int> GetFleet_area() const {
            return fleet_area;
        }

        void SetFleet_area(std::vector<int> fleet_area) {
            this->fleet_area = fleet_area;
        }

        std::string GetFleetnameread() const {
            return fleetnameread;
        }

        void SetFleetnameread(std::string fleetnameread) {
            this->fleetnameread = fleetnameread;
        }

        int GetFloop() const {
            return floop;
        }

        void SetFloop(int floop) {
            this->floop = floop;
        }

        std::vector<REAL_T> GetFrac_ages() const {
            return frac_ages;
        }

        void SetFrac_ages(std::vector<REAL_T> frac_ages) {
            this->frac_ages = frac_ages;
        }

        REAL_T GetFracfemale() const {
            return fracfemale;
        }

        void SetFracfemale(REAL_T fracfemale) {
            this->fracfemale = fracfemale;
        }

        std::vector<REAL_T> GetFunc_conv() const {
            return func_conv;
        }

        void SetFunc_conv(std::vector<REAL_T> func_conv) {
            this->func_conv = func_conv;
        }

        std::vector<REAL_T> GetFunc_eval() const {
            return func_eval;
        }

        void SetFunc_eval(std::vector<REAL_T> func_eval) {
            this->func_eval = func_eval;
        }

        int GetG() const {
            return g;
        }

        void SetG(int g) {
            this->g = g;
        }

        std::valarray<std::valarray<int> > GetGen_a() const {
            return gen_a;
        }

        void SetGen_a(std::valarray<std::valarray<int> > gen_a) {
            this->gen_a = gen_a;
        }

        std::valarray<std::valarray<int> > GetGen_l() const {
            return gen_l;
        }

        void SetGen_l(std::valarray<std::valarray<int> > gen_l) {
            this->gen_l = gen_l;
        }

        std::valarray<std::valarray<int> > GetGen_ms() const {
            return gen_ms;
        }

        void SetGen_ms(std::valarray<std::valarray<int> > gen_ms) {
            this->gen_ms = gen_ms;
        }

        int GetGender() const {
            return gender;
        }

        void SetGender(int gender) {
            this->gender = gender;
        }

        int GetGg() const {
            return gg;
        }

        void SetGg(int gg) {
            this->gg = gg;
        }

        int GetGmorph() const {
            return gmorph;
        }

        void SetGmorph(int gmorph) {
            this->gmorph = gmorph;
        }

        int GetGp() const {
            return gp;
        }

        void SetGp(int gp) {
            this->gp = gp;
        }

        int GetGp2() const {
            return gp2;
        }

        void SetGp2(int gp2) {
            this->gp2 = gp2;
        }

        std::vector<int> GetHave_catch() const {
            return have_catch;
        }

        void SetHave_catch(std::vector<int> have_catch) {
            this->have_catch = have_catch;
        }

        std::valarray<std::valarray<REAL_T> > GetHave_data() const {
            return have_data;
        }

        void SetHave_data(std::valarray<std::valarray<REAL_T> > have_data) {
            this->have_data = have_data;
        }

        std::vector<REAL_T> GetHeader_a() const {
            return header_a;
        }

        void SetHeader_a(std::vector<REAL_T> header_a) {
            this->header_a = header_a;
        }

        std::vector<REAL_T> GetHeader_l() const {
            return header_l;
        }

        void SetHeader_l(std::vector<REAL_T> header_l) {
            this->header_l = header_l;
        }

        std::vector<REAL_T> GetHeader_ms() const {
            return header_ms;
        }

        void SetHeader_ms(std::vector<REAL_T> header_ms) {
            this->header_ms = header_ms;
        }

        int GetI() const {
            return i;
        }

        void SetI(int i) {
            this->i = i;
        }

        int GetIbin() const {
            return ibin;
        }

        void SetIbin(int ibin) {
            this->ibin = ibin;
        }

        int GetIbinsave() const {
            return ibinsave;
        }

        void SetIbinsave(int ibinsave) {
            this->ibinsave = ibinsave;
        }

        int GetIcycle() const {
            return icycle;
        }

        void SetIcycle(int icycle) {
            this->icycle = icycle;
        }

        int GetIn_superperiod() const {
            return in_superperiod;
        }

        void SetIn_superperiod(int in_superperiod) {
            this->in_superperiod = in_superperiod;
        }

        std::valarray<std::valarray<REAL_T> > GetIndexdata() const {
            return indexdata;
        }

        void SetIndexdata(std::valarray<std::valarray<REAL_T> > indexdata) {
            this->indexdata = indexdata;
        }

        std::vector<REAL_T> GetInit_F_CV() const {
            return init_F_CV;
        }

        void SetInit_F_CV(std::vector<REAL_T> init_F_CV) {
            this->init_F_CV = init_F_CV;
        }

        std::vector<REAL_T> GetInit_F_HI() const {
            return init_F_HI;
        }

        void SetInit_F_HI(std::vector<REAL_T> init_F_HI) {
            this->init_F_HI = init_F_HI;
        }

        std::vector<REAL_T> GetInit_F_LO() const {
            return init_F_LO;
        }

        void SetInit_F_LO(std::vector<REAL_T> init_F_LO) {
            this->init_F_LO = init_F_LO;
        }

        std::vector<int> GetInit_F_PH() const {
            return init_F_PH;
        }

        void SetInit_F_PH(std::vector<int> init_F_PH) {
            this->init_F_PH = init_F_PH;
        }

        std::vector<REAL_T> GetInit_F_PR() const {
            return init_F_PR;
        }

        void SetInit_F_PR(std::vector<REAL_T> init_F_PR) {
            this->init_F_PR = init_F_PR;
        }

        std::vector<REAL_T> GetInit_F_PRtype() const {
            return init_F_PRtype;
        }

        void SetInit_F_PRtype(std::vector<REAL_T> init_F_PRtype) {
            this->init_F_PRtype = init_F_PRtype;
        }

        std::vector<REAL_T> GetInit_F_RD() const {
            return init_F_RD;
        }

        void SetInit_F_RD(std::vector<REAL_T> init_F_RD) {
            this->init_F_RD = init_F_RD;
        }

        std::valarray<std::valarray<REAL_T> > GetInit_F_parm_1() const {
            return init_F_parm_1;
        }

        void SetInit_F_parm_1(std::valarray<std::valarray<REAL_T> > init_F_parm_1) {
            this->init_F_parm_1 = init_F_parm_1;
        }

        std::vector<REAL_T> GetInit_equ_lambda() const {
            return init_equ_lambda;
        }

        void SetInit_equ_lambda(std::vector<REAL_T> init_equ_lambda) {
            this->init_equ_lambda = init_equ_lambda;
        }

        int GetIobs() const {
            return iobs;
        }

        void SetIobs(int iobs) {
            this->iobs = iobs;
        }

        std::vector<int> GetIshadow() const {
            return ishadow;
        }

        void SetIshadow(std::vector<int> ishadow) {
            this->ishadow = ishadow;
        }

        int GetJ() const {
            return j;
        }

        void SetJ(int j) {
            this->j = j;
        }

        int GetJ1() const {
            return j1;
        }

        void SetJ1(int j1) {
            this->j1 = j1;
        }

        int GetJ2() const {
            return j2;
        }

        void SetJ2(int j2) {
            this->j2 = j2;
        }

        REAL_T GetJitter() const {
            return jitter;
        }

        void SetJitter(REAL_T jitter) {
            this->jitter = jitter;
        }

        std::vector<REAL_T> GetJunkvec() const {
            return junkvec;
        }

        void SetJunkvec(std::vector<REAL_T> junkvec) {
            this->junkvec = junkvec;
        }

        std::vector<REAL_T> GetJunkvec2() const {
            return junkvec2;
        }

        void SetJunkvec2(std::vector<REAL_T> junkvec2) {
            this->junkvec2 = junkvec2;
        }

        int GetK() const {
            return k;
        }

        void SetK(int k) {
            this->k = k;
        }

        int GetK1() const {
            return k1;
        }

        void SetK1(int k1) {
            this->k1 = k1;
        }

        int GetK2() const {
            return k2;
        }

        void SetK2(int k2) {
            this->k2 = k2;
        }

        int GetK3() const {
            return k3;
        }

        void SetK3(int k3) {
            this->k3 = k3;
        }

        std::vector<REAL_T> GetLen_bins() const {
            return len_bins;
        }

        void SetLen_bins(std::vector<REAL_T> len_bins) {
            this->len_bins = len_bins;
        }

        std::vector<REAL_T> GetLen_bins2() const {
            return len_bins2;
        }

        void SetLen_bins2(std::vector<REAL_T> len_bins2) {
            this->len_bins2 = len_bins2;
        }

        std::vector<REAL_T> GetLen_bins_dat() const {
            return len_bins_dat;
        }

        void SetLen_bins_dat(std::vector<REAL_T> len_bins_dat) {
            this->len_bins_dat = len_bins_dat;
        }

        std::vector<REAL_T> GetLen_bins_dat2() const {
            return len_bins_dat2;
        }

        void SetLen_bins_dat2(std::vector<REAL_T> len_bins_dat2) {
            this->len_bins_dat2 = len_bins_dat2;
        }

        std::vector<REAL_T> GetLen_bins_dat_m() const {
            return len_bins_dat_m;
        }

        void SetLen_bins_dat_m(std::vector<REAL_T> len_bins_dat_m) {
            this->len_bins_dat_m = len_bins_dat_m;
        }

        std::vector<REAL_T> GetLen_bins_m() const {
            return len_bins_m;
        }

        void SetLen_bins_m(std::vector<REAL_T> len_bins_m) {
            this->len_bins_m = len_bins_m;
        }

        std::vector<REAL_T> GetLen_bins_m2() const {
            return len_bins_m2;
        }

        void SetLen_bins_m2(std::vector<REAL_T> len_bins_m2) {
            this->len_bins_m2 = len_bins_m2;
        }

        std::vector<REAL_T> GetLen_bins_rd() const {
            return len_bins_rd;
        }

        void SetLen_bins_rd(std::vector<REAL_T> len_bins_rd) {
            this->len_bins_rd = len_bins_rd;
        }

        std::vector<REAL_T> GetLen_bins_sq() const {
            return len_bins_sq;
        }

        void SetLen_bins_sq(std::vector<REAL_T> len_bins_sq) {
            this->len_bins_sq = len_bins_sq;
        }

        std::valarray<std::valarray<REAL_T> > GetLendata() const {
            return lendata;
        }

        void SetLendata(std::valarray<std::valarray<REAL_T> > lendata) {
            this->lendata = lendata;
        }

        std::valarray<std::valarray<REAL_T> > GetLength_lambda() const {
            return length_lambda;
        }

        void SetLength_lambda(std::valarray<std::valarray<REAL_T> > length_lambda) {
            this->length_lambda = length_lambda;
        }

        std::vector<REAL_T> GetLin_grow1() const {
            return lin_grow1;
        }

        void SetLin_grow1(std::vector<REAL_T> lin_grow1) {
            this->lin_grow1 = lin_grow1;
        }

        std::vector<REAL_T> GetLin_grow2() const {
            return lin_grow2;
        }

        void SetLin_grow2(std::vector<REAL_T> lin_grow2) {
            this->lin_grow2 = lin_grow2;
        }

        int GetLoop() const {
            return loop;
        }

        void SetLoop(int loop) {
            this->loop = loop;
        }

        std::valarray<std::valarray<REAL_T> > GetMake_len_bin() const {
            return make_len_bin;
        }

        void SetMake_len_bin(std::valarray<std::valarray<REAL_T> > make_len_bin) {
            this->make_len_bin = make_len_bin;
        }

        int GetMakefishsel_yr() const {
            return makefishsel_yr;
        }

        void SetMakefishsel_yr(int makefishsel_yr) {
            this->makefishsel_yr = makefishsel_yr;
        }

        std::vector<REAL_T> GetMale_offset() const {
            return male_offset;
        }

        void SetMale_offset(std::vector<REAL_T> male_offset) {
            this->male_offset = male_offset;
        }

        REAL_T GetMaxL() const {
            return maxL;
        }

        void SetMaxL(REAL_T maxL) {
            this->maxL = maxL;
        }

        REAL_T GetMaxLread() const {
            return maxLread;
        }

        void SetMaxLread(REAL_T maxLread) {
            this->maxLread = maxLread;
        }

        REAL_T GetMax_harvest_rate() const {
            return max_harvest_rate;
        }

        void SetMax_harvest_rate(REAL_T max_harvest_rate) {
            this->max_harvest_rate = max_harvest_rate;
        }

        REAL_T GetMax_lambda_phase() const {
            return max_lambda_phase;
        }

        void SetMax_lambda_phase(REAL_T max_lambda_phase) {
            this->max_lambda_phase = max_lambda_phase;
        }

        int GetMax_phase() const {
            return max_phase;
        }

        void SetMax_phase(int max_phase) {
            this->max_phase = max_phase;
        }

        std::vector<REAL_T> GetMc_temp() const {
            return mc_temp;
        }

        void SetMc_temp(std::vector<REAL_T> mc_temp) {
            this->mc_temp = mc_temp;
        }

        int GetMceval_counter() const {
            return mceval_counter;
        }

        void SetMceval_counter(int mceval_counter) {
            this->mceval_counter = mceval_counter;
        }

        int GetMceval_header() const {
            return mceval_header;
        }

        void SetMceval_header(int mceval_header) {
            this->mceval_header = mceval_header;
        }

        REAL_T GetMcmcFlag() const {
            return mcmcFlag;
        }

        void SetMcmcFlag(REAL_T mcmcFlag) {
            this->mcmcFlag = mcmcFlag;
        }

        int GetMcmc_counter() const {
            return mcmc_counter;
        }

        void SetMcmc_counter(int mcmc_counter) {
            this->mcmc_counter = mcmc_counter;
        }

        std::vector<int> GetMgp_type() const {
            return mgp_type;
        }

        void SetMgp_type(std::vector<int> mgp_type) {
            this->mgp_type = mgp_type;
        }

        REAL_T GetMigr_firstage() const {
            return migr_firstage;
        }

        void SetMigr_firstage(REAL_T migr_firstage) {
            this->migr_firstage = migr_firstage;
        }

        std::valarray<std::valarray<REAL_T> > GetMigr_start() const {
            return migr_start;
        }

        void SetMigr_start(std::valarray<std::valarray<REAL_T> > migr_start) {
            this->migr_start = migr_start;
        }

        REAL_T GetMinL() const {
            return minL;
        }

        void SetMinL(REAL_T minL) {
            this->minL = minL;
        }

        REAL_T GetMinL_m() const {
            return minL_m;
        }

        void SetMinL_m(REAL_T minL_m) {
            this->minL_m = minL_m;
        }

        REAL_T GetMinLread() const {
            return minLread;
        }

        void SetMinLread(REAL_T minLread) {
            this->minLread = minLread;
        }

        REAL_T GetMin_comp() const {
            return min_comp;
        }

        void SetMin_comp(REAL_T min_comp) {
            this->min_comp = min_comp;
        }

        REAL_T GetMin_tail() const {
            return min_tail;
        }

        void SetMin_tail(REAL_T min_tail) {
            this->min_tail = min_tail;
        }

        std::valarray<std::valarray<int> > GetMkt_a() const {
            return mkt_a;
        }

        void SetMkt_a(std::valarray<std::valarray<int> > mkt_a) {
            this->mkt_a = mkt_a;
        }

        std::valarray<std::valarray<int> > GetMkt_l() const {
            return mkt_l;
        }

        void SetMkt_l(std::valarray<std::valarray<int> > mkt_l) {
            this->mkt_l = mkt_l;
        }

        std::valarray<std::valarray<int> > GetMkt_ms() const {
            return mkt_ms;
        }

        void SetMkt_ms(std::valarray<std::valarray<int> > mkt_ms) {
            this->mkt_ms = mkt_ms;
        }

        std::valarray<std::valarray<REAL_T> > GetMnwt_lambda() const {
            return mnwt_lambda;
        }

        void SetMnwt_lambda(std::valarray<std::valarray<REAL_T> > mnwt_lambda) {
            this->mnwt_lambda = mnwt_lambda;
        }

        std::valarray<std::valarray<REAL_T> > GetMnwtdata() const {
            return mnwtdata;
        }

        void SetMnwtdata(std::valarray<std::valarray<REAL_T> > mnwtdata) {
            this->mnwtdata = mnwtdata;
        }

        std::valarray<std::valarray<REAL_T> > GetMnwtdata1() const {
            return mnwtdata1;
        }

        void SetMnwtdata1(std::valarray<std::valarray<REAL_T> > mnwtdata1) {
            this->mnwtdata1 = mnwtdata1;
        }

        std::vector<int> GetMore_Fcast_input() const {
            return more_Fcast_input;
        }

        void SetMore_Fcast_input(std::vector<int> more_Fcast_input) {
            this->more_Fcast_input = more_Fcast_input;
        }

        std::valarray<std::valarray<REAL_T> > GetMove_def() const {
            return move_def;
        }

        void SetMove_def(std::valarray<std::valarray<REAL_T> > move_def) {
            this->move_def = move_def;
        }

        std::valarray<std::valarray<REAL_T> > GetMove_def2() const {
            return move_def2;
        }

        void SetMove_def2(std::valarray<std::valarray<REAL_T> > move_def2) {
            this->move_def2 = move_def2;
        }

        std::vector<REAL_T> GetMove_pattern() const {
            return move_pattern;
        }

        void SetMove_pattern(std::vector<REAL_T> move_pattern) {
            this->move_pattern = move_pattern;
        }

        int GetN_abins() const {
            return n_abins;
        }

        void SetN_abins(int n_abins) {
            this->n_abins = n_abins;
        }

        int GetN_abins1() const {
            return n_abins1;
        }

        void SetN_abins1(int n_abins1) {
            this->n_abins1 = n_abins1;
        }

        int GetN_abins2() const {
            return n_abins2;
        }

        void SetN_abins2(int n_abins2) {
            this->n_abins2 = n_abins2;
        }

        int GetNages() const {
            return nages;
        }

        void SetNages(int nages) {
            this->nages = nages;
        }

        REAL_T GetNatM_amax() const {
            return natM_amax;
        }

        void SetNatM_amax(REAL_T natM_amax) {
            this->natM_amax = natM_amax;
        }

        REAL_T GetNatM_amin() const {
            return natM_amin;
        }

        void SetNatM_amin(REAL_T natM_amin) {
            this->natM_amin = natM_amin;
        }

        int GetNatM_type() const {
            return natM_type;
        }

        void SetNatM_type(int natM_type) {
            this->natM_type = natM_type;
        }

        REAL_T GetNeglog19() const {
            return neglog19;
        }

        void SetNeglog19(REAL_T neglog19) {
            this->neglog19 = neglog19;
        }

        int GetNiter() const {
            return niter;
        }

        void SetNiter(int niter) {
            this->niter = niter;
        }

        int GetNlen_bin() const {
            return nlen_bin;
        }

        void SetNlen_bin(int nlen_bin) {
            this->nlen_bin = nlen_bin;
        }

        int GetNlen_bin2() const {
            return nlen_bin2;
        }

        void SetNlen_bin2(int nlen_bin2) {
            this->nlen_bin2 = nlen_bin2;
        }

        int GetNlen_binP() const {
            return nlen_binP;
        }

        void SetNlen_binP(int nlen_binP) {
            this->nlen_binP = nlen_binP;
        }

        int GetNlength() const {
            return nlength;
        }

        void SetNlength(int nlength) {
            this->nlength = nlength;
        }

        int GetNlength1() const {
            return nlength1;
        }

        void SetNlength1(int nlength1) {
            this->nlength1 = nlength1;
        }

        int GetNlength2() const {
            return nlength2;
        }

        void SetNlength2(int nlength2) {
            this->nlength2 = nlength2;
        }

        int GetNobs_cr() const {
            return nobs_cr;
        }

        void SetNobs_cr(int nobs_cr) {
            this->nobs_cr = nobs_cr;
        }

        int GetNobs_cr_rd() const {
            return nobs_cr_rd;
        }

        void SetNobs_cr_rd(int nobs_cr_rd) {
            this->nobs_cr_rd = nobs_cr_rd;
        }

        int GetNobs_disc() const {
            return nobs_disc;
        }

        void SetNobs_disc(int nobs_disc) {
            this->nobs_disc = nobs_disc;
        }

        int GetNobs_disc_rd() const {
            return nobs_disc_rd;
        }

        void SetNobs_disc_rd(int nobs_disc_rd) {
            this->nobs_disc_rd = nobs_disc_rd;
        }

        int GetNobs_mnwt() const {
            return nobs_mnwt;
        }

        void SetNobs_mnwt(int nobs_mnwt) {
            this->nobs_mnwt = nobs_mnwt;
        }

        int GetNobs_mnwt_rd() const {
            return nobs_mnwt_rd;
        }

        void SetNobs_mnwt_rd(int nobs_mnwt_rd) {
            this->nobs_mnwt_rd = nobs_mnwt_rd;
        }

        int GetNobs_ms_rd() const {
            return nobs_ms_rd;
        }

        void SetNobs_ms_rd(int nobs_ms_rd) {
            this->nobs_ms_rd = nobs_ms_rd;
        }

        int GetNobs_ms_tot() const {
            return nobs_ms_tot;
        }

        void SetNobs_ms_tot(int nobs_ms_tot) {
            this->nobs_ms_tot = nobs_ms_tot;
        }

        int GetNobsa_rd() const {
            return nobsa_rd;
        }

        void SetNobsa_rd(int nobsa_rd) {
            this->nobsa_rd = nobsa_rd;
        }

        int GetNobsl_rd() const {
            return nobsl_rd;
        }

        void SetNobsl_rd(int nobsl_rd) {
            this->nobsl_rd = nobsl_rd;
        }

        std::valarray<std::valarray<REAL_T> > GetNsamp_a() const {
            return nsamp_a;
        }

        void SetNsamp_a(std::valarray<std::valarray<REAL_T> > nsamp_a) {
            this->nsamp_a = nsamp_a;
        }

        std::valarray<std::valarray<REAL_T> > GetNsamp_a_read() const {
            return nsamp_a_read;
        }

        void SetNsamp_a_read(std::valarray<std::valarray<REAL_T> > nsamp_a_read) {
            this->nsamp_a_read = nsamp_a_read;
        }

        std::valarray<std::valarray<REAL_T> > GetNsamp_l() const {
            return nsamp_l;
        }

        void SetNsamp_l(std::valarray<std::valarray<REAL_T> > nsamp_l) {
            this->nsamp_l = nsamp_l;
        }

        std::valarray<std::valarray<REAL_T> > GetNsamp_l_read() const {
            return nsamp_l_read;
        }

        void SetNsamp_l_read(std::valarray<std::valarray<REAL_T> > nsamp_l_read) {
            this->nsamp_l_read = nsamp_l_read;
        }

        int GetNseas() const {
            return nseas;
        }

        void SetNseas(int nseas) {
            this->nseas = nseas;
        }

        std::vector<int> GetNyr_cr() const {
            return nyr_cr;
        }

        void SetNyr_cr(std::vector<int> nyr_cr) {
            this->nyr_cr = nyr_cr;
        }

        std::vector<int> GetNyr_disc() const {
            return nyr_disc;
        }

        void SetNyr_disc(std::vector<int> nyr_disc) {
            this->nyr_disc = nyr_disc;
        }

        std::vector<REAL_T> GetObs_a() const {
            return obs_a;
        }

        void SetObs_a(std::vector<REAL_T> obs_a) {
            this->obs_a = obs_a;
        }

        std::valarray<std::valarray<REAL_T> > GetObs_cr() const {
            return obs_cr;
        }

        void SetObs_cr(std::valarray<std::valarray<REAL_T> > obs_cr) {
            this->obs_cr = obs_cr;
        }

        std::valarray<std::valarray<REAL_T> > GetObs_disc() const {
            return obs_disc;
        }

        void SetObs_disc(std::valarray<std::valarray<REAL_T> > obs_disc) {
            this->obs_disc = obs_disc;
        }

        std::vector<REAL_T> GetObs_equ_catch() const {
            return obs_equ_catch;
        }

        void SetObs_equ_catch(std::vector<REAL_T> obs_equ_catch) {
            this->obs_equ_catch = obs_equ_catch;
        }

        std::vector<REAL_T> GetObs_l() const {
            return obs_l;
        }

        void SetObs_l(std::vector<REAL_T> obs_l) {
            this->obs_l = obs_l;
        }

        std::valarray<std::valarray<REAL_T> > GetObs_l_all() const {
            return obs_l_all;
        }

        void SetObs_l_all(std::valarray<std::valarray<REAL_T> > obs_l_all) {
            this->obs_l_all = obs_l_all;
        }

        std::vector<REAL_T> GetObs_ms() const {
            return obs_ms;
        }

        void SetObs_ms(std::vector<REAL_T> obs_ms) {
            this->obs_ms = obs_ms;
        }

        std::vector<REAL_T> GetObs_ms_n() const {
            return obs_ms_n;
        }

        void SetObs_ms_n(std::vector<REAL_T> obs_ms_n) {
            this->obs_ms_n = obs_ms_n;
        }

        std::vector<REAL_T> GetObs_ms_n_read() const {
            return obs_ms_n_read;
        }

        void SetObs_ms_n_read(std::vector<REAL_T> obs_ms_n_read) {
            this->obs_ms_n_read = obs_ms_n_read;
        }

        int GetP() const {
            return p;
        }

        void SetP(int p) {
            this->p = p;
        }

        int GetP1() const {
            return p1;
        }

        void SetP1(int p1) {
            this->p1 = p1;
        }

        int GetP2() const {
            return p2;
        }

        void SetP2(int p2) {
            this->p2 = p2;
        }

        std::vector<REAL_T> GetParm_dev_lambda() const {
            return parm_dev_lambda;
        }

        void SetParm_dev_lambda(std::vector<REAL_T> parm_dev_lambda) {
            this->parm_dev_lambda = parm_dev_lambda;
        }

        std::vector<REAL_T> GetParm_prior_lambda() const {
            return parm_prior_lambda;
        }

        void SetParm_prior_lambda(std::vector<REAL_T> parm_prior_lambda) {
            this->parm_prior_lambda = parm_prior_lambda;
        }

        std::valarray<std::valarray<int> > GetPfleetname() const {
            return pfleetname;
        }

        void SetPfleetname(std::valarray<std::valarray<int> > pfleetname) {
            this->pfleetname = pfleetname;
        }

        REAL_T GetPi() const {
            return pi;
        }

        void SetPi(REAL_T pi) {
            this->pi = pi;
        }

        int GetPop() const {
            return pop;
        }

        void SetPop(int pop) {
            this->pop = pop;
        }

        int GetProf_junk() const {
            return prof_junk;
        }

        void SetProf_junk(int prof_junk) {
            this->prof_junk = prof_junk;
        }

        std::vector<REAL_T> GetProf_var() const {
            return prof_var;
        }

        void SetProf_var(std::vector<REAL_T> prof_var) {
            this->prof_var = prof_var;
        }

        int GetProf_var_cnt() const {
            return prof_var_cnt;
        }

        void SetProf_var_cnt(int prof_var_cnt) {
            this->prof_var_cnt = prof_var_cnt;
        }

        std::vector<REAL_T> GetR_ages() const {
            return r_ages;
        }

        void SetR_ages(std::vector<REAL_T> r_ages) {
            this->r_ages = r_ages;
        }

        std::vector<REAL_T> GetR_years() const {
            return r_years;
        }

        void SetR_years(std::vector<REAL_T> r_years) {
            this->r_years = r_years;
        }

        int GetReadparfile() const {
            return readparfile;
        }

        void SetReadparfile(int readparfile) {
            this->readparfile = readparfile;
        }

        REAL_T GetRecdev_HI() const {
            return recdev_HI;
        }

        void SetRecdev_HI(REAL_T recdev_HI) {
            this->recdev_HI = recdev_HI;
        }

        REAL_T GetRecdev_LO() const {
            return recdev_LO;
        }

        void SetRecdev_LO(REAL_T recdev_LO) {
            this->recdev_LO = recdev_LO;
        }

        int GetRecdev_PH() const {
            return recdev_PH;
        }

        void SetRecdev_PH(int recdev_PH) {
            this->recdev_PH = recdev_PH;
        }

        int GetRecdev_PH_rd() const {
            return recdev_PH_rd;
        }

        void SetRecdev_PH_rd(int recdev_PH_rd) {
            this->recdev_PH_rd = recdev_PH_rd;
        }

        std::vector<REAL_T> GetRecdev_adj() const {
            return recdev_adj;
        }

        void SetRecdev_adj(std::vector<REAL_T> recdev_adj) {
            this->recdev_adj = recdev_adj;
        }

        int GetRecdev_adv() const {
            return recdev_adv;
        }

        void SetRecdev_adv(int recdev_adv) {
            this->recdev_adv = recdev_adv;
        }

        int GetRecdev_cycle() const {
            return recdev_cycle;
        }

        void SetRecdev_cycle(int recdev_cycle) {
            this->recdev_cycle = recdev_cycle;
        }

        std::vector<REAL_T> GetRecdev_cycle_HI() const {
            return recdev_cycle_HI;
        }

        void SetRecdev_cycle_HI(std::vector<REAL_T> recdev_cycle_HI) {
            this->recdev_cycle_HI = recdev_cycle_HI;
        }

        std::vector<REAL_T> GetRecdev_cycle_LO() const {
            return recdev_cycle_LO;
        }

        void SetRecdev_cycle_LO(std::vector<REAL_T> recdev_cycle_LO) {
            this->recdev_cycle_LO = recdev_cycle_LO;
        }

        std::vector<int> GetRecdev_cycle_PH() const {
            return recdev_cycle_PH;
        }

        void SetRecdev_cycle_PH(std::vector<int> recdev_cycle_PH) {
            this->recdev_cycle_PH = recdev_cycle_PH;
        }

        std::valarray<std::valarray<REAL_T> > GetRecdev_cycle_parm_RD() const {
            return recdev_cycle_parm_RD;
        }

        void SetRecdev_cycle_parm_RD(std::valarray<std::valarray<REAL_T> > recdev_cycle_parm_RD) {
            this->recdev_cycle_parm_RD = recdev_cycle_parm_RD;
        }

        int GetRecdev_do_early() const {
            return recdev_do_early;
        }

        void SetRecdev_do_early(int recdev_do_early) {
            this->recdev_do_early = recdev_do_early;
        }

        std::vector<int> GetRecdev_doit() const {
            return recdev_doit;
        }

        void SetRecdev_doit(std::vector<int> recdev_doit) {
            this->recdev_doit = recdev_doit;
        }

        int GetRecdev_early_PH() const {
            return recdev_early_PH;
        }

        void SetRecdev_early_PH(int recdev_early_PH) {
            this->recdev_early_PH = recdev_early_PH;
        }

        int GetRecdev_early_end() const {
            return recdev_early_end;
        }

        void SetRecdev_early_end(int recdev_early_end) {
            this->recdev_early_end = recdev_early_end;
        }

        int GetRecdev_early_start() const {
            return recdev_early_start;
        }

        void SetRecdev_early_start(int recdev_early_start) {
            this->recdev_early_start = recdev_early_start;
        }

        int GetRecdev_early_start_rd() const {
            return recdev_early_start_rd;
        }

        void SetRecdev_early_start_rd(int recdev_early_start_rd) {
            this->recdev_early_start_rd = recdev_early_start_rd;
        }

        int GetRecdev_end() const {
            return recdev_end;
        }

        void SetRecdev_end(int recdev_end) {
            this->recdev_end = recdev_end;
        }

        int GetRecdev_first() const {
            return recdev_first;
        }

        void SetRecdev_first(int recdev_first) {
            this->recdev_first = recdev_first;
        }

        std::valarray<std::valarray<REAL_T> > GetRecdev_input() const {
            return recdev_input;
        }

        void SetRecdev_input(std::valarray<std::valarray<REAL_T> > recdev_input) {
            this->recdev_input = recdev_input;
        }

        std::vector<REAL_T> GetRecdev_options() const {
            return recdev_options;
        }

        void SetRecdev_options(std::vector<REAL_T> recdev_options) {
            this->recdev_options = recdev_options;
        }

        std::vector<REAL_T> GetRecdev_options_rd() const {
            return recdev_options_rd;
        }

        void SetRecdev_options_rd(std::vector<REAL_T> recdev_options_rd) {
            this->recdev_options_rd = recdev_options_rd;
        }

        int GetRecdev_read() const {
            return recdev_read;
        }

        void SetRecdev_read(int recdev_read) {
            this->recdev_read = recdev_read;
        }

        int GetRecdev_start() const {
            return recdev_start;
        }

        void SetRecdev_start(int recdev_start) {
            this->recdev_start = recdev_start;
        }

        std::vector<int> GetRecr_dist_input() const {
            return recr_dist_input;
        }

        void SetRecr_dist_input(std::vector<int> recr_dist_input) {
            this->recr_dist_input = recr_dist_input;
        }

        int GetRecr_dist_inx() const {
            return recr_dist_inx;
        }

        void SetRecr_dist_inx(int recr_dist_inx) {
            this->recr_dist_inx = recr_dist_inx;
        }

        int GetRecr_dist_parms() const {
            return recr_dist_parms;
        }

        void SetRecr_dist_parms(int recr_dist_parms) {
            this->recr_dist_parms = recr_dist_parms;
        }

        std::vector<REAL_T> GetRecr_dist_pattern() const {
            return recr_dist_pattern;
        }

        void SetRecr_dist_pattern(std::vector<REAL_T> recr_dist_pattern) {
            this->recr_dist_pattern = recr_dist_pattern;
        }

        std::valarray<std::valarray<int> > GetRecr_dist_pattern_1() const {
            return recr_dist_pattern_1;
        }

        void SetRecr_dist_pattern_1(std::valarray<std::valarray<int> > recr_dist_pattern_1) {
            this->recr_dist_pattern_1 = recr_dist_pattern_1;
        }

        std::valarray<std::valarray<int> > GetRecr_dist_pattern_2() const {
            return recr_dist_pattern_2;
        }

        void SetRecr_dist_pattern_2(std::valarray<std::valarray<int> > recr_dist_pattern_2) {
            this->recr_dist_pattern_2 = recr_dist_pattern_2;
        }

        int GetRecr_dist_read() const {
            return recr_dist_read;
        }

        void SetRecr_dist_read(int recr_dist_read) {
            this->recr_dist_read = recr_dist_read;
        }

        std::vector<REAL_T> GetRecrdev_lambda() const {
            return recrdev_lambda;
        }

        void SetRecrdev_lambda(std::vector<REAL_T> recrdev_lambda) {
            this->recrdev_lambda = recrdev_lambda;
        }

        int GetReport_phase() const {
            return report_phase;
        }

        void SetReport_phase(int report_phase) {
            this->report_phase = report_phase;
        }

        int GetReportdetail() const {
            return reportdetail;
        }

        void SetReportdetail(int reportdetail) {
            this->reportdetail = reportdetail;
        }

        int GetRetro_yr() const {
            return retro_yr;
        }

        void SetRetro_yr(int retro_yr) {
            this->retro_yr = retro_yr;
        }

        int GetRundetail() const {
            return rundetail;
        }

        void SetRundetail(int rundetail) {
            this->rundetail = rundetail;
        }

        int GetRunnumber() const {
            return runnumber;
        }

        void SetRunnumber(int runnumber) {
            this->runnumber = runnumber;
        }

        int GetS() const {
            return s;
        }

        void SetS(int s) {
            this->s = s;
        }

        int GetS2() const {
            return s2;
        }

        void SetS2(int s2) {
            this->s2 = s2;
        }

        int GetS_off() const {
            return s_off;
        }

        void SetS_off(int s_off) {
            this->s_off = s_off;
        }

        std::valarray<std::valarray<REAL_T> > GetSave_G_parm() const {
            return save_G_parm;
        }

        void SetSave_G_parm(std::valarray<std::valarray<REAL_T> > save_G_parm) {
            this->save_G_parm = save_G_parm;
        }

        int GetSave_for_report() const {
            return save_for_report;
        }

        void SetSave_for_report(int save_for_report) {
            this->save_for_report = save_for_report;
        }

        int GetSave_gparm() const {
            return save_gparm;
        }

        void SetSave_gparm(int save_gparm) {
            this->save_gparm = save_gparm;
        }

        int GetSave_gparm_print() const {
            return save_gparm_print;
        }

        void SetSave_gparm_print(int save_gparm_print) {
            this->save_gparm_print = save_gparm_print;
        }

        std::valarray<std::valarray<REAL_T> > GetSave_seas_parm() const {
            return save_seas_parm;
        }

        void SetSave_seas_parm(std::valarray<std::valarray<REAL_T> > save_seas_parm) {
            this->save_seas_parm = save_seas_parm;
        }

        REAL_T GetSd_betweenmorph() const {
            return sd_betweenmorph;
        }

        void SetSd_betweenmorph(REAL_T sd_betweenmorph) {
            this->sd_betweenmorph = sd_betweenmorph;
        }

        std::valarray<std::valarray<REAL_T> > GetSd_disc() const {
            return sd_disc;
        }

        void SetSd_disc(std::valarray<std::valarray<REAL_T> > sd_disc) {
            this->sd_disc = sd_disc;
        }

        REAL_T GetSd_offset() const {
            return sd_offset;
        }

        void SetSd_offset(REAL_T sd_offset) {
            this->sd_offset = sd_offset;
        }

        REAL_T GetSd_ratio() const {
            return sd_ratio;
        }

        void SetSd_ratio(REAL_T sd_ratio) {
            this->sd_ratio = sd_ratio;
        }

        REAL_T GetSd_withinmorph() const {
            return sd_withinmorph;
        }

        void SetSd_withinmorph(REAL_T sd_withinmorph) {
            this->sd_withinmorph = sd_withinmorph;
        }

        std::valarray<std::valarray<REAL_T> > GetSe_cr_obs() const {
            return se_cr_obs;
        }

        void SetSe_cr_obs(std::valarray<std::valarray<REAL_T> > se_cr_obs) {
            this->se_cr_obs = se_cr_obs;
        }

        std::vector<REAL_T> GetSeasdur() const {
            return seasdur;
        }

        void SetSeasdur(std::vector<REAL_T> seasdur) {
            this->seasdur = seasdur;
        }

        std::vector<REAL_T> GetSeasdur2() const {
            return seasdur2;
        }

        void SetSeasdur2(std::vector<REAL_T> seasdur2) {
            this->seasdur2 = seasdur2;
        }

        std::valarray<std::valarray<REAL_T> > GetSelparm_1() const {
            return selparm_1;
        }

        void SetSelparm_1(std::valarray<std::valarray<REAL_T> > selparm_1) {
            this->selparm_1 = selparm_1;
        }

        std::vector<REAL_T> GetSelparm_CV() const {
            return selparm_CV;
        }

        void SetSelparm_CV(std::vector<REAL_T> selparm_CV) {
            this->selparm_CV = selparm_CV;
        }

        std::vector<REAL_T> GetSelparm_HI() const {
            return selparm_HI;
        }

        void SetSelparm_HI(std::vector<REAL_T> selparm_HI) {
            this->selparm_HI = selparm_HI;
        }

        std::vector<REAL_T> GetSelparm_LO() const {
            return selparm_LO;
        }

        void SetSelparm_LO(std::vector<REAL_T> selparm_LO) {
            this->selparm_LO = selparm_LO;
        }

        std::vector<int> GetSelparm_PH() const {
            return selparm_PH;
        }

        void SetSelparm_PH(std::vector<int> selparm_PH) {
            this->selparm_PH = selparm_PH;
        }

        std::vector<REAL_T> GetSelparm_PR() const {
            return selparm_PR;
        }

        void SetSelparm_PR(std::vector<REAL_T> selparm_PR) {
            this->selparm_PR = selparm_PR;
        }

        std::vector<REAL_T> GetSelparm_PRtype() const {
            return selparm_PRtype;
        }

        void SetSelparm_PRtype(std::vector<REAL_T> selparm_PRtype) {
            this->selparm_PRtype = selparm_PRtype;
        }

        std::vector<REAL_T> GetSelparm_RD() const {
            return selparm_RD;
        }

        void SetSelparm_RD(std::vector<REAL_T> selparm_RD) {
            this->selparm_RD = selparm_RD;
        }

        int GetSelparm_adjust_method() const {
            return selparm_adjust_method;
        }

        void SetSelparm_adjust_method(int selparm_adjust_method) {
            this->selparm_adjust_method = selparm_adjust_method;
        }

        std::valarray<std::valarray<REAL_T> > GetSelparm_blk_1() const {
            return selparm_blk_1;
        }

        void SetSelparm_blk_1(std::valarray<std::valarray<REAL_T> > selparm_blk_1) {
            this->selparm_blk_1 = selparm_blk_1;
        }

        REAL_T GetSelparm_dev_PH() const {
            return selparm_dev_PH;
        }

        void SetSelparm_dev_PH(REAL_T selparm_dev_PH) {
            this->selparm_dev_PH = selparm_dev_PH;
        }

        std::vector<int> GetSelparm_dev_maxyr() const {
            return selparm_dev_maxyr;
        }

        void SetSelparm_dev_maxyr(std::vector<int> selparm_dev_maxyr) {
            this->selparm_dev_maxyr = selparm_dev_maxyr;
        }

        std::vector<int> GetSelparm_dev_minyr() const {
            return selparm_dev_minyr;
        }

        void SetSelparm_dev_minyr(std::vector<int> selparm_dev_minyr) {
            this->selparm_dev_minyr = selparm_dev_minyr;
        }

        std::vector<int> GetSelparm_dev_select() const {
            return selparm_dev_select;
        }

        void SetSelparm_dev_select(std::vector<int> selparm_dev_select) {
            this->selparm_dev_select = selparm_dev_select;
        }

        std::vector<REAL_T> GetSelparm_dev_stddev() const {
            return selparm_dev_stddev;
        }

        void SetSelparm_dev_stddev(std::vector<REAL_T> selparm_dev_stddev) {
            this->selparm_dev_stddev = selparm_dev_stddev;
        }

        std::vector<int> GetSelparm_dev_type() const {
            return selparm_dev_type;
        }

        void SetSelparm_dev_type(std::vector<int> selparm_dev_type) {
            this->selparm_dev_type = selparm_dev_type;
        }

        std::vector<int> GetSelparm_env() const {
            return selparm_env;
        }

        void SetSelparm_env(std::vector<int> selparm_env) {
            this->selparm_env = selparm_env;
        }

        std::valarray<std::valarray<REAL_T> > GetSelparm_env_1() const {
            return selparm_env_1;
        }

        void SetSelparm_env_1(std::valarray<std::valarray<REAL_T> > selparm_env_1) {
            this->selparm_env_1 = selparm_env_1;
        }

        std::vector<int> GetSelparm_envtype() const {
            return selparm_envtype;
        }

        void SetSelparm_envtype(std::vector<int> selparm_envtype) {
            this->selparm_envtype = selparm_envtype;
        }

        std::vector<int> GetSelparm_envuse() const {
            return selparm_envuse;
        }

        void SetSelparm_envuse(std::vector<int> selparm_envuse) {
            this->selparm_envuse = selparm_envuse;
        }

        std::valarray<std::valarray<REAL_T> > GetSelparm_trend_1() const {
            return selparm_trend_1;
        }

        void SetSelparm_trend_1(std::valarray<std::valarray<REAL_T> > selparm_trend_1) {
            this->selparm_trend_1 = selparm_trend_1;
        }

        std::vector<int> GetSelparm_trend_point() const {
            return selparm_trend_point;
        }

        void SetSelparm_trend_point(std::vector<int> selparm_trend_point) {
            this->selparm_trend_point = selparm_trend_point;
        }

        std::vector<int> GetSelparm_trend_rev() const {
            return selparm_trend_rev;
        }

        void SetSelparm_trend_rev(std::vector<int> selparm_trend_rev) {
            this->selparm_trend_rev = selparm_trend_rev;
        }

        std::vector<int> GetSelparm_trend_rev_1() const {
            return selparm_trend_rev_1;
        }

        void SetSelparm_trend_rev_1(std::vector<int> selparm_trend_rev_1) {
            this->selparm_trend_rev_1 = selparm_trend_rev_1;
        }

        std::valarray<std::valarray<int> > GetSeltype() const {
            return seltype;
        }

        void SetSeltype(std::valarray<std::valarray<int> > seltype) {
            this->seltype = seltype;
        }

        std::vector<int> GetSeltype_Nparam() const {
            return seltype_Nparam;
        }

        void SetSeltype_Nparam(std::vector<int> seltype_Nparam) {
            this->seltype_Nparam = seltype_Nparam;
        }

        std::vector<REAL_T> GetShadow() const {
            return shadow;
        }

        void SetShadow(std::vector<REAL_T> shadow) {
            this->shadow = shadow;
        }

        int GetShow_MSY() const {
            return show_MSY;
        }

        void SetShow_MSY(int show_MSY) {
            this->show_MSY = show_MSY;
        }

        std::valarray<std::valarray<REAL_T> > GetSizeage_lambda() const {
            return sizeage_lambda;
        }

        void SetSizeage_lambda(std::valarray<std::valarray<REAL_T> > sizeage_lambda) {
            this->sizeage_lambda = sizeage_lambda;
        }

        std::valarray<std::valarray<REAL_T> > GetSizeagedata() const {
            return sizeagedata;
        }

        void SetSizeagedata(std::valarray<std::valarray<REAL_T> > sizeagedata) {
            this->sizeagedata = sizeagedata;
        }

        int GetSmid() const {
            return smid;
        }

        void SetSmid(int smid) {
            this->smid = smid;
        }

        int GetSpawn_seas() const {
            return spawn_seas;
        }

        void SetSpawn_seas(int spawn_seas) {
            this->spawn_seas = spawn_seas;
        }

        REAL_T GetStartbin() const {
            return startbin;
        }

        void SetStartbin(REAL_T startbin) {
            this->startbin = startbin;
        }

        int GetStyr() const {
            return styr;
        }

        void SetStyr(int styr) {
            this->styr = styr;
        }

        std::vector<REAL_T> GetSubmorphdist() const {
            return submorphdist;
        }

        void SetSubmorphdist(std::vector<REAL_T> submorphdist) {
            this->submorphdist = submorphdist;
        }

        REAL_T GetSumseas() const {
            return sumseas;
        }

        void SetSumseas(REAL_T sumseas) {
            this->sumseas = sumseas;
        }

        std::valarray<std::valarray<int> > GetSuprper_a1() const {
            return suprper_a1;
        }

        void SetSuprper_a1(std::valarray<std::valarray<int> > suprper_a1) {
            this->suprper_a1 = suprper_a1;
        }

        std::valarray<std::valarray<int> > GetSuprper_a2() const {
            return suprper_a2;
        }

        void SetSuprper_a2(std::valarray<std::valarray<int> > suprper_a2) {
            this->suprper_a2 = suprper_a2;
        }

        std::valarray<std::valarray<int> > GetSuprper_cr1() const {
            return suprper_cr1;
        }

        void SetSuprper_cr1(std::valarray<std::valarray<int> > suprper_cr1) {
            this->suprper_cr1 = suprper_cr1;
        }

        std::valarray<std::valarray<int> > GetSuprper_cr2() const {
            return suprper_cr2;
        }

        void SetSuprper_cr2(std::valarray<std::valarray<int> > suprper_cr2) {
            this->suprper_cr2 = suprper_cr2;
        }

        std::valarray<std::valarray<int> > GetSuprper_disc1() const {
            return suprper_disc1;
        }

        void SetSuprper_disc1(std::valarray<std::valarray<int> > suprper_disc1) {
            this->suprper_disc1 = suprper_disc1;
        }

        std::valarray<std::valarray<int> > GetSuprper_disc2() const {
            return suprper_disc2;
        }

        void SetSuprper_disc2(std::valarray<std::valarray<int> > suprper_disc2) {
            this->suprper_disc2 = suprper_disc2;
        }

        std::valarray<std::valarray<int> > GetSuprper_end_SzFreq() const {
            return suprper_end_SzFreq;
        }

        void SetSuprper_end_SzFreq(std::valarray<std::valarray<int> > suprper_end_SzFreq) {
            this->suprper_end_SzFreq = suprper_end_SzFreq;
        }

        std::valarray<std::valarray<int> > GetSuprper_l1() const {
            return suprper_l1;
        }

        void SetSuprper_l1(std::valarray<std::valarray<int> > suprper_l1) {
            this->suprper_l1 = suprper_l1;
        }

        std::valarray<std::valarray<int> > GetSuprper_l2() const {
            return suprper_l2;
        }

        void SetSuprper_l2(std::valarray<std::valarray<int> > suprper_l2) {
            this->suprper_l2 = suprper_l2;
        }

        std::valarray<std::valarray<int> > GetSuprper_ms1() const {
            return suprper_ms1;
        }

        void SetSuprper_ms1(std::valarray<std::valarray<int> > suprper_ms1) {
            this->suprper_ms1 = suprper_ms1;
        }

        std::valarray<std::valarray<int> > GetSuprper_ms2() const {
            return suprper_ms2;
        }

        void SetSuprper_ms2(std::valarray<std::valarray<int> > suprper_ms2) {
            this->suprper_ms2 = suprper_ms2;
        }

        std::valarray<std::valarray<int> > GetSuprper_start_SzFreq() const {
            return suprper_start_SzFreq;
        }

        void SetSuprper_start_SzFreq(std::valarray<std::valarray<int> > suprper_start_SzFreq) {
            this->suprper_start_SzFreq = suprper_start_SzFreq;
        }

        std::valarray<std::valarray<REAL_T> > GetSurv_lambda() const {
            return surv_lambda;
        }

        void SetSurv_lambda(std::valarray<std::valarray<REAL_T> > surv_lambda) {
            this->surv_lambda = surv_lambda;
        }

        std::vector<REAL_T> GetSurveytime() const {
            return surveytime;
        }

        void SetSurveytime(std::vector<REAL_T> surveytime) {
            this->surveytime = surveytime;
        }

        std::vector<int> GetSx() const {
            return sx;
        }

        void SetSx(std::vector<int> sx) {
            this->sx = sx;
        }

        int GetT() const {
            return t;
        }

        void SetT(int t) {
            this->t = t;
        }

        int GetT_base() const {
            return t_base;
        }

        void SetT_base(int t_base) {
            this->t_base = t_base;
        }

        std::vector<REAL_T> GetTails_a() const {
            return tails_a;
        }

        void SetTails_a(std::vector<REAL_T> tails_a) {
            this->tails_a = tails_a;
        }

        std::vector<REAL_T> GetTails_l() const {
            return tails_l;
        }

        void SetTails_l(std::vector<REAL_T> tails_l) {
            this->tails_l = tails_l;
        }

        std::vector<int> GetTails_w() const {
            return tails_w;
        }

        void SetTails_w(std::vector<int> tails_w) {
            this->tails_w = tails_w;
        }

        REAL_T GetTemp() const {
            return temp;
        }

        void SetTemp(REAL_T temp) {
            this->temp = temp;
        }

        REAL_T GetTemp1() const {
            return temp1;
        }

        void SetTemp1(REAL_T temp1) {
            this->temp1 = temp1;
        }

        std::vector<REAL_T> GetTempvec4() const {
            return tempvec4;
        }

        void SetTempvec4(std::vector<REAL_T> tempvec4) {
            this->tempvec4 = tempvec4;
        }

        std::vector<REAL_T> GetTempvec5() const {
            return tempvec5;
        }

        void SetTempvec5(std::vector<REAL_T> tempvec5) {
            this->tempvec5 = tempvec5;
        }

        int GetThin_intvl() const {
            return thin_intvl;
        }

        void SetThin_intvl(int thin_intvl) {
            this->thin_intvl = thin_intvl;
        }

        std::valarray<std::valarray<int> > GetTime_vary_MG() const {
            return time_vary_MG;
        }

        void SetTime_vary_MG(std::valarray<std::valarray<int> > time_vary_MG) {
            this->time_vary_MG = time_vary_MG;
        }

        std::valarray<std::valarray<int> > GetTime_vary_makefishsel() const {
            return time_vary_makefishsel;
        }

        void SetTime_vary_makefishsel(std::valarray<std::valarray<int> > time_vary_makefishsel) {
            this->time_vary_makefishsel = time_vary_makefishsel;
        }

        std::valarray<std::valarray<int> > GetTime_vary_sel() const {
            return time_vary_sel;
        }

        void SetTime_vary_sel(std::valarray<std::valarray<int> > time_vary_sel) {
            this->time_vary_sel = time_vary_sel;
        }

        int GetTloop() const {
            return tloop;
        }

        void SetTloop(int tloop) {
            this->tloop = tloop;
        }

        REAL_T GetTopbin() const {
            return topbin;
        }

        void SetTopbin(REAL_T topbin) {
            this->topbin = topbin;
        }

        std::vector<REAL_T> GetTotcat() const {
            return totcat;
        }

        void SetTotcat(std::vector<REAL_T> totcat) {
            this->totcat = totcat;
        }

        std::valarray<std::valarray<REAL_T> > GetTotcatch_byarea() const {
            return totcatch_byarea;
        }

        void SetTotcatch_byarea(std::valarray<std::valarray<REAL_T> > totcatch_byarea) {
            this->totcatch_byarea = totcatch_byarea;
        }

        std::valarray<std::valarray<int> > GetUse_Lbin_filter() const {
            return use_Lbin_filter;
        }

        void SetUse_Lbin_filter(std::valarray<std::valarray<int> > use_Lbin_filter) {
            this->use_Lbin_filter = use_Lbin_filter;
        }

        std::vector<int> GetUse_morph() const {
            return use_morph;
        }

        void SetUse_morph(std::vector<int> use_morph) {
            this->use_morph = use_morph;
        }

        std::valarray<std::valarray<int> > GetUse_ms() const {
            return use_ms;
        }

        void SetUse_ms(std::valarray<std::valarray<int> > use_ms) {
            this->use_ms = use_ms;
        }

        std::valarray<std::valarray<REAL_T> > GetVar_adjust() const {
            return var_adjust;
        }

        void SetVar_adjust(std::valarray<std::valarray<REAL_T> > var_adjust) {
            this->var_adjust = var_adjust;
        }

        std::valarray<std::valarray<REAL_T> > GetVar_adjust1() const {
            return var_adjust1;
        }

        void SetVar_adjust1(std::valarray<std::valarray<REAL_T> > var_adjust1) {
            this->var_adjust1 = var_adjust1;
        }

        std::valarray<std::valarray<REAL_T> > GetVul_bio() const {
            return vul_bio;
        }

        void SetVul_bio(std::valarray<std::valarray<REAL_T> > vul_bio) {
            this->vul_bio = vul_bio;
        }

        int GetY() const {
            return y;
        }

        void SetY(int y) {
            this->y = y;
        }

        int GetY2() const {
            return y2;
        }

        void SetY2(int y2) {
            this->y2 = y2;
        }

        std::vector<int> GetYears() const {
            return years;
        }

        void SetYears(std::vector<int> years) {
            this->years = years;
        }

        std::valarray<std::valarray<int> > GetYr_a_t() const {
            return yr_a_t;
        }

        void SetYr_a_t(std::valarray<std::valarray<int> > yr_a_t) {
            this->yr_a_t = yr_a_t;
        }

        std::valarray<std::valarray<int> > GetYr_a_y() const {
            return yr_a_y;
        }

        void SetYr_a_y(std::valarray<std::valarray<int> > yr_a_y) {
            this->yr_a_y = yr_a_y;
        }

        std::valarray<std::valarray<int> > GetYr_cr() const {
            return yr_cr;
        }

        void SetYr_cr(std::valarray<std::valarray<int> > yr_cr) {
            this->yr_cr = yr_cr;
        }

        std::valarray<std::valarray<REAL_T> > GetYr_cr2() const {
            return yr_cr2;
        }

        void SetYr_cr2(std::valarray<std::valarray<REAL_T> > yr_cr2) {
            this->yr_cr2 = yr_cr2;
        }

        std::valarray<std::valarray<int> > GetYr_cr_s() const {
            return yr_cr_s;
        }

        void SetYr_cr_s(std::valarray<std::valarray<int> > yr_cr_s) {
            this->yr_cr_s = yr_cr_s;
        }

        std::valarray<std::valarray<int> > GetYr_cr_use() const {
            return yr_cr_use;
        }

        void SetYr_cr_use(std::valarray<std::valarray<int> > yr_cr_use) {
            this->yr_cr_use = yr_cr_use;
        }

        std::valarray<std::valarray<int> > GetYr_cr_y() const {
            return yr_cr_y;
        }

        void SetYr_cr_y(std::valarray<std::valarray<int> > yr_cr_y) {
            this->yr_cr_y = yr_cr_y;
        }

        std::valarray<std::valarray<int> > GetYr_disc() const {
            return yr_disc;
        }

        void SetYr_disc(std::valarray<std::valarray<int> > yr_disc) {
            this->yr_disc = yr_disc;
        }

        std::valarray<std::valarray<REAL_T> > GetYr_disc2() const {
            return yr_disc2;
        }

        void SetYr_disc2(std::valarray<std::valarray<REAL_T> > yr_disc2) {
            this->yr_disc2 = yr_disc2;
        }

        std::valarray<std::valarray<int> > GetYr_disc_s() const {
            return yr_disc_s;
        }

        void SetYr_disc_s(std::valarray<std::valarray<int> > yr_disc_s) {
            this->yr_disc_s = yr_disc_s;
        }

        std::valarray<std::valarray<int> > GetYr_disc_use() const {
            return yr_disc_use;
        }

        void SetYr_disc_use(std::valarray<std::valarray<int> > yr_disc_use) {
            this->yr_disc_use = yr_disc_use;
        }

        std::valarray<std::valarray<int> > GetYr_disc_y() const {
            return yr_disc_y;
        }

        void SetYr_disc_y(std::valarray<std::valarray<int> > yr_disc_y) {
            this->yr_disc_y = yr_disc_y;
        }

        std::valarray<std::valarray<int> > GetYr_l_t() const {
            return yr_l_t;
        }

        void SetYr_l_t(std::valarray<std::valarray<int> > yr_l_t) {
            this->yr_l_t = yr_l_t;
        }

        std::valarray<std::valarray<int> > GetYr_l_y() const {
            return yr_l_y;
        }

        void SetYr_l_y(std::valarray<std::valarray<int> > yr_l_y) {
            this->yr_l_y = yr_l_y;
        }

        std::vector<REAL_T> GetYr_mnwt2() const {
            return yr_mnwt2;
        }

        void SetYr_mnwt2(std::vector<REAL_T> yr_mnwt2) {
            this->yr_mnwt2 = yr_mnwt2;
        }

        std::valarray<std::valarray<int> > GetYr_ms_t() const {
            return yr_ms_t;
        }

        void SetYr_ms_t(std::valarray<std::valarray<int> > yr_ms_t) {
            this->yr_ms_t = yr_ms_t;
        }

        std::valarray<std::valarray<int> > GetYr_ms_y() const {
            return yr_ms_y;
        }

        void SetYr_ms_y(std::valarray<std::valarray<int> > yr_ms_y) {
            this->yr_ms_y = yr_ms_y;
        }

        int GetYz() const {
            return yz;
        }

        void SetYz(int yz) {
            this->yz = yz;
        }

        int GetZ() const {
            return z;
        }

        void SetZ(int z) {
            this->z = z;
        }

        int GetZ1() const {
            return z1;
        }

        void SetZ1(int z1) {
            this->z1 = z1;
        }

        int GetZ2() const {
            return z2;
        }

        void SetZ2(int z2) {
            this->z2 = z2;
        }





    };

    template<class REAL_T, class EVAL_T = REAL_T>
    class StockSynthesisModel : public ss::ModelBase<REAL_T, EVAL_T> {
        std::vector<int> integer_control_flags;
        std::vector<REAL_T> double_control_flags;
        EVAL_T dummy_parm;
        std::vector<EVAL_T> MGparm;
        std::valarray<std::valarray<EVAL_T> > MGparm_trend;
        std::valarray<std::valarray<EVAL_T> > MGparm_block_val;
        std::valarray<std::valarray<EVAL_T> > MGparm_dev;
        std::valarray<std::valarray<EVAL_T> > MGparm_dev_rwalk;
        std::vector<EVAL_T> L_inf;
        std::vector<EVAL_T> Lmax_temp;
        std::vector<EVAL_T> VBK;
        std::vector<EVAL_T> Richards;
        std::vector<EVAL_T> Lmin;
        std::vector<EVAL_T> Lmin_last;
        std::valarray<std::valarray<EVAL_T> > natMparms;
        std::vector<EVAL_T> natM;
        std::vector<EVAL_T> surv1;
        std::vector<EVAL_T> surv2;
        std::vector<EVAL_T> CVLmin;
        std::vector<EVAL_T> CVLmax;
        std::vector<EVAL_T> CV_const;
        std::valarray<std::valarray<EVAL_T> > mgp_save;
        std::vector<EVAL_T> mgp_adj;
        std::valarray<std::valarray<EVAL_T> > Cohort_Growth;
        std::vector<EVAL_T> Cohort_Lmin;
        std::vector<EVAL_T> VBK_seas;
        std::valarray<std::valarray<EVAL_T> > wtlen_seas;
        std::vector<EVAL_T> wtlen_p;
        std::vector<EVAL_T> wt_len;
        std::valarray<std::valarray<EVAL_T> > wt_len2;
        std::valarray<std::valarray<EVAL_T> > wt_len2_sq;
        std::valarray<std::valarray<EVAL_T> > wt_len_low;
        std::valarray<std::valarray<EVAL_T> > wt_len_fd;
        std::vector<EVAL_T> mat_len;
        std::vector<EVAL_T> fec_len;
        std::vector<EVAL_T> mat_len_wt;
        std::vector<EVAL_T> mat_age;
        std::valarray<std::valarray<EVAL_T> > Hermaphro_val;
        std::vector<EVAL_T> age_age;
        std::vector<EVAL_T> age_err;
        std::vector<EVAL_T> ALK;
        std::valarray<std::valarray<EVAL_T> > exp_AL;
        std::vector<EVAL_T> Sd_Size_within;
        std::vector<EVAL_T> Sd_Size_between;
        std::vector<EVAL_T> Ave_Size;
        std::valarray<std::valarray<EVAL_T> > CV_G;
        std::vector<EVAL_T> Save_Wt_Age;
        std::vector<EVAL_T> Wt_Age_beg;
        std::vector<EVAL_T> Wt_Age_mid;
        std::vector<EVAL_T> migrrate;
        std::vector<EVAL_T> recr_dist;
        std::vector<EVAL_T> SR_parm;
        EVAL_T two_sigmaRsq;
        EVAL_T half_sigmaRsq;
        EVAL_T sigmaR;
        EVAL_T rho;
        std::vector<EVAL_T> biasadj;
        std::vector<EVAL_T> biasadj_full;
        EVAL_T sd_offset_rec;
        std::vector<EVAL_T> recdev_cycle_parm;
        std::vector<EVAL_T> recdev_early;
        std::vector<EVAL_T> recdev1;
        std::vector<EVAL_T> recdev2;
        std::vector<EVAL_T> Fcast_recruitments;
        std::vector<EVAL_T> recdev;
        std::vector<EVAL_T> Fcast_impl_error;
        EVAL_T SPB_current;
        EVAL_T SPB_vir_LH;
        EVAL_T Recr_virgin;
        EVAL_T SPB_virgin;
        EVAL_T SPR_unf;
        EVAL_T SPR_trial;
        std::vector<EVAL_T> SPB_pop_gp;
        std::vector<EVAL_T> SPB_yr;
        std::vector<EVAL_T> MaleSPB;
        std::valarray<std::valarray<EVAL_T> > SPB_equil_pop_gp;
        std::valarray<std::valarray<EVAL_T> > MaleSPB_equil_pop_gp;
        EVAL_T SPB_equil;
        EVAL_T Recruits;
        std::valarray<std::valarray<EVAL_T> > Recr;
        std::vector<EVAL_T> exp_rec;
        std::vector<EVAL_T> pred_rec;
        std::vector<EVAL_T> use_rec;
        std::valarray<std::valarray<EVAL_T> > Nmid;
        std::valarray<std::valarray<EVAL_T> > Nsurv;
        std::vector<EVAL_T> natage_temp;
        std::vector<EVAL_T> agetemp;
        std::vector<EVAL_T> Save_PopLen;
        std::vector<EVAL_T> Save_PopAge;
        EVAL_T ave_age;
        std::vector<EVAL_T> init_F;
        std::vector<EVAL_T> est_equ_catch;
        std::vector<EVAL_T> natage;
        std::vector<EVAL_T> catage;
        std::vector<EVAL_T> equ_catage;
        std::vector<EVAL_T> equ_numbers;
        std::vector<EVAL_T> equ_Z;
        std::valarray<std::valarray<EVAL_T> > catage_tot;
        std::valarray<std::valarray<EVAL_T> > Hrate;
        std::vector<EVAL_T> catch_fleet;
        std::valarray<std::valarray<EVAL_T> > equ_catch_fleet;
        std::valarray<std::valarray<EVAL_T> > fec;
        std::valarray<std::valarray<EVAL_T> > virg_fec;
        EVAL_T fish_bio;
        EVAL_T fish_bio_r;
        EVAL_T fish_bio_e;
        EVAL_T fish_num_e;
        EVAL_T fish_num;
        EVAL_T fish_num_r;
        EVAL_T vbio;
        EVAL_T totbio;
        EVAL_T smrybio;
        EVAL_T smrynum;
        EVAL_T harvest_rate;
        EVAL_T maxpossF;
        std::vector<EVAL_T> Get_EquilCalc;
        std::vector<EVAL_T> Z_rate;
        std::valarray<std::valarray<EVAL_T> > Zrate2;
        std::vector<EVAL_T> F_rate;
        std::vector<EVAL_T> Nmigr;
        EVAL_T Nsurvive;
        EVAL_T YPR_tgt_enc;
        EVAL_T YPR_tgt_dead;
        EVAL_T YPR_tgt_N_dead;
        EVAL_T YPR_tgt_ret;
        EVAL_T YPR_spr;
        EVAL_T Vbio_spr;
        EVAL_T Vbio1_spr;
        EVAL_T SPR_actual;
        EVAL_T YPR_Btgt_enc;
        EVAL_T YPR_Btgt_dead;
        EVAL_T YPR_Btgt_N_dead;
        EVAL_T YPR_Btgt_ret;
        EVAL_T YPR_Btgt;
        EVAL_T Vbio_Btgt;
        EVAL_T Vbio1_Btgt;
        EVAL_T Btgt;
        EVAL_T Btgttgt;
        EVAL_T SPR_Btgt;
        EVAL_T Btgt_Rec;
        EVAL_T Bspr;
        EVAL_T Bspr_rec;
        EVAL_T YPR;
        EVAL_T MSY;
        EVAL_T Bmsy;
        EVAL_T Recr_msy;
        EVAL_T YPR_msy_enc;
        EVAL_T YPR_msy_dead;
        EVAL_T YPR_msy_N_dead;
        EVAL_T YPR_msy_ret;
        EVAL_T YPR_enc;
        EVAL_T YPR_dead;
        EVAL_T YPR_N_dead;
        EVAL_T YPR_ret;
        EVAL_T MSY_Fmult;
        EVAL_T SPR_Fmult;
        EVAL_T Btgt_Fmult;
        EVAL_T caa;
        EVAL_T Fmult;
        EVAL_T Fcast_Fmult;
        EVAL_T Fchange;
        EVAL_T last_calc;
        std::valarray<std::valarray<EVAL_T> > Fcast_RelF_Use;
        std::valarray<std::valarray<EVAL_T> > Bmark_RelF_Use;
        EVAL_T alpha;
        EVAL_T beta;
        EVAL_T MSY_SPR;
        EVAL_T GenTime;
        std::vector<EVAL_T> cumF;
        std::vector<EVAL_T> maxF;
        EVAL_T Yield;
        EVAL_T Adj4010;
        std::vector<EVAL_T> Q_parm;
        std::valarray<std::valarray<EVAL_T> > log_q_cr;
        std::valarray<std::valarray<EVAL_T> > q_cr;
        std::valarray<std::valarray<EVAL_T> > se_cr_use;
        std::valarray<std::valarray<EVAL_T> > exp_cr;
        std::vector<EVAL_T> surv_like;
        std::valarray<std::valarray<EVAL_T> > Q_dev_like;
        std::vector<EVAL_T> disc_like;
        std::vector<EVAL_T> mnwt_like;
        std::valarray<std::valarray<EVAL_T> > exp_disc;
        std::vector<EVAL_T> retain;
        std::vector<EVAL_T> retain_M;
        std::vector<EVAL_T> discmort;
        std::vector<EVAL_T> discmort_M;
        std::vector<EVAL_T> exp_mnwt;
        std::valarray<std::valarray<EVAL_T> > Morphcomp_exp;
        std::vector<EVAL_T> SzFreqTrans;
        std::valarray<std::valarray<EVAL_T> > TG_alive;
        std::valarray<std::valarray<EVAL_T> > TG_alive_temp;
        std::vector<EVAL_T> TG_recap_exp;
        std::vector<EVAL_T> TG_like1;
        std::vector<EVAL_T> TG_like2;
        EVAL_T overdisp;
        std::vector<EVAL_T> TG_parm;
        std::vector<EVAL_T> selparm;
        std::valarray<std::valarray<EVAL_T> > selparm_trend;
        std::valarray<std::valarray<EVAL_T> > selparm_block_val;
        std::valarray<std::valarray<EVAL_T> > selparm_dev;
        std::valarray<std::valarray<EVAL_T> > selparm_dev_rwalk;
        std::vector<EVAL_T> sel_l;
        std::vector<EVAL_T> sel_l_r;
        std::vector<EVAL_T> discmort2;
        std::vector<EVAL_T> sel_a;
        std::vector<EVAL_T> sel;
        std::vector<EVAL_T> fish_body_wt;
        std::vector<EVAL_T> sel_al_1;
        std::vector<EVAL_T> sel_al_2;
        std::vector<EVAL_T> sel_al_3;
        std::vector<EVAL_T> sel_al_4;
        std::vector<EVAL_T> deadfish;
        std::vector<EVAL_T> deadfish_B;
        std::vector<EVAL_T> save_sel_fec;
        std::vector<EVAL_T> Sel_for_tag;
        std::vector<EVAL_T> TG_report;
        std::vector<EVAL_T> TG_rep_decay;
        std::vector<EVAL_T> save_sp_len;
        std::vector<EVAL_T> exp_l;
        std::valarray<std::valarray<EVAL_T> > neff_l;
        std::vector<EVAL_T> tempvec_l;
        std::vector<EVAL_T> exp_l_temp;
        std::vector<EVAL_T> exp_l_temp_ret;
        std::vector<EVAL_T> exp_l_temp_dat;
        std::vector<EVAL_T> offset_l;
        std::vector<EVAL_T> length_like;
        std::valarray<std::valarray<EVAL_T> > SzFreq_exp;
        std::vector<EVAL_T> SzFreq_like;
        std::vector<EVAL_T> exp_a;
        std::vector<EVAL_T> exp_a_temp;
        std::vector<EVAL_T> tempvec;
        std::valarray<std::valarray<EVAL_T> > neff_a;
        std::vector<EVAL_T> offset_a;
        std::vector<EVAL_T> age_like;
        std::vector<EVAL_T> sizeage_like;
        std::vector<EVAL_T> exp_ms;
        std::vector<EVAL_T> exp_ms_sq;
        EVAL_T Morphcomp_like;
        EVAL_T equ_catch_like;
        std::vector<EVAL_T> catch_like;
        EVAL_T recr_like;
        EVAL_T Fcast_recr_like;
        EVAL_T parm_like;
        EVAL_T parm_dev_like;
        EVAL_T CrashPen;
        EVAL_T SoftBoundPen;
        EVAL_T Equ_penalty;
        EVAL_T F_ballpark_like;
        EVAL_T F_ballpark_lambda;
        EVAL_T R1;
        EVAL_T R1_exp;
        EVAL_T t1;
        EVAL_T t2;
        EVAL_T temp;
        EVAL_T temp1;
        EVAL_T temp2;
        EVAL_T temp3;
        EVAL_T temp4;
        EVAL_T join1;
        EVAL_T join2;
        EVAL_T join3;
        EVAL_T upselex;
        EVAL_T downselex;
        EVAL_T peak;
        EVAL_T peak2;
        EVAL_T point1;
        EVAL_T point2;
        EVAL_T point3;
        EVAL_T point4;
        EVAL_T timing;
        EVAL_T equ_Recr;
        EVAL_T equ_F_std;
        std::valarray<std::valarray<EVAL_T> > smry;
        std::valarray<std::valarray<EVAL_T> > env_data;
        std::valarray<std::valarray<EVAL_T> > TG_save;
        std::vector<EVAL_T> SPB_std;
        std::vector<EVAL_T> recr_std;
        std::vector<EVAL_T> SPR_std;
        std::vector<EVAL_T> F_std;
        std::vector<EVAL_T> depletion;
        std::vector<EVAL_T> Mgmt_quant;
        std::vector<EVAL_T> Extra_Std;
        std::vector<EVAL_T> MGparm_Like;
        std::vector<EVAL_T> init_F_Like;
        std::vector<EVAL_T> Q_parm_Like;
        std::vector<EVAL_T> selparm_Like;
        std::vector<EVAL_T> SR_parm_Like;
        std::vector<EVAL_T> recdev_cycle_Like;
        std::vector<EVAL_T> TG_parm_Like;
        EVAL_T prior_function_value;
        EVAL_T likelihood_function_value;
        EVAL_T last_objfun;
        std::vector<EVAL_T> phase_output;


    public:

        void get_MGsetup(void) {
        }

        void get_growth1(void) {
        }

        void get_growth2(void) {
        }

        void get_natmort(void) {
        }

        void get_recr_distribution(void) {
        }

        void get_wtlen(void) {
        }

        void get_migration(void) {
        }

        void get_saveGparm(void) {
        }

        void get_selectivity(void) {
        }

        void get_age_selectivity(void) {
        }

        void get_size_selectivity(void) {
        }

        void get_initial_conditions(void) {
        }

        void get_time_series(void) {
        }

        void evaluate_the_objective_function(void) {
        }

        void Process_STDquant(void) {
        }

        EVAL_T Check_Parm(const REAL_T& Pmin, const REAL_T& Pmax, const REAL_T& jitter, const EVAL_T& Pval) {
        }

        void Report_Parm(const int NParm, const int AC, const int Activ, const int PH, const REAL_T& Pmin, const REAL_T& Pmax, const int PR_T, const REAL_T& PR, const REAL_T& CV, const REAL_T& RD, const EVAL_T& Pval, const EVAL_T& Like) {
        }

        EVAL_T Get_Prior(const int T, const REAL_T& Pmin, const REAL_T& Pmax, const REAL_T& Pr, const REAL_T& Psd, const EVAL_T& Pval) {
        }

        void Do_Equil_Calc(void) {
        }

        void Make_AgeLength_Key(void) {
        }

        void Make_FishSelex(void) {
        }

        void get_posteriors(void) {
        }

        void Get_Benchmarks(void) {
        }

        void Get_Forecast(void) {
        }

        void write_summaryoutput(void) {
        }

        void write_rebuilder_output(void) {
        }

        void write_nudata(void) {
        }

        void write_nucontrol(void) {
        }

        void write_bigoutput(void) {
        }

        void write_Bzero_output(void) {
        }

        EVAL_T Join_Fxn(const EVAL_T& MinPoss, const EVAL_T& MaxPoss, const EVAL_T& Inflec, const EVAL_T& Xvar, const EVAL_T& Y1, const EVAL_T& Y2) {
        }

        EVAL_T Spawn_Recr(const EVAL_T& SPB_current) {
        }

        std::vector<EVAL_T> Equil_Spawn_Recr_Fxn(const EVAL_T& SPB_current) {
        }

        void get_age_age(const int Keynum) {
        }

        void Evaluate(EVAL_T &f) {

        }




        //getters and setters

        std::vector<EVAL_T> GetALK() const {
            return ALK;
        }

        void SetALK(std::vector<EVAL_T> ALK) {
            this->ALK = ALK;
        }

        EVAL_T GetAdj4010() const {
            return Adj4010;
        }

        void SetAdj4010(EVAL_T Adj4010) {
            this->Adj4010 = Adj4010;
        }

        std::vector<EVAL_T> GetAve_Size() const {
            return Ave_Size;
        }

        void SetAve_Size(std::vector<EVAL_T> Ave_Size) {
            this->Ave_Size = Ave_Size;
        }

        std::valarray<std::valarray<EVAL_T> > GetBmark_RelF_Use() const {
            return Bmark_RelF_Use;
        }

        void SetBmark_RelF_Use(std::valarray<std::valarray<EVAL_T> > Bmark_RelF_Use) {
            this->Bmark_RelF_Use = Bmark_RelF_Use;
        }

        EVAL_T GetBmsy() const {
            return Bmsy;
        }

        void SetBmsy(EVAL_T Bmsy) {
            this->Bmsy = Bmsy;
        }

        EVAL_T GetBspr() const {
            return Bspr;
        }

        void SetBspr(EVAL_T Bspr) {
            this->Bspr = Bspr;
        }

        EVAL_T GetBspr_rec() const {
            return Bspr_rec;
        }

        void SetBspr_rec(EVAL_T Bspr_rec) {
            this->Bspr_rec = Bspr_rec;
        }

        EVAL_T GetBtgt() const {
            return Btgt;
        }

        void SetBtgt(EVAL_T Btgt) {
            this->Btgt = Btgt;
        }

        EVAL_T GetBtgt_Fmult() const {
            return Btgt_Fmult;
        }

        void SetBtgt_Fmult(EVAL_T Btgt_Fmult) {
            this->Btgt_Fmult = Btgt_Fmult;
        }

        EVAL_T GetBtgt_Rec() const {
            return Btgt_Rec;
        }

        void SetBtgt_Rec(EVAL_T Btgt_Rec) {
            this->Btgt_Rec = Btgt_Rec;
        }

        EVAL_T GetBtgttgt() const {
            return Btgttgt;
        }

        void SetBtgttgt(EVAL_T Btgttgt) {
            this->Btgttgt = Btgttgt;
        }

        std::vector<EVAL_T> GetCVLmax() const {
            return CVLmax;
        }

        void SetCVLmax(std::vector<EVAL_T> CVLmax) {
            this->CVLmax = CVLmax;
        }

        std::vector<EVAL_T> GetCVLmin() const {
            return CVLmin;
        }

        void SetCVLmin(std::vector<EVAL_T> CVLmin) {
            this->CVLmin = CVLmin;
        }

        std::valarray<std::valarray<EVAL_T> > GetCV_G() const {
            return CV_G;
        }

        void SetCV_G(std::valarray<std::valarray<EVAL_T> > CV_G) {
            this->CV_G = CV_G;
        }

        std::vector<EVAL_T> GetCV_const() const {
            return CV_const;
        }

        void SetCV_const(std::vector<EVAL_T> CV_const) {
            this->CV_const = CV_const;
        }

        std::valarray<std::valarray<EVAL_T> > GetCohort_Growth() const {
            return Cohort_Growth;
        }

        void SetCohort_Growth(std::valarray<std::valarray<EVAL_T> > Cohort_Growth) {
            this->Cohort_Growth = Cohort_Growth;
        }

        std::vector<EVAL_T> GetCohort_Lmin() const {
            return Cohort_Lmin;
        }

        void SetCohort_Lmin(std::vector<EVAL_T> Cohort_Lmin) {
            this->Cohort_Lmin = Cohort_Lmin;
        }

        EVAL_T GetCrashPen() const {
            return CrashPen;
        }

        void SetCrashPen(EVAL_T CrashPen) {
            this->CrashPen = CrashPen;
        }

        EVAL_T GetEqu_penalty() const {
            return Equ_penalty;
        }

        void SetEqu_penalty(EVAL_T Equ_penalty) {
            this->Equ_penalty = Equ_penalty;
        }

        std::vector<EVAL_T> GetExtra_Std() const {
            return Extra_Std;
        }

        void SetExtra_Std(std::vector<EVAL_T> Extra_Std) {
            this->Extra_Std = Extra_Std;
        }

        EVAL_T GetF_ballpark_lambda() const {
            return F_ballpark_lambda;
        }

        void SetF_ballpark_lambda(EVAL_T F_ballpark_lambda) {
            this->F_ballpark_lambda = F_ballpark_lambda;
        }

        EVAL_T GetF_ballpark_like() const {
            return F_ballpark_like;
        }

        void SetF_ballpark_like(EVAL_T F_ballpark_like) {
            this->F_ballpark_like = F_ballpark_like;
        }

        std::vector<EVAL_T> GetF_rate() const {
            return F_rate;
        }

        void SetF_rate(std::vector<EVAL_T> F_rate) {
            this->F_rate = F_rate;
        }

        std::vector<EVAL_T> GetF_std() const {
            return F_std;
        }

        void SetF_std(std::vector<EVAL_T> F_std) {
            this->F_std = F_std;
        }

        EVAL_T GetFcast_Fmult() const {
            return Fcast_Fmult;
        }

        void SetFcast_Fmult(EVAL_T Fcast_Fmult) {
            this->Fcast_Fmult = Fcast_Fmult;
        }

        std::valarray<std::valarray<EVAL_T> > GetFcast_RelF_Use() const {
            return Fcast_RelF_Use;
        }

        void SetFcast_RelF_Use(std::valarray<std::valarray<EVAL_T> > Fcast_RelF_Use) {
            this->Fcast_RelF_Use = Fcast_RelF_Use;
        }

        std::vector<EVAL_T> GetFcast_impl_error() const {
            return Fcast_impl_error;
        }

        void SetFcast_impl_error(std::vector<EVAL_T> Fcast_impl_error) {
            this->Fcast_impl_error = Fcast_impl_error;
        }

        EVAL_T GetFcast_recr_like() const {
            return Fcast_recr_like;
        }

        void SetFcast_recr_like(EVAL_T Fcast_recr_like) {
            this->Fcast_recr_like = Fcast_recr_like;
        }

        std::vector<EVAL_T> GetFcast_recruitments() const {
            return Fcast_recruitments;
        }

        void SetFcast_recruitments(std::vector<EVAL_T> Fcast_recruitments) {
            this->Fcast_recruitments = Fcast_recruitments;
        }

        EVAL_T GetFchange() const {
            return Fchange;
        }

        void SetFchange(EVAL_T Fchange) {
            this->Fchange = Fchange;
        }

        EVAL_T GetFmult() const {
            return Fmult;
        }

        void SetFmult(EVAL_T Fmult) {
            this->Fmult = Fmult;
        }

        EVAL_T GetGenTime() const {
            return GenTime;
        }

        void SetGenTime(EVAL_T GenTime) {
            this->GenTime = GenTime;
        }

        std::vector<EVAL_T> GetGet_EquilCalc() const {
            return Get_EquilCalc;
        }

        void SetGet_EquilCalc(std::vector<EVAL_T> Get_EquilCalc) {
            this->Get_EquilCalc = Get_EquilCalc;
        }

        std::valarray<std::valarray<EVAL_T> > GetHermaphro_val() const {
            return Hermaphro_val;
        }

        void SetHermaphro_val(std::valarray<std::valarray<EVAL_T> > Hermaphro_val) {
            this->Hermaphro_val = Hermaphro_val;
        }

        std::valarray<std::valarray<EVAL_T> > GetHrate() const {
            return Hrate;
        }

        void SetHrate(std::valarray<std::valarray<EVAL_T> > Hrate) {
            this->Hrate = Hrate;
        }

        std::vector<EVAL_T> GetL_inf() const {
            return L_inf;
        }

        void SetL_inf(std::vector<EVAL_T> L_inf) {
            this->L_inf = L_inf;
        }

        std::vector<EVAL_T> GetLmax_temp() const {
            return Lmax_temp;
        }

        void SetLmax_temp(std::vector<EVAL_T> Lmax_temp) {
            this->Lmax_temp = Lmax_temp;
        }

        std::vector<EVAL_T> GetLmin() const {
            return Lmin;
        }

        void SetLmin(std::vector<EVAL_T> Lmin) {
            this->Lmin = Lmin;
        }

        std::vector<EVAL_T> GetLmin_last() const {
            return Lmin_last;
        }

        void SetLmin_last(std::vector<EVAL_T> Lmin_last) {
            this->Lmin_last = Lmin_last;
        }

        std::vector<EVAL_T> GetMGparm() const {
            return MGparm;
        }

        void SetMGparm(std::vector<EVAL_T> MGparm) {
            this->MGparm = MGparm;
        }

        std::vector<EVAL_T> GetMGparm_Like() const {
            return MGparm_Like;
        }

        void SetMGparm_Like(std::vector<EVAL_T> MGparm_Like) {
            this->MGparm_Like = MGparm_Like;
        }

        std::valarray<std::valarray<EVAL_T> > GetMGparm_block_val() const {
            return MGparm_block_val;
        }

        void SetMGparm_block_val(std::valarray<std::valarray<EVAL_T> > MGparm_block_val) {
            this->MGparm_block_val = MGparm_block_val;
        }

        std::valarray<std::valarray<EVAL_T> > GetMGparm_dev() const {
            return MGparm_dev;
        }

        void SetMGparm_dev(std::valarray<std::valarray<EVAL_T> > MGparm_dev) {
            this->MGparm_dev = MGparm_dev;
        }

        std::valarray<std::valarray<EVAL_T> > GetMGparm_dev_rwalk() const {
            return MGparm_dev_rwalk;
        }

        void SetMGparm_dev_rwalk(std::valarray<std::valarray<EVAL_T> > MGparm_dev_rwalk) {
            this->MGparm_dev_rwalk = MGparm_dev_rwalk;
        }

        std::valarray<std::valarray<EVAL_T> > GetMGparm_trend() const {
            return MGparm_trend;
        }

        void SetMGparm_trend(std::valarray<std::valarray<EVAL_T> > MGparm_trend) {
            this->MGparm_trend = MGparm_trend;
        }

        EVAL_T GetMSY() const {
            return MSY;
        }

        void SetMSY(EVAL_T MSY) {
            this->MSY = MSY;
        }

        EVAL_T GetMSY_Fmult() const {
            return MSY_Fmult;
        }

        void SetMSY_Fmult(EVAL_T MSY_Fmult) {
            this->MSY_Fmult = MSY_Fmult;
        }

        EVAL_T GetMSY_SPR() const {
            return MSY_SPR;
        }

        void SetMSY_SPR(EVAL_T MSY_SPR) {
            this->MSY_SPR = MSY_SPR;
        }

        std::vector<EVAL_T> GetMaleSPB() const {
            return MaleSPB;
        }

        void SetMaleSPB(std::vector<EVAL_T> MaleSPB) {
            this->MaleSPB = MaleSPB;
        }

        std::valarray<std::valarray<EVAL_T> > GetMaleSPB_equil_pop_gp() const {
            return MaleSPB_equil_pop_gp;
        }

        void SetMaleSPB_equil_pop_gp(std::valarray<std::valarray<EVAL_T> > MaleSPB_equil_pop_gp) {
            this->MaleSPB_equil_pop_gp = MaleSPB_equil_pop_gp;
        }

        std::vector<EVAL_T> GetMgmt_quant() const {
            return Mgmt_quant;
        }

        void SetMgmt_quant(std::vector<EVAL_T> Mgmt_quant) {
            this->Mgmt_quant = Mgmt_quant;
        }

        std::valarray<std::valarray<EVAL_T> > GetMorphcomp_exp() const {
            return Morphcomp_exp;
        }

        void SetMorphcomp_exp(std::valarray<std::valarray<EVAL_T> > Morphcomp_exp) {
            this->Morphcomp_exp = Morphcomp_exp;
        }

        EVAL_T GetMorphcomp_like() const {
            return Morphcomp_like;
        }

        void SetMorphcomp_like(EVAL_T Morphcomp_like) {
            this->Morphcomp_like = Morphcomp_like;
        }

        std::valarray<std::valarray<EVAL_T> > GetNmid() const {
            return Nmid;
        }

        void SetNmid(std::valarray<std::valarray<EVAL_T> > Nmid) {
            this->Nmid = Nmid;
        }

        std::vector<EVAL_T> GetNmigr() const {
            return Nmigr;
        }

        void SetNmigr(std::vector<EVAL_T> Nmigr) {
            this->Nmigr = Nmigr;
        }

        std::valarray<std::valarray<EVAL_T> > GetNsurv() const {
            return Nsurv;
        }

        void SetNsurv(std::valarray<std::valarray<EVAL_T> > Nsurv) {
            this->Nsurv = Nsurv;
        }

        EVAL_T GetNsurvive() const {
            return Nsurvive;
        }

        void SetNsurvive(EVAL_T Nsurvive) {
            this->Nsurvive = Nsurvive;
        }

        std::valarray<std::valarray<EVAL_T> > GetQ_dev_like() const {
            return Q_dev_like;
        }

        void SetQ_dev_like(std::valarray<std::valarray<EVAL_T> > Q_dev_like) {
            this->Q_dev_like = Q_dev_like;
        }

        std::vector<EVAL_T> GetQ_parm() const {
            return Q_parm;
        }

        void SetQ_parm(std::vector<EVAL_T> Q_parm) {
            this->Q_parm = Q_parm;
        }

        std::vector<EVAL_T> GetQ_parm_Like() const {
            return Q_parm_Like;
        }

        void SetQ_parm_Like(std::vector<EVAL_T> Q_parm_Like) {
            this->Q_parm_Like = Q_parm_Like;
        }

        EVAL_T GetR1() const {
            return R1;
        }

        void SetR1(EVAL_T R1) {
            this->R1 = R1;
        }

        EVAL_T GetR1_exp() const {
            return R1_exp;
        }

        void SetR1_exp(EVAL_T R1_exp) {
            this->R1_exp = R1_exp;
        }

        std::valarray<std::valarray<EVAL_T> > GetRecr() const {
            return Recr;
        }

        void SetRecr(std::valarray<std::valarray<EVAL_T> > Recr) {
            this->Recr = Recr;
        }

        EVAL_T GetRecr_msy() const {
            return Recr_msy;
        }

        void SetRecr_msy(EVAL_T Recr_msy) {
            this->Recr_msy = Recr_msy;
        }

        EVAL_T GetRecr_virgin() const {
            return Recr_virgin;
        }

        void SetRecr_virgin(EVAL_T Recr_virgin) {
            this->Recr_virgin = Recr_virgin;
        }

        EVAL_T GetRecruits() const {
            return Recruits;
        }

        void SetRecruits(EVAL_T Recruits) {
            this->Recruits = Recruits;
        }

        std::vector<EVAL_T> GetRichards() const {
            return Richards;
        }

        void SetRichards(std::vector<EVAL_T> Richards) {
            this->Richards = Richards;
        }

        EVAL_T GetSPB_current() const {
            return SPB_current;
        }

        void SetSPB_current(EVAL_T SPB_current) {
            this->SPB_current = SPB_current;
        }

        EVAL_T GetSPB_equil() const {
            return SPB_equil;
        }

        void SetSPB_equil(EVAL_T SPB_equil) {
            this->SPB_equil = SPB_equil;
        }

        std::valarray<std::valarray<EVAL_T> > GetSPB_equil_pop_gp() const {
            return SPB_equil_pop_gp;
        }

        void SetSPB_equil_pop_gp(std::valarray<std::valarray<EVAL_T> > SPB_equil_pop_gp) {
            this->SPB_equil_pop_gp = SPB_equil_pop_gp;
        }

        std::vector<EVAL_T> GetSPB_pop_gp() const {
            return SPB_pop_gp;
        }

        void SetSPB_pop_gp(std::vector<EVAL_T> SPB_pop_gp) {
            this->SPB_pop_gp = SPB_pop_gp;
        }

        std::vector<EVAL_T> GetSPB_std() const {
            return SPB_std;
        }

        void SetSPB_std(std::vector<EVAL_T> SPB_std) {
            this->SPB_std = SPB_std;
        }

        EVAL_T GetSPB_vir_LH() const {
            return SPB_vir_LH;
        }

        void SetSPB_vir_LH(EVAL_T SPB_vir_LH) {
            this->SPB_vir_LH = SPB_vir_LH;
        }

        EVAL_T GetSPB_virgin() const {
            return SPB_virgin;
        }

        void SetSPB_virgin(EVAL_T SPB_virgin) {
            this->SPB_virgin = SPB_virgin;
        }

        std::vector<EVAL_T> GetSPB_yr() const {
            return SPB_yr;
        }

        void SetSPB_yr(std::vector<EVAL_T> SPB_yr) {
            this->SPB_yr = SPB_yr;
        }

        EVAL_T GetSPR_Btgt() const {
            return SPR_Btgt;
        }

        void SetSPR_Btgt(EVAL_T SPR_Btgt) {
            this->SPR_Btgt = SPR_Btgt;
        }

        EVAL_T GetSPR_Fmult() const {
            return SPR_Fmult;
        }

        void SetSPR_Fmult(EVAL_T SPR_Fmult) {
            this->SPR_Fmult = SPR_Fmult;
        }

        EVAL_T GetSPR_actual() const {
            return SPR_actual;
        }

        void SetSPR_actual(EVAL_T SPR_actual) {
            this->SPR_actual = SPR_actual;
        }

        std::vector<EVAL_T> GetSPR_std() const {
            return SPR_std;
        }

        void SetSPR_std(std::vector<EVAL_T> SPR_std) {
            this->SPR_std = SPR_std;
        }

        EVAL_T GetSPR_trial() const {
            return SPR_trial;
        }

        void SetSPR_trial(EVAL_T SPR_trial) {
            this->SPR_trial = SPR_trial;
        }

        EVAL_T GetSPR_unf() const {
            return SPR_unf;
        }

        void SetSPR_unf(EVAL_T SPR_unf) {
            this->SPR_unf = SPR_unf;
        }

        std::vector<EVAL_T> GetSR_parm() const {
            return SR_parm;
        }

        void SetSR_parm(std::vector<EVAL_T> SR_parm) {
            this->SR_parm = SR_parm;
        }

        std::vector<EVAL_T> GetSR_parm_Like() const {
            return SR_parm_Like;
        }

        void SetSR_parm_Like(std::vector<EVAL_T> SR_parm_Like) {
            this->SR_parm_Like = SR_parm_Like;
        }

        std::vector<EVAL_T> GetSave_PopAge() const {
            return Save_PopAge;
        }

        void SetSave_PopAge(std::vector<EVAL_T> Save_PopAge) {
            this->Save_PopAge = Save_PopAge;
        }

        std::vector<EVAL_T> GetSave_PopLen() const {
            return Save_PopLen;
        }

        void SetSave_PopLen(std::vector<EVAL_T> Save_PopLen) {
            this->Save_PopLen = Save_PopLen;
        }

        std::vector<EVAL_T> GetSave_Wt_Age() const {
            return Save_Wt_Age;
        }

        void SetSave_Wt_Age(std::vector<EVAL_T> Save_Wt_Age) {
            this->Save_Wt_Age = Save_Wt_Age;
        }

        std::vector<EVAL_T> GetSd_Size_between() const {
            return Sd_Size_between;
        }

        void SetSd_Size_between(std::vector<EVAL_T> Sd_Size_between) {
            this->Sd_Size_between = Sd_Size_between;
        }

        std::vector<EVAL_T> GetSd_Size_within() const {
            return Sd_Size_within;
        }

        void SetSd_Size_within(std::vector<EVAL_T> Sd_Size_within) {
            this->Sd_Size_within = Sd_Size_within;
        }

        std::vector<EVAL_T> GetSel_for_tag() const {
            return Sel_for_tag;
        }

        void SetSel_for_tag(std::vector<EVAL_T> Sel_for_tag) {
            this->Sel_for_tag = Sel_for_tag;
        }

        EVAL_T GetSoftBoundPen() const {
            return SoftBoundPen;
        }

        void SetSoftBoundPen(EVAL_T SoftBoundPen) {
            this->SoftBoundPen = SoftBoundPen;
        }

        std::vector<EVAL_T> GetSzFreqTrans() const {
            return SzFreqTrans;
        }

        void SetSzFreqTrans(std::vector<EVAL_T> SzFreqTrans) {
            this->SzFreqTrans = SzFreqTrans;
        }

        std::valarray<std::valarray<EVAL_T> > GetSzFreq_exp() const {
            return SzFreq_exp;
        }

        void SetSzFreq_exp(std::valarray<std::valarray<EVAL_T> > SzFreq_exp) {
            this->SzFreq_exp = SzFreq_exp;
        }

        std::vector<EVAL_T> GetSzFreq_like() const {
            return SzFreq_like;
        }

        void SetSzFreq_like(std::vector<EVAL_T> SzFreq_like) {
            this->SzFreq_like = SzFreq_like;
        }

        std::valarray<std::valarray<EVAL_T> > GetTG_alive() const {
            return TG_alive;
        }

        void SetTG_alive(std::valarray<std::valarray<EVAL_T> > TG_alive) {
            this->TG_alive = TG_alive;
        }

        std::valarray<std::valarray<EVAL_T> > GetTG_alive_temp() const {
            return TG_alive_temp;
        }

        void SetTG_alive_temp(std::valarray<std::valarray<EVAL_T> > TG_alive_temp) {
            this->TG_alive_temp = TG_alive_temp;
        }

        std::vector<EVAL_T> GetTG_like1() const {
            return TG_like1;
        }

        void SetTG_like1(std::vector<EVAL_T> TG_like1) {
            this->TG_like1 = TG_like1;
        }

        std::vector<EVAL_T> GetTG_like2() const {
            return TG_like2;
        }

        void SetTG_like2(std::vector<EVAL_T> TG_like2) {
            this->TG_like2 = TG_like2;
        }

        std::vector<EVAL_T> GetTG_parm() const {
            return TG_parm;
        }

        void SetTG_parm(std::vector<EVAL_T> TG_parm) {
            this->TG_parm = TG_parm;
        }

        std::vector<EVAL_T> GetTG_parm_Like() const {
            return TG_parm_Like;
        }

        void SetTG_parm_Like(std::vector<EVAL_T> TG_parm_Like) {
            this->TG_parm_Like = TG_parm_Like;
        }

        std::vector<EVAL_T> GetTG_recap_exp() const {
            return TG_recap_exp;
        }

        void SetTG_recap_exp(std::vector<EVAL_T> TG_recap_exp) {
            this->TG_recap_exp = TG_recap_exp;
        }

        std::vector<EVAL_T> GetTG_rep_decay() const {
            return TG_rep_decay;
        }

        void SetTG_rep_decay(std::vector<EVAL_T> TG_rep_decay) {
            this->TG_rep_decay = TG_rep_decay;
        }

        std::vector<EVAL_T> GetTG_report() const {
            return TG_report;
        }

        void SetTG_report(std::vector<EVAL_T> TG_report) {
            this->TG_report = TG_report;
        }

        std::valarray<std::valarray<EVAL_T> > GetTG_save() const {
            return TG_save;
        }

        void SetTG_save(std::valarray<std::valarray<EVAL_T> > TG_save) {
            this->TG_save = TG_save;
        }

        std::vector<EVAL_T> GetVBK() const {
            return VBK;
        }

        void SetVBK(std::vector<EVAL_T> VBK) {
            this->VBK = VBK;
        }

        std::vector<EVAL_T> GetVBK_seas() const {
            return VBK_seas;
        }

        void SetVBK_seas(std::vector<EVAL_T> VBK_seas) {
            this->VBK_seas = VBK_seas;
        }

        EVAL_T GetVbio1_Btgt() const {
            return Vbio1_Btgt;
        }

        void SetVbio1_Btgt(EVAL_T Vbio1_Btgt) {
            this->Vbio1_Btgt = Vbio1_Btgt;
        }

        EVAL_T GetVbio1_spr() const {
            return Vbio1_spr;
        }

        void SetVbio1_spr(EVAL_T Vbio1_spr) {
            this->Vbio1_spr = Vbio1_spr;
        }

        EVAL_T GetVbio_Btgt() const {
            return Vbio_Btgt;
        }

        void SetVbio_Btgt(EVAL_T Vbio_Btgt) {
            this->Vbio_Btgt = Vbio_Btgt;
        }

        EVAL_T GetVbio_spr() const {
            return Vbio_spr;
        }

        void SetVbio_spr(EVAL_T Vbio_spr) {
            this->Vbio_spr = Vbio_spr;
        }

        std::vector<EVAL_T> GetWt_Age_beg() const {
            return Wt_Age_beg;
        }

        void SetWt_Age_beg(std::vector<EVAL_T> Wt_Age_beg) {
            this->Wt_Age_beg = Wt_Age_beg;
        }

        std::vector<EVAL_T> GetWt_Age_mid() const {
            return Wt_Age_mid;
        }

        void SetWt_Age_mid(std::vector<EVAL_T> Wt_Age_mid) {
            this->Wt_Age_mid = Wt_Age_mid;
        }

        EVAL_T GetYPR() const {
            return YPR;
        }

        void SetYPR(EVAL_T YPR) {
            this->YPR = YPR;
        }

        EVAL_T GetYPR_Btgt() const {
            return YPR_Btgt;
        }

        void SetYPR_Btgt(EVAL_T YPR_Btgt) {
            this->YPR_Btgt = YPR_Btgt;
        }

        EVAL_T GetYPR_Btgt_N_dead() const {
            return YPR_Btgt_N_dead;
        }

        void SetYPR_Btgt_N_dead(EVAL_T YPR_Btgt_N_dead) {
            this->YPR_Btgt_N_dead = YPR_Btgt_N_dead;
        }

        EVAL_T GetYPR_Btgt_dead() const {
            return YPR_Btgt_dead;
        }

        void SetYPR_Btgt_dead(EVAL_T YPR_Btgt_dead) {
            this->YPR_Btgt_dead = YPR_Btgt_dead;
        }

        EVAL_T GetYPR_Btgt_enc() const {
            return YPR_Btgt_enc;
        }

        void SetYPR_Btgt_enc(EVAL_T YPR_Btgt_enc) {
            this->YPR_Btgt_enc = YPR_Btgt_enc;
        }

        EVAL_T GetYPR_Btgt_ret() const {
            return YPR_Btgt_ret;
        }

        void SetYPR_Btgt_ret(EVAL_T YPR_Btgt_ret) {
            this->YPR_Btgt_ret = YPR_Btgt_ret;
        }

        EVAL_T GetYPR_N_dead() const {
            return YPR_N_dead;
        }

        void SetYPR_N_dead(EVAL_T YPR_N_dead) {
            this->YPR_N_dead = YPR_N_dead;
        }

        EVAL_T GetYPR_dead() const {
            return YPR_dead;
        }

        void SetYPR_dead(EVAL_T YPR_dead) {
            this->YPR_dead = YPR_dead;
        }

        EVAL_T GetYPR_enc() const {
            return YPR_enc;
        }

        void SetYPR_enc(EVAL_T YPR_enc) {
            this->YPR_enc = YPR_enc;
        }

        EVAL_T GetYPR_msy_N_dead() const {
            return YPR_msy_N_dead;
        }

        void SetYPR_msy_N_dead(EVAL_T YPR_msy_N_dead) {
            this->YPR_msy_N_dead = YPR_msy_N_dead;
        }

        EVAL_T GetYPR_msy_dead() const {
            return YPR_msy_dead;
        }

        void SetYPR_msy_dead(EVAL_T YPR_msy_dead) {
            this->YPR_msy_dead = YPR_msy_dead;
        }

        EVAL_T GetYPR_msy_enc() const {
            return YPR_msy_enc;
        }

        void SetYPR_msy_enc(EVAL_T YPR_msy_enc) {
            this->YPR_msy_enc = YPR_msy_enc;
        }

        EVAL_T GetYPR_msy_ret() const {
            return YPR_msy_ret;
        }

        void SetYPR_msy_ret(EVAL_T YPR_msy_ret) {
            this->YPR_msy_ret = YPR_msy_ret;
        }

        EVAL_T GetYPR_ret() const {
            return YPR_ret;
        }

        void SetYPR_ret(EVAL_T YPR_ret) {
            this->YPR_ret = YPR_ret;
        }

        EVAL_T GetYPR_spr() const {
            return YPR_spr;
        }

        void SetYPR_spr(EVAL_T YPR_spr) {
            this->YPR_spr = YPR_spr;
        }

        EVAL_T GetYPR_tgt_N_dead() const {
            return YPR_tgt_N_dead;
        }

        void SetYPR_tgt_N_dead(EVAL_T YPR_tgt_N_dead) {
            this->YPR_tgt_N_dead = YPR_tgt_N_dead;
        }

        EVAL_T GetYPR_tgt_dead() const {
            return YPR_tgt_dead;
        }

        void SetYPR_tgt_dead(EVAL_T YPR_tgt_dead) {
            this->YPR_tgt_dead = YPR_tgt_dead;
        }

        EVAL_T GetYPR_tgt_enc() const {
            return YPR_tgt_enc;
        }

        void SetYPR_tgt_enc(EVAL_T YPR_tgt_enc) {
            this->YPR_tgt_enc = YPR_tgt_enc;
        }

        EVAL_T GetYPR_tgt_ret() const {
            return YPR_tgt_ret;
        }

        void SetYPR_tgt_ret(EVAL_T YPR_tgt_ret) {
            this->YPR_tgt_ret = YPR_tgt_ret;
        }

        EVAL_T GetYield() const {
            return Yield;
        }

        void SetYield(EVAL_T Yield) {
            this->Yield = Yield;
        }

        std::vector<EVAL_T> GetZ_rate() const {
            return Z_rate;
        }

        void SetZ_rate(std::vector<EVAL_T> Z_rate) {
            this->Z_rate = Z_rate;
        }

        std::valarray<std::valarray<EVAL_T> > GetZrate2() const {
            return Zrate2;
        }

        void SetZrate2(std::valarray<std::valarray<EVAL_T> > Zrate2) {
            this->Zrate2 = Zrate2;
        }

        std::vector<EVAL_T> GetAge_age() const {
            return age_age;
        }

        void SetAge_age(std::vector<EVAL_T> age_age) {
            this->age_age = age_age;
        }

        std::vector<EVAL_T> GetAge_err() const {
            return age_err;
        }

        void SetAge_err(std::vector<EVAL_T> age_err) {
            this->age_err = age_err;
        }

        std::vector<EVAL_T> GetAge_like() const {
            return age_like;
        }

        void SetAge_like(std::vector<EVAL_T> age_like) {
            this->age_like = age_like;
        }

        std::vector<EVAL_T> GetAgetemp() const {
            return agetemp;
        }

        void SetAgetemp(std::vector<EVAL_T> agetemp) {
            this->agetemp = agetemp;
        }

        EVAL_T GetAlpha() const {
            return alpha;
        }

        void SetAlpha(EVAL_T alpha) {
            this->alpha = alpha;
        }

        EVAL_T GetAve_age() const {
            return ave_age;
        }

        void SetAve_age(EVAL_T ave_age) {
            this->ave_age = ave_age;
        }

        EVAL_T GetBeta() const {
            return beta;
        }

        void SetBeta(EVAL_T beta) {
            this->beta = beta;
        }

        std::vector<EVAL_T> GetBiasadj() const {
            return biasadj;
        }

        void SetBiasadj(std::vector<EVAL_T> biasadj) {
            this->biasadj = biasadj;
        }

        std::vector<EVAL_T> GetBiasadj_full() const {
            return biasadj_full;
        }

        void SetBiasadj_full(std::vector<EVAL_T> biasadj_full) {
            this->biasadj_full = biasadj_full;
        }

        EVAL_T GetCaa() const {
            return caa;
        }

        void SetCaa(EVAL_T caa) {
            this->caa = caa;
        }

        std::vector<EVAL_T> GetCatage() const {
            return catage;
        }

        void SetCatage(std::vector<EVAL_T> catage) {
            this->catage = catage;
        }

        std::valarray<std::valarray<EVAL_T> > GetCatage_tot() const {
            return catage_tot;
        }

        void SetCatage_tot(std::valarray<std::valarray<EVAL_T> > catage_tot) {
            this->catage_tot = catage_tot;
        }

        std::vector<EVAL_T> GetCatch_fleet() const {
            return catch_fleet;
        }

        void SetCatch_fleet(std::vector<EVAL_T> catch_fleet) {
            this->catch_fleet = catch_fleet;
        }

        std::vector<EVAL_T> GetCatch_like() const {
            return catch_like;
        }

        void SetCatch_like(std::vector<EVAL_T> catch_like) {
            this->catch_like = catch_like;
        }

        std::vector<EVAL_T> GetCumF() const {
            return cumF;
        }

        void SetCumF(std::vector<EVAL_T> cumF) {
            this->cumF = cumF;
        }

        std::vector<EVAL_T> GetDeadfish() const {
            return deadfish;
        }

        void SetDeadfish(std::vector<EVAL_T> deadfish) {
            this->deadfish = deadfish;
        }

        std::vector<EVAL_T> GetDeadfish_B() const {
            return deadfish_B;
        }

        void SetDeadfish_B(std::vector<EVAL_T> deadfish_B) {
            this->deadfish_B = deadfish_B;
        }

        std::vector<EVAL_T> GetDepletion() const {
            return depletion;
        }

        void SetDepletion(std::vector<EVAL_T> depletion) {
            this->depletion = depletion;
        }

        std::vector<EVAL_T> GetDisc_like() const {
            return disc_like;
        }

        void SetDisc_like(std::vector<EVAL_T> disc_like) {
            this->disc_like = disc_like;
        }

        std::vector<EVAL_T> GetDiscmort() const {
            return discmort;
        }

        void SetDiscmort(std::vector<EVAL_T> discmort) {
            this->discmort = discmort;
        }

        std::vector<EVAL_T> GetDiscmort2() const {
            return discmort2;
        }

        void SetDiscmort2(std::vector<EVAL_T> discmort2) {
            this->discmort2 = discmort2;
        }

        std::vector<EVAL_T> GetDiscmort_M() const {
            return discmort_M;
        }

        void SetDiscmort_M(std::vector<EVAL_T> discmort_M) {
            this->discmort_M = discmort_M;
        }

        std::vector<REAL_T> GetDouble_control_flags() const {
            return double_control_flags;
        }

        void SetDouble_control_flags(std::vector<REAL_T> double_control_flags) {
            this->double_control_flags = double_control_flags;
        }

        EVAL_T GetDownselex() const {
            return downselex;
        }

        void SetDownselex(EVAL_T downselex) {
            this->downselex = downselex;
        }

        EVAL_T GetDummy_parm() const {
            return dummy_parm;
        }

        void SetDummy_parm(EVAL_T dummy_parm) {
            this->dummy_parm = dummy_parm;
        }

        std::valarray<std::valarray<EVAL_T> > GetEnv_data() const {
            return env_data;
        }

        void SetEnv_data(std::valarray<std::valarray<EVAL_T> > env_data) {
            this->env_data = env_data;
        }

        EVAL_T GetEqu_F_std() const {
            return equ_F_std;
        }

        void SetEqu_F_std(EVAL_T equ_F_std) {
            this->equ_F_std = equ_F_std;
        }

        EVAL_T GetEqu_Recr() const {
            return equ_Recr;
        }

        void SetEqu_Recr(EVAL_T equ_Recr) {
            this->equ_Recr = equ_Recr;
        }

        std::vector<EVAL_T> GetEqu_Z() const {
            return equ_Z;
        }

        void SetEqu_Z(std::vector<EVAL_T> equ_Z) {
            this->equ_Z = equ_Z;
        }

        std::vector<EVAL_T> GetEqu_catage() const {
            return equ_catage;
        }

        void SetEqu_catage(std::vector<EVAL_T> equ_catage) {
            this->equ_catage = equ_catage;
        }

        std::valarray<std::valarray<EVAL_T> > GetEqu_catch_fleet() const {
            return equ_catch_fleet;
        }

        void SetEqu_catch_fleet(std::valarray<std::valarray<EVAL_T> > equ_catch_fleet) {
            this->equ_catch_fleet = equ_catch_fleet;
        }

        EVAL_T GetEqu_catch_like() const {
            return equ_catch_like;
        }

        void SetEqu_catch_like(EVAL_T equ_catch_like) {
            this->equ_catch_like = equ_catch_like;
        }

        std::vector<EVAL_T> GetEqu_numbers() const {
            return equ_numbers;
        }

        void SetEqu_numbers(std::vector<EVAL_T> equ_numbers) {
            this->equ_numbers = equ_numbers;
        }

        std::vector<EVAL_T> GetEst_equ_catch() const {
            return est_equ_catch;
        }

        void SetEst_equ_catch(std::vector<EVAL_T> est_equ_catch) {
            this->est_equ_catch = est_equ_catch;
        }

        std::valarray<std::valarray<EVAL_T> > GetExp_AL() const {
            return exp_AL;
        }

        void SetExp_AL(std::valarray<std::valarray<EVAL_T> > exp_AL) {
            this->exp_AL = exp_AL;
        }

        std::vector<EVAL_T> GetExp_a() const {
            return exp_a;
        }

        void SetExp_a(std::vector<EVAL_T> exp_a) {
            this->exp_a = exp_a;
        }

        std::vector<EVAL_T> GetExp_a_temp() const {
            return exp_a_temp;
        }

        void SetExp_a_temp(std::vector<EVAL_T> exp_a_temp) {
            this->exp_a_temp = exp_a_temp;
        }

        std::valarray<std::valarray<EVAL_T> > GetExp_cr() const {
            return exp_cr;
        }

        void SetExp_cr(std::valarray<std::valarray<EVAL_T> > exp_cr) {
            this->exp_cr = exp_cr;
        }

        std::valarray<std::valarray<EVAL_T> > GetExp_disc() const {
            return exp_disc;
        }

        void SetExp_disc(std::valarray<std::valarray<EVAL_T> > exp_disc) {
            this->exp_disc = exp_disc;
        }

        std::vector<EVAL_T> GetExp_l() const {
            return exp_l;
        }

        void SetExp_l(std::vector<EVAL_T> exp_l) {
            this->exp_l = exp_l;
        }

        std::vector<EVAL_T> GetExp_l_temp() const {
            return exp_l_temp;
        }

        void SetExp_l_temp(std::vector<EVAL_T> exp_l_temp) {
            this->exp_l_temp = exp_l_temp;
        }

        std::vector<EVAL_T> GetExp_l_temp_dat() const {
            return exp_l_temp_dat;
        }

        void SetExp_l_temp_dat(std::vector<EVAL_T> exp_l_temp_dat) {
            this->exp_l_temp_dat = exp_l_temp_dat;
        }

        std::vector<EVAL_T> GetExp_l_temp_ret() const {
            return exp_l_temp_ret;
        }

        void SetExp_l_temp_ret(std::vector<EVAL_T> exp_l_temp_ret) {
            this->exp_l_temp_ret = exp_l_temp_ret;
        }

        std::vector<EVAL_T> GetExp_mnwt() const {
            return exp_mnwt;
        }

        void SetExp_mnwt(std::vector<EVAL_T> exp_mnwt) {
            this->exp_mnwt = exp_mnwt;
        }

        std::vector<EVAL_T> GetExp_ms() const {
            return exp_ms;
        }

        void SetExp_ms(std::vector<EVAL_T> exp_ms) {
            this->exp_ms = exp_ms;
        }

        std::vector<EVAL_T> GetExp_ms_sq() const {
            return exp_ms_sq;
        }

        void SetExp_ms_sq(std::vector<EVAL_T> exp_ms_sq) {
            this->exp_ms_sq = exp_ms_sq;
        }

        std::vector<EVAL_T> GetExp_rec() const {
            return exp_rec;
        }

        void SetExp_rec(std::vector<EVAL_T> exp_rec) {
            this->exp_rec = exp_rec;
        }

        std::valarray<std::valarray<EVAL_T> > GetFec() const {
            return fec;
        }

        void SetFec(std::valarray<std::valarray<EVAL_T> > fec) {
            this->fec = fec;
        }

        std::vector<EVAL_T> GetFec_len() const {
            return fec_len;
        }

        void SetFec_len(std::vector<EVAL_T> fec_len) {
            this->fec_len = fec_len;
        }

        EVAL_T GetFish_bio() const {
            return fish_bio;
        }

        void SetFish_bio(EVAL_T fish_bio) {
            this->fish_bio = fish_bio;
        }

        EVAL_T GetFish_bio_e() const {
            return fish_bio_e;
        }

        void SetFish_bio_e(EVAL_T fish_bio_e) {
            this->fish_bio_e = fish_bio_e;
        }

        EVAL_T GetFish_bio_r() const {
            return fish_bio_r;
        }

        void SetFish_bio_r(EVAL_T fish_bio_r) {
            this->fish_bio_r = fish_bio_r;
        }

        std::vector<EVAL_T> GetFish_body_wt() const {
            return fish_body_wt;
        }

        void SetFish_body_wt(std::vector<EVAL_T> fish_body_wt) {
            this->fish_body_wt = fish_body_wt;
        }

        EVAL_T GetFish_num() const {
            return fish_num;
        }

        void SetFish_num(EVAL_T fish_num) {
            this->fish_num = fish_num;
        }

        EVAL_T GetFish_num_e() const {
            return fish_num_e;
        }

        void SetFish_num_e(EVAL_T fish_num_e) {
            this->fish_num_e = fish_num_e;
        }

        EVAL_T GetFish_num_r() const {
            return fish_num_r;
        }

        void SetFish_num_r(EVAL_T fish_num_r) {
            this->fish_num_r = fish_num_r;
        }

        EVAL_T GetHalf_sigmaRsq() const {
            return half_sigmaRsq;
        }

        void SetHalf_sigmaRsq(EVAL_T half_sigmaRsq) {
            this->half_sigmaRsq = half_sigmaRsq;
        }

        EVAL_T GetHarvest_rate() const {
            return harvest_rate;
        }

        void SetHarvest_rate(EVAL_T harvest_rate) {
            this->harvest_rate = harvest_rate;
        }

        std::vector<EVAL_T> GetInit_F() const {
            return init_F;
        }

        void SetInit_F(std::vector<EVAL_T> init_F) {
            this->init_F = init_F;
        }

        std::vector<EVAL_T> GetInit_F_Like() const {
            return init_F_Like;
        }

        void SetInit_F_Like(std::vector<EVAL_T> init_F_Like) {
            this->init_F_Like = init_F_Like;
        }

        std::vector<int> GetInteger_control_flags() const {
            return integer_control_flags;
        }

        void SetInteger_control_flags(std::vector<int> integer_control_flags) {
            this->integer_control_flags = integer_control_flags;
        }

        EVAL_T GetJoin1() const {
            return join1;
        }

        void SetJoin1(EVAL_T join1) {
            this->join1 = join1;
        }

        EVAL_T GetJoin2() const {
            return join2;
        }

        void SetJoin2(EVAL_T join2) {
            this->join2 = join2;
        }

        EVAL_T GetJoin3() const {
            return join3;
        }

        void SetJoin3(EVAL_T join3) {
            this->join3 = join3;
        }

        EVAL_T GetLast_calc() const {
            return last_calc;
        }

        void SetLast_calc(EVAL_T last_calc) {
            this->last_calc = last_calc;
        }

        EVAL_T GetLast_objfun() const {
            return last_objfun;
        }

        void SetLast_objfun(EVAL_T last_objfun) {
            this->last_objfun = last_objfun;
        }

        std::vector<EVAL_T> GetLength_like() const {
            return length_like;
        }

        void SetLength_like(std::vector<EVAL_T> length_like) {
            this->length_like = length_like;
        }

        EVAL_T GetLikelihood_function_value() const {
            return likelihood_function_value;
        }

        void SetLikelihood_function_value(EVAL_T likelihood_function_value) {
            this->likelihood_function_value = likelihood_function_value;
        }

        std::valarray<std::valarray<EVAL_T> > GetLog_q_cr() const {
            return log_q_cr;
        }

        void SetLog_q_cr(std::valarray<std::valarray<EVAL_T> > log_q_cr) {
            this->log_q_cr = log_q_cr;
        }

        std::vector<EVAL_T> GetMat_age() const {
            return mat_age;
        }

        void SetMat_age(std::vector<EVAL_T> mat_age) {
            this->mat_age = mat_age;
        }

        std::vector<EVAL_T> GetMat_len() const {
            return mat_len;
        }

        void SetMat_len(std::vector<EVAL_T> mat_len) {
            this->mat_len = mat_len;
        }

        std::vector<EVAL_T> GetMat_len_wt() const {
            return mat_len_wt;
        }

        void SetMat_len_wt(std::vector<EVAL_T> mat_len_wt) {
            this->mat_len_wt = mat_len_wt;
        }

        std::vector<EVAL_T> GetMaxF() const {
            return maxF;
        }

        void SetMaxF(std::vector<EVAL_T> maxF) {
            this->maxF = maxF;
        }

        EVAL_T GetMaxpossF() const {
            return maxpossF;
        }

        void SetMaxpossF(EVAL_T maxpossF) {
            this->maxpossF = maxpossF;
        }

        std::vector<EVAL_T> GetMgp_adj() const {
            return mgp_adj;
        }

        void SetMgp_adj(std::vector<EVAL_T> mgp_adj) {
            this->mgp_adj = mgp_adj;
        }

        std::valarray<std::valarray<EVAL_T> > GetMgp_save() const {
            return mgp_save;
        }

        void SetMgp_save(std::valarray<std::valarray<EVAL_T> > mgp_save) {
            this->mgp_save = mgp_save;
        }

        std::vector<EVAL_T> GetMigrrate() const {
            return migrrate;
        }

        void SetMigrrate(std::vector<EVAL_T> migrrate) {
            this->migrrate = migrrate;
        }

        std::vector<EVAL_T> GetMnwt_like() const {
            return mnwt_like;
        }

        void SetMnwt_like(std::vector<EVAL_T> mnwt_like) {
            this->mnwt_like = mnwt_like;
        }

        std::vector<EVAL_T> GetNatM() const {
            return natM;
        }

        void SetNatM(std::vector<EVAL_T> natM) {
            this->natM = natM;
        }

        std::valarray<std::valarray<EVAL_T> > GetNatMparms() const {
            return natMparms;
        }

        void SetNatMparms(std::valarray<std::valarray<EVAL_T> > natMparms) {
            this->natMparms = natMparms;
        }

        std::vector<EVAL_T> GetNatage() const {
            return natage;
        }

        void SetNatage(std::vector<EVAL_T> natage) {
            this->natage = natage;
        }

        std::vector<EVAL_T> GetNatage_temp() const {
            return natage_temp;
        }

        void SetNatage_temp(std::vector<EVAL_T> natage_temp) {
            this->natage_temp = natage_temp;
        }

        std::valarray<std::valarray<EVAL_T> > GetNeff_a() const {
            return neff_a;
        }

        void SetNeff_a(std::valarray<std::valarray<EVAL_T> > neff_a) {
            this->neff_a = neff_a;
        }

        std::valarray<std::valarray<EVAL_T> > GetNeff_l() const {
            return neff_l;
        }

        void SetNeff_l(std::valarray<std::valarray<EVAL_T> > neff_l) {
            this->neff_l = neff_l;
        }

        std::vector<EVAL_T> GetOffset_a() const {
            return offset_a;
        }

        void SetOffset_a(std::vector<EVAL_T> offset_a) {
            this->offset_a = offset_a;
        }

        std::vector<EVAL_T> GetOffset_l() const {
            return offset_l;
        }

        void SetOffset_l(std::vector<EVAL_T> offset_l) {
            this->offset_l = offset_l;
        }

        EVAL_T GetOverdisp() const {
            return overdisp;
        }

        void SetOverdisp(EVAL_T overdisp) {
            this->overdisp = overdisp;
        }

        EVAL_T GetParm_dev_like() const {
            return parm_dev_like;
        }

        void SetParm_dev_like(EVAL_T parm_dev_like) {
            this->parm_dev_like = parm_dev_like;
        }

        EVAL_T GetParm_like() const {
            return parm_like;
        }

        void SetParm_like(EVAL_T parm_like) {
            this->parm_like = parm_like;
        }

        EVAL_T GetPeak() const {
            return peak;
        }

        void SetPeak(EVAL_T peak) {
            this->peak = peak;
        }

        EVAL_T GetPeak2() const {
            return peak2;
        }

        void SetPeak2(EVAL_T peak2) {
            this->peak2 = peak2;
        }

        std::vector<EVAL_T> GetPhase_output() const {
            return phase_output;
        }

        void SetPhase_output(std::vector<EVAL_T> phase_output) {
            this->phase_output = phase_output;
        }

        EVAL_T GetPoint1() const {
            return point1;
        }

        void SetPoint1(EVAL_T point1) {
            this->point1 = point1;
        }

        EVAL_T GetPoint2() const {
            return point2;
        }

        void SetPoint2(EVAL_T point2) {
            this->point2 = point2;
        }

        EVAL_T GetPoint3() const {
            return point3;
        }

        void SetPoint3(EVAL_T point3) {
            this->point3 = point3;
        }

        EVAL_T GetPoint4() const {
            return point4;
        }

        void SetPoint4(EVAL_T point4) {
            this->point4 = point4;
        }

        std::vector<EVAL_T> GetPred_rec() const {
            return pred_rec;
        }

        void SetPred_rec(std::vector<EVAL_T> pred_rec) {
            this->pred_rec = pred_rec;
        }

        EVAL_T GetPrior_function_value() const {
            return prior_function_value;
        }

        void SetPrior_function_value(EVAL_T prior_function_value) {
            this->prior_function_value = prior_function_value;
        }

        std::valarray<std::valarray<EVAL_T> > GetQ_cr() const {
            return q_cr;
        }

        void SetQ_cr(std::valarray<std::valarray<EVAL_T> > q_cr) {
            this->q_cr = q_cr;
        }

        std::vector<EVAL_T> GetRecdev() const {
            return recdev;
        }

        void SetRecdev(std::vector<EVAL_T> recdev) {
            this->recdev = recdev;
        }

        std::vector<EVAL_T> GetRecdev1() const {
            return recdev1;
        }

        void SetRecdev1(std::vector<EVAL_T> recdev1) {
            this->recdev1 = recdev1;
        }

        std::vector<EVAL_T> GetRecdev2() const {
            return recdev2;
        }

        void SetRecdev2(std::vector<EVAL_T> recdev2) {
            this->recdev2 = recdev2;
        }

        std::vector<EVAL_T> GetRecdev_cycle_Like() const {
            return recdev_cycle_Like;
        }

        void SetRecdev_cycle_Like(std::vector<EVAL_T> recdev_cycle_Like) {
            this->recdev_cycle_Like = recdev_cycle_Like;
        }

        std::vector<EVAL_T> GetRecdev_cycle_parm() const {
            return recdev_cycle_parm;
        }

        void SetRecdev_cycle_parm(std::vector<EVAL_T> recdev_cycle_parm) {
            this->recdev_cycle_parm = recdev_cycle_parm;
        }

        std::vector<EVAL_T> GetRecdev_early() const {
            return recdev_early;
        }

        void SetRecdev_early(std::vector<EVAL_T> recdev_early) {
            this->recdev_early = recdev_early;
        }

        std::vector<EVAL_T> GetRecr_dist() const {
            return recr_dist;
        }

        void SetRecr_dist(std::vector<EVAL_T> recr_dist) {
            this->recr_dist = recr_dist;
        }

        EVAL_T GetRecr_like() const {
            return recr_like;
        }

        void SetRecr_like(EVAL_T recr_like) {
            this->recr_like = recr_like;
        }

        std::vector<EVAL_T> GetRecr_std() const {
            return recr_std;
        }

        void SetRecr_std(std::vector<EVAL_T> recr_std) {
            this->recr_std = recr_std;
        }

        std::vector<EVAL_T> GetRetain() const {
            return retain;
        }

        void SetRetain(std::vector<EVAL_T> retain) {
            this->retain = retain;
        }

        std::vector<EVAL_T> GetRetain_M() const {
            return retain_M;
        }

        void SetRetain_M(std::vector<EVAL_T> retain_M) {
            this->retain_M = retain_M;
        }

        EVAL_T GetRho() const {
            return rho;
        }

        void SetRho(EVAL_T rho) {
            this->rho = rho;
        }

        std::vector<EVAL_T> GetSave_sel_fec() const {
            return save_sel_fec;
        }

        void SetSave_sel_fec(std::vector<EVAL_T> save_sel_fec) {
            this->save_sel_fec = save_sel_fec;
        }

        std::vector<EVAL_T> GetSave_sp_len() const {
            return save_sp_len;
        }

        void SetSave_sp_len(std::vector<EVAL_T> save_sp_len) {
            this->save_sp_len = save_sp_len;
        }

        EVAL_T GetSd_offset_rec() const {
            return sd_offset_rec;
        }

        void SetSd_offset_rec(EVAL_T sd_offset_rec) {
            this->sd_offset_rec = sd_offset_rec;
        }

        std::valarray<std::valarray<EVAL_T> > GetSe_cr_use() const {
            return se_cr_use;
        }

        void SetSe_cr_use(std::valarray<std::valarray<EVAL_T> > se_cr_use) {
            this->se_cr_use = se_cr_use;
        }

        std::vector<EVAL_T> GetSel() const {
            return sel;
        }

        void SetSel(std::vector<EVAL_T> sel) {
            this->sel = sel;
        }

        std::vector<EVAL_T> GetSel_a() const {
            return sel_a;
        }

        void SetSel_a(std::vector<EVAL_T> sel_a) {
            this->sel_a = sel_a;
        }

        std::vector<EVAL_T> GetSel_al_1() const {
            return sel_al_1;
        }

        void SetSel_al_1(std::vector<EVAL_T> sel_al_1) {
            this->sel_al_1 = sel_al_1;
        }

        std::vector<EVAL_T> GetSel_al_2() const {
            return sel_al_2;
        }

        void SetSel_al_2(std::vector<EVAL_T> sel_al_2) {
            this->sel_al_2 = sel_al_2;
        }

        std::vector<EVAL_T> GetSel_al_3() const {
            return sel_al_3;
        }

        void SetSel_al_3(std::vector<EVAL_T> sel_al_3) {
            this->sel_al_3 = sel_al_3;
        }

        std::vector<EVAL_T> GetSel_al_4() const {
            return sel_al_4;
        }

        void SetSel_al_4(std::vector<EVAL_T> sel_al_4) {
            this->sel_al_4 = sel_al_4;
        }

        std::vector<EVAL_T> GetSel_l() const {
            return sel_l;
        }

        void SetSel_l(std::vector<EVAL_T> sel_l) {
            this->sel_l = sel_l;
        }

        std::vector<EVAL_T> GetSel_l_r() const {
            return sel_l_r;
        }

        void SetSel_l_r(std::vector<EVAL_T> sel_l_r) {
            this->sel_l_r = sel_l_r;
        }

        std::vector<EVAL_T> GetSelparm() const {
            return selparm;
        }

        void SetSelparm(std::vector<EVAL_T> selparm) {
            this->selparm = selparm;
        }

        std::vector<EVAL_T> GetSelparm_Like() const {
            return selparm_Like;
        }

        void SetSelparm_Like(std::vector<EVAL_T> selparm_Like) {
            this->selparm_Like = selparm_Like;
        }

        std::valarray<std::valarray<EVAL_T> > GetSelparm_block_val() const {
            return selparm_block_val;
        }

        void SetSelparm_block_val(std::valarray<std::valarray<EVAL_T> > selparm_block_val) {
            this->selparm_block_val = selparm_block_val;
        }

        std::valarray<std::valarray<EVAL_T> > GetSelparm_dev() const {
            return selparm_dev;
        }

        void SetSelparm_dev(std::valarray<std::valarray<EVAL_T> > selparm_dev) {
            this->selparm_dev = selparm_dev;
        }

        std::valarray<std::valarray<EVAL_T> > GetSelparm_dev_rwalk() const {
            return selparm_dev_rwalk;
        }

        void SetSelparm_dev_rwalk(std::valarray<std::valarray<EVAL_T> > selparm_dev_rwalk) {
            this->selparm_dev_rwalk = selparm_dev_rwalk;
        }

        std::valarray<std::valarray<EVAL_T> > GetSelparm_trend() const {
            return selparm_trend;
        }

        void SetSelparm_trend(std::valarray<std::valarray<EVAL_T> > selparm_trend) {
            this->selparm_trend = selparm_trend;
        }

        EVAL_T GetSigmaR() const {
            return sigmaR;
        }

        void SetSigmaR(EVAL_T sigmaR) {
            this->sigmaR = sigmaR;
        }

        std::vector<EVAL_T> GetSizeage_like() const {
            return sizeage_like;
        }

        void SetSizeage_like(std::vector<EVAL_T> sizeage_like) {
            this->sizeage_like = sizeage_like;
        }

        std::valarray<std::valarray<EVAL_T> > GetSmry() const {
            return smry;
        }

        void SetSmry(std::valarray<std::valarray<EVAL_T> > smry) {
            this->smry = smry;
        }

        EVAL_T GetSmrybio() const {
            return smrybio;
        }

        void SetSmrybio(EVAL_T smrybio) {
            this->smrybio = smrybio;
        }

        EVAL_T GetSmrynum() const {
            return smrynum;
        }

        void SetSmrynum(EVAL_T smrynum) {
            this->smrynum = smrynum;
        }

        std::vector<EVAL_T> GetSurv1() const {
            return surv1;
        }

        void SetSurv1(std::vector<EVAL_T> surv1) {
            this->surv1 = surv1;
        }

        std::vector<EVAL_T> GetSurv2() const {
            return surv2;
        }

        void SetSurv2(std::vector<EVAL_T> surv2) {
            this->surv2 = surv2;
        }

        std::vector<EVAL_T> GetSurv_like() const {
            return surv_like;
        }

        void SetSurv_like(std::vector<EVAL_T> surv_like) {
            this->surv_like = surv_like;
        }

        EVAL_T GetT1() const {
            return t1;
        }

        void SetT1(EVAL_T t1) {
            this->t1 = t1;
        }

        EVAL_T GetT2() const {
            return t2;
        }

        void SetT2(EVAL_T t2) {
            this->t2 = t2;
        }

        EVAL_T GetTemp() const {
            return temp;
        }

        void SetTemp(EVAL_T temp) {
            this->temp = temp;
        }

        EVAL_T GetTemp1() const {
            return temp1;
        }

        void SetTemp1(EVAL_T temp1) {
            this->temp1 = temp1;
        }

        EVAL_T GetTemp2() const {
            return temp2;
        }

        void SetTemp2(EVAL_T temp2) {
            this->temp2 = temp2;
        }

        EVAL_T GetTemp3() const {
            return temp3;
        }

        void SetTemp3(EVAL_T temp3) {
            this->temp3 = temp3;
        }

        EVAL_T GetTemp4() const {
            return temp4;
        }

        void SetTemp4(EVAL_T temp4) {
            this->temp4 = temp4;
        }

        std::vector<EVAL_T> GetTempvec() const {
            return tempvec;
        }

        void SetTempvec(std::vector<EVAL_T> tempvec) {
            this->tempvec = tempvec;
        }

        std::vector<EVAL_T> GetTempvec_l() const {
            return tempvec_l;
        }

        void SetTempvec_l(std::vector<EVAL_T> tempvec_l) {
            this->tempvec_l = tempvec_l;
        }

        EVAL_T GetTiming() const {
            return timing;
        }

        void SetTiming(EVAL_T timing) {
            this->timing = timing;
        }

        EVAL_T GetTotbio() const {
            return totbio;
        }

        void SetTotbio(EVAL_T totbio) {
            this->totbio = totbio;
        }

        EVAL_T GetTwo_sigmaRsq() const {
            return two_sigmaRsq;
        }

        void SetTwo_sigmaRsq(EVAL_T two_sigmaRsq) {
            this->two_sigmaRsq = two_sigmaRsq;
        }

        EVAL_T GetUpselex() const {
            return upselex;
        }

        void SetUpselex(EVAL_T upselex) {
            this->upselex = upselex;
        }

        std::vector<EVAL_T> GetUse_rec() const {
            return use_rec;
        }

        void SetUse_rec(std::vector<EVAL_T> use_rec) {
            this->use_rec = use_rec;
        }

        EVAL_T GetVbio() const {
            return vbio;
        }

        void SetVbio(EVAL_T vbio) {
            this->vbio = vbio;
        }

        std::valarray<std::valarray<EVAL_T> > GetVirg_fec() const {
            return virg_fec;
        }

        void SetVirg_fec(std::valarray<std::valarray<EVAL_T> > virg_fec) {
            this->virg_fec = virg_fec;
        }

        std::vector<EVAL_T> GetWt_len() const {
            return wt_len;
        }

        void SetWt_len(std::vector<EVAL_T> wt_len) {
            this->wt_len = wt_len;
        }

        std::valarray<std::valarray<EVAL_T> > GetWt_len2() const {
            return wt_len2;
        }

        void SetWt_len2(std::valarray<std::valarray<EVAL_T> > wt_len2) {
            this->wt_len2 = wt_len2;
        }

        std::valarray<std::valarray<EVAL_T> > GetWt_len2_sq() const {
            return wt_len2_sq;
        }

        void SetWt_len2_sq(std::valarray<std::valarray<EVAL_T> > wt_len2_sq) {
            this->wt_len2_sq = wt_len2_sq;
        }

        std::valarray<std::valarray<EVAL_T> > GetWt_len_fd() const {
            return wt_len_fd;
        }

        void SetWt_len_fd(std::valarray<std::valarray<EVAL_T> > wt_len_fd) {
            this->wt_len_fd = wt_len_fd;
        }

        std::valarray<std::valarray<EVAL_T> > GetWt_len_low() const {
            return wt_len_low;
        }

        void SetWt_len_low(std::valarray<std::valarray<EVAL_T> > wt_len_low) {
            this->wt_len_low = wt_len_low;
        }

        std::vector<EVAL_T> GetWtlen_p() const {
            return wtlen_p;
        }

        void SetWtlen_p(std::vector<EVAL_T> wtlen_p) {
            this->wtlen_p = wtlen_p;
        }

        std::valarray<std::valarray<EVAL_T> > GetWtlen_seas() const {
            return wtlen_seas;
        }

        void SetWtlen_seas(std::valarray<std::valarray<EVAL_T> > wtlen_seas) {
            this->wtlen_seas = wtlen_seas;
        }






    };


}



#endif	/* STOCKSYNTHESISMODEL_HPP */

