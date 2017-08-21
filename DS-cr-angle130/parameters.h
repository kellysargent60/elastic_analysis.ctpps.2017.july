#include <string>
#include <vector>
#include <map>
#include <cmath>

double timestamp0 = 1499464800;

vector<string> distilledNtuples;

vector<AlignmentSource> alignmentSources;
Analysis anal;
Environment env;

string unsmearing_file;
string unsmearing_object;

string luminosity_data_file;

void Init_base()
{
	// list of (sub-)directories with distilled ntuples
	distilledNtuples.push_back(".");

	// selection of bunches
	keepAllBunches = true;

	// alignment settings
	/*
	AlignmentSource alSrc;
	alSrc.SetAlignmentA(atConstant);
	alSrc.SetAlignmentB(atConstant);
	alSrc.SetAlignmentC(atConstant);

	alSrc.cnst.a_L_1_F = 0E-3; alSrc.cnst.b_L_1_F = 0.; alSrc.cnst.c_L_1_F = 0E-3;
	alSrc.cnst.a_L_1_N = 0E-3; alSrc.cnst.b_L_1_N = 0.; alSrc.cnst.c_L_1_N = 0E-3;
                                                                       
	alSrc.cnst.a_R_1_N = 0E-3; alSrc.cnst.b_R_1_N = 0.; alSrc.cnst.c_R_1_N = 0E-3;
	alSrc.cnst.a_R_1_F = 0E-3; alSrc.cnst.b_R_1_F = 0.; alSrc.cnst.c_R_1_F = 0E-3;
	
	alignmentSources.push_back(alSrc);
	*/

	// environment settings
	env.InitNominal();

	// binning
	// TODO
	anal.t_min = 0.; anal.t_max = 0.5;
	anal.t_min_full = 0.; anal.t_max_full = 0.8;
	
	// approximate (time independent) resolutions
	// TODO
	anal.si_th_y_1arm = 1.7E-6 / sqrt(2.);
	anal.si_th_y_1arm_unc = 0.E-6 / sqrt(2.);

	anal.si_th_y_2arm = anal.si_th_y_1arm / sqrt(2.);
	anal.si_th_y_2arm_unc = 0E-6;

	anal.si_th_x_1arm_L = 8E-6;
	anal.si_th_x_1arm_R = 8E-6;
	anal.si_th_x_1arm_unc = 0E-6;

	anal.si_th_x_2arm = 0E-6;
	anal.si_th_x_2arm_unc = 0E-6;

	// analysis settings
	anal.th_x_lcut = -1.;	
	anal.th_x_hcut = +1.;
	
#if 0
	anal.use_time_dependent_resolutions = false;

	anal.use_3outof4_efficiency_fits = false;
	anal.use_pileup_efficiency_fits = false;
	anal.inefficiency_3outof4 = 0.;				// diagonal dependent
	anal.inefficiency_shower_near = 0.03;
	anal.inefficiency_pile_up = 0.;				// diagonal dependent
	anal.inefficiency_trigger = 0.;

	anal.bckg_corr = 1.;
	
	anal.L_int_eff = 0.;	// mb^-1, diagonal dependent
	
	anal.eff_th_y_min = 200E-6; // TODO
	
	anal.t_min_fit = 0.027; // TODO
#endif

	anal.alignment_t0 = 0.;			// beginning of the first time-slice
	anal.alignment_ts = 120.*60.;	// time-slice in s
	
	// TODO
	anal.alignmentYRanges["L_2_F"] = Analysis::AlignmentYRange(-40., +10, -10., +40);
	anal.alignmentYRanges["L_1_F"] = Analysis::AlignmentYRange(-40., +10, -10., +40);
	anal.alignmentYRanges["R_1_F"] = Analysis::AlignmentYRange(-40., +10, -10., +40);
	anal.alignmentYRanges["R_2_F"] = Analysis::AlignmentYRange(-40., +10, -10., +40);

#if 0
	// TODO
	unsmearing_file = "";	// diagonal dependent
	//unsmearing_object = "cf,<binning>/exp3/corr_final";
	//unsmearing_object = "cf,<binning>/exp3+exp4/corr_final";
	unsmearing_object = "ff";

	// TODO
	luminosity_data_file = "../fill_3549_lumiCalc2.py_V04-02-04_lumibylsXing.csv";
#endif
}

//----------------------------------------------------------------------------------------------------

void Init_45b_56t()
{
	// analysis settings
	anal.cut1_a = 1.; anal.cut1_c =  0.0005; anal.cut1_si = 5E-4;
	anal.cut2_a = 1.; anal.cut2_c =  5.1E-4; anal.cut2_si = 6.5E-5;

	anal.cut3_c = 0.46; anal.cut3_si = 0.041;
	anal.cut4_c = -0.043; anal.cut4_si = 0.026;
	
	anal.cut5_a = -0.14; anal.cut5_c = -0.55; anal.cut5_si = 0.038;
	anal.cut6_a = -0.0; anal.cut6_c = 0.14; anal.cut6_si = 10.0;

	anal.cut7_a = 0.0; anal.cut7_c = 0.51; anal.cut7_si = 10.0;

	anal.th_y_lcut_L = 13.5E-6; anal.th_y_lcut_R = 14E-6; anal.th_y_lcut = 14.5E-6;
	anal.th_y_hcut_L = 102E-6; anal.th_y_hcut_R = 102E-6; anal.th_y_hcut = 100E-6;

#if 0
	// TODO
	//unsmearing_file = "unfolding_fit_45b_56t_old.root";

	anal.inefficiency_3outof4 = 0.0; // TODO
	anal.inefficiency_pile_up = 0.0; // TODO

	anal.L_int_eff = 79.42E3;	// TODO	
#endif
}

//----------------------------------------------------------------------------------------------------

void Init_45t_56b()
{
	// analysis settings
	anal.cut1_a = 1.; anal.cut1_c =  6.1E-4; anal.cut1_si = 4.9E-4;
	anal.cut2_a = 1.; anal.cut2_c =  4.7E-4; anal.cut2_si = 2.9E-5;

	anal.cut3_c = 0.43; anal.cut3_si = 0.038;
	anal.cut4_c = -0.048; anal.cut4_si = 0.026;
	
	anal.cut5_a = -0.0; anal.cut5_c = -0.5; anal.cut5_si = 10.0;
	anal.cut6_a = -0.0; anal.cut6_c = -0.04; anal.cut6_si = 10.0;

	anal.cut7_a = 153.; anal.cut7_c = 0.; anal.cut7_si = 0.011;

	anal.th_y_lcut_L = 13.5E-6; anal.th_y_lcut_R = 14E-6; anal.th_y_lcut = 14.5E-6;
	anal.th_y_hcut_L = 102E-6; anal.th_y_hcut_R = 102E-6; anal.th_y_hcut = 100E-6;

#if 0

	// TODO
	//unsmearing_file = "unfolding_fit_45b_56t_old.root";

	anal.inefficiency_3outof4 = 0.0; // TODO
	anal.inefficiency_pile_up = 0.0; // TODO

	anal.L_int_eff = 79.42E3;	// TODO	
#endif
}
