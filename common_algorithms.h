#ifndef _common_algorithms_h_
#define _common_algorithms_h_

#include "TMath.h"

#include "common_definitions.h"

//----------------------------------------------------------------------------------------------------

Kinematics DoReconstruction(const HitData &h, const Environment & env)
{
	Kinematics k;

	// ----- single-arm kinematics reconstruction -----

	double D_x_L = - env.v_x_L_1_F * env.L_x_L_2_F + env.v_x_L_2_F * env.L_x_L_1_F;
	k.th_x_L = (env.v_x_L_1_F * h.L_2_F.x - env.v_x_L_2_F * h.L_1_F.x) / D_x_L;
	k.vtx_x_L = (- h.L_1_F.x * env.L_x_L_2_F + h.L_2_F.x * env.L_x_L_1_F) / D_x_L;

	double D_x_R = + env.v_x_R_1_F * env.L_x_R_2_F - env.v_x_R_2_F * env.L_x_R_1_F;
	k.th_x_R = (env.v_x_R_1_F * h.R_2_F.x - env.v_x_R_2_F * h.R_1_F.x) / D_x_R;
	k.vtx_x_R = (+ h.R_1_F.x * env.L_x_R_2_F - h.R_2_F.x * env.L_x_R_1_F) / D_x_R;

	k.th_y_L_1_F = - h.L_1_F.y / env.L_y_L_1_F;
  	k.th_y_L_2_F = - h.L_2_F.y / env.L_y_L_2_F;
  	k.th_y_L = (k.th_y_L_1_F + k.th_y_L_2_F) / 2.;
  	
	k.th_y_R_1_F = + h.R_1_F.y / env.L_y_R_1_F;
  	k.th_y_R_2_F = + h.R_2_F.y / env.L_y_R_2_F;
  	k.th_y_R = (k.th_y_R_1_F + k.th_y_R_2_F) / 2.;
	
	double D_y_L = - env.v_y_L_1_F * env.L_y_L_2_F + env.v_y_L_2_F * env.L_y_L_1_F;
	//k.th_y_L = (env.v_y_L_1_F * h.L_2_F.y - env.v_y_L_2_F * h.L_1_F.y) / D_y_L;
	k.vtx_y_L = (- h.L_1_F.y * env.L_y_L_2_F + h.L_2_F.y * env.L_y_L_1_F) / D_y_L;

	double D_y_R = + env.v_y_R_1_F * env.L_y_R_2_F - env.v_y_R_2_F * env.L_y_R_1_F;
	//k.th_y_R = (env.v_y_R_1_F * h.R_2_F.y - env.v_y_R_2_F * h.R_1_F.y) / D_y_R;
	k.vtx_y_R = (+ h.R_1_F.y * env.L_y_R_2_F - h.R_2_F.y * env.L_y_R_1_F) / D_y_R;

	// ----- double-arm kinematics reconstruction -----
	
	k.th_x = (k.th_x_L + k.th_x_R) / 2.;
	k.th_y = (k.th_y_L + k.th_y_R) / 2.;
	
	k.vtx_x = (k.vtx_x_L + k.vtx_x_R) / 2.;
	k.vtx_y = (k.vtx_y_L + k.vtx_y_R) / 2.;

	// ----- theta reconstruction -----
	double th_sq = k.th_x*k.th_x + k.th_y*k.th_y;
	k.th = sqrt(th_sq);
	k.phi = atan2(k.th_y, k.th_x);

	// t reconstruction
	k.t_x = env.p*env.p * k.th_x * k.th_x;
	k.t_y = env.p*env.p * k.th_y * k.th_y;
	k.t = k.t_x + k.t_y;

	return k;
}

//----------------------------------------------------------------------------------------------------

void BuildBinning(const Analysis &anal, const string &type, double* &binEdges, unsigned int &bins)
{
	bool verbose = false;

	if (verbose)
		printf(">> BuildBinning(%s)\n", type.c_str());

	std::vector<double> be;

	// low-|t| part
	double t = anal.t_min;

	// central part
	if (type.compare("ub") == 0)
	{
		double w = 2E-3;
		for (; t <= anal.t_max; t += w)
			be.push_back(t);
	}

	if (type.compare("ub-90m") == 0)
	{
		unsigned int N_bins_cen = 200;
		double t_min_90 = 0.;
		double t_max_90 = 1.4;
		double w = (t_max_90 - t_min_90) / N_bins_cen;
		for (unsigned int i = 0; i < N_bins_cen; i++)
		{
			t = t_min_90 + w * i;
			if (t > anal.t_max)
				break;
			be.push_back(t);
		}
	}
	
	if (type.compare("eb") == 0)
	{
		double w = 0.2E-3;
		for (; t < 0.6E-3; t += w)
			be.push_back(t);
		
		w = 0.2E-3;
		for (; t < 0.02; t += w)
			be.push_back(t);

		for (; t < anal.t_max; t += w)
		{
			be.push_back(t);
			w *= 1.05;
		}
	}

	if (type.compare("ob-0-1") == 0)
	{
		double t = anal.t_min;
		while (t < anal.t_max)
		{
			be.push_back(t);

			double w = 0.01;
			if (t >= 0.00) w = 2E-4;
			if (t >= 5E-4) w = (t - 0.00) / (0.01 - 0.00) * (0.0037 - 0.0001) + 0.0001;
			if (t >= 0.01) w = (t - 0.01) / (0.06 - 0.01) * (0.0100 - 0.0037) + 0.0037;
			if (t >= 0.06) w = 0.01;

			t += w;
		}
	}

	if (type.compare("ob-1-4") == 0)
	{
		double t = anal.t_min;
		while (t < anal.t_max)
		{
			be.push_back(t);

			double w = 0.01;
			if (t >= 0.00) w = 2E-4;
			if (t >= 5E-4) w = 2.68E-3 * sqrt(t);
			if (t >= 0.06) w = 1.37E-4 * exp(26.82 * t);
			if (t >= 0.16) w = 0.01;

			t += w;
		}
	}

	if (type.compare("ob-3-4") == 0)
	{
		double t = anal.t_min;
		while (t < anal.t_max)
		{
			be.push_back(t);

			double w = 0.01;
			if (t >= 0.00) w = 2E-4;
			if (t >= 5E-4) w = 2.68E-3 * sqrt(t);
			if (t >= 2E-3) w = 2.68E-3 * (1. + (t - 2E-3) / 2E-3 * 2.) * sqrt(t);
			if (t >= 4E-3) w = 2.68E-3*3. * sqrt(t);
			if (t >= 0.115) w = 1.37E-4 * exp(26.82 * t);
			if (t >= 0.16) w = 0.01;

			t += w;
		}
	}

	if (type.compare("ob-5-4") == 0)
	{
		double t = anal.t_min;
		while (t < anal.t_max)
		{
			be.push_back(t);

			double w = 0.01;
			if (t >= 0.00) w = 2E-4;
			if (t >= 5E-4) w = 2.68E-3 * sqrt(t);
			if (t >= 2E-3) w = 2.68E-3 * (1. + (t - 2E-3) / 2E-3 * 4.) * sqrt(t);
			if (t >= 4E-3) w = 2.68E-3*5. * sqrt(t);
			if (t >= 0.14) w = 1.37E-4 * exp(26.82 * t);
			if (t >= 0.16) w = 0.01;

			t += w;
		}
	}

	if (type.compare("ob-0.2-10") == 0)
	{
		double t = anal.t_min;
		while (t < anal.t_max)
		{
			be.push_back(t);
			
			double w = 0.01;
			if (t >= 0.00) w = 2E-4;
			if (t >= 5E-4) w = 2.68E-3 * 0.2 * sqrt(t);
			if (t >= 0.077) w = 1E-4 * exp(-6.56*t + 150.*t*t);
			if (t >= 0.19) w = 0.01;

			t += w;
		}
	}

	if (type.compare("ob-0.5-10") == 0)
	{
		double t = anal.t_min;
		while (t < anal.t_max)
		{
			be.push_back(t);
			
			double w = 0.01;
			if (t >= 0.00) w = 2E-4;
			if (t >= 5E-4) w = 2.68E-3 * 0.5 * sqrt(t);
			if (t >= 0.126) w = 1E-4 * exp(-6.56*t + 150.*t*t);
			if (t >= 0.19) w = 0.01;

			t += w;
		}
	}

	if (type.compare("ob-1-10") == 0)
	{
		double t = anal.t_min;
		while (t < anal.t_max)
		{
			be.push_back(t);
			
			double w = 0.01;
			if (t >= 0.00) w = 2E-4;
			if (t >= 5E-4) w = 2.68E-3 * sqrt(t);
			if (t >= 0.15) w = 1E-4 * exp(-6.56*t + 150.*t*t);
			if (t >= 0.19) w = 0.01;

			t += w;
		}
	}

	// high-|t| part
	unsigned int N_bins_high = 10;
	double w = (anal.t_max_full - anal.t_max) / N_bins_high;
	for (unsigned int i = 0; i <= N_bins_high; i++)
		be.push_back(anal.t_max + w * i);

	bins = be.size() - 1;
	binEdges = new double[be.size()];
	for (unsigned int i = 0; i < be.size(); i++)
	{
		binEdges[i] = be[i];
		if (verbose)
			printf("\tbi = %u: %.2E\n", i, binEdges[i]);
	}
}

//----------------------------------------------------------------------------------------------------

bool CalculateAcceptanceCorrections(double th_y_sign,
		const Kinematics &k, const Analysis &anal,
		double &phi_corr, double &div_corr)
{
	// ---------- smearing component ----------

	/*
	if ((k.th_x_L < anal.th_x_lcut_L) || (k.th_x_R < anal.th_x_lcut_R) || (k.th_x_L > anal.th_x_hcut_L) || (k.th_x_R > anal.th_x_hcut_R))
		return true;
	*/

	if ((th_y_sign * k.th_y_L < anal.th_y_lcut_L) || (th_y_sign * k.th_y_R < anal.th_y_lcut_R)
		|| (th_y_sign * k.th_y_L > anal.th_y_hcut_L) || (th_y_sign * k.th_y_R > anal.th_y_hcut_R))
		return true;
	
	/*
	double LB_x_L = anal.th_x_lcut_L - k.th_x, UB_x_L = anal.th_x_hcut_L - k.th_x;
	double LB_x_R = anal.th_x_lcut_R - k.th_x, UB_x_R = anal.th_x_hcut_R - k.th_x;
	double F_x_L = (UB_x_L > LB_x_L) ? ( TMath::Erf(UB_x_L / anal.si_th_x_1arm_L / sqrt(2.)) - TMath::Erf(LB_x_L / anal.si_th_x_1arm_L / sqrt(2.)) ) / 2. : 0.;
	double F_x_R = (UB_x_R > LB_x_R) ? ( TMath::Erf(UB_x_R / anal.si_th_x_1arm_R / sqrt(2.)) - TMath::Erf(LB_x_R / anal.si_th_x_1arm_R / sqrt(2.)) ) / 2. : 0.;
	double F_x = F_x_L * F_x_R;
	*/
	double F_x = 1.;

	double th_y_abs = th_y_sign * k.th_y;

	double UB_y = min(anal.th_y_hcut_R - th_y_abs, th_y_abs - anal.th_y_lcut_L);
	double LB_y = max(anal.th_y_lcut_R - th_y_abs, th_y_abs - anal.th_y_hcut_L);
	double F_y = (UB_y > LB_y) ? ( TMath::Erf(UB_y / anal.si_th_y_1arm) - TMath::Erf(LB_y / anal.si_th_y_1arm) ) / 2. : 0.;

	//printf(">> F_x_L = %E, F_x_R = %E, F_y = %E\n", F_x_L, F_x_R, F_y);

	div_corr = 1./(F_x * F_y);

	// ---------- phi component ----------
	
	// apply safety margins to avoid excessive smearing component
	//double th_x_lcut = max(anal.th_x_lcut_L, anal.th_x_lcut_R) + 3.0E-6;
	//double th_x_hcut = min(anal.th_x_hcut_L, anal.th_x_hcut_R) - 3.0E-6;
	double th_x_lcut = anal.th_x_lcut;
	double th_x_hcut = anal.th_x_hcut;

	//double th_y_lcut = max(anal.th_y_lcut_L, anal.th_y_lcut_R) + 0.2E-6;
	//double th_y_hcut = min(anal.th_y_hcut_L, anal.th_y_hcut_R) - 1.0E-6;

	double th_y_lcut = anal.th_y_lcut;
	double th_y_hcut = anal.th_y_hcut;

	//printf("-------------------------\n");
	//printf("k.th_x=%E, th_x_lcut=%E, th_x_hcut=%E\n", k.th_x, th_x_lcut, th_x_hcut);
	//printf("k.th_y=%E, th_y_lcut=%E, th_y_hcut=%E\n", k.th_y, th_y_lcut, th_y_hcut);

	if (k.th_x <= th_x_lcut || k.th_x >= th_x_hcut || th_y_abs <= th_y_lcut || th_y_abs >= th_y_hcut)
		return true;

	// get all intersections
	set<double> phis;

	if (k.th > th_y_lcut)
	{
		double phi = asin(th_y_lcut / k.th);
		double ta_x = k.th * cos(phi);
		if (th_x_lcut < ta_x && ta_x < th_x_hcut)
			phis.insert(phi);
		if (th_x_lcut < -ta_x && -ta_x < th_x_hcut)
			phis.insert(M_PI - phi);
	}
	
	if (k.th > th_y_hcut)
	{
		double phi = asin(th_y_hcut / k.th);
		double ta_x = k.th * cos(phi);
		if (th_x_lcut < ta_x && ta_x < th_x_hcut)
			phis.insert(phi);
		if (th_x_lcut < -ta_x && -ta_x < th_x_hcut)
			phis.insert(M_PI - phi);
	}

	if (k.th > fabs(th_x_hcut))
	{
		double phi = acos(fabs(th_x_hcut) / k.th);
		double ta_y = k.th * sin(phi);
		if (th_y_lcut < ta_y && ta_y < th_y_hcut)
			phis.insert(phi);
	}

	if (k.th > fabs(th_x_lcut))
	{
		double phi = acos(fabs(th_x_lcut) / k.th);
		double ta_y = k.th * sin(phi);
		if (th_y_lcut < ta_y && ta_y < th_y_hcut)
			phis.insert(M_PI - phi);
	}

	// the number of intersections must be even
	if ((phis.size() % 2) == 1)
	{
		printf("ERROR: odd number of intersections in acceptance calculation\n");
	}

	// no intersection => no acceptances
	if (phis.size() == 0)
		return true;

	// calculate arc-length in within acceptance
	double phiSum = 0.;
	for (set<double>::iterator it = phis.begin(); it != phis.end(); ++it)
	{
		double phi_start = *it;
		++it;
		double phi_end = *it;

		phiSum += phi_end - phi_start;
	}
	
	phi_corr = 2. * M_PI / phiSum;

	return false;
}

//----------------------------------------------------------------------------------------------------

bool SkipRun(unsigned int /*run*/)
{
	return false;
}

//----------------------------------------------------------------------------------------------------

// map: run number (8372) --> list of triggered bunches
typedef std::map<unsigned int, std::vector<unsigned int> > BunchMap;

bool keepAllBunches;
BunchMap bunchMap;

bool SkipBunch(unsigned int run, unsigned bunch)
{
	if (keepAllBunches)
		return false;

	const std::vector<unsigned int> &bunches = bunchMap[run];
	
	return (find(bunches.begin(), bunches.end(), bunch) == bunches.end());
}

//----------------------------------------------------------------------------------------------------

// returns the beam for which the bunch is non-colliding
// for colliding bunches returns zero
unsigned int NonCollidingBunch(unsigned int /*run*/, unsigned /*bunch*/)
{
	/*
	if (run == 8318) {
		if (bunch == 994)
			return 1;
		if (bunch == 991)
			return 2;
	}
	
	if (run >= 8333 && run <= 8341)
	{
		if (bunch == 900)
			return 1;
		if (bunch == 991)
			return 2;
	}

	if (run >= 8367 && run <= 8372)
	{
		if (bunch == 3104 || bunch == 3130 || bunch == 3156 || bunch == 3078)
			return 1;
		if (bunch == 3143 || bunch == 3169 || bunch == 3195 || bunch == 3117)
			return 2;
	}
	*/

	return 0;
}

//----------------------------------------------------------------------------------------------------

bool IsZeroBias(unsigned int trigger, unsigned int /*run*/, unsigned int /*event*/)
{
	return ((trigger & 512) != 0);
}

//----------------------------------------------------------------------------------------------------

HitData ProtonTransport(const Kinematics & k, const Environment & env)
{
	HitData h;

	h.L_2_F.x = -env.L_x_L_2_F*k.th_x_L + env.v_x_L_2_F*k.vtx_x;
	h.L_2_F.y = -env.L_y_L_2_F*k.th_y_L + env.v_y_L_2_F*k.vtx_y;
                                                               
	h.L_1_F.x = -env.L_x_L_1_F*k.th_x_L + env.v_x_L_1_F*k.vtx_x;
	h.L_1_F.y = -env.L_y_L_1_F*k.th_y_L + env.v_y_L_1_F*k.vtx_y;
	                                                           
	h.R_1_F.x = +env.L_x_R_1_F*k.th_x_R + env.v_x_R_1_F*k.vtx_x;
	h.R_1_F.y = +env.L_y_R_1_F*k.th_y_R + env.v_y_R_1_F*k.vtx_y;
	                                                           
	h.R_2_F.x = +env.L_x_R_2_F*k.th_x_R + env.v_x_R_2_F*k.vtx_x;
	h.R_2_F.y = +env.L_y_R_2_F*k.th_y_R + env.v_y_R_2_F*k.vtx_y;

	return h;
}

#endif
