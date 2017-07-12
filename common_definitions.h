#ifndef _common_definitions_h_
#define _common_definitions_h_

#include <string>
#include <vector>
#include <set>
#include <map>
#include <cmath>

#include "TGraph.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSymEigen.h"
#include "TRandom2.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

enum DiagonalType { dUnknown, d45b_56t, d45t_56b, dCombined, ad45b_56b, ad45t_56t };

DiagonalType diagonal = dUnknown;

double th_y_sign = 0.;

//----------------------------------------------------------------------------------------------------

struct AlignmentData
{
	//	a: xy coupling in rad
	//	b: x shift in mm
	//	c: y shift in mm
	double a_L_2_F, b_L_2_F, c_L_2_F;
	double a_L_1_F, b_L_1_F, c_L_1_F;
	double a_R_1_F, b_R_1_F, c_R_1_F;
	double a_R_2_F, b_R_2_F, c_R_2_F;

	AlignmentData()
	{
		a_L_2_F = b_L_2_F = c_L_2_F = 0.;
		a_L_1_F = b_L_1_F = c_L_1_F = 0.;
		a_R_1_F = b_R_1_F = c_R_1_F = 0.;
		a_R_2_F = b_R_2_F = c_R_2_F = 0.;
	}

	/*
	AlignmentData Interpolate(double s_N, double s_F, double s_NH, double s_FH) const
	{
		AlignmentData r;

		r.a_L_F = a_L_N + (a_L_F - a_L_N)/(s_F - s_N) * (s_FH - s_N); r.a_L_N = a_L_N + (a_L_F - a_L_N)/(s_F - s_N) * (s_NH - s_N);
		r.a_R_F = a_R_N + (a_R_F - a_R_N)/(s_F - s_N) * (s_FH - s_N); r.a_R_N = a_R_N + (a_R_F - a_R_N)/(s_F - s_N) * (s_NH - s_N);
                                                                          
		r.b_L_F = b_L_N + (b_L_F - b_L_N)/(s_F - s_N) * (s_FH - s_N); r.b_L_N = b_L_N + (b_L_F - b_L_N)/(s_F - s_N) * (s_NH - s_N);
		r.b_R_F = b_R_N + (b_R_F - b_R_N)/(s_F - s_N) * (s_FH - s_N); r.b_R_N = b_R_N + (b_R_F - b_R_N)/(s_F - s_N) * (s_NH - s_N);
                                                                          
		r.c_L_F = c_L_N + (c_L_F - c_L_N)/(s_F - s_N) * (s_FH - s_N); r.c_L_N = c_L_N + (c_L_F - c_L_N)/(s_F - s_N) * (s_NH - s_N);
		r.c_R_F = c_R_N + (c_R_F - c_R_N)/(s_F - s_N) * (s_FH - s_N); r.c_R_N = c_R_N + (c_R_F - c_R_N)/(s_F - s_N) * (s_NH - s_N);

		return r;
	}
	*/
};

//----------------------------------------------------------------------------------------------------

enum AlignmentType { atNone, atConstant, atTimeDependent };

struct AlignmentSource
{
	struct GraphSet
	{
		TGraph *L_2_F, *L_1_F, *R_1_F, *R_2_F;
		GraphSet() : L_2_F(NULL), L_1_F(NULL), R_1_F(NULL), R_2_F(NULL) {}
	} gs_a, gs_b, gs_c;

	AlignmentData cnst;

	AlignmentType type_a, type_b, type_c;
	string src_a, src_b, src_c;

	AlignmentSource() : type_a(atNone), type_b(atNone), type_c(atNone)
	{
	}

	void SetAlignmentA(AlignmentType t, const string &fn = "")
	{
		type_a = t;
		src_a = fn;
	}

	void SetAlignmentB(AlignmentType t, const string &fn = "")
	{
		type_b = t;
		src_b = fn;
	}

	void SetAlignmentC(AlignmentType t, const string &fn = "")
	{
		type_c = t;
		src_c = fn;
	}

	void InitOne(const string label, AlignmentType t, const string &fn, GraphSet &gs, const string &obj)
	{
		printf(">> AlignmentSource::InitOne > alignment `%s': type %u\n", label.c_str(), t);

		if (t == atTimeDependent)
		{
			TFile *alF = new TFile(fn.c_str());

			if (alF->IsZombie())
			{
				printf("\tERROR: cannot open file with alignment graphs.\n");
				delete alF;
				return;
			}
			
			TGraph *g_L_2_F = (TGraph *) alF->Get(( string("L_2_F/") + obj).c_str() );
			TGraph *g_L_1_F = (TGraph *) alF->Get(( string("L_1_F/") + obj).c_str() );
			TGraph *g_R_1_F = (TGraph *) alF->Get(( string("R_1_F/") + obj).c_str() );
			TGraph *g_R_2_F = (TGraph *) alF->Get(( string("R_2_F/") + obj).c_str() );

			if (g_L_2_F && g_L_1_F && g_R_1_F && g_R_2_F)
				printf("\talignment graphs successfully loaded\n");
			else {
				printf("\tERROR: unable to load some alignment graphs\n");
				delete alF;
				return;
			}

			gs.L_2_F = new TGraph(*g_L_2_F);
			gs.L_1_F = new TGraph(*g_L_1_F);
			gs.R_1_F = new TGraph(*g_R_1_F);
			gs.R_2_F = new TGraph(*g_R_2_F);

			delete alF;
		}
	}

	void Init()
	{
		printf(">> AlignmentSource::Init\n");
		InitOne("a", type_a, src_a, gs_a, "a_fit");
		InitOne("b", type_b, src_b, gs_b, "b_fit");
		InitOne("c", type_c, src_c, gs_c, "c_fit");
	}

	AlignmentData Eval(double timestamp) const
	{
		AlignmentData d;

		if (type_a == atNone)
		{
			d.a_L_2_F = 0.; d.a_L_1_F = 0.; d.a_R_1_F = 0.; d.a_R_2_F = 0.;
		}

		if (type_a == atConstant)
		{
			d.a_L_2_F = cnst.a_L_2_F; d.a_L_1_F = cnst.a_L_1_F; d.a_R_1_F = cnst.a_R_1_F; d.a_R_2_F = cnst.a_R_2_F;
		}

		if (type_a == atTimeDependent)
		{
			d.a_L_2_F = gs_a.L_2_F->Eval(timestamp)*1E-3; d.a_L_1_F = gs_a.L_1_F->Eval(timestamp)*1E-3; d.a_R_1_F = gs_a.R_1_F->Eval(timestamp)*1E-3; d.a_R_2_F = gs_a.R_2_F->Eval(timestamp)*1E-3;
		}

		if (type_b == atNone)
		{
			d.b_L_2_F = 0.; d.b_L_1_F = 0.; d.b_R_1_F = 0.; d.b_R_2_F = 0.;
		}

		if (type_b == atConstant)
		{
			d.b_L_2_F = cnst.b_L_2_F; d.b_L_1_F = cnst.b_L_1_F; d.b_R_1_F = cnst.b_R_1_F; d.b_R_2_F = cnst.b_R_2_F;
		}

		if (type_b == atTimeDependent)
		{
			d.b_L_2_F = gs_b.L_2_F->Eval(timestamp)*1E-3; d.b_L_1_F = gs_b.L_1_F->Eval(timestamp)*1E-3; d.b_R_1_F = gs_b.R_1_F->Eval(timestamp)*1E-3; d.b_R_2_F = gs_b.R_2_F->Eval(timestamp)*1E-3;
		}
		
		if (type_c == atNone)
		{
			d.c_L_2_F = 0.; d.c_L_1_F = 0.; d.c_R_1_F = 0.; d.c_R_2_F = 0.;
		}

		if (type_c == atConstant)
		{
			d.c_L_2_F = cnst.c_L_2_F; d.c_L_1_F = cnst.c_L_1_F; d.c_R_1_F = cnst.c_R_1_F; d.c_R_2_F = cnst.c_R_2_F;
		}

		if (type_c == atTimeDependent)
		{
			d.c_L_2_F = gs_c.L_2_F->Eval(timestamp)*1E-3; d.c_L_1_F = gs_c.L_1_F->Eval(timestamp)*1E-3; d.c_R_1_F = gs_c.R_1_F->Eval(timestamp)*1E-3; d.c_R_2_F = gs_c.R_2_F->Eval(timestamp)*1E-3;
		}

		return d;
	}
};

//----------------------------------------------------------------------------------------------------

struct UnitHitData
{
	// validity flag
	unsigned int v;

	// hit position in mm
	double x, y;

	UnitHitData() : v(0), x(0.), y(0.) {}

	void operator += (const UnitHitData &add)
	{
		x += add.x;
		y += add.y;
	}
};

//----------------------------------------------------------------------------------------------------

struct HitData
{
	UnitHitData L_1_F, L_2_F /*, L_2_N, L_2_F*/;
	UnitHitData R_1_F, R_2_F /*, R_2_N, R_2_F*/;


	void operator += (const HitData &add)
	{
		L_1_F += add.L_1_F;
		L_2_F += add.L_2_F;

		R_1_F += add.R_1_F;
		R_2_F += add.R_2_F;
	}

	HitData ApplyAlignment(const AlignmentData &al) const
	{
		HitData r;

		r.L_2_F.x = L_2_F.x - al.a_L_2_F * L_2_F.y - al.b_L_2_F; r.L_2_F.y = L_2_F.y - al.c_L_2_F;
		r.L_1_F.x = L_1_F.x - al.a_L_1_F * L_1_F.y - al.b_L_1_F; r.L_1_F.y = L_1_F.y - al.c_L_1_F;
		r.R_1_F.x = R_1_F.x - al.a_R_1_F * R_1_F.y - al.b_R_1_F; r.R_1_F.y = R_1_F.y - al.c_R_1_F;
		r.R_2_F.x = R_2_F.x - al.a_R_2_F * R_2_F.y - al.b_R_2_F; r.R_2_F.y = R_2_F.y - al.c_R_2_F;

		return r;
	}

	// TODO: remove hard-coded z positions
	/*
	HitData ApplyInterpolatedAlignment(const AlignmentData &a, double sN, double sF) const
	{
		AlignmentData a_int = a.Interpolate(214.628, 220.000, sN, sF);

		return ApplyAlignment(a_int);
	}
	*/
};

//----------------------------------------------------------------------------------------------------

struct EventRed
{
	unsigned int timestamp;
	unsigned int run_num, bunch_num, event_num;
	//unsigned int trigger_num;
	//unsigned int trigger_bits;

	// vertical RPs
	HitData h;

	//HitData hH;	// horizontal RPs
};

//----------------------------------------------------------------------------------------------------

struct Environment
{
	// beam momentum (GeV)
	double p, p_L, p_R;

	// beam momentum uncertainty
	double si_de_p;
	
	// beam divergence
	double si_th_x_L, si_th_y_L;		// rad
	double si_th_x_R, si_th_y_R;		// rad

	double si_th_y_RL_assym_unc;		// uncertainty of the L-R assymetry

	// vertex smearing
	double si_vtx_x, si_vtx_y; 		// mm

	// pitch-induced error
	double si_de_P_L, si_de_P_R;	// mm

	// optics
	double v_x_L_1_F, v_x_L_2_F, v_x_R_1_F, v_x_R_2_F;	// 1
	double v_y_L_1_F, v_y_L_2_F, v_y_R_1_F, v_y_R_2_F;	// 1
	double L_x_L_1_F, L_x_L_2_F, L_x_R_1_F, L_x_R_2_F;	// mm
	double L_y_L_1_F, L_y_L_2_F, L_y_R_1_F, L_y_R_2_F;	// mm

	// optics: x-y coupling (x = L_x * th_x + v_x * x^* + la_x * th_y)
	/*
	double la_x_L_F, la_x_L_N, la_x_R_N, la_x_R_F;	// mm
	double la_y_L_F, la_y_L_N, la_y_R_N, la_y_R_F;	// mm
	*/

	// optics perturbation covariance matrices
	// order of elements:
	//		left arm:  v_x_L_N, L_x_L_N, v_y_L_N, L_y_L_N, v_x_L_F, L_x_L_F, v_y_L_F, L_y_L_F
	//		right arm: v_x_R_N, L_x_R_N, v_y_R_N, L_y_R_N, v_x_R_F, L_x_R_F, v_y_R_F, L_y_R_F
	// units: v's in 1, L's in m
	TMatrixDSym opt_cov;

	// optics perturbation generator matrices
	TMatrixD opt_per_gen;

	// alignment uncertainties
	double si_de_x, si_de_y_R, si_de_y_D, si_tilt;

	// misalignments (mm)
	double de_x_L_N, de_y_L_N, tilt_L_N;
	double de_x_L_F, de_y_L_F, tilt_L_F;
	double de_x_R_N, de_y_R_N, tilt_R_N;
	double de_x_R_F, de_y_R_F, tilt_R_F;

	Environment() : opt_cov(16), opt_per_gen(16, 16)
	{
	}

	void InitNominal();
	void UseMatchedOptics();

	void PrintOpticsUncertainties() const;

	void Print() const
	{
		printf("p=%E, p_L=%E, p_R=%E\n", p, p_L, p_R);
		printf("\n");
		printf("si_th_x_L=%E, si_th_y_L=%E\n", si_th_x_L, si_th_y_L);
		printf("si_th_x_R=%E, si_th_y_R=%E\n", si_th_x_R, si_th_y_R);
		printf("si_vtx_x=%E, si_vtx_y=%E\n", si_vtx_x, si_vtx_y);
		printf("si_de_P_L=%E, si_de_P_R=%E\n", si_de_P_L, si_de_P_R);
		printf("\n");

		printf("L_x_L_1_F = %E, v_x_L_1_F = %E, L_y_L_1_F = %E, v_y_L_1_F = %E\n", L_x_L_1_F, v_x_L_1_F, L_y_L_1_F, v_y_L_1_F);
		printf("L_x_L_2_F = %E, v_x_L_2_F = %E, L_y_L_2_F = %E, v_y_L_2_F = %E\n", L_x_L_2_F, v_x_L_2_F, L_y_L_2_F, v_y_L_2_F);
		printf("L_x_R_1_F = %E, v_x_R_1_F = %E, L_y_R_1_F = %E, v_y_R_1_F = %E\n", L_x_R_1_F, v_x_R_1_F, L_y_R_1_F, v_y_R_1_F);
		printf("L_x_R_2_F = %E, v_x_R_2_F = %E, L_y_R_2_F = %E, v_y_R_2_F = %E\n", L_x_R_2_F, v_x_R_2_F, L_y_R_2_F, v_y_R_2_F);

		printf("\n");
		printf("si_de_x=%E, si_de_y_R=%E, si_de_y_D=%E, si_tilt=%E\n", si_de_x, si_de_y_R, si_de_y_D, si_tilt);
		printf("\n");
		printf("de_x_L_N=%E, de_y_L_N=%E, tilt_L_N=%E\n", de_x_L_N, de_y_L_N, tilt_L_N);
		printf("de_x_L_F=%E, de_y_L_F=%E, tilt_L_F=%E\n", de_x_L_F, de_y_L_F, tilt_L_F);
		printf("de_x_R_N=%E, de_y_R_N=%E, tilt_R_N=%E\n", de_x_R_N, de_y_R_N, tilt_R_N);
		printf("de_x_R_F=%E, de_y_R_F=%E, tilt_R_F=%E\n", de_x_R_F, de_y_R_F, tilt_R_F);
		printf("\n");
		printf("si_th_y_RL_assym_unc=%E\n", si_th_y_RL_assym_unc);

		PrintOpticsUncertainties();
	}

	void ApplyRandomOpticsPerturbations(TVectorD &de);
	
	void ApplyRandomOpticsPerturbations()
	{
		TVectorD de(16);
		ApplyRandomOpticsPerturbations(de);
	}

	/// modes counted from 0 to 15
	void ApplyOpticsPerturbationMode(int mode, double coef);

	/// modes counted from 0 to 7
	void ApplyEffectiveLengthPerturbationMode(int mode, double coef);
};

//----------------------------------------------------------------------------------------------------

void Environment::ApplyRandomOpticsPerturbations(TVectorD & /*de*/)
{
	/*
	TVectorD r(16);

	for (unsigned int i = 0; i < 16; i++)
		r(i) = gRandom->Gaus();
	de = opt_per_gen * r;
	
	v_x_L_N += de(0) * 1E0;
	L_x_L_N += de(1) * 1E3;
	v_y_L_N += de(2) * 1E0;
	L_y_L_N += de(3) * 1E3;
	v_x_L_F += de(4) * 1E0;
	L_x_L_F += de(5) * 1E3;
	v_y_L_F += de(6) * 1E0;
	L_y_L_F += de(7) * 1E3;

	v_x_R_N += de(8) * 1E0;
	L_x_R_N += de(9) * 1E3;
	v_y_R_N += de(10) * 1E0;
	L_y_R_N += de(11) * 1E3;
	v_x_R_F += de(12) * 1E0;
	L_x_R_F += de(13) * 1E3;
	v_y_R_F += de(14) * 1E0;
	L_y_R_F += de(15) * 1E3;
	*/
}

//----------------------------------------------------------------------------------------------------

void Environment::ApplyOpticsPerturbationMode(int /*mode*/, double /*coef*/)
{
	/*
	printf(">> Environment::ApplyOpticsPerturbationMode\n");

	// prepare correlation matrix
	TMatrixDSym cor(opt_cov);
	TMatrixDSym Sigma(opt_cov);
	for (int i = 0; i < opt_cov.GetNrows(); i++)
		for (int j = 0; j < opt_cov.GetNcols(); j++)
		{
			cor(i, j) /= sqrt( opt_cov(i, i) * opt_cov(j, j) );
			Sigma(i, j) = (i == j) ? sqrt( opt_cov(i, i) ) : 0.;
		}

	// eigen decomposition
	TMatrixDSymEigen eig_decomp(cor);
	TVectorD eig_values(eig_decomp.GetEigenValues());

	// construct mode
	TVectorD vm(opt_cov.GetNrows());
	for (int i = 0; i < opt_cov.GetNrows(); i++)
	{
		double l = eig_values(i);
		double sl = (l > 0.) ? sqrt(l) : 0.;
		vm(i) = (i == mode) ? sl * coef : 0.;
	}

	vm = Sigma * eig_decomp.GetEigenVectors() * vm;

	printf("\tleft arm: mode %u, coefficient %+.3f\n", mode, coef);
	vm.Print();

	v_x_L_N += vm(0) * 1E0;
	L_x_L_N += vm(1) * 1E3;
	v_y_L_N += vm(2) * 1E0;
	L_y_L_N += vm(3) * 1E3;
	v_x_L_F += vm(4) * 1E0;
	L_x_L_F += vm(5) * 1E3;
	v_y_L_F += vm(6) * 1E0;
	L_y_L_F += vm(7) * 1E3;

	v_x_R_N += vm(8) * 1E0;
	L_x_R_N += vm(9) * 1E3;
	v_y_R_N += vm(10) * 1E0;
	L_y_R_N += vm(11) * 1E3;
	v_x_R_F += vm(12) * 1E0;
	L_x_R_F += vm(13) * 1E3;
	v_y_R_F += vm(14) * 1E0;
	L_y_R_F += vm(15) * 1E3;
	*/
}

//----------------------------------------------------------------------------------------------------

void Environment::ApplyEffectiveLengthPerturbationMode(int /*mode*/, double /*coef*/)
{
	/*

	printf(">> Environment::ApplyEffectiveLengthPerturbationMode\n");

	// prepare reduced covariance matrix
	TMatrixDSym cov_red(8);
	for (unsigned int i = 0; i < 8; i++)
		for (unsigned int j = 0; j < 8; j++)
			cov_red(i, j) = opt_cov(2*i+1, 2*j+1);

	// eigen decomposition
	TMatrixDSymEigen eig_decomp(cov_red);
	TVectorD eig_values(eig_decomp.GetEigenValues());

	// construct mode
	TVectorD vm(cov_red.GetNrows());
	for (int i = 0; i < cov_red.GetNrows(); i++)
	{
		double l = eig_values(i);
		double sl = (l > 0.) ? sqrt(l) : 0.;
		vm(i) = (i == mode) ? sl * coef : 0.;
	}

	vm = eig_decomp.GetEigenVectors() * vm;

	printf("\tmode %u, coefficient %+.3f\n", mode, coef);
	//vm.Print();

	L_x_L_N += vm(0) * 1E3;
	L_y_L_N += vm(1) * 1E3;
	L_x_L_F += vm(2) * 1E3;
	L_y_L_F += vm(3) * 1E3;
	L_x_R_N += vm(4) * 1E3;
	L_y_R_N += vm(5) * 1E3;
	L_x_R_F += vm(6) * 1E3;
	L_y_R_F += vm(7) * 1E3;

	*/
}

//----------------------------------------------------------------------------------------------------

void Environment::PrintOpticsUncertainties() const
{
	printf("optics uncertainties: left arm\n");
	printf("\tv_x_N: %.4f\n", sqrt(opt_cov(0, 0)));
	printf("\tL_x_N: %.4f m\n", sqrt(opt_cov(1, 1)));
	printf("\tv_y_N: %.4f\n", sqrt(opt_cov(2, 2)));
	printf("\tL_y_N: %.4f m\n", sqrt(opt_cov(3, 3)));
	printf("\tv_x_F: %.4f\n", sqrt(opt_cov(4, 4)));
	printf("\tL_x_F: %.4f m\n", sqrt(opt_cov(5, 5)));
	printf("\tv_y_F: %.4f\n", sqrt(opt_cov(6, 6)));
	printf("\tL_y_F: %.4f m\n", sqrt(opt_cov(7, 7)));

	printf("optics uncertainties: right arm\n");
	printf("\tv_x_N: %.4f\n", sqrt(opt_cov(8, 8)));
	printf("\tL_x_N: %.4f m\n", sqrt(opt_cov(9, 9)));
	printf("\tv_y_N: %.4f\n", sqrt(opt_cov(10, 10)));
	printf("\tL_y_N: %.4f m\n", sqrt(opt_cov(11, 11)));
	printf("\tv_x_F: %.4f\n", sqrt(opt_cov(12, 12)));
	printf("\tL_x_F: %.4f m\n", sqrt(opt_cov(13, 13)));
	printf("\tv_y_F: %.4f\n", sqrt(opt_cov(14, 14)));
	printf("\tL_y_F: %.4f m\n", sqrt(opt_cov(15, 15)));
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct Kinematics
{
	double th_x_L_1_F, th_x_L_2_F, th_x_R_1_F, th_x_R_2_F, th_x_L, th_x_R, th_x;	//	rad
	double th_y_L_1_F, th_y_L_2_F, th_y_R_1_F, th_y_R_2_F, th_y_L, th_y_R, th_y;	//	rad

	double vtx_x_L_1_F, vtx_x_L_2_F, vtx_x_R_1_F, vtx_x_R_2_F, vtx_x_L, vtx_x_R, vtx_x;	// in mm
	double vtx_y_L_1_F, vtx_y_L_2_F, vtx_y_R_1_F, vtx_y_R_2_F, vtx_y_L, vtx_y_R, vtx_y;	// in mm

	double th;				// in rad
	double phi;				// in rad
	double t_x, t_y, t;		// in GeV^2

	Kinematics() : th_y(0.) {}

	void ThetasToTPhi(const Environment &env)
	{
		th = sqrt(th_x*th_x + th_y*th_y);
		t_x = env.p*env.p * th_x * th_x;
		t_y = env.p*env.p * th_y * th_y;
		t = t_x + t_y;
		phi = atan2(th_y, th_x);
	}

	void TPhiToThetas(const Environment &env)
	{
		th = sqrt(t) / env.p;
		th_x_L = th_x_R = th_x = th * cos(phi);
		th_y_L = th_y_R = th_y = th * sin(phi);

		t_x = t * cos(phi) * cos(phi);
		t_y = t * sin(phi) * sin(phi);
	}
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct CutData
{
	double cqa[9];	///< array of quantities qa
	double cqb[9];	///< array of quantities qb
	double cv[9];	///< array of cut quantities v = a*qa + b*qb + c
	bool ct[9];		///< array of flags whether |v| < n_si * si
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct Analysis
{
	// binning, |t| in GeV^2
	double t_min, t_max;
	double t_min_full, t_max_full;
	double t_min_fit;

	// elastic selection cuts
	double n_si;

	double cut1_a, cut1_c, cut1_si;
	double cut2_a, cut2_c, cut2_si;
	double cut3_a, cut3_c, cut3_si;
	double cut4_a, cut4_c, cut4_si;
	double cut5_a, cut5_c, cut5_si;
	double cut6_a, cut6_c, cut6_si;
	double cut7_a, cut7_c, cut7_si;
	double cut8_a, cut8_c, cut8_si;

	std::vector< std::pair<double, double> > timeIntervals;

	unsigned int N_cuts;	// number of cuts - indexed from 1!
	string cqaN[9], cqbN[9];
	double cca[9], ccb[9], ccc[9], csi[9];
	std::vector<unsigned int> cuts;	// list of active cuts

	// analysis cuts (rad)
	double th_y_lcut_L, th_y_lcut_R, th_y_lcut;
	double th_y_hcut_L, th_y_hcut_R, th_y_hcut;
	
	double th_x_lcut;
	double th_x_hcut;

	// (un)-smearing parameters
	double si_th_x_1arm_L;
	double si_th_x_1arm_R;
	double si_th_x_1arm_unc;
	double si_th_x_2arm;
	double si_th_x_2arm_unc;

	double si_th_y_1arm;
	double si_th_y_1arm_unc;
	double si_th_y_2arm;
	double si_th_y_2arm_unc;

	// efficiency parameters
	bool use_3outof4_efficiency_fits;		// whether to use time-dependent fits of 3-out-of-4 efficiency
	bool use_pileup_efficiency_fits;		// whether to use time-dependent fits of pile-up efficiency

	double inefficiency_3outof4;			// inefficiency from 3-out-of-4 method, used only if use_3outof4_efficiency_fits=false
	double inefficiency_shower_near;		// inefficiency due to shower in near RP
	double inefficiency_pile_up;			// inefficiency due to pile-up, used only if use_pileup_efficiency_fits=false
	double inefficiency_trigger;			// trigger inefficiency
	double inefficiency_DAQ;				// DAQ inefficiency

	// normalisation correction to subtract background
	double bckg_corr;

	// (delivered) luminosity
	double L_int;	// mb^-1

	// 3-out-of-4 efficiency uncertainty (only used in MC simulation)
	double eff_3outof4_fixed_point, eff_3outof4_slope, eff_3outof4_slope_unc;

	// normalisation correction and its uncertainty (only used in MC simulation)
	double norm_corr, norm_corr_unc;

	double alignment_t0;	// beginning of the first time-slice
	double alignment_ts;	// time-slice in s
	
	double eff_th_y_min;

	// y ranges for alignment
	struct AlignmentYRange
	{
		double bot_min, bot_max, top_min, top_max;
		AlignmentYRange(double bmi=0., double bma=0., double tmi=0., double tma=0.) :
			bot_min(bmi), bot_max(bma), top_min(tmi), top_max(tma) {}
	};
	map<std::string, AlignmentYRange> alignmentYRanges;

	void BuildCuts();
	bool EvaluateCuts(const HitData &, const Kinematics &, CutData &) const;

	bool SkipTime(unsigned int timestamp) const
	{
		if (timeIntervals.size() == 0)
			return false;

		bool selected = false;
		for (unsigned int i = 0; i < timeIntervals.size(); i++)
		{
			if (timestamp >= timeIntervals[i].first && timestamp <= timeIntervals[i].second)
			{
				selected = true;
				break;
			}
		}

		return !selected;
	}


	void Print() const
	{
		printf("t_min=%E, t_max=%E, t_min_full=%E, t_max_full=%E\n", t_min, t_max, t_min_full, t_max_full);
		printf("t_min_fit=%E\n", t_min_fit);

		printf("\n");
		printf("%lu time intervals:\n", timeIntervals.size());
		for (std::vector< std::pair<double, double> >::const_iterator it = timeIntervals.begin(); it != timeIntervals.end(); ++it)
			printf("\tfrom %.1f to %.1f\n", it->first, it->second);

		printf("\n");
		printf("n_si=%E\n", n_si);

		printf("\n");
		printf("cut1_a=%E, cut1_c=%E, cut1_si=%E\n", cut1_a, cut1_c, cut1_si);
		printf("cut2_a=%E, cut2_c=%E, cut2_si=%E\n", cut2_a, cut2_c, cut2_si);
		printf("cut3_a=%E, cut3_c=%E, cut3_si=%E\n", cut3_a, cut3_c, cut3_si);
		printf("cut4_a=%E, cut4_c=%E, cut4_si=%E\n", cut4_a, cut4_c, cut4_si);
		printf("cut5_a=%E, cut5_c=%E, cut5_si=%E\n", cut5_a, cut5_c, cut5_si);
		printf("cut6_a=%E, cut6_c=%E, cut6_si=%E\n", cut6_a, cut6_c, cut6_si);
		printf("cut7_a=%E, cut7_c=%E, cut7_si=%E\n", cut7_a, cut7_c, cut7_si);
		printf("cut8_a=%E, cut8_c=%E, cut8_si=%E\n", cut8_a, cut8_c, cut8_si);

		printf("\n");
		printf("cut parameters:\n");
		for (unsigned int i = 1; i <= N_cuts; i++)
		{
			printf("%u| cqaN=%s, cqbN=%s | cca=%E, ccb=%E, ccc=%E, csi=%E\n", i,
				cqaN[i].c_str(), cqbN[i].c_str(), cca[i], ccb[i], ccc[i], csi[i]);
		}

		printf("\n");
		printf("%lu enabled cuts: ", cuts.size());
		for (unsigned int i = 0; i < cuts.size(); i++)
			printf((i == 0) ? "%i" : ", %i", cuts[i]);

		printf("\n");
		printf("th_x_lcut=%E\n", th_x_lcut);
		printf("th_x_hcut=%E\n", th_x_hcut);
		printf("th_y_lcut_L=%E, th_y_lcut_R=%E, th_y_lcut=%E\n", th_y_lcut_L, th_y_lcut_R, th_y_lcut);
		printf("th_y_hcut_L=%E, th_y_hcut_R=%E, th_y_hcut=%E\n", th_y_hcut_L, th_y_hcut_R, th_y_hcut);

		printf("\n");
		printf("si_th_x_1arm_L=%E, si_th_x_1arm_R=%E, si_th_x_1arm_unc=%E\n", si_th_x_1arm_L, si_th_x_1arm_R, si_th_x_1arm_unc);
		printf("si_th_x_2arm=%E, si_th_x_2arm_unc=%E\n", si_th_x_2arm, si_th_x_2arm_unc);
		printf("si_th_y_1arm=%E, si_th_y_1arm_unc=%E\n", si_th_y_1arm, si_th_y_1arm_unc);
		printf("si_th_y_2arm=%E, si_th_y_2arm_unc=%E\n", si_th_y_2arm, si_th_y_2arm_unc);
	
		printf("\n");
		printf("use_3outof4_efficiency_fits = %i\n", use_3outof4_efficiency_fits);
		printf("use_pileup_efficiency_fits= %i\n", use_pileup_efficiency_fits);
		printf("inefficiency_3outof4 = %.3f\n", inefficiency_3outof4);
		printf("inefficiency_shower_near = %.3f\n", inefficiency_shower_near);
		printf("inefficiency_pile_up = %.3f\n", inefficiency_pile_up);
		printf("inefficiency_trigger = %.3f\n", inefficiency_trigger);
		printf("inefficiency_DAQ = %.3f\n", inefficiency_DAQ);
		printf("bckg_corr = %.3f\n", bckg_corr);
		printf("L_int=%E\n", L_int);
		printf("eff_3outof4_fixed_point=%E, eff_3outof4_slope=%E, eff_3outof4_slope_unc=%E\n", eff_3outof4_fixed_point, eff_3outof4_slope, eff_3outof4_slope_unc);
		printf("norm_corr=%E, norm_corr_unc=%E\n", norm_corr, norm_corr_unc);
	}
};

//----------------------------------------------------------------------------------------------------

#include "common_cuts.h"
#include "common_parameters.h"

//----------------------------------------------------------------------------------------------------

#endif
