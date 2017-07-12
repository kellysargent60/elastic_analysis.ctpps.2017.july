#include <vector>
#include <string>
#include <cmath>

#include "TVectorD.h"

using namespace std;


struct Stat
{
	unsigned int dim;
	double S1;
	vector<double> Sv, Svv, Svvv, Svvvv;

	vector< vector<double> > Sxy, Sxxy, Sxyy, Sxxyy;

	vector<string> labels;

	void Init(unsigned int _dim = 1)
	{
		dim = _dim;

		S1 = 0.;
		for (unsigned int i = 0; i < dim; i++)
		{
			Sv.push_back(0);
			Svv.push_back(0);
			Svvv.push_back(0);
			Svvvv.push_back(0);

			vector<double> temp;
			for (unsigned int j = 0; j < dim; j++)
			{
				temp.push_back(0);
			}
			Sxy.push_back(temp);
			Sxxy.push_back(temp);
			Sxyy.push_back(temp);
			Sxxyy.push_back(temp);
		}
	}

	Stat() {}

	Stat(unsigned int _dim)
	{
		Init(_dim);
	}

	void SetLabels(const vector<string> &_l)
	{
		labels.resize(dim);
		for (unsigned int i = 0; i < dim; i++)
			labels[i] = _l[i];
	}

	template <class T>
	void Fill(const T &v)
	{
		S1 += 1.;
		for (unsigned int i = 0; i < dim; i++)
		{
			Sv[i] += v[i];
			Svv[i] += v[i]*v[i];
			Svvv[i] += v[i]*v[i]*v[i];
			Svvvv[i] += v[i]*v[i]*v[i]*v[i];

			for (unsigned int j = 0; j < dim; j++)
			{
				Sxy[i][j] += v[i] * v[j];
				Sxxy[i][j] += v[i]*v[i] * v[j];
				Sxyy[i][j] += v[i] * v[j]*v[j];
				Sxxyy[i][j] += v[i]*v[i] * v[j]*v[j];
			}
		}
	}

	void Fill(double v1, double v2 = 0., double v3 = 0., double v4 = 0., double v5 = 0.)
	{
		vector<double> v(5);
		v[0] = v1;
		v[1] = v2;
		v[2] = v3;
		v[3] = v4;
		v[4] = v5;

		Fill(v);
	}

	string QLabel(unsigned int i) const
	{
		if (labels.empty())
		{
			char buf[10];
			sprintf(buf, "qu.%3i", i);
			return buf;
		} else
			return labels[i];
	}
	
	//--------------------
	// 1D getters
	//--------------------

	double GetEntries()
	{
		return S1;
	}

	double GetMean(unsigned int i) const
	{
		double mu = Sv[i] / S1;
		return mu;
	}

	double GetStdDev(unsigned int i) const
	{
		double v = (Svv[i] - Sv[i]*Sv[i] / S1) / (S1 - 1.);
		double s = (v > 0.) ? sqrt(v) : 0.;
		return s;
	}

	double GetMeanUnc(unsigned int i) const
	{
		double mu_unc = GetStdDev(i) / sqrt(S1);
		return mu_unc;
	}

	double GetStdDevUnc(unsigned int i) const
	{
		double mu = GetMean(i);
		double s = GetStdDev(i);
		double v = s*s;

		double sum = Svvvv[i] - 4.*mu*Svvv[i] + 6.*mu*mu*Svv[i] - 4.*mu*mu*mu*Sv[i] + mu*mu*mu*mu*S1;
		double E4 = (S1 > 1.) ? sum / (S1 - 1.) : 0.;

		double v_var = (S1 > 3.) ? (E4 - (S1 - 3.)/(S1 - 1.)*v*v) / S1 : 0.;
		double s_var = v_var / 4. / v;
		double s_s = (s_var > 0.) ? sqrt(s_var) : 0.;
		return s_s;
	}

	// approximation of GetStdDevUnc valid for Gaussian distributions
	double GetStdDevUncGauss(unsigned int i) const
	{
		double s = GetStdDev(i);
		double s_s = (S1 > 0.) ? s / sqrt(2. * S1) : 0.;
		return s_s;
	}
	
	//--------------------
	// 2D getters
	//--------------------

	double GetCovariance(unsigned int i, unsigned int j) const
	{
		double C = (S1 > 1.) ? (Sxy[i][j] - Sv[i]*Sv[j] / S1) / (S1 - 1.) : 0.;
		return C;
	}

	double GetCorrelation(unsigned int i, unsigned int j) const
	{
		double C = GetCovariance(i, j);
		double den = GetStdDev(i) * GetStdDev(j);
		double rho = (den > 0.) ? C / den : 0.;
		return rho;
	}

	double GetCovarianceUnc(unsigned int i, unsigned int j) const
	{
		double mx = GetMean(i);
		double my = GetMean(j);
		double sx = GetStdDev(i);
		double sy = GetStdDev(j);
		double C = GetCovariance(i, j);

		double sum =
			Sxxyy[i][j] 
			-2.*Sxyy[i][j]*mx - 2.*Sxxy[i][j]*my
			+ 4.*Sxy[i][j]*mx*my
			+ Svv[i]*my*my + Svv[j]*mx*mx
			- 2.*Sv[i]*mx*my*my - 2.*Sv[j]*mx*mx*my
			+ mx*mx*my*my;
		double D = (S1 > 1.) ? sum / (S1 - 1.) : 0.;

		double C_var = (S1 > 2.) ? (D + sx*sx*sy*sy/(S1 - 1.) - (S1-2.)/(S1-1.)*C*C) / S1 : 0.;
		double C_s = (C_var > 0.) ? sqrt(C_var) : 0.;

		return C_s;
	}

	double GetCorrelationUnc(unsigned int i, unsigned int j) const
	{
		// WARNING: the calculation below assumes no correlation between C, si_i and si_j, which
		// might not be correct - in that case it gives an upper bound for the uncertainty

		double C = GetCovariance(i, j), C_unc = GetCovarianceUnc(i, j);
		double si_i = GetStdDev(i), si_i_unc = GetStdDevUnc(i);
		double si_j = GetStdDev(j), si_j_unc = GetStdDevUnc(j);
		double rho = C / (si_i * si_j);
		double sum =
			(C != 0. && si_i != 0. && si_j != 0.) ? pow(C_unc / C, 2.) + pow(si_i_unc / si_i, 2.) + pow(si_j_unc / si_j, 2.) : 0.;
		double rho_unc = fabs(rho) * sqrt(sum);
		return rho_unc;
	}

	TMatrixDSym GetCovarianceMatrix() const
	{
		TMatrixDSym m(dim);

		for (unsigned int i = 0; i < dim; i++)
			for (unsigned int j = 0; j < dim; j++)
				m(i, j) = GetCovariance(i, j);

		return m;
	}

	//--------------------
	// print methods
	//--------------------

	void PrintStat() const
	{
		printf("entries: %.3E\n", S1);
	}

	void PrintMeanAndStdDev() const
	{
		for (unsigned int i = 0; i < dim; i++)
		{
			double mu = GetMean(i);
			double mu_unc = GetMeanUnc(i);
			double s = GetStdDev(i);
			double s_unc = GetStdDevUnc(i);
			printf("%s: mean %+.3E +- %.3E, std. dev. = %.3E +- %.3E\n", QLabel(i).c_str(), mu, mu_unc, s, s_unc);
		}
	}
	
	void PrintCovariance() const
	{
		printf("      ");
		for (unsigned int i = 0; i < dim; i++)
			printf("   %6s", QLabel(i).c_str());
		printf("\n");

		for (unsigned int i = 0; i < dim; i++)
		{
			printf("%6s", QLabel(i).c_str());
			for (unsigned int j = 0; j < dim; j++)
				printf("   %+.3f", GetCovariance(i, j));
			printf("\n");
		}
	}

	void PrintCorrelation()
	{
		printf("      ");
		for (unsigned int i = 0; i < dim; i++)
			printf("   %6s", QLabel(i).c_str());
		printf("\n");

		for (unsigned int i = 0; i < dim; i++)
		{
			printf("%6s", QLabel(i).c_str());
			for (unsigned int j = 0; j < dim; j++)
				printf("   %+.3f", GetCorrelation(i, j));
			printf("\n");
		}
	}
};


