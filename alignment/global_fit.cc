#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"

#include <vector>

using namespace std;

//----------------------------------------------------------------------------------------------------

void FitAndWriteOne(TGraphErrors *g, const string label)
{
	// make fit
	TF1 *ff = new TF1("ff", "[0]");
	g->Fit(ff);

	g->Write(label.c_str());

	// determine graph range
	double x_min = 1E100, x_max = -1E100;
	for (signed int i = 0; i < g->GetN(); i++)
	{
		double x, y;
		g->GetPoint(i, x, y);
		double x_unc = g->GetErrorX(i);

		x_min = min(x - x_unc, x_min);
		x_max = max(x + x_unc, x_max);
	}

	// make fit graph
	double x[2] = { x_min, x_max };
	double y[2] = { ff->Eval(x_min), ff->Eval(x_max) };
	TGraph *gf = new TGraph(2, x, y);

	gf->Write((label + "_fit").c_str());

	delete ff;
}

//----------------------------------------------------------------------------------------------------

struct UnitData
{
	TGraphErrors *g_a, *g_b;

	UnitData()
	{
		g_a = new TGraphErrors();
		g_b = new TGraphErrors();
	}

	void FitAndWrite()
	{
		FitAndWriteOne(g_a, "a");
		FitAndWriteOne(g_b, "b");
	}
};

//----------------------------------------------------------------------------------------------------

void Append(TGraphErrors *g_dest, const TGraphErrors *g_src)
{
	for (signed int i = 0; i < g_src->GetN(); i++)
	{
		double x, y;
		g_src->GetPoint(i, x, y);
		double x_unc = g_src->GetErrorX(i);
		double y_unc = g_src->GetErrorY(i);

		int idx = g_dest->GetN();
		g_dest->SetPoint(idx, x, y);
		g_dest->SetPointError(idx, x_unc, y_unc);
	}
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// input
	vector<string> datasets = {
		"DS-cr-angle120",
		"DS-cr-angle130",
		"DS-cr-angle140",
	};

	vector<string> units = {
		"L_2_F",
		"L_1_F",
		"R_1_F",
		"R_2_F",
	};

	// prepare data structures
	vector<UnitData> unitData(units.size());

	// process input
	for (const auto &dataset : datasets)
	{
		string f_in_name = "../" + dataset + "/alignment_fit.root";
		TFile *f_in = TFile::Open(f_in_name.c_str());

		for (unsigned int ui = 0; ui < units.size(); ++ui)
		{
			printf("\n*%s, %s\n", dataset.c_str(), units[ui].c_str());


			TGraphErrors *g_a = (TGraphErrors *) f_in->Get((units[ui] + "/a").c_str());
			TGraphErrors *g_b = (TGraphErrors *) f_in->Get((units[ui] + "/b").c_str());

			Append(unitData[ui].g_a, g_a);
			Append(unitData[ui].g_b, g_b);
		}

		delete f_in;
	}

	// save output
	TFile *f_out = TFile::Open("global_fit.root", "recreate");

	for (unsigned int ui = 0; ui < units.size(); ++ui)
	{
		gDirectory = f_out->mkdir(units[ui].c_str());

		unitData[ui].FitAndWrite();
	}

	delete f_out;

	// clean up
	return 0;
}
