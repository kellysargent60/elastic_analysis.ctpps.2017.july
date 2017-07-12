#include "common_definitions.h"
#include "common_algorithms.h"
#include "parameters.h"
#include "common.h"

#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"

#include <vector>

using namespace std;

//----------------------------------------------------------------------------------------------------

struct GraphSet
{
	TGraph *a, *b, *c;

	GraphSet() : a(NULL), b(NULL), c(NULL)
	{
	}
};

//----------------------------------------------------------------------------------------------------

TGraph* WriteFitGraph(TGraph *gi, TF1 *f)
{
	double x_min, ex_min, x_max, ex_max, dummy;
	gi->GetPoint(0, x_min, dummy);
	ex_min = gi->GetErrorX(0);
	x_min -= 1.*ex_min;

	gi->GetPoint(gi->GetN()-1, x_max, dummy);
	ex_max = gi->GetErrorX(gi->GetN()-1);
	x_max += 1.*ex_max;

	double x[2] = { x_min, x_max };
	double y[2] = { f->Eval(x_min), f->Eval(x_max) };
	TGraph *go = new TGraph(2, x, y);
	char buf[100];
	sprintf(buf, "%s_fit", gi->GetName());
	go->SetName(buf);
	go->SetLineColor(4);
	go->Write();

	return go;
}

//----------------------------------------------------------------------------------------------------

TGraphErrors *g_idx = NULL;

//----------------------------------------------------------------------------------------------------

void FitOneUnit(TDirectory *inDir, GraphSet &gs)
{
	TGraphErrors *g_a = (TGraphErrors *) inDir->Get("a");
	TGraphErrors *g_b = (TGraphErrors *) inDir->Get("b");
	TGraphErrors *g_c = (TGraphErrors *) inDir->Get("c");

	if (g_idx == NULL)
		g_idx = g_a;

	TF1 *pol0 = new TF1("pol0", "[0]");
	//TF1 *pol1 = new TF1("pol1", "[0] + [1]*x");

	TF1 *ff;

	if (g_a)
	{
		g_a->SetTitle(";time   (s);a   (mrad)");
		printf(">> a\n");
		ff = pol0;
		g_a->Fit(ff);
		g_a->Write();
		gs.a = WriteFitGraph(g_a, ff);
	}
	
	if (g_b)
	{
		g_b->SetTitle(";time   (s);b   (#murad)");
		printf(">> b\n");
		ff = pol0;
		g_b->Fit(ff);
		g_b->Write();
		WriteFitGraph(g_b, ff);
		gs.b = g_b; 
	}
	
	if (g_c)
	{
		g_c->SetTitle(";time   (s);c   (#murad)");
		printf(">> c\n");
		ff = pol0;
		g_c->Fit(ff);
		g_c->Write();
		gs.c = WriteFitGraph(g_c, ff);
	}
}

//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	// init diagonal
	Init(argv[1]);
	if (diagonal != dCombined)
		return rcIncompatibleDiagonal;

	// files
	TFile *inF = new TFile("alignment.root");
	TFile *outF = new TFile("alignment_fit.root", "recreate");

	char buf[1000];

	vector<string> units;
	units.push_back("L_1_N");
	units.push_back("L_1_F");
	units.push_back("R_1_N");
	units.push_back("R_1_F");

	vector<GraphSet> graphs;

	// make fits
	for (unsigned int ui = 0; ui < units.size(); ui++)
	{
		sprintf(buf, "global/%s", units[ui].c_str());
		TDirectory *inDir = (TDirectory *) inF->Get(buf);

		gDirectory = outF->mkdir(units[ui].c_str());

		printf("\n-------------------- %s --------------------\n", units[ui].c_str());
		GraphSet gs;
		FitOneUnit(inDir, gs);
		graphs.push_back(gs);
	}

	// print out samples
	FILE *fo = fopen("alignment_fit.out", "w");
	for (int pi = 0; pi < g_idx->GetN(); ++pi)
	{
		double center, bla;
		g_idx->GetPoint(pi, center, bla);
		double w = g_idx->GetErrorX(pi);
		unsigned long from = floor(timestamp0 + center - w);
		unsigned long to = ceil(timestamp0 + center + w); to--;
		fprintf(fo, "%lu,%lu", from, to);
		for (unsigned int ui = 0; ui < units.size(); ui++)
		{
			double a_v = (graphs[ui].a) ? graphs[ui].a->Eval(center) : 0.;
			double b_v = (graphs[ui].b) ? graphs[ui].b->Eval(center) : 0.;
			double c_v = (graphs[ui].c) ? graphs[ui].c->Eval(center) : 0.;
			fprintf(fo, ",%E,%E,%E", a_v, b_v, c_v);
		}
		fprintf(fo, "\n");
	}
	fclose(fo);

	delete outF;
	return 0;
}
