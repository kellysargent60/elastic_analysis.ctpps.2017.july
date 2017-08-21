#include "TFile.h"
#include "TH1D.h"

#include <vector>

using namespace std;

//----------------------------------------------------------------------------------------------------

int main()
{
	vector<string> datasets = {
		"DS-cr-angle120",
		"DS-cr-angle130",
		"DS-cr-angle140",
	};

	vector<string> diagonals = {
		"45b_56t",
		"45t_56b",
	};

	double entries_tot = 0.;

	for (const auto &dataset : datasets)
	{
		for (const auto &diagonal : diagonals)
		{
			printf("\n*%s, %s\n", dataset.c_str(), diagonal.c_str());

			string f_in_name = "../" + dataset + "/distributions_" + diagonal + ".root";
			TFile *f_in = TFile::Open(f_in_name.c_str());

			string h_path = "elastic cuts/cut 1/h_cq1";
			TH1D *h_in = (TH1D *) f_in->Get(h_path.c_str());

			double entries = h_in->GetEntries();

			printf("    %.0f\n", entries);

			entries_tot += entries;
		}
	}

	printf("\n*total\n");
	printf("    %.2E\n", entries_tot);

	return 0;
}
