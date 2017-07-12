#include <string>
#include <vector>

//----------------------------------------------------------------------------------------------------

std::vector<std::string> input_files;

void InitInputFiles()
{
	std::string prefix = "root://eostotem.cern.ch//eos/totem/data/ctpps/reconstruction/2017/alignment_run_July/version2/";
	input_files = {
		prefix + "run_298593.0_re_reco.root",
		prefix + "run_298593.11_re_reco.root",
		prefix + "run_298593.12_re_reco.root",
		prefix + "run_298593.2_re_reco.root",
		prefix + "run_298593.3_re_reco.root",
		prefix + "run_298593.5_re_reco.root",
		prefix + "run_298593.6_re_reco.root",
		prefix + "run_298593.7_re_reco.root",
		prefix + "run_298593.9_re_reco.root",
	};
}
