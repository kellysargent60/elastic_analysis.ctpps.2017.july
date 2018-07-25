#include <string>
#include <vector>

//----------------------------------------------------------------------------------------------------

std::vector<std::string> input_files;

void InitInputFiles()
{
	std::string prefix;

	prefix="root://eostotem.cern.ch//eos/totem/data/ctpps/reconstruction/2017/alignment_run_July/version6/";
	input_files.push_back(prefix + "run_298593.0_re_reco.root");
	input_files.push_back(prefix + "run_298593.10_re_reco.root");
	input_files.push_back(prefix + "run_298593.11_re_reco.root");
	input_files.push_back(prefix + "run_298593.12_re_reco.root");
	input_files.push_back(prefix + "run_298593.13_re_reco.root");
	input_files.push_back(prefix + "run_298593.14_re_reco.root");
	input_files.push_back(prefix + "run_298593.15_re_reco.root");
	input_files.push_back(prefix + "run_298593.16_re_reco.root");
	input_files.push_back(prefix + "run_298593.17_re_reco.root");
	input_files.push_back(prefix + "run_298593.18_re_reco.root");
	input_files.push_back(prefix + "run_298593.19_re_reco.root");
	input_files.push_back(prefix + "run_298593.1_re_reco.root");
	input_files.push_back(prefix + "run_298593.2_re_reco.root");
	input_files.push_back(prefix + "run_298593.3_re_reco.root");
	input_files.push_back(prefix + "run_298593.4_re_reco.root");
	input_files.push_back(prefix + "run_298593.5_re_reco.root");
	input_files.push_back(prefix + "run_298593.6_re_reco.root");
	input_files.push_back(prefix + "run_298593.7_re_reco.root");
	input_files.push_back(prefix + "run_298593.8_re_reco.root");
	input_files.push_back(prefix + "run_298593.9_re_reco.root");
}
