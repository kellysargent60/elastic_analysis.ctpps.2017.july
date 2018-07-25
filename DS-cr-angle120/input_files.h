#include <string>
#include <vector>

//----------------------------------------------------------------------------------------------------

std::vector<std::string> input_files;

void InitInputFiles()
{
	std::string prefix;

	prefix="root://eostotem.cern.ch//eos/totem/data/ctpps/reconstruction/2017/alignment_run_July/version6/";
	input_files.push_back(prefix + "run_298597.0_re_reco.root");
	input_files.push_back(prefix + "run_298597.10_re_reco.root");
	input_files.push_back(prefix + "run_298597.11_re_reco.root");
	input_files.push_back(prefix + "run_298597.12_re_reco.root");
	input_files.push_back(prefix + "run_298597.13_re_reco.root");
	input_files.push_back(prefix + "run_298597.14_re_reco.root");
	input_files.push_back(prefix + "run_298597.15_re_reco.root");
	input_files.push_back(prefix + "run_298597.16_re_reco.root");
	input_files.push_back(prefix + "run_298597.17_re_reco.root");
	input_files.push_back(prefix + "run_298597.18_re_reco.root");
	input_files.push_back(prefix + "run_298597.19_re_reco.root");
	input_files.push_back(prefix + "run_298597.1_re_reco.root");
	input_files.push_back(prefix + "run_298597.20_re_reco.root");
	input_files.push_back(prefix + "run_298597.21_re_reco.root");
	input_files.push_back(prefix + "run_298597.22_re_reco.root");
	input_files.push_back(prefix + "run_298597.23_re_reco.root");
	input_files.push_back(prefix + "run_298597.2_re_reco.root");
	input_files.push_back(prefix + "run_298597.3_re_reco.root");
	input_files.push_back(prefix + "run_298597.4_re_reco.root");
	input_files.push_back(prefix + "run_298597.5_re_reco.root");
	input_files.push_back(prefix + "run_298597.6_re_reco.root");
	input_files.push_back(prefix + "run_298597.7_re_reco.root");
	input_files.push_back(prefix + "run_298597.8_re_reco.root");
	input_files.push_back(prefix + "run_298597.9_re_reco.root");
}
