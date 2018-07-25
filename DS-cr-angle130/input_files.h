#include <string>
#include <vector>

//----------------------------------------------------------------------------------------------------

std::vector<std::string> input_files;

void InitInputFiles()
{
	std::string prefix;

	prefix="root://eostotem.cern.ch//eos/totem/data/ctpps/reconstruction/2017/alignment_run_July/version6/";
	input_files.push_back(prefix + "run_298594.0_re_reco.root");
	input_files.push_back(prefix + "run_298594.10_re_reco.root");
	input_files.push_back(prefix + "run_298594.11_re_reco.root");
	input_files.push_back(prefix + "run_298594.12_re_reco.root");
	input_files.push_back(prefix + "run_298594.13_re_reco.root");
	input_files.push_back(prefix + "run_298594.14_re_reco.root");
	input_files.push_back(prefix + "run_298594.15_re_reco.root");
	input_files.push_back(prefix + "run_298594.16_re_reco.root");
	input_files.push_back(prefix + "run_298594.17_re_reco.root");
	input_files.push_back(prefix + "run_298594.18_re_reco.root");
	input_files.push_back(prefix + "run_298594.19_re_reco.root");
	input_files.push_back(prefix + "run_298594.1_re_reco.root");
	input_files.push_back(prefix + "run_298594.2_re_reco.root");
	input_files.push_back(prefix + "run_298594.3_re_reco.root");
	input_files.push_back(prefix + "run_298594.4_re_reco.root");
	input_files.push_back(prefix + "run_298594.5_re_reco.root");
	input_files.push_back(prefix + "run_298594.6_re_reco.root");
	input_files.push_back(prefix + "run_298594.7_re_reco.root");
	input_files.push_back(prefix + "run_298594.8_re_reco.root");
	input_files.push_back(prefix + "run_298594.9_re_reco.root");

	prefix="root://eostotem.cern.ch//eos/totem/data/ctpps/reconstruction/2017/alignment_run_July/version6/";
	input_files.push_back(prefix + "run_298596.0_re_reco.root");
	input_files.push_back(prefix + "run_298596.1_re_reco.root");
	input_files.push_back(prefix + "run_298596.2_re_reco.root");
	input_files.push_back(prefix + "run_298596.3_re_reco.root");
	input_files.push_back(prefix + "run_298596.4_re_reco.root");
	input_files.push_back(prefix + "run_298596.5_re_reco.root");
	input_files.push_back(prefix + "run_298596.6_re_reco.root");
	input_files.push_back(prefix + "run_298596.7_re_reco.root");
	input_files.push_back(prefix + "run_298596.8_re_reco.root");
	input_files.push_back(prefix + "run_298596.9_re_reco.root");
}
