#include <string>
#include <vector>

//----------------------------------------------------------------------------------------------------

std::vector<std::string> input_files;

void InitInputFiles()
{
        std::string prefix;

        prefix="root://eostotem.cern.ch//eos/totem/data/ctpps/reconstruction/2017/alignment_run_July/version4-pilot/";
        input_files.push_back(prefix + "run_298593.0_re_reco.root");
        input_files.push_back(prefix + "run_298593.10_re_reco.root");
        input_files.push_back(prefix + "run_298593.11_re_reco.root");
        input_files.push_back(prefix + "run_298593.12_re_reco.root");
        input_files.push_back(prefix + "run_298593.1_re_reco.root");
        input_files.push_back(prefix + "run_298593.2_re_reco.root");
        input_files.push_back(prefix + "run_298593.3_re_reco.root");
        input_files.push_back(prefix + "run_298593.4_re_reco.root");
        input_files.push_back(prefix + "run_298593.5_re_reco.root");
        input_files.push_back(prefix + "run_298593.6_re_reco.root");
        input_files.push_back(prefix + "run_298593.7_re_reco.root");
        input_files.push_back(prefix + "run_298593.8_re_reco.root");
        input_files.push_back(prefix + "run_298593.9_re_reco.root");

        prefix="root://eostotem.cern.ch//eos/totem/data/ctpps/reconstruction/2017/alignment_run_July/version4-rest/";
        input_files.push_back(prefix + "run_298593.0_re_reco.root");
        input_files.push_back(prefix + "run_298593.100_re_reco.root");
        input_files.push_back(prefix + "run_298593.101_re_reco.root");
        input_files.push_back(prefix + "run_298593.102_re_reco.root");
        input_files.push_back(prefix + "run_298593.103_re_reco.root");
        input_files.push_back(prefix + "run_298593.104_re_reco.root");
        input_files.push_back(prefix + "run_298593.108_re_reco.root");
        input_files.push_back(prefix + "run_298593.109_re_reco.root");
        input_files.push_back(prefix + "run_298593.10_re_reco.root");
        input_files.push_back(prefix + "run_298593.110_re_reco.root");
        input_files.push_back(prefix + "run_298593.13_re_reco.root");
        input_files.push_back(prefix + "run_298593.14_re_reco.root");
        input_files.push_back(prefix + "run_298593.16_re_reco.root");
        input_files.push_back(prefix + "run_298593.17_re_reco.root");
        input_files.push_back(prefix + "run_298593.19_re_reco.root");
        input_files.push_back(prefix + "run_298593.1_re_reco.root");
        input_files.push_back(prefix + "run_298593.22_re_reco.root");
        input_files.push_back(prefix + "run_298593.23_re_reco.root");
        input_files.push_back(prefix + "run_298593.24_re_reco.root");
        input_files.push_back(prefix + "run_298593.25_re_reco.root");
        input_files.push_back(prefix + "run_298593.26_re_reco.root");
        input_files.push_back(prefix + "run_298593.27_re_reco.root");
        input_files.push_back(prefix + "run_298593.28_re_reco.root");
        input_files.push_back(prefix + "run_298593.29_re_reco.root");
        input_files.push_back(prefix + "run_298593.2_re_reco.root");
        input_files.push_back(prefix + "run_298593.30_re_reco.root");
        input_files.push_back(prefix + "run_298593.31_re_reco.root");
        input_files.push_back(prefix + "run_298593.32_re_reco.root");
        input_files.push_back(prefix + "run_298593.33_re_reco.root");
        input_files.push_back(prefix + "run_298593.34_re_reco.root");
        input_files.push_back(prefix + "run_298593.35_re_reco.root");
        input_files.push_back(prefix + "run_298593.36_re_reco.root");
        input_files.push_back(prefix + "run_298593.38_re_reco.root");
        input_files.push_back(prefix + "run_298593.40_re_reco.root");
        input_files.push_back(prefix + "run_298593.41_re_reco.root");
        input_files.push_back(prefix + "run_298593.43_re_reco.root");
        input_files.push_back(prefix + "run_298593.45_re_reco.root");
        input_files.push_back(prefix + "run_298593.47_re_reco.root");
        input_files.push_back(prefix + "run_298593.48_re_reco.root");
        input_files.push_back(prefix + "run_298593.49_re_reco.root");
        input_files.push_back(prefix + "run_298593.4_re_reco.root");
        input_files.push_back(prefix + "run_298593.51_re_reco.root");
        input_files.push_back(prefix + "run_298593.52_re_reco.root");
        input_files.push_back(prefix + "run_298593.53_re_reco.root");
        input_files.push_back(prefix + "run_298593.55_re_reco.root");
        input_files.push_back(prefix + "run_298593.56_re_reco.root");
        input_files.push_back(prefix + "run_298593.57_re_reco.root");
        input_files.push_back(prefix + "run_298593.58_re_reco.root");
        input_files.push_back(prefix + "run_298593.60_re_reco.root");
        input_files.push_back(prefix + "run_298593.61_re_reco.root");
        input_files.push_back(prefix + "run_298593.62_re_reco.root");
        input_files.push_back(prefix + "run_298593.63_re_reco.root");
        input_files.push_back(prefix + "run_298593.64_re_reco.root");
        input_files.push_back(prefix + "run_298593.65_re_reco.root");
        input_files.push_back(prefix + "run_298593.66_re_reco.root");
        input_files.push_back(prefix + "run_298593.67_re_reco.root");
        input_files.push_back(prefix + "run_298593.68_re_reco.root");
        input_files.push_back(prefix + "run_298593.69_re_reco.root");
        input_files.push_back(prefix + "run_298593.6_re_reco.root");
        input_files.push_back(prefix + "run_298593.70_re_reco.root");
        input_files.push_back(prefix + "run_298593.71_re_reco.root");
        input_files.push_back(prefix + "run_298593.72_re_reco.root");
        input_files.push_back(prefix + "run_298593.73_re_reco.root");
        input_files.push_back(prefix + "run_298593.75_re_reco.root");
        input_files.push_back(prefix + "run_298593.76_re_reco.root");
        input_files.push_back(prefix + "run_298593.77_re_reco.root");
        input_files.push_back(prefix + "run_298593.79_re_reco.root");
        input_files.push_back(prefix + "run_298593.80_re_reco.root");
        input_files.push_back(prefix + "run_298593.82_re_reco.root");
        input_files.push_back(prefix + "run_298593.83_re_reco.root");
        input_files.push_back(prefix + "run_298593.85_re_reco.root");
        input_files.push_back(prefix + "run_298593.86_re_reco.root");
        input_files.push_back(prefix + "run_298593.87_re_reco.root");
        input_files.push_back(prefix + "run_298593.88_re_reco.root");
        input_files.push_back(prefix + "run_298593.89_re_reco.root");
        input_files.push_back(prefix + "run_298593.8_re_reco.root");
        input_files.push_back(prefix + "run_298593.90_re_reco.root");
        input_files.push_back(prefix + "run_298593.91_re_reco.root");
        input_files.push_back(prefix + "run_298593.93_re_reco.root");
        input_files.push_back(prefix + "run_298593.94_re_reco.root");
        input_files.push_back(prefix + "run_298593.95_re_reco.root");
        input_files.push_back(prefix + "run_298593.96_re_reco.root");
        input_files.push_back(prefix + "run_298593.97_re_reco.root");
        input_files.push_back(prefix + "run_298593.98_re_reco.root");
        input_files.push_back(prefix + "run_298593.99_re_reco.root");
        input_files.push_back(prefix + "run_298593.9_re_reco.root");
}
