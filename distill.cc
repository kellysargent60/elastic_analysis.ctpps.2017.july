#include "input_files.h"

#include "common_definitions.h"
#include "common_algorithms.h"
#include "parameters.h"
#include "common.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TChain.h"
#include "TH1D.h"

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 2)
		return 1;

	Init(argv[1]);
	if (diagonal == dCombined)
		return rcIncompatibleDiagonal;

	// prepare input
	InitInputFiles();
	printf("* input files: \n");
	for (unsigned int i = 0; i < input_files.size(); i++)
		printf("		%s\n", input_files[i].c_str());

	fwlite::ChainEvent event(input_files);

	printf("* events in input chain: %llu\n", event.size());

	// define RP ids
	unsigned int rpId_L_2_F = 0, rpId_L_1_F = 0, rpId_R_1_F = 0, rpId_R_2_F = 0;

	if (diagonal == ad45b_56b)
	{
		rpId_L_2_F = 25; rpId_L_1_F = 5; rpId_R_1_F = 105; rpId_R_2_F = 125;
	}

	if (diagonal == d45b_56t)
	{
		rpId_L_2_F = 25; rpId_L_1_F = 5; rpId_R_1_F = 104; rpId_R_2_F = 124;
	}

	if (diagonal == d45t_56b)
	{
		rpId_L_2_F = 24; rpId_L_1_F = 4; rpId_R_1_F = 105; rpId_R_2_F = 125;
	}

	if (diagonal == ad45t_56t)
	{
		rpId_L_2_F = 24; rpId_L_1_F = 4; rpId_R_1_F = 104; rpId_R_2_F = 124;
	}

	printf("* RP ids:\n");
	printf("    rpId_L_2_F = %u\n", rpId_L_2_F);
	printf("    rpId_L_1_F = %u\n", rpId_L_1_F);
	printf("    rpId_R_1_F = %u\n", rpId_R_1_F);
	printf("    rpId_R_2_F = %u\n", rpId_R_2_F);

	// prepare ouput file
	TFile *f_out = new TFile((string("distill_") + argv[1] + ".root").c_str(), "recreate");

	// set up output tree
	EventRed ev;
	TTree *t_out = new TTree("distilled", "bla");

	t_out->Branch("v_L_2_F", &ev.h.L_2_F.v); t_out->Branch("x_L_2_F", &ev.h.L_2_F.x); t_out->Branch("y_L_2_F", &ev.h.L_2_F.y);
	t_out->Branch("v_L_1_F", &ev.h.L_1_F.v); t_out->Branch("x_L_1_F", &ev.h.L_1_F.x); t_out->Branch("y_L_1_F", &ev.h.L_1_F.y);

	t_out->Branch("v_R_1_F", &ev.h.R_1_F.v); t_out->Branch("x_R_1_F", &ev.h.R_1_F.x); t_out->Branch("y_R_1_F", &ev.h.R_1_F.y);
	t_out->Branch("v_R_2_F", &ev.h.R_2_F.v); t_out->Branch("x_R_2_F", &ev.h.R_2_F.x); t_out->Branch("y_R_2_F", &ev.h.R_2_F.y);

	t_out->Branch("timestamp", &ev.timestamp);
	t_out->Branch("run_num", &ev.run_num);
	t_out->Branch("bunch_num", &ev.bunch_num);
	t_out->Branch("event_num", &ev.event_num);
	//t_out->Branch("trigger_num", &ev.trigger_num);
	//t_out->Branch("trigger_bits", &ev.trigger_bits);

	// loop over the chain entries
	unsigned int ev_count = 0;
	unsigned int ev_accepted = 0;
	for (event.toBegin(); ! event.atEnd(); ++event)
	{
		ev_count++;

		// TODO: remove
		//if (ev_count > 100000)
		//	break;

		// fill meta data
		ev.timestamp = event.time().unixTime() - timestamp0;
		ev.run_num = event.id().run();
		ev.bunch_num = event.bunchCrossing();
		ev.event_num = event.id().event();
		//ev.trigger_num = 0;
		//ev.trigger_bits = 0;

		// default output track data
		ev.h.L_2_F.v = false; ev.h.L_2_F.x = 0.; ev.h.L_2_F.y = 0.;
		ev.h.L_1_F.v = false; ev.h.L_1_F.x = 0.; ev.h.L_1_F.y = 0.;
		ev.h.R_1_F.v = false; ev.h.R_1_F.x = 0.; ev.h.R_1_F.y = 0.;
		ev.h.R_2_F.v = false; ev.h.R_2_F.x = 0.; ev.h.R_2_F.y = 0.;

		// load track data
		fwlite::Handle< vector<CTPPSLocalTrackLite> > tracks;
		tracks.getByLabel(event, "ctppsLocalTrackLiteProducer");

		// process track data
		for (const auto &tr : *tracks)
		{
			CTPPSDetId rpId(tr.getRPId());
			unsigned int rpDecId = 100*rpId.arm() + 10*rpId.station() + 1*rpId.rp();

			if (rpDecId == rpId_L_2_F)
			{
				ev.h.L_2_F.v = true; ev.h.L_2_F.x = tr.getX(); ev.h.L_2_F.y = tr.getY();
			}

			if (rpDecId == rpId_L_1_F)
			{
				ev.h.L_1_F.v = true; ev.h.L_1_F.x = tr.getX(); ev.h.L_1_F.y = tr.getY();
			}

			if (rpDecId == rpId_R_1_F)
			{
				ev.h.R_1_F.v = true; ev.h.R_1_F.x = tr.getX(); ev.h.R_1_F.y = tr.getY();
			}

			if (rpDecId == rpId_R_2_F)
			{
				ev.h.R_2_F.v = true; ev.h.R_2_F.x = tr.getX(); ev.h.R_2_F.y = tr.getY();
			}
		}

		// accept event ?
		bool accept = (ev.h.L_2_F.v && ev.h.L_1_F.v && ev.h.R_1_F.v && ev.h.R_2_F.v);

		if (!accept)
			continue;

		ev_accepted++;

		// save record
		t_out->Fill();
	}

	printf("* events processed: %i\n", ev_count);
	printf("* events accepted: %i\n", ev_accepted);

	// save output tree
	gDirectory = f_out;
	t_out->Write();

	// clean up
	delete f_out;
	return 0;
}
