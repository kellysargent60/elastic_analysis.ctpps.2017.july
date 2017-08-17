import root;
import pad_layout;
include "../run_info.asy";

string topDir = "../../";

TH2_palette = Gradient(blue, heavygreen, yellow, red);

int rebin = 50;

string datasets[] = { "DS1" };

string dgns[] = { "45b_56t", "45t_56b" };
string dgn_labs[] = { "45 bot -- 56 top", "45 top -- 56 bot" };

xSizeDef = 12cm;
xTicksDef = LeftTicks(Step=1, step=0.5);

//----------------------------------------------------------------------------------------------------

for (int dgni : dgns.keys)
{
	NewPad("time$\ung{s}$", "rate$\ung{Hz}$");
	DrawRunBands(0, +3);
	for (int di : datasets.keys)
	{
		RootObject obj = RootGetObject(topDir+datasets[di]+"/distributions_"+dgns[dgni]+".root", "metadata/h_timestamp_dgn");
		obj.vExec("Rebin", rebin);
		draw(scale(1, 1./rebin)*swToHours, obj, "vl", black);
	
		RootObject obj = RootGetObject(topDir+datasets[di]+"/distributions_"+dgns[dgni]+".root", "metadata/h_timestamp_sel");
		obj.vExec("Rebin", rebin);
		draw(scale(1, 1./rebin)*swToHours, obj, "vl", red);
	}

	limits((time_min, 0), (time_max, +3), Crop);
	layer();
	AttachLegend(dgn_labs[dgni], E, E);
}
