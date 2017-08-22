import root;
import pad_layout;
include "../run_info.asy";

string topDir = "../../";

string datasets[] = {
	"DS-cr-angle120",
	"DS-cr-angle130",
	"DS-cr-angle140",
};

string units[] = { "L_2_F", "L_1_F", "R_1_F", "R_2_F" };
string unit_labels[] = { "45-220-fr", "45-210-fr", "56-210-fr", "56-220-fr" };

xSizeDef = 12cm;
//xTicksDef = LeftTicks(Step=1, step=0.5);
drawGridDef = true;

TGraph_errorBar = None;

bool drawFits = true;

real samplingTime = 20 * 3600;

//----------------------------------------------------------------------------------------------------
NewRow();

for (int ui : units.keys)
{
	NewPad("time $\ung{h}$", "tilt $\ung{mrad}$", axesAbove=false);
	//currentpad.yTicks = RightTicks(5., 1.);
	DrawRunBands(-60, +60);

	for (int di : datasets.keys)
	{
		draw(swToHours, RootGetObject(topDir+datasets[di]+"/alignment.root", "global/"+units[ui]+"/a_p"), "p,l,eb", blue, mCi+1pt+blue);
		draw(swToHours, RootGetObject(topDir+datasets[di]+"/alignment.root", "global/"+units[ui]+"/a_g"), "p,l,eb", heavygreen, mCi+1pt+heavygreen);
		
		//draw(swToHours, RootGetObject(topDir+datasets[di]+"/alignment.root", "global/"+units[ui]+"/a"), "p,l,eb", blue, mCi+1pt+blue);
	}
	
	if (drawFits)
	{
		real unc = 15;
		RootGetObject(topDir+"alignment/global_fit.root", units[ui]+"/a_fit");
		draw(shift(0, +unc)*swToHours, robj, "l", red+dashed);
		draw(shift(0,    0)*swToHours, robj, "l", red+2pt);
		draw(shift(0, -unc)*swToHours, robj, "l", red+dashed);

		real corr = robj.rExec("Eval", samplingTime);
		AddToLegend(format("%.0f", corr), red);
	}

	//limits((time_min, 0), (time_max, +40), Crop);
	AttachLegend(unit_labels[ui], SE, SE);
}

//----------------------------------------------------------------------------------------------------
NewRow();

for (int ui : units.keys)
{
	NewPad("time $\ung{h}$", "horizontal position $\ung{\mu m}$", axesAbove=false);
	//currentpad.yTicks = RightTicks(20., 10.);
	DrawRunBands(-100, +1000);

	/*
	TGraph_reducePoints = 30;
	draw(unixToHours * shift(0, sh_x[ui]), RootGetObject("bpm.root", "LHC.BOFSU:POSITIONS_H::"+bpms[ui]), black);
	TGraph_reducePoints = 1; 
	*/

	for (int di : datasets.keys)
	{
		// TODO
		draw(swToHours, RootGetObject(topDir+datasets[di]+"/alignment.root", "global/"+units[ui]+"/b_p"), "p,l,eb", blue, mCi+1pt+blue);
		draw(swToHours, RootGetObject(topDir+datasets[di]+"/alignment.root", "global/"+units[ui]+"/b_g"), "p,l,eb", heavygreen, mCi+1pt+heavygreen);

		//draw(swToHours, RootGetObject(topDir+datasets[di]+"/alignment.root", "global/"+units[ui]+"/b"), "p,l,eb", blue+1pt, mCi+1pt+blue);
	
	}

	if (drawFits)
	{
		real unc = 80;
		RootGetObject(topDir+"/alignment/global_fit.root", units[ui]+"/b_fit");
		draw(shift(0, +unc)*swToHours, robj, "l", red+dashed);
		draw(shift(0,    0)*swToHours, robj, "l", red+2pt);
		draw(shift(0, -unc)*swToHours, robj, "l", red+dashed);

		real corr = robj.rExec("Eval", samplingTime);
		AddToLegend(format("%.0f", corr), red);
	}

	//limits((time_min, -100), (time_max, +100), Crop);
	AttachLegend(unit_labels[ui], SE, SE);
}

//----------------------------------------------------------------------------------------------------
NewRow();

/*
for (int ui : units.keys)
{
	NewPad("time $\ung{h}$", "vertical position $\ung{\mu m}$", axesAbove=false);
	currentpad.yTicks = RightTicks(100., 20.);
	//DrawRunBands(-500, +500);

	for (int di : datasets.keys)
	{
		pen p = StdPen(di+1);
		draw(swToHours, RootGetObject(topDir+datasets[di]+"/alignment.root", "global/"+units[ui]+"/c_min_diff"), "p,l,eb", blue, mCi+1pt+blue);
		//draw(swToHours, RootGetObject(topDir+datasets[di]+"/alignment.root", "global/"+units[ui]+"/c_prob"), "p,l,eb", blue, mCi+1pt+blue); // same as c_min_diff
		draw(swToHours, RootGetObject(topDir+datasets[di]+"/alignment.root", "global/"+units[ui]+"/c_mean_diff_sq"), "p,l,eb", magenta, mCi+1pt+magenta);
		draw(swToHours, RootGetObject(topDir+datasets[di]+"/alignment.root", "global/"+units[ui]+"/c_hist_chi_sq"), "p,l,eb", heavygreen, mCi+1pt+heavygreen);
		
		//draw(swToHours, RootGetObject(topDir+datasets[di]+"/alignment.root", "global/"+units[ui]+"/c"), "p,l,eb", blue+1pt, mCi+1pt+blue);

		if (drawFits)
		{
			real unc = 50;
			RootGetObject(topDir+datasets[di]+"/alignment_fit.root", ""+units[ui]+"/c_fit");
			draw(swToHours*shift(0, +unc), robj, "l", red+dashed);
			draw(swToHours*shift(0,    0), robj, "l", red+2pt);
			draw(swToHours*shift(0, -unc), robj, "l", red+dashed);
		}
	}

	//limits((time_min, -500), (time_max, +500), Crop);
	AttachLegend(unit_labels[ui], SE, SE);
}
*/

//----------------------------------------------------------------------------------------------------

GShipout();
