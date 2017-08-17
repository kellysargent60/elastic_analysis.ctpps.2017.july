import root;
import pad_layout;
include "../run_info.asy";

string topDir = "../../";

string dataset = "DS-cr-angle140";

string units[] = { "L_2_F", "L_1_F", "R_1_F", "R_2_F" };
string unit_labels[] = { "left, 220, far", "left, 210, far", "right, 210, far", "right, 220, far" };

xSizeDef = 10cm;
drawGridDef = true;

TGraph_errorBar = None;

string period = "period 6";

//----------------------------------------------------------------------------------------------------
NewRow();

for (int ui : units.keys)
{
	NewPad("$y\ung{mm}$", "$\hbox{mean } x\ung{mm}$");

	draw(RootGetObject(topDir+dataset+"/alignment.root", period + "/unit "+units[ui]+"/horizontal/horizontal profile/p"), "d0,eb", blue);
	draw(RootGetObject(topDir+dataset+"/alignment.root", period + "/unit "+units[ui]+"/horizontal/horizontal profile/p|ff"), "l", red+1pt);
	
	limits((-7, 2), (+7, 4), Crop);
	AttachLegend(unit_labels[ui], NE, NE);
}

//----------------------------------------------------------------------------------------------------
NewRow();

/*

for (int ui : units.keys)
{
	NewPad("$y\ung{mm}$", "");
	scale(Linear, Log);

	draw(RootGetObject(topDir+dataset+"/alignment.root", period + "/unit "+units[ui]+"/vertical/y_hist|y_hist"), "d0,vl", blue);
	draw(RootGetObject(topDir+dataset+"/alignment.root", period + "/unit "+units[ui]+"/vertical/y_hist|y_hist_range"), "d0,vl", red);

	limits((-7, 1), (+7, 2e2), Crop);
	AttachLegend(unit_labels[ui], NE, NE);
}


//----------------------------------------------------------------------------------------------------

NewRow();

for (int ui : units.keys)
{
	NewPad("$\de y\ung{mm}$", "");

	if (units[ui] == "R_1_N")
		continue;

	draw(RootGetObject(topDir+dataset+"/alignment.root", period + "/unit "+units[ui]+"/vertical/g_max_diff"), "l,p", heavygreen, mCi+1pt+heavygreen);

	limits((-2.5, 0), (+2.5, 1), Crop);
	AttachLegend(unit_labels[ui], NE, NE);
}

*/
