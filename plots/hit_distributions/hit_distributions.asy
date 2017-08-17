import root;
import pad_layout;

string topDir = "../../";

TH2_palette = Gradient(blue, heavygreen, yellow, red);

string dataset = "DS-cr-angle140";
string file_45b = topDir + dataset + "/distributions_45b_56t.root";
string file_45t = topDir + dataset + "/distributions_45t_56b.root";
string rps[] = { "L_2_F", "L_1_F", "R_1_F", "R_2_F" };
string rp_labels[] = { "left, 220, far", "left, 210, far", "right, 210, far", "right, 220, far" };

real sh_top[] = { 6.6, 6.0, 5.7,     7.2, 7.5, 8.2 };
real sh_bot[] = { -6.8, -6.4, -6.2, -6.3, -6.4, -7.1 };

drawGridDef = true;

TH2_x_min = -3;
TH2_x_max = +3;
TH2_y_min = -3;
TH2_y_max = +7;

//----------------------------------------------------------------------------------------------------

real edge = 3.6101;
real cutEdge = 2.22721 / sqrt(2);
int strips = 11;
real margin_v_e = 0.2;
real margin_v_b = 0.4;
real margin_u = 0.1;

path det_shape = (cutEdge, 0)--(edge, 0)--(edge, edge)--(0, edge)--(0, cutEdge)--cycle;
det_shape = scale(10) * rotate(45) * shift(-cutEdge/2, -cutEdge/2) * det_shape;

path hor_det_shape = shift(0, -cutEdge/sqrt(2)*10) * det_shape;

//----------------------------------------------------------------------------------------------------

NewRow();

//NewPad(false);
//label("{\SetFontSizesXX after selection}");

for (int ri : rps.keys)
{
	NewPad("$x\ung{mm}$", "$y\ung{mm}$", 6cm, 6cm);
	scale(Linear, Linear, Log);

	//TH2_z_max = log10(6e2);

	//draw(RootGetObject(file_45b, "hit distributions/vertical, aligned, before selection/h_y_"+rps[ri]+"_vs_x_"+rps[ri]+"_al_nosel"), "p,bar");
	//draw(RootGetObject(file_45t, "hit distributions/vertical, aligned, before selection/h_y_"+rps[ri]+"_vs_x_"+rps[ri]+"_al_nosel"), "p");

	RootGetObject(file_45b, "hit distributions/vertical, aligned, after selection/h_y_"+rps[ri]+"_vs_x_"+rps[ri]+"_al_sel");
	//if (abs(robj.rExec("GetMean", 2)) > 0.1)
		draw(robj, "p,bar");
	draw(RootGetObject(file_45t, "hit distributions/vertical, aligned, after selection/h_y_"+rps[ri]+"_vs_x_"+rps[ri]+"_al_sel"), "p");
	
	//draw(RootGetObject(file_45t, "hit distributions/vertical, not aligned, after selection/h_y_"+rps[ri]+"_vs_x_"+rps[ri]+"_noal_sel"), "p");
	
	/*
	draw(shift(0, sh_top[ri])*det_shape);
	draw(shift(0, sh_bot[ri])*scale(1, -1)*det_shape);
	*/

	//draw(shift(6.2, 0)*rotate(-90)*det_shape);

	limits((-3, -3), (+3, +7), Crop);
	//limits((-30, -9), (+30, +9), Crop);

	AttachLegend(replace(rp_labels[ri], "_", "\_"));
}
