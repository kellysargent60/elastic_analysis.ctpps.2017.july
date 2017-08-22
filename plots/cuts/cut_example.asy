import root;
import pad_layout;

string topDir = "../../";

TH2_palette = Gradient(blue, heavygreen, yellow, red);

string dataset = "DS-cr-angle140";
string diagonal = "45b_56t";

int cuts[] = { 1, 2, 3, 4, 5, 6, 7 };

int cuts[] = { 1, 2, 3, 4, 5, 6, 7 };

real scale_x[] = { 1e6, 1e6, 1e6, 1e6, 1e0, 1e0, 1e6, 1e6 };
real scale_y[] = { 1e6, 1e6, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0 };

string label_x[] = { "$\th_x^{*R}\ung{\mu rad}$", "$\th_y^{*R}\ung{\mu rad}$", "$\th_x^{*R}\ung{\mu rad}$", "$\th_x^{*L}\ung{\mu rad}$", "$y^{R,N}\ung{mm}$", "$y^{L,N}\ung{mm}$", "$\th_x^*\ung{\mu rad}$", "$\th_y^*\ung{\mu rad}$" };
string label_y[] = { "$\th_x^{*L}\ung{\mu rad}$", "$\th_y^{*L}\ung{\mu rad}$", "$x^{*R}\ung{mm}$", "$x^{*L}\ung{mm}$", "$y^{R,F} - y^{R,N}\ung{mm}$", "$y^{L,F} - y^{L,N}\ung{mm}$", "$\De^{R-L} x^*\ung{mm}$", "$\De^{R-L} y^*\ung{mm}$" };
string label_cut[] = { "$\De^{R-L} \th_x^{*}\ung{\mu rad}$", "$\De^{R-L} \th_y^{*}\ung{\mu rad}$", "$x^{*R}\ung{mm}$", "$x^{*L}\ung{mm}$", "$cq5$", "$cq6$", "$cq7$", "$cq8$" };

real lim_x_low[] = { -2000, -1000, -1000, -1000, -3, -3, -1000, -1000 };
real lim_x_high[] = { +2000, 1000, +1000, +1000, +7, +7, +1000, +1000 };

real lim_y_low[] = { -2000, -1000, -1, -1, -1.5, -1.5, -1.0, -4 };
real lim_y_high[] = { +2000, 1000, +1, +1, +1.5, +1.5, +1.0, +4 };

for (int ci : cuts.keys)
{
	if (ci == 4)
		NewRow();

	int cut = cuts[ci];
	int idx = cut - 1;

	string f = topDir + dataset+"/distributions_" + diagonal + ".root";
	
	NewPad(label_x[idx], label_y[idx]);
	scale(Linear, Linear, Log);
	string objC = format("elastic cuts/cut %i", cut) + format("/plot_after_cq%i", cut);
	draw(scale(scale_x[idx], scale_y[idx]), RootGetObject(f, objC+"#0"), "p,d0,bar");
	draw(scale(scale_x[idx], scale_y[idx]), RootGetObject(f, objC+"#1"));
	draw(scale(scale_x[idx], scale_y[idx]), RootGetObject(f, objC+"#2"));
	limits((lim_x_low[idx], lim_y_low[idx]), (lim_x_high[idx], lim_y_high[idx]), Crop);
	//limits((-lim_x_high[idx], -lim_y_high[idx]), (-lim_x_low[idx], -lim_y_low[idx]), Crop);

	AttachLegend(format("cut %i", cut));
}

GShipout(hSkip=1mm, vSkip=1mm);
