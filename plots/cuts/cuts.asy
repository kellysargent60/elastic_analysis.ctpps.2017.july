import root;
import pad_layout;

string topDir = "../../";

TH2_palette = Gradient(blue, heavygreen, yellow, red);

string datasets[] = { "DS1" };
//string dgns[] = { "45b_56t", "45t_56b" };
string dgns[] = { "45t_56b" };

//int cuts[] = { 1, 2, 3, 4, 5, 6, 7, 8 };
int cuts[] = { 2, 5, 6 };

real scale_x[] = { 1e6, 1e6, 1e6, 1e6, 1e0, 1e0, 1e6, 1e6 };
real scale_y[] = { 1e6, 1e6, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0 };

string label_x[] = { "$\th_x^{*R}\ung{\mu rad}$", "$\th_y^{*R}\ung{\mu rad}$", "$\th_x^{*R}\ung{\mu rad}$", "$\th_x^{*L}\ung{\mu rad}$", "$y^{R,N}\ung{mm}$", "$y^{L,N}\ung{mm}$", "$\th_x^*\ung{\mu rad}$", "$\th_y^*\ung{\mu rad}$" };
string label_y[] = { "$\th_x^{*L}\ung{\mu rad}$", "$\th_y^{*L}\ung{\mu rad}$", "$x^{*R}\ung{mm}$", "$x^{*L}\ung{mm}$", "$y^{R,F} - y^{R,N}\ung{mm}$", "$y^{L,F} - y^{L,N}\ung{mm}$", "$\De^{R-L} x^*\ung{mm}$", "$\De^{R-L} y^*\ung{mm}$" };
string label_cut[] = { "$\De^{R-L} \th_x^{*}\ung{\mu rad}$", "$\De^{R-L} \th_y^{*}\ung{\mu rad}$", "$x^{*R}\ung{mm}$", "$x^{*L}\ung{mm}$", "$cq5$", "$cq6$", "$cq7$", "$cq8$" };

real lim_x_low[] = { -200, -150, -1000, -1000, -15, -15, -200, -600 };
real lim_x_high[] = { +200, 150, +1000, +1000, +15, +15, +200, +600 };

real lim_y_low[] = { -200, -150, -0.8, -0.8, -0.5, -0.5, -0.2, -4 };
real lim_y_high[] = { +200, 150, +0.8, +0.8, +0.5, +0.5, +0.2, +4 };

for (int ci : cuts.keys)
{
	int cut = cuts[ci];
	int idx = cut - 1;

	NewPad(false);
	
	NewPad(false);
	label("\SetFontSizesXX before cuts");
	
	NewPad(false);
	label("\SetFontSizesXX after cuts");

	for (int dsi : datasets.keys)
	{
		string dataset = datasets[dsi];

		for (int dgi : dgns.keys)
		{
			string dgn = dgns[dgi];
			string f = topDir + dataset+"/distributions_" + dgn + ".root";
	
			NewRow();

			NewPad(false);
			label("\vbox{\SetFontSizesXX\hbox{"+format("cut %i", cut)+"}\hbox{"+dataset+"}\hbox{"+replace(dgn, "_", "--")+"}}");
	
			NewPad(label_x[idx], label_y[idx]);
			scale(Linear, Linear, Log);
			string objC = format("elastic cuts/cut %i", cut) + format("/plot_before_cq%i", cut);
			draw(scale(scale_x[idx], scale_y[idx]), RootGetObject(f, objC+"#0"), "p,d0,bar");
			draw(scale(scale_x[idx], scale_y[idx]), RootGetObject(f, objC+"#1"));
			draw(scale(scale_x[idx], scale_y[idx]), RootGetObject(f, objC+"#2"));
			limits((-lim_x_high[idx], -lim_y_high[idx]), (-lim_x_low[idx], -lim_y_low[idx]), Crop);
			
			NewPad(label_x[idx], label_y[idx]);
			scale(Linear, Linear, Log);
			string objC = format("elastic cuts/cut %i", cut) + format("/plot_after_cq%i", cut);
			draw(scale(scale_x[idx], scale_y[idx]), RootGetObject(f, objC+"#0"), "p,d0,bar");
			draw(scale(scale_x[idx], scale_y[idx]), RootGetObject(f, objC+"#1"));
			draw(scale(scale_x[idx], scale_y[idx]), RootGetObject(f, objC+"#2"));
			limits((-lim_x_high[idx], -lim_y_high[idx]), (-lim_x_low[idx], -lim_y_low[idx]), Crop);
			
			NewPad(label_cut[idx]);
			string objH = format("elastic cuts/cut %i", cut) + format("/h_cq%i", cut);
			draw(scale(scale_y[idx], 1.), RootGetObject(f, objH+""), "vl,eb,lR");
			//draw(scale(scale_x[idx], scale_y[idx]), RootGetObject(f, objH+"|gaus"));
			AttachLegend();
		}
	}

	GShipout(format("cut_%i", cut));
}
