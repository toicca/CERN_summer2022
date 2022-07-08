void tdrDraw(TH1* h, string opt,
	     int marker=kFullCircle, int mcolor = kBlack,
	     int lstyle=kSolid, int lcolor=-1,
	     int fstyle=1001, int fcolor=kYellow+1);

void tdrDraw(TGraph* g, string opt,
	     int marker=kFullCircle, int mcolor = kBlack,
	     int lstyle=kSolid, int lcolor=-1,
	     int fstyle=1001, int fcolor=kYellow+1,
	     double msize=1);

TLegend *tdrLeg(double x1, double y1, double x2, double y2);

TH1D *tdrHist(string name="h", string ylabel="Response",
	      double y1=0, double y2=1,
	      string xlabel="p_{T} (GeV)",
	      double x1=15, double x2=3500);

void tdrGrid(bool gridOn);

void fixOverlay();

void setTDRStyle();

void CMS_lumi( TPad* pad, int iPeriod, int iPosX );

TCanvas* tdrCanvas(const char* canvName, TH1D *h,
		   int iPeriod = 2, int iPos = 11,
		   bool square = kRectangular);

TCanvas* tdrDiCanvas(const char* canvName, TH1D *hup, TH1D *hdw,
		   int iPeriod = 2, int iPos = 11);