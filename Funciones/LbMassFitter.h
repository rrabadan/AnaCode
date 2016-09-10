#ifndef __LbMassFitter__
#define __LbMassFitter__
#include <iostream>
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooChebychev.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooMinuit.h"
#include "RooWorkspace.h"
#include "RooAbsArg.h"
#include "RooAbsRealLValue.h"

using namespace std;

class LbMassFitter {
		//int width, height;
		RooWorkspace *ws;
	public:
		//void fit(RooWorkspace* ws, int cut);
		//void set_values (int,int);
	    //int area() {return width*height;}
		void fitMC(RooRealVar& mass, RooDataSet* data);
		void fitMCmatch(RooRealVar& mass, RooDataSet* data);
		void fitdatos(RooRealVar& mass, RooDataSet* data, RooWorkspace* ws);
		RooRealVar* fitBin(RooRealVar& mass, RooDataSet* data, const char* binVar, Int_t& nbins, Int_t& bin, Double_t& xmin, Double_t& width);
		Double_t fitBinBDT(RooRealVar& mass, RooDataSet* data, const char* binVar, Int_t& nbins, Int_t& bin, Double_t& xmin, Double_t& width);
	private:
		void doFit(RooWorkspace* ws);
		double significance(const double& nsig, const double& nbkg){return nsig/sqrt(nsig+nbkg);};
};

#endif
