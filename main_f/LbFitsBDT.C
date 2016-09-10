#include "../TreeClasses/streeLb.h"
#include "../Funciones/LbMassFitter.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TVector3.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooChebychev.h"
#include "RooAbsReal.h"
#include <iostream>
#include <vector>
#include <string>
using namespace RooFit;

int main( int argc ,char *argv[] ){
	//if ( argc < 2){
	//	std::cout << "1: MCmatch, 0: no MCmatch" << std::endl;
	//	return 1;
	//}
	const char* binname = argv[1];
	Double_t    min     = atof(argv[2]);
	Double_t    max     = atof(argv[3]);
	Int_t       bins    = atoi(argv[4]);	
	
	//std::cout << binname << std::endl;
	TFile *input(0);
	TString fname;
	//fname =".datos/datosLb.root";
	fname ="MVA/TMVApp.root";
	std::cout << "using " << fname << std::endl;
	input = TFile::Open( fname );

	TTree* tree = (TTree*)input->Get("OniattTree");
	//streeLb *t = new streeLb(tree);
	
	Double_t Lbmass, binvar; 
	
	tree->SetBranchAddress( "Lb_mass", &Lbmass);
	tree->SetBranchAddress( binname,   &binvar);

	Double_t st = (max - min )/bins;

	RooRealVar mass("mass","m(J/#psi Kp)",5.6197,5.494,5.746,"GeV");
	RooDataSet *datasets[bins];


	for (Int_t b = 0; b < bins; ++b){
		std::string s = "bin"+std::to_string(b);
		const char* name = s.c_str();
		datasets[b] = new RooDataSet(name,"",RooArgSet(mass));
	}
	
	std::vector<double> bdt, signif;

	std::cout << "--- Processing: " << tree->GetEntries() << " events" << std::endl;
	
	Long64_t nevents = tree->GetEntries();
	for (Long64_t entry = 0 ; entry < nevents; ++entry){
		//if (entry > 200000) break;
		if (entry%100000 == 0) std::cout << "--- ... Processing event: " << entry << std::endl;
		tree->GetEntry(entry);
		if ( Lbmass < 5.494 ) continue;
		if ( Lbmass > 5.746 ) continue;
		mass.setVal(Lbmass);
		for (Int_t j = 0; j < bins; ++j){
			if ( binvar > min + st*j ) datasets[j]->add(RooArgSet(mass));
		}
	}
	
	LbMassFitter* LbFitter = new LbMassFitter();
	//else LbFitter->fitMCmatch(mass,dataset);
	for (Int_t b = 0; b < bins; ++b){
		Double_t significance = LbFitter->fitBinBDT(mass,datasets[b], binname, bins, b, min, st);
		bdt.push_back(min + st*b );
		signif.push_back(significance);
	}

	Int_t nn = bdt.size();
	Double_t arrBDT[nn];
    Double_t arrsignif[nn];

	for ( Int_t ii = 0; ii < bdt.size(); ++ii){
		arrBDT[ii]    = bdt.at(ii);
		arrsignif[ii] = signif.at(ii);
	}
	TGraph *gr = new TGraph(nn,arrBDT,arrsignif);
	gr->SetLineColor(8);
	gr->SetTitle("Optimal Cut");
	gr->GetXaxis()->SetTitle("Cut Value applied on BDT output");
	gr->GetYaxis()->SetTitle("Significance");
	gr->Draw("AC*");
	TCanvas *c = new TCanvas("c", "c", 1000,600);
	c->cd();
	gr->Draw();
	//gStyle->SetEndErrorSize(3);
	//gStyle->SetErrorX(1.);
	//h->SetMarkerStyle(20);
	//h->Draw();
	c->SaveAs("sigBDT.png");
	return 0;
}
