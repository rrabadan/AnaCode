#include "LbMassFitter.h"
#include "RooFitResult.h"
#include "RooChebychev.h"
#include "RooPlot.h"
#include "RooAbsReal.h"
#include <string> 
#include <sstream>

using namespace RooFit;

RooRealVar* LbMassFitter::fitBin(RooRealVar& mass, RooDataSet* data, const char* binVar, Int_t& nbins, Int_t& bin, Double_t& xmin, Double_t& width){
	ostringstream ssname, sstitle;
	
	ssname << "Lbfit_"; ssname << binVar; ssname << data->GetName();
	std::string name = "someplots/"+ssname.str()+".png";
	
	sstitle << "fit "; sstitle << binVar << " " << xmin + bin*width << "-" << xmin + (bin+1)*width; 
	std::string title = sstitle.str();

	const char *namechar  = name.c_str();
	const char *titlechar = title.c_str();
	//std::cout << namechar << std::endl;
	//std::cout << titlechar << std::endl;
	
	ws = new RooWorkspace("workspace");
	fitdatos(mass, data, ws);

	RooRealVar* Nsig = ws->var("Nsig");
	
	RooAbsPdf* model = ws->pdf("model");
    TCanvas *c = new TCanvas("c", "c", 900,700);
    RooPlot* massframe = mass.frame(Title(titlechar), Bins(80));
	data->plotOn(massframe);
	model->plotOn(massframe);
	model->plotOn(massframe,Components("sigmod"),LineStyle(kDashed),LineColor(kMagenta));
	model->plotOn(massframe,Components("bkgmod"),LineStyle(kDashed),LineColor(kGray));
	model->paramOn(massframe, Layout(0.60,0.97,0.88));
	c->cd(); gPad->SetLeftMargin(0.15) ; massframe->GetYaxis()->SetTitleOffset(1.6) ;
	massframe->Draw();
	//c->SetTitle(titlechar);
	c->Update();
	c->SaveAs(namechar);
	delete c;
	delete ws;

	return Nsig;
}

Double_t LbMassFitter::fitBinBDT(RooRealVar& mass, RooDataSet* data, const char* binVar, Int_t& nbins, Int_t& bin, Double_t& xmin, Double_t& width){
    ostringstream ssname, sstitle;
   
    ssname << "Lbfit_"; ssname << binVar; ssname << data->GetName();
    std::string name = "someplots/"+ssname.str()+".png";
   
	sstitle << "fit "; sstitle << binVar << " > " << xmin + bin*width;
	std::string title = sstitle.str();

	const char *namechar  = name.c_str();
	const char *titlechar = title.c_str();
	//std::cout << namechar << std::endl;
	//std::cout << titlechar << std::endl;
	
	ws = new RooWorkspace("workspace");
	fitdatos(mass, data, ws);

	RooRealVar* Nsig = ws->var("Nsig");
	RooRealVar* Nbkg = ws->var("Nbkg");

	RooAbsPdf* model = ws->pdf("model");
	TCanvas *c = new TCanvas("c", "c", 900,700);
	RooPlot* massframe = mass.frame(Title(titlechar), Bins(80));
	data->plotOn(massframe);
	model->plotOn(massframe);
	model->plotOn(massframe,Components("sigmod"),LineStyle(kDashed),LineColor(kMagenta));
	model->plotOn(massframe,Components("bkgmod"),LineStyle(kDashed),LineColor(kGray));
	model->paramOn(massframe, Layout(0.60,0.97,0.88));
	c->cd(); gPad->SetLeftMargin(0.15) ; massframe->GetYaxis()->SetTitleOffset(1.6) ;
	massframe->Draw();
	//c->SetTitle(titlechar);
	c->Update();
	c->SaveAs(namechar);

	delete c;
    delete ws;

    return Nsig->getVal()/sqrt(Nsig->getVal()+Nbkg->getVal());
}

void LbMassFitter::fitdatos(RooRealVar& mass, RooDataSet* data, RooWorkspace* ws){

	data->SetName("data");
	bool plot = 0;
	if (!ws){ ws = new RooWorkspace("workspace"); plot = 1; }
	
	ws->import(mass);
	ws->import(*data);

	ws->factory("Gaussian::sig1(mass,mean[5.6197],s1[0.00856])");
	ws->factory("Gaussian::sig2(mass,mean,s2[0.02112])");
	ws->factory("Gaussian::sig3(mass,mean,s3[0.0903])");
	ws->factory("SUM::sigmod(frac[0.6]*sig1,frac2[0.385]*sig2,sig3)");
	//if (!ws) 
	ws->factory("Chebychev::bkgmod(mass,{B1[0.,-10.,10.],B2[0.,-10.,10.],B3[0.,-10.,10.]})");
	//else ws->factory("Chebychev::bkgmod(mass,{B1[0.,-10.,10.]})");	
	ws->factory("SUM::model( Nsig[0,15000]*sigmod, Nbkg[0,1500000]*bkgmod )") ;

	doFit(ws);

	if (plot){
		RooAbsPdf* model = ws->pdf("model");
    	TCanvas *c = new TCanvas("c", "c", 900,700);
		RooPlot* massframe = mass.frame(Title(" "), Bins(80));
    	data->plotOn(massframe);
    	model->plotOn(massframe);
    	model->plotOn(massframe,Components("sigmod"),LineStyle(kDashed),LineColor(kMagenta));
    	model->plotOn(massframe,Components("bkgmod"),LineStyle(kDashed),LineColor(kGray));
    	model->paramOn(massframe, Layout(0.60,0.97,0.88));
    	c->cd(); gPad->SetLeftMargin(0.15) ; massframe->GetYaxis()->SetTitleOffset(1.6) ;
    	massframe->Draw();
    	c->Update();
    	c->SaveAs("someplots/Datafit.png");
    	mass.setRange("R3",5.54,5.70);
		RooRealVar* Nsig = ws->var("Nsig");
		RooRealVar* Nbkg = ws->var("Nbkg");
		RooAbsPdf* bkgmod = ws->pdf("bkgmod");
		RooAbsReal* Int = bkgmod->createIntegral(mass,NormSet(mass),Range("R3"));
		std::cout << "Ns: " << Nsig->getVal() << std::endl;
		std::cout << "Nb: " << Nbkg->getVal()*Int->getVal() << std::endl;
		delete c;
		delete ws;
	}

}

void LbMassFitter::fitMC(RooRealVar& mass, RooDataSet* data){
	
	mass.Print();
	data->Print();
	data->SetName("data");

	ws = new RooWorkspace("workspace");
	ws->import(mass);
	ws->import(*data);

	ws->factory("Gaussian::sig1(mass,mean[5.6197,5.38,5.86],s1[0.01065,0.0,1.0])");
	ws->factory("Gaussian::sig2(mass,mean,s2[0.03606,0.0,1.0])");
	ws->factory("SUM::sigmod(frac[0.8723,0.0,1.0]*sig1,sig2)");
	ws->factory("Chebychev::bkgmod(mass,{B1[0.,-10.,10.]})");
	ws->factory("SUM::model( Nsig[0,2000000]*sigmod, Nbkg[0,2000000]*bkgmod )") ;
	ws->Print();
	
	doFit(ws);

	RooAbsPdf* model = ws->pdf("model");
	TCanvas *c = new TCanvas("c", "c", 900,700);
	RooPlot* massframe = mass.frame(Title(" "), Bins(80));	
	data->plotOn(massframe);
	model->plotOn(massframe); 
	model->plotOn(massframe,Components("sig1"),LineStyle(kDashed),LineColor(kRed));
	model->plotOn(massframe,Components("sig2"),LineStyle(kDashed),LineColor(kMagenta));
	model->plotOn(massframe,Components("bkgmod"),LineStyle(kDashed),LineColor(kGray));
	model->paramOn(massframe, Layout(0.60,0.97,0.88));
	c->cd(); gPad->SetLeftMargin(0.15) ; massframe->GetYaxis()->SetTitleOffset(1.6) ;
	massframe->Draw();
	c->Update();
	c->SaveAs("someplots/MCfit.png");
	delete c;
	
	delete ws;

	 /*ostringstream ss;
	 ss << cut;
	 std::string name = "Lbfit_"+ ss.str()+".png";
	 const char *namechar = name.c_str();
	 ws->factory("SUM::model( Nsig[0,2000000]*sigmod, Nbkg[0,2000000]*bkgmod )") ;
	 */
}

void LbMassFitter::fitMCmatch(RooRealVar& mass, RooDataSet* data){

    mass.Print();
    data->Print();
    data->SetName("data");
	//mass.setRange(5.46,5.78);
    ws = new RooWorkspace("workspace");
    ws->import(mass);
    ws->import(*data);

    ws->factory("Gaussian::sig1(mass,mean[5.6197,5.38,5.86],s1[0.01065,0.0,1.0])");
    ws->factory("Gaussian::sig2(mass,mean,s2[0.03606,0.0,1.0])");
	ws->factory("Gaussian::sig3(mass,mean,s3[0.05,0.0,1.0])");
    //ws->factory("SUM::model(frac[0.8723,0.0,1.0]*sig1,sig2)");
	ws->factory("SUM::model(frac[0.8723,0.0,1.0]*sig1,frac2[0.1,0.0,1.0]*sig2,sig3)");
	ws->Print();
    
	doFit(ws);

	RooAbsPdf* model = ws->pdf("model");
	TCanvas *c = new TCanvas("c", "c", 900,700);
	RooPlot* massframe = mass.frame(Title(" "), Bins(80));
	data->plotOn(massframe);
	model->plotOn(massframe);
	model->plotOn(massframe,Components("sig1"),LineStyle(kDashed),LineColor(kRed));
	model->plotOn(massframe,Components("sig2"),LineStyle(kDashed),LineColor(kMagenta));
	model->plotOn(massframe,Components("sig3"),LineStyle(kDashed),LineColor(kGray));
	model->paramOn(massframe, Layout(0.60,0.97,0.88));
	c->cd(); c->cd(); gPad->SetLeftMargin(0.15) ; massframe->GetYaxis()->SetTitleOffset(1.6) ;
	massframe->Draw();
	c->Update();
	c->SaveAs("someplots/MCmatchfit.png");
	delete c;

    delete ws;

}

void LbMassFitter::doFit(RooWorkspace* ws){
	
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    RooMsgService::instance().setSilentMode(kTRUE);
    RooMsgService::instance().setStreamStatus(1,false);
    RooMsgService::instance().getStream(1).removeTopic(Integration) ;
    RooMsgService::instance().getStream(1).removeTopic(Minimization) ;
    RooMsgService::instance().getStream(1).removeTopic(Fitting) ;
    RooMsgService::instance().getStream(1).removeTopic(NumIntegration) ;
    RooMsgService::instance().getStream(1).removeTopic(Optimization) ;
    RooMsgService::instance().getStream(1).removeTopic(ObjectHandling) ;
    RooMsgService::instance().getStream(1).removeTopic(Eval) ;
    RooMsgService::instance().Print() ;

	//RooRealVar *mass  = ws->var("mass");
	RooAbsData *data  = ws->data("data");
	RooAbsPdf  *model = ws->pdf("model");
	model->getComponents()->Print();
	
	RooFitResult* r = model->fitTo(*data, RooFit::Extended(kTRUE), RooFit::Hesse(kTRUE)/*, Minos(kTRUE)*/, RooFit::Save(), RooFit::NumCPU(2));
	r->Print("v");
	/*TCanvas *c = new TCanvas("c", "c", 900,700);
	RooPlot* massframe = mass->frame(Title("JpsiKp mass fit"), Bins(100));
	data->plotOn(massframe);
	model->plotOn(massframe);
	model->paramOn(massframe, Layout(0.60,0.97,0.88));
	c->cd();
	massframe->Draw();
	c->Update();
	c->SaveAs("MCfit.png");
	delete c;*/
}
