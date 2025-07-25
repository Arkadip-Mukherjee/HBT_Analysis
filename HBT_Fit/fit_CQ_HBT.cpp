
// TO RUN THIS PROGRAM, USE : g++ -o run fit_CQ_HBT.cpp `root-config --cflags --glibs`

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TH3F.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF3.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TMath.h>
#include <TLegend.h>
#include <TString.h>
#include <TError.h>
#include <Math/MinimizerOptions.h>
#include <TROOT.h>
#include <TVirtualFitter.h>

using namespace std;

// Alter these parameters as per your requirements
const int q_num_bins = 31;
const double qmin = -0.05;
const double qmax = 0.05;
const int kT_num_bins = 4;
const double kT_min = 0.15;
const double kT_max = 0.55;
double kT_bin_width = (kT_max - kT_min) / kT_num_bins ;

// Initial guesses for HBT Radii
const double lambda_ini = 1.0;
const double Rout_ini = 4.0;
const double Rside_ini = 4.0;
const double Rlong_ini = 5.0;

// Limits for fitting HBT radii
const double Rout_min = 2.0; 
const double Rout_max = 7.0; 
const double Rside_min = 2.0; 
const double Rside_max = 7.0; 
const double Rlong_min = 4.0; 
const double Rlong_max = 14.0; 

const double hbar_c = 0.197326980;

// Defining the 3D model function: C(q_out, q_side, q_long)
Double_t model_func_3D(Double_t *x, Double_t *par) {
	Double_t qoutsq = x[0]*x[0];
	Double_t qsidesq = x[1]*x[1];
	Double_t qlongsq = x[2]*x[2];
    
	Double_t lambda = TMath::Abs(par[1]);
  // Here, I used a complete ellipsoid. If you need only the diagonal terms (Rout, Rside & Rlong only), remove the three cross terms from expoterm below 
	Double_t expoterm = TMath::Exp((-par[2]*qoutsq -par[3]*qsidesq -par[4]*qlongsq
								 -2*par[5]*x[0]*x[1] -2*par[6]*x[1]*x[2] -2*par[7]*x[0]*x[2]) / (hbar_c*hbar_c));
    
	return (par[0]*(1 + lambda*expoterm));
}

// Function to read data from file
void extract_data(const std::string& filename, double* q_out_vals, double* q_side_vals, 
		  double* q_long_vals, double* C_vals, double* err_vals) {
	std::ifstream file(filename);
	if (!file) {
	        std::cerr << "Error: Unable to open file " << filename << std::endl;
	        return;
	}
	double qo, qs, ql, Cq, err;
	int ii = 0;
	std::string line;
	while (std::getline(file, line)){
		file >> qo >> qs >> ql >> Cq >> err;
		q_out_vals[ii]  = qo;
		q_side_vals[ii] = qs;
		q_long_vals[ii] = ql;
		C_vals[ii] = Cq;
    err_vals[ii] = err;
	} 
	file.close();
}

int main() {

	gErrorIgnoreLevel = kWarning;

//	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Simplex");
//	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
	ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-6);
//	ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(3);

	std::vector<double>* Rout_storeval = new std::vector<double>;
	std::vector<double>* Rside_storeval = new std::vector<double>;
	std::vector<double>* Rlong_storeval = new std::vector<double>;
	std::vector<double>* Ros_storeval = new std::vector<double>;
	std::vector<double>* Rsl_storeval = new std::vector<double>;
	std::vector<double>* Rol_storeval = new std::vector<double>;
	
	double kT_bins[kT_num_bins];
	for (int ii = 0; ii < kT_num_bins; ii++){
		kT_bins[ii] = kT_min + (ii + 0.5)*kT_bin_width;
	} 

	std::ofstream OutputFile ("HBT_Radii_Output.dat");
	
	for (int kTbin = 0; kTbin < kT_num_bins; kTbin++){

		// Create a 3D histogram to hold combined data
		TH3D* h3 = new TH3D("h3_bin", "Bin", q_num_bins, qmin, qmax, q_num_bins, qmin, qmax, q_num_bins, qmin, qmax);
		h3->Reset("ICE");
		h3->GetXaxis()->SetTitle("q_{out} [GeV/c]");
		h3->GetYaxis()->SetTitle("q_{side} [GeV/c]");
		h3->GetZaxis()->SetTitle("q_{long} [GeV/c]");

		std::string filename = "Cq3D_kT_" + std::to_string(kTbin) + ".data";
		std::ifstream file2(filename);
		double qout_val, qside_val, qlong_val, Cq_val, err_val;
		
		std::vector<double>* qout_store = new std::vector<double>;
		std::vector<double>* qside_store = new std::vector<double>;
		std::vector<double>* qlong_store = new std::vector<double>;
		std::vector<double>* Cq_store = new std::vector<double>;

		while (file2 >> qout_val >> qside_val >> qlong_val >> Cq_val >> err_val) {
		
			qout_store->push_back(qout_val);
			qside_store->push_back(qside_val);
			qlong_store->push_back(qlong_val);
			Cq_store->push_back(Cq_val);
			
	   	h3->Fill(qout_val, qside_val, qlong_val, Cq_val);
    	h3->SetBinError(h3->FindBin(qout_val, qside_val, qlong_val), err_val);
		}
		
		file2.close();
		

		// Define the 3D fitting function
		TF3* fitFunc = new TF3("fitFunc", model_func_3D, qmin, qmax, qmin, qmax, qmin, qmax, 8);

		fitFunc->SetRange(qmin, qmin, qmin, qmax, qmax, qmax);

//		fitFunc->SetParLimits(1, 0.0, 1.0);
//		fitFunc->SetParLimits(2, Rout_min, Rout_max); 
//		fitFunc->SetParLimits(3, Rside_min, Rside_max); 
//		fitFunc->SetParLimits(4, Rlong_min, Rlong_max); 

		fitFunc->FixParameter(0, 1.0);
		fitFunc->SetParameter(1, lambda_ini);   
		fitFunc->SetParameter(2, Rout_ini);   
		fitFunc->SetParameter(3, Rside_ini);   
		fitFunc->SetParameter(4, Rlong_ini);   

		// Perform the fit
		h3->Fit(fitFunc, "RB");
		
		TVirtualFitter *fitter = TVirtualFitter::GetFitter();
		Double_t chi2 = fitFunc->GetChisquare();
		Int_t ndf = fitFunc->GetNDF();
		Double_t chi2_ndf = chi2/ndf;


		Rout_storeval->push_back(fitFunc->GetParameter(2));
		Rside_storeval->push_back(fitFunc->GetParameter(3));
		Rlong_storeval->push_back(fitFunc->GetParameter(4));
		Ros_storeval->push_back(fitFunc->GetParameter(5));
		Rsl_storeval->push_back(fitFunc->GetParameter(6));
		Rol_storeval->push_back(fitFunc->GetParameter(7));


		OutputFile << std::scientific << kT_bins[kTbin] << "  " << chi2_ndf << "  "
									  << fitFunc->GetParameter(2) << "  " << fitFunc->GetParError(2) << "  "
									  << fitFunc->GetParameter(3) << "  " << fitFunc->GetParError(3) << "  "
									  << fitFunc->GetParameter(4) << "  " << fitFunc->GetParError(4) << "  ";
									  << (fitFunc->GetParameter(2) / fitFunc->GetParameter(3)) << "\n";
									  << fitFunc->GetParameter(5) << "  " << fitFunc->GetParError(5) << "  "
									  << fitFunc->GetParameter(6) << "  " << fitFunc->GetParError(6) << "  "
									  << fitFunc->GetParameter(7) << "  " << fitFunc->GetParError(7) << "\n";
		
		delete fitFunc;
		delete h3;
	}
	
	OutputFile.close();
	
	TCanvas *c1 = new TCanvas ("c1", "Radii_vs_kT_plots", 1200, 1200);
	c1->Divide(3, 2, 0.00001, 0.00001);
	
	for (int i = 1; i <= 6; i++) {
		c1->cd(i);
		gPad->SetMargin(0.2, 0.03, 0.2, 0.08);  // (left, right, bottom, top)
	}
	
	c1->cd(1);
	TGraph *gr1 = new TGraph(kT_num_bins, kT_bins, Rout_storeval->data());
	gr1->SetTitle("");
	gr1->GetXaxis()->SetTitle("k_{T} [GeV/c]");
	gr1->GetYaxis()->SetTitle("R_{out} [fm]");
	gr1->GetXaxis()->SetTitleSize(0.08);  
	gr1->GetYaxis()->SetTitleSize(0.08);  
	gr1->GetXaxis()->SetLabelSize(0.04); 
	gr1->GetYaxis()->SetLabelSize(0.04);
	gr1->SetMarkerColor(kRed);
	gr1->SetMarkerStyle(0);
	gr1->SetLineColor(kRed);
	gr1->SetLineWidth(3);
	gr1->GetXaxis()->SetLimits(kT_min, kT_max);
//	gr1->SetMinimum(3.5);
//	gr1->SetMaximum(8.5);
	gr1->Draw("APC");
	
	c1->cd(2);
	TGraph *gr2 = new TGraph(kT_num_bins, kT_bins, Rside_storeval->data());
	gr2->SetTitle("");
	gr2->GetXaxis()->SetTitle("k_{T} [GeV/c]");
	gr2->GetYaxis()->SetTitle("R_{side} [fm]");
	gr2->GetXaxis()->SetTitleSize(0.08);  
	gr2->GetYaxis()->SetTitleSize(0.08);  
	gr2->GetXaxis()->SetLabelSize(0.04); 
	gr2->GetYaxis()->SetLabelSize(0.04);
	gr2->SetMarkerColor(kBlue);
	gr2->SetMarkerStyle(0);
	gr2->SetLineColor(kBlue);
	gr2->SetLineWidth(3);
	gr2->GetXaxis()->SetLimits(kT_min, kT_max);
//	gr2->SetMinimum(3.5);
//	gr2->SetMaximum(8.5);
	gr2->Draw("APC");
	
	c1->cd(3);
	TGraph *gr3 = new TGraph(kT_num_bins, kT_bins, Rlong_storeval->data());
	gr3->SetTitle("");
	gr3->GetXaxis()->SetTitle("k_{T} [GeV/c]");
	gr3->GetYaxis()->SetTitle("R_{long} [fm]");
	gr3->GetXaxis()->SetTitleSize(0.08);  
	gr3->GetYaxis()->SetTitleSize(0.08);  
	gr3->GetXaxis()->SetLabelSize(0.04); 
	gr3->GetYaxis()->SetLabelSize(0.04);
	gr3->SetMarkerColor(kGreen);
	gr3->SetMarkerStyle(0);
	gr3->SetLineColor(kGreen);
	gr3->SetLineWidth(3);
	gr3->GetXaxis()->SetLimits(kT_min, kT_max);
//	gr3->SetMinimum(3.5);
//	gr3->SetMaximum(8.5);
	gr3->Draw("APC");
	
	c1->cd(4);
	TGraph *gr4 = new TGraph(kT_num_bins, kT_bins, Ros_storeval->data());
	gr4->SetTitle("");
	gr4->GetXaxis()->SetTitle("k_{T} [GeV/c]");
	gr4->GetYaxis()->SetTitle("R_{os} [fm]");
	gr4->GetXaxis()->SetTitleSize(0.08);  
	gr4->GetYaxis()->SetTitleSize(0.08);  
	gr4->GetXaxis()->SetLabelSize(0.04); 
	gr4->GetYaxis()->SetLabelSize(0.04);
	gr4->SetMarkerColor(kRed);
	gr4->SetMarkerStyle(0);
	gr4->SetLineColor(kRed);
	gr4->SetLineWidth(3);
	gr4->GetXaxis()->SetLimits(kT_min, kT_max);
//	gr4->SetMinimum(3.5);
//	gr4->SetMaximum(8.5);
	gr4->Draw("APC");
	
	c1->cd(5);
	TGraph *gr5 = new TGraph(kT_num_bins, kT_bins, Rsl_storeval->data());
	gr5->SetTitle("");
	gr5->GetXaxis()->SetTitle("k_{T} [GeV/c]");
	gr5->GetYaxis()->SetTitle("R_{sl} [fm]");
	gr5->GetXaxis()->SetTitleSize(0.08);  
	gr5->GetYaxis()->SetTitleSize(0.08);  
	gr5->GetXaxis()->SetLabelSize(0.04); 
	gr5->GetYaxis()->SetLabelSize(0.04);
	gr5->SetMarkerColor(kRed);
	gr5->SetMarkerStyle(0);
	gr5->SetLineColor(kRed);
	gr5->SetLineWidth(3);
	gr5->GetXaxis()->SetLimits(kT_min, kT_max);
//	gr5->SetMinimum(3.5);
//	gr5->SetMaximum(8.5);
	gr5->Draw("APC");
	
	c1->cd(6);
	TGraph *gr6 = new TGraph(kT_num_bins, kT_bins, Rol_storeval->data());
	gr6->SetTitle("");
	gr6->GetXaxis()->SetTitle("k_{T} [GeV/c]");
	gr6->GetYaxis()->SetTitle("R_{ol} [fm]");
	gr6->GetXaxis()->SetTitleSize(0.08);  
	gr6->GetYaxis()->SetTitleSize(0.08);  
	gr6->GetXaxis()->SetLabelSize(0.04); 
	gr6->GetYaxis()->SetLabelSize(0.04);
	gr6->SetMarkerColor(kRed);
	gr6->SetMarkerStyle(0);
	gr6->SetLineColor(kRed);
	gr6->SetLineWidth(3);
	gr6->GetXaxis()->SetLimits(kT_min, kT_max);
//	gr6->SetMinimum(3.5);
//	gr6->SetMaximum(8.5);
	gr6->Draw("APC");
	
	c1->SaveAs("_HBT_Plot_.pdf");
	
	delete gr1;
	delete gr2;
	delete gr3;
	delete gr4;
	delete gr5;
	delete gr6;
	delete c1;
	
	delete Rout_storeval;
	delete Rside_storeval;
	delete Rlong_storeval;
	delete Ros_storeval;
	delete Rsl_storeval;
	delete Rol_storeval;

	return 0;
}
