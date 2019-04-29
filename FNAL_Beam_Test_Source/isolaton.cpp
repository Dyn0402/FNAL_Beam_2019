/*
 * isolaton.cpp
 *
 *  Created on: Apr 8, 2019
 *      Author: dylan
 */

#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TF1.h"
#include "TSpectrum.h"

#include "../FNAL_Beam_Test_Source/config.h"

using namespace std;


void get_hCal_iso(TTree *tree, vector<int> events, vector<vector<float>> ped_fits) {
	//Set branches from tree.
	float hCal[config::hCal_channels];
	TBranch *b_hCal = tree->GetBranch("HCalRawADC");
	b_hCal->SetAddress(hCal);

	//Set hCal histograms.
	TH1D *hCal_hists[config::hCal_channels];
	TH1D *hCal_iso_hists[config::hCal_channels];
	string name;
	for(int i=0; i<config::hCal_channels; i++) {
		name = config::hCal_iso_hist_name + to_string(i);
		hCal_iso_hists[i] = new TH1D(name.data(), name.data(), config::bins/config::iso_hCal_rebinning, config::low_bin, config::high_bin);
		name = config::hCal_iso_hist_name + "_Raw_" + to_string(i);
		hCal_hists[i] = new TH1D(name.data(), name.data(), config::bins/config::iso_hCal_rebinning, config::low_bin, config::high_bin);
	}

	//Fill channel histograms.
	for(auto event_index:events) {
		b_hCal->GetEntry(event_index);

		//Fill raw hCal histograms and then determine and fill isolated hCal histograms.
		bool adj_sig = false;
		for(int i = 0; i<config::hCal_channels; i++) {
			//Fill raw hCal hists.
			hCal_hists[i]->Fill(hCal[i]-ped_fits[0][i+config::hCal_start_channel]);

			//Check for signal.
			if(hCal[i] > ped_fits[0][i+config::hCal_start_channel] + ped_fits[1][i+config::hCal_start_channel] * config::iso_hCal_signal_sigma) {
				//Signal in hCal[i]. Check adjacent channels.
				for(auto j:config::adj_channels[i]) {
					if(hCal[j] > ped_fits[0][j+config::hCal_start_channel] + ped_fits[1][j+config::hCal_start_channel] * config::iso_hCal_adjsignal_sigma) {
						adj_sig = true;
						break;
					}
				}
				if(!adj_sig) {
					hCal_iso_hists[i]->Fill(hCal[i]-ped_fits[0][i+config::hCal_start_channel]);
				}
				adj_sig = false;
			}
		}
	}

	//Set hCal fits TCanvas.
	TCanvas *hCal_canvas = new TCanvas("hCal_Canvas" , "hCal MIP Peak Fits", config::canvas_x, config::canvas_y);
	hCal_canvas->cd();
	hCal_canvas->Divide(4,4);
	gStyle->SetOptFit(11);
	gStyle->SetOptStat(0000);
	gStyle->SetStatW(0.3);
	gStyle->SetStatH(0.3);


	//Fit hCal histograms.
	TF1 *gauss = new TF1("hCal Iso Gauss Fit", "gaus");
	vector<float> hCal_peaks;
	//Hard coded fit ranges.
//	vector<vector<int>> fit_range = {{600,1250}, {450,1100}, {600,1500}, {500,1150}, {500,1100}, {1100,2000}, {900,1700}, {900, 2000}, {750,1650}, {800,1600}, {400,1000}, {800,1900}, {700, 1450}, {550,1350}, {500, 1400}, {900, 1700}};
	vector<vector<int>> fit_range = {{150,400}, //0
			{180,425}, //1
			{150,370}, //2
			{190,370}, //3
			{70,200}, //4
			{190,415}, //5
			{190,425}, //6
			{225,500}, //7
			{250,500}, //8
			{180,370}, //9
			{220,550}, //10
			{150,450}, //11
			{210,450}, //12
			{200,450}, //13
			{200,500}, //14
			{225,500}};//15
	for(int i=0; i<config::hCal_channels; i++) {
		gauss->SetRange(fit_range[i][0], fit_range[i][1]);
		hCal_canvas->cd(i+1);
		gPad->SetLogy();
		hCal_iso_hists[i]->Fit(gauss, "RQ");
		hCal_iso_hists[i]->Draw();
		hCal_iso_hists[i]->Write();
		hCal_hists[i]->Write();
		hCal_peaks.push_back(gauss->GetParameter(1));
	}
	hCal_canvas->Update();
	hCal_canvas->Write();

	//Set Peak Histogram Canvas
	TCanvas *hCal_peaks_canvas = new TCanvas("hCal_peaks_canvas", "hCal MIP Peak Means", config::canvas_x, config::canvas_y);

	//Set Peak Histogram
	TH1D *hCal_peaks_hist = new TH1D("hCal_peaks", "hCal MIP Peaks", 10, 600, 1500); //Hard coded binning and range.

	for(auto peak:hCal_peaks) {
		hCal_peaks_hist->Fill(peak);
	}

	hCal_peaks_hist->Draw();
	hCal_peaks_canvas->Update();
	hCal_peaks_canvas->Write();
	hCal_peaks_hist->Write();


	//Get mean of peaks.
	double mean = 0;
	for(int i=0; i<config::hCal_channels; i++) {mean+=hCal_peaks[i];}
	mean/=config::hCal_channels;


	//Delete objects from heap
	delete hCal_canvas;
	delete gauss;
	delete hCal_peaks_canvas;
	delete hCal_peaks_hist;
	for(int i=0; i<config::hCal_channels; i++) {
		cout << hCal_peaks[i] << "\t" << hCal_peaks[i] / mean << "\t" << ped_fits[0][i] << endl;
		delete hCal_iso_hists[i];
	}
}




void get_eCal_iso(TTree *tree, vector<int> hadron_events, vector<vector<float>> ped_fits) {
	//Set branches from tree.
	float eCal[config::eCal_channels];
	TBranch *b_event_num = tree->GetBranch("eventno");
	TBranch *b_eCal = tree->GetBranch("ECalRawADC");
	b_eCal->SetAddress(eCal);

	//Set eCal histograms.
	TH1D *eCal_hists[config::eCal_channels];
	string name;
	for(int i=0; i<config::eCal_channels; i++) {
		name = config::eCal_iso_hist_name + to_string(i);
		eCal_hists[i] = new TH1D(name.data(), name.data(), config::bins, config::low_bin, config::high_bin);
	}

	//Fill channel histograms.
	int event_index = 0;
	int hadron_event_index = 0;
	int next_hadron_event = hadron_events[hadron_event_index];
	while(b_event_num->GetEntry(event_index)) {
		//Check filters
		if(event_index == next_hadron_event) {
			next_hadron_event = hadron_events[++hadron_event_index];

			//Get branch entries
			b_eCal->GetEntry(event_index);

			//Fill eCal histograms with isolated signals.
			bool adj_sig = false;
			for(int i = 0; i<config::eCal_channels; i++) {
				//Check for signal.
				if(eCal[i] > ped_fits[0][i+config::eCal_start_channel] + ped_fits[1][i+config::eCal_start_channel] * config::iso_eCal_signal_sigma) {
					//Signal in eCal[i]. Check adjacent channels.
					for(auto j:config::adj_channels[i]) {
						if(eCal[j] > ped_fits[0][j+config::eCal_start_channel] + ped_fits[1][j+config::eCal_start_channel] * config::iso_eCal_adjsignal_sigma + config::iso_eCal_adj_const) {
							adj_sig = true;
							break;
						}
					}
					if(!adj_sig) {
						eCal_hists[i]->Fill(eCal[i]);
					}
					adj_sig = false;
				}
			}
		}
		event_index++;
	}

	//Fit eCal histograms.
	double amp, mean, rms;
	TF1 *gauss = new TF1("eCal Iso Gauss Fit", "gaus");
	TF1 *bkg_gauss = new TF1("eCal Iso Background Gauss fit", "gaus(0) + pol1(3)");
	//Not checked. May or may not work. Fix hard coding and clean up either way.
	for(int i=0; i<config::eCal_channels; i++) {
		TSpectrum *spec = new TSpectrum(2); //Hard-coded, fix
		spec->Search(eCal_hists[i], 30, "nobackground", 0.02);
		Double_t *x_peak = spec->GetPositionX();
		gauss->SetParameters(100, x_peak[0] + 400, 30);
		gauss->SetRange(x_peak[0]+200, x_peak[0]+1000);
		eCal_hists[i]->Fit(gauss, "RQ");
		amp = gauss->GetParameter(0);
		mean = gauss->GetParameter(1);
		rms = gauss->GetParameter(2);
		bkg_gauss->SetParameters(amp, mean, rms, 2, -0.001);
		bkg_gauss->SetRange(mean - 3*rms, mean + 3*rms);
		eCal_hists[i]->Fit(bkg_gauss, "RQ");
		eCal_hists[i]->Write();
		delete spec;
	}

	//Delete objects from heap
	for(int i=0; i<config::eCal_channels; i++) {
		delete eCal_hists[i];
	}
}


void fit_eCal_iso(TFile *file) {
	//Get eCal histograms from file.
	TH1D *eCal_hists[16];
//	for(int i=0; i<16; i++) {
//		TH1D *hist = (TH1D*)file->Get((config::eCal_iso_hist_name + to_string(i)).data());
//		eCal_hists[i] = hist->Clone(("eCal Fit" + to_string(i)).data());
//	}

	//Fit eCal histograms.
	Double_t *x_peak, *y_peak;
	TF1 *gauss = new TF1("eCal Iso Gauss Fit", "gaus");
	for(int i=0; i<config::eCal_channels; i++) {
		TSpectrum *spec = new TSpectrum(4); //Hard-coded, fix
		spec->Search(eCal_hists[i], 30, "nobackground", 0.02);
		x_peak = spec->GetPositionX();
		y_peak = spec->GetPositionY();
		gauss->SetParameters(y_peak[1], x_peak[1], 30);
		gauss->SetRange(x_peak[1]-60, x_peak[1]+60);
		eCal_hists[i]->Fit("eCal Iso Gauss Fit", "RQ"); //Seg violation here.
//		eCal_hists[i]->Write();
		delete spec;
	}
}


void get_hCal_iso2(TTree *tree, vector<int> hadron_events, vector<vector<float>> ped_fits) {
	//Set branches from tree.
	float hCal[config::hCal_channels];
	TBranch *b_event_num = tree->GetBranch("eventno");
	TBranch *b_hCal = tree->GetBranch("HCalRawADC");
	b_hCal->SetAddress(hCal);

	//Set hCal histograms.
	TH1D *hCal_hists[config::hCal_channels];
	TH1D *hCal_iso_hists[config::hCal_channels];
	string name;
	for(int i=0; i<config::hCal_channels; i++) {
		name = config::hCal_iso_hist_name + to_string(i);
		hCal_iso_hists[i] = new TH1D(name.data(), name.data(), config::bins/config::iso_hCal_rebinning, config::low_bin, config::high_bin);
		name = config::hCal_iso_hist_name + "_Raw_" + to_string(i);
		hCal_hists[i] = new TH1D(name.data(), name.data(), config::bins/config::iso_hCal_rebinning, config::low_bin, config::high_bin);
	}

	//Fill channel histograms.
	int event_index = 0;
	int hadron_event_index = 0;
	int next_hadron_event = hadron_events[hadron_event_index];
	while(b_event_num->GetEntry(event_index)) {
		//Check filters
		if(event_index == next_hadron_event) {
			next_hadron_event = hadron_events[++hadron_event_index];

			//Get branch entries
			b_hCal->GetEntry(event_index);

			//Fill raw hCal histograms and then determine and fill isolated hCal histograms.
			bool adj_sig = false;
			for(int i = 0; i<config::hCal_channels; i++) {
				//Fill raw hCal hists.
				hCal_hists[i]->Fill(hCal[i]-ped_fits[0][i+config::hCal_start_channel]);

				//Check for signal.
				if(hCal[i] > ped_fits[0][i+config::hCal_start_channel] + ped_fits[1][i+config::hCal_start_channel] * config::iso_hCal_signal_sigma) {
					//Signal in hCal[i]. Check adjacent channels.
					for(int j=0; j<config::hCal_channels; j++) {
						if(j!=i) {
							if(hCal[j] > ped_fits[0][j+config::hCal_start_channel] + ped_fits[1][j+config::hCal_start_channel] * config::iso_hCal_adjsignal_sigma) {
								adj_sig = true;
								break;
							}
						}
					}
					if(!adj_sig) {
						hCal_iso_hists[i]->Fill(hCal[i]-ped_fits[0][i+config::hCal_start_channel]);
					}
					adj_sig = false;
				}
			}
		}
		event_index++;
	}

	//Set hCal fits TCanvas.
	TCanvas *hCal_canvas = new TCanvas("hCal_Canvas" , "hCal MIP Peak Fits", config::canvas_x, config::canvas_y);
	hCal_canvas->cd();
	hCal_canvas->Divide(4,4);
	gStyle->SetOptFit(11);
	gStyle->SetOptStat(0000);
//	gStyle->SetStatX(0.5);
//	gStyle->SetStatY(0.9);
	gStyle->SetStatW(0.3);
	gStyle->SetStatH(0.3);


	//Fit hCal histograms.
//	double amp, mean, rms;
	TF1 *gauss = new TF1("hCal Iso Gauss Fit", "gaus");
//	TF1 *bkg_gauss = new TF1("hCal Iso Background Gauss Fit", "gaus(0)");
	vector<float> hCal_peaks;
	//Hard coded fit ranges.
	vector<vector<int>> fit_range = {{600,1250}, {450,1100}, {600,1500}, {500,1150}, {500,1100}, {1100,2000}, {900,1700}, {900, 2000}, {750,1650}, {800,1600}, {400,1000}, {800,1900}, {700, 1450}, {550,1350}, {500, 1400}, {900, 1700}};
	for(int i=0; i<config::hCal_channels; i++) {
		gauss->SetRange(fit_range[i][0], fit_range[i][1]);
		hCal_canvas->cd(i+1);
		gPad->SetLogy();
		hCal_iso_hists[i]->Fit(gauss, "RQ");
		hCal_iso_hists[i]->Draw();
		hCal_iso_hists[i]->Write();
		hCal_hists[i]->Write();
		hCal_peaks.push_back(gauss->GetParameter(1));
	}
	hCal_canvas->Update();
	hCal_canvas->Write();

	//Set Peak Histogram Canvas
	TCanvas *hCal_peaks_canvas = new TCanvas("hCal_peaks_canvas", "hCal MIP Peak Means", config::canvas_x, config::canvas_y);

	//Set Peak Histogram
	TH1D *hCal_peaks_hist = new TH1D("hCal_peaks", "hCal MIP Peaks", 10, 600, 1500); //Hard coded binning and range.

	for(auto peak:hCal_peaks) {
		hCal_peaks_hist->Fill(peak);
	}

	hCal_peaks_hist->Draw();
	hCal_peaks_canvas->Update();
	hCal_peaks_canvas->Write();
	hCal_peaks_hist->Write();



	//Delete objects from heap
	delete hCal_canvas;
	delete gauss;
	delete hCal_peaks_canvas;
	delete hCal_peaks_hist;
	for(int i=0; i<config::hCal_channels; i++) {
		cout << hCal_peaks[i] << " | " << ped_fits[0][i] << endl;
		delete hCal_iso_hists[i];
	}
}


void get_hCal_iso_dynamic(TTree *tree, vector<int> events, vector<vector<float>> ped_fits) {
	//Set branches from tree.
	float hCal[config::hCal_channels];
	TBranch *b_hCal = tree->GetBranch("HCalRawADC");
	b_hCal->SetAddress(hCal);

	//Set hCal histograms.
	TH1D *hCal_hists[config::hCal_channels];
	TH1D *hCal_iso_hists[config::hCal_channels];
	string name;
	for(int i=0; i<config::hCal_channels; i++) {
		name = config::hCal_iso_hist_name + to_string(i);
		hCal_iso_hists[i] = new TH1D(name.data(), name.data(), config::bins/config::iso_hCal_rebinning, config::low_bin, config::high_bin);
		name = config::hCal_iso_hist_name + "_Raw_" + to_string(i);
		hCal_hists[i] = new TH1D(name.data(), name.data(), config::bins/config::iso_hCal_rebinning, config::low_bin, config::high_bin);
	}

	//Fill channel histograms.
	for(auto event_index:events) {
		b_hCal->GetEntry(event_index);

		//Fill raw hCal histograms and then determine and fill isolated hCal histograms.
		bool adj_sig = false;
		for(int i = 0; i<config::hCal_channels; i++) {
			//Fill raw hCal hists.
			hCal_hists[i]->Fill(hCal[i]-ped_fits[0][i+config::hCal_start_channel]);

			//Check for signal.
			if(hCal[i] > ped_fits[0][i+config::hCal_start_channel] + ped_fits[1][i+config::hCal_start_channel] * config::iso_hCal_signal_sigma) {
				//Signal in hCal[i]. Check adjacent channels.
				for(auto j:config::adj_channels[i]) {
					if(hCal[j] > ped_fits[0][j+config::hCal_start_channel] + ped_fits[1][j+config::hCal_start_channel] * config::iso_hCal_adjsignal_sigma) {
						adj_sig = true;
						break;
					}
				}
				if(!adj_sig) {
					hCal_iso_hists[i]->Fill(hCal[i]-ped_fits[0][i+config::hCal_start_channel]);
				}
				adj_sig = false;
			}
		}
	}

	//Set hCal fits TCanvas.
	TCanvas *hCal_canvas = new TCanvas("hCal_Canvas" , "hCal MIP Peak Fits", config::canvas_x, config::canvas_y);
	hCal_canvas->cd();
	hCal_canvas->Divide(4,4);
	gStyle->SetOptFit(11);
	gStyle->SetOptStat(0000);
	gStyle->SetStatW(0.3);
	gStyle->SetStatH(0.3);
	hCal_canvas->Update();


	//Fit hCal histograms.
	TF1 *gauss = new TF1("hCal Iso Gauss Fit", "gaus");
	vector<float> hCal_peaks;
	TSpectrum *specs[config::hCal_channels];
	Double_t *x_peaks, *y_peaks;
	double x_peak, y_peak;
	int num_peaks;
	for(int i=0; i<config::hCal_channels; i++) {
		specs[i] = new TSpectrum(4); //Fix hardcode.
		num_peaks = specs[i]->Search(hCal_iso_hists[i], 10.0, "nobackground", 0.01); //Fix hardcode.
		x_peaks = specs[i]->GetPositionX();
		y_peaks = specs[i]->GetPositionY();
		x_peak = x_peaks[0];
		y_peak = y_peaks[0];
		for(int i=1; i<num_peaks; i++) {
			if(x_peaks[i] > x_peak) {
				x_peak = x_peaks[i];
				y_peak = y_peaks[i];
			}
		}
//		cout << "Channel " << i << " num_peaks=" << num_peaks << " x_peak=" << x_peak << endl;
		gauss->SetRange(x_peak * 0.8, x_peak * 1.2);
		hCal_canvas->cd(i+1);
		gPad->SetLogy();
		gauss->SetParameters(y_peak, x_peak, x_peak * 0.2); //Hard-code fix.
		hCal_iso_hists[i]->Fit(gauss, "RQ");
		hCal_iso_hists[i]->Draw();
		hCal_canvas->Update();
//		specs[i]->Draw("same");
		hCal_peaks.push_back(gauss->GetParameter(1));
	}

	//Set Peak Histogram Canvas
	TCanvas *hCal_peaks_canvas = new TCanvas("hCal_peaks_canvas", "hCal MIP Peak Means", config::canvas_x, config::canvas_y);

	//Set Peak Histogram
	TH1D *hCal_peaks_hist = new TH1D("hCal_peaks", "hCal MIP Peaks", 10, 600, 1500); //Hard coded binning and range.

	for(auto peak:hCal_peaks) {
		hCal_peaks_hist->Fill(peak);
	}
	hCal_peaks_hist->Draw();


	//Write objects to file.
	for(int i=0; i<config::hCal_channels; i++) {
		hCal_iso_hists[i]->Write();
		hCal_hists[i]->Write();
	}
	hCal_canvas->Update();
	hCal_canvas->Write();
	hCal_peaks_canvas->Update();
	hCal_peaks_canvas->Write();
	hCal_peaks_hist->Write();


	//Get mean of peaks.
	double mean = 0;
	for(int i=0; i<config::hCal_channels; i++) {mean+=hCal_peaks[i];}
	mean/=config::hCal_channels;


	//Delete objects from heap
	delete hCal_canvas;
	delete gauss;
	delete hCal_peaks_canvas;
	delete hCal_peaks_hist;
	for(int i=0; i<config::hCal_channels; i++) {
		cout << hCal_peaks[i] << "\t" << hCal_peaks[i] / mean << "\t" << ped_fits[0][i] << endl;
		delete hCal_iso_hists[i];
		delete specs[i];
	}
}



void check_hcalibration(TTree *tree, vector<int> events, vector<vector<float>> ped_fits) {
	//Set branches from tree.
	float hCal[config::hCal_channels];
	TBranch *b_hCal = tree->GetBranch("HCalRawADC");
	b_hCal->SetAddress(hCal);

	//Set hCal histograms.
	TH1D *hCal_iso_hists[config::hCal_channels];
	string name;
	for(int i=0; i<config::hCal_channels; i++) {
		name = config::hCal_iso_hist_name + to_string(i);
		hCal_iso_hists[i] = new TH1D(name.data(), name.data(), config::bins/config::iso_hCal_rebinning, config::low_bin, config::high_bin);
	}

	//Fill channel histograms.
	for(auto event_index:events) {
		b_hCal->GetEntry(event_index);

		//Fill raw hCal histograms and then determine and fill isolated hCal histograms.
		bool adj_sig = false;
		for(int i = 0; i<config::hCal_channels; i++) {
			//Check for signal.
			if(hCal[i] > ped_fits[0][i+config::hCal_start_channel] + ped_fits[1][i+config::hCal_start_channel] * config::iso_hCal_signal_sigma) {
				//Signal in hCal[i]. Check adjacent channels.
				for(auto j:config::adj_channels[i]) {
					if(hCal[j] > ped_fits[0][j+config::hCal_start_channel] + ped_fits[1][j+config::hCal_start_channel] * config::iso_hCal_adjsignal_sigma) {
						adj_sig = true;
						break;
					}
				}
				if(!adj_sig) {
					hCal_iso_hists[i]->Fill((hCal[i]-ped_fits[0][i+config::hCal_start_channel])/config::hCal_weights[i]);
				}
				adj_sig = false;
			}
		}
	}

	//Set hCal fits TCanvas.
	TCanvas *hCal_canvas = new TCanvas("hCal_corrected_Canvas" , "hCal Corrected MIP Peak Fits", config::canvas_x, config::canvas_y);
	hCal_canvas->cd();
	hCal_canvas->Divide(4,4);
	gStyle->SetOptFit(11);
	gStyle->SetOptStat(0000);
	gStyle->SetStatW(0.3);
	gStyle->SetStatH(0.3);
	hCal_canvas->Update();


	//Fit hCal histograms.
	TF1 *gauss = new TF1("hCal Corrected Iso Gauss Fit", "gaus");
	vector<float> hCal_peaks;
	//Hard coded fit ranges.
//	vector<vector<int>> fit_range = {{600,1250}, {450,1100}, {600,1500}, {500,1150}, {500,1100}, {1100,2000}, {900,1700}, {900, 2000}, {750,1650}, {800,1600}, {400,1000}, {800,1900}, {700, 1450}, {550,1350}, {500, 1400}, {900, 1700}};
	vector<int> fit_range = {150,500};
	for(int i=0; i<config::hCal_channels; i++) {
			gauss->SetRange(fit_range[0], fit_range[1]);
			hCal_canvas->cd(i+1);
			gPad->SetLogy();
			hCal_iso_hists[i]->Fit(gauss, "RQ");
			hCal_iso_hists[i]->Draw();
			hCal_iso_hists[i]->Write();
			hCal_peaks.push_back(gauss->GetParameter(1));
		}

	//Set Peak Histogram Canvas
	TCanvas *hCal_peaks_canvas = new TCanvas("hCal_corrected_peaks_canvas", "hCal Corrected MIP Peak Means", config::canvas_x, config::canvas_y);

	//Set Peak Histogram
	TH1D *hCal_peaks_hist = new TH1D("hCal_corrected_peaks", "hCal Corrected MIP Peaks", 10, 300, 320); //Hard coded binning and range.

	for(auto peak:hCal_peaks) {
		hCal_peaks_hist->Fill(peak);
	}
	hCal_peaks_hist->Draw();


	//Write objects to file.
	for(int i=0; i<config::hCal_channels; i++) {
		hCal_iso_hists[i]->Write();
	}
	hCal_canvas->Update();
	hCal_canvas->Write();
	hCal_peaks_canvas->Update();
	hCal_peaks_canvas->Write();
	hCal_peaks_hist->Write();


	//Get mean of peaks.
	double mean = 0;
	for(int i=0; i<config::hCal_channels; i++) {mean+=hCal_peaks[i];}
	mean/=config::hCal_channels;


	//Delete objects from heap
	delete hCal_canvas;
	delete gauss;
	delete hCal_peaks_canvas;
	delete hCal_peaks_hist;
	for(int i=0; i<config::hCal_channels; i++) {
		cout << hCal_peaks[i] << "\t" << hCal_peaks[i] / mean << "\t" << ped_fits[0][i] << endl;
		delete hCal_iso_hists[i];
	}
}

