/*
 * split_on_sc1.cpp
 *
 *  Created on: Apr 7, 2019
 *      Author: dylan
 */

#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TLine.h>

#include "../FNAL_Beam_Test_Source/config.h"

using namespace std;


vector<double> get_sc1_cuts(TTree *tree, vector<int> noise_events) {
	//Set branches from tree.
	float sc1;
	TBranch *b_event_num = tree->GetBranch("eventno");
	TBranch *b_sc1 = tree->GetBranch("Sc1RawADC");
	b_sc1->SetAddress(&sc1);

	//Set and fill sc1 histogram.
	TH1D *sc1_hist = new TH1D(config::sc1_hist_name.data(), config::sc1_hist_name.data(), config::bins, config::low_bin, config::high_bin);

	int event_index = 0;
	int noise_event_index = 0;
	int next_noise_event = noise_events[noise_event_index];
	while(b_event_num->GetEntry(event_index)) {
		if(event_index == next_noise_event) {
			next_noise_event = noise_events[++noise_event_index];
		} else {
			b_sc1->GetEntry(event_index);
			sc1_hist->Fill(sc1);
		}
		event_index++;
	}

	//Fit sc1 low peak and write fit to root file.
	double mean, rms;
	int num_peaks;
	Double_t *x_peaks, *y_peaks;
	double x_peak, y_peak;
	TF1 *gauss = new TF1(config::sc1_fit_name.data(), "gaus");
	TSpectrum *spec = new TSpectrum(2); //Hard-coded, fix
	num_peaks = spec->Search(sc1_hist, config::sc1_spec_sigma, "nobackground", config::sc1_spec_thresh);
	x_peaks = spec->GetPositionX();
	y_peaks = spec->GetPositionY();
	x_peak = x_peaks[0];
	y_peak = y_peaks[0];
	for(int i=1; i<num_peaks; i++) {
		if(x_peaks[i] < x_peak) {
			x_peak = x_peaks[i];
			y_peak = y_peaks[i];
		}
	}
	gauss->SetParameters(y_peak, x_peak, config::sc1_spec_sigma);
	gauss->SetRange(x_peak-config::sc1_spec_sigma*config::sc1_fit_range, x_peak+config::sc1_spec_sigma*config::sc1_fit_range);
	sc1_hist->Fit(config::sc1_fit_name.data(), "RQ");
	mean = gauss->GetParameter(1);
	rms = gauss->GetParameter(2);

	sc1_hist->Write();


	//Delete objects from heap
	delete sc1_hist;
	delete gauss;
	delete spec;

	//Determine cuts based on fit.
	vector<double> cuts = {0.0,0.0};
	cuts[0] = mean + config::sc1_low_sigmas*rms;
	cuts[1] = mean + config::sc1_high_sigmas*rms;

	return(cuts);
}


vector<int> get_bkg_events_sc1(TTree *tree, vector<int> noise_events) {
	//Set branches from tree.
	float sc1;
	TBranch *b_event_num = tree->GetBranch("eventno");
	TBranch *b_sc1 = tree->GetBranch("Sc1RawADC");
	b_sc1->SetAddress(&sc1);

	//Set and fill sc1 histogram.
	TH1D *sc1_hist = new TH1D(config::sc1_bkg_hist_name.data(), config::sc1_bkg_hist_name.data(), config::bins, config::low_bin, config::high_bin);

	int event_index = 0;
	int noise_event_index = 0;
	int next_noise_event = noise_events[noise_event_index];
	while(b_event_num->GetEntry(event_index)) {
		if(event_index == next_noise_event) {
			next_noise_event = noise_events[++noise_event_index];
		} else {
			b_sc1->GetEntry(event_index);
			sc1_hist->Fill(sc1);
		}
		event_index++;
	}

	//Fit sc1 low peak.
	double mean, rms;
	int num_peaks;
	Double_t *x_peaks, *y_peaks;
	double x_peak, y_peak;
	TF1 *gauss = new TF1(config::sc1_bkg_fit_name.data(), "gaus");
	TSpectrum *spec = new TSpectrum(2); //Hard-coded, fix
	num_peaks = spec->Search(sc1_hist, config::sc1_bkg_spec_sigma, "nobackground", config::sc1_bkg_spec_thresh);
	x_peaks = spec->GetPositionX();
	y_peaks = spec->GetPositionY();
	x_peak = x_peaks[0];
	y_peak = y_peaks[0];
	for(int i=1; i<num_peaks; i++) {
		if(x_peaks[i] < x_peak) {
			x_peak = x_peaks[i];
			y_peak = y_peaks[i];
		}
	}
	gauss->SetParameters(y_peak, x_peak, config::sc1_bkg_spec_sigma);
	gauss->SetRange(x_peak-config::sc1_bkg_spec_sigma*config::sc1_fit_range, x_peak+config::sc1_bkg_spec_sigma*config::sc1_fit_range);
	sc1_hist->Fit(gauss, "RQ");
	mean = gauss->GetParameter(1);
	rms = gauss->GetParameter(2);

	sc1_hist->Write();

	//Determine background cut based on fit.
	double cut = mean + rms*config::sc1_bkg_sigmas;

	//Determine background events
	vector<int> bkg_events;
	event_index = 0;
	noise_event_index = 0;
	next_noise_event = noise_events[noise_event_index];
	while(b_event_num->GetEntry(event_index)) {
		if(event_index == next_noise_event) {
			next_noise_event = noise_events[++noise_event_index];
		} else {
			b_sc1->GetEntry(event_index);
			if(sc1 < cut) {
				bkg_events.push_back(event_index);
			}
		}
		event_index++;
	}


	//Set and save canvas showing fits and cuts.
	TCanvas *sc1_canvas = new TCanvas(config::sc1_bkg_hist_name.data(), config::sc1_bkg_hist_name.data(), config::canvas_x, config::canvas_y);

	TLine *cut_line = new TLine(cut, 0, cut, gauss->GetParameter(0));
	cut_line->SetLineColor(kRed);
	sc1_hist->Draw();
	cut_line->Draw();
	sc1_canvas->Write();


	//Delete objects from heap
	delete sc1_canvas;
	delete cut_line;
	delete sc1_hist;
	delete gauss;
	delete spec;


	return(bkg_events);
}



vector<int> get_beam_events_sc1(TTree *tree, vector<int> noise_events) {
	//Set branches from tree.
	float sc1;
	TBranch *b_event_num = tree->GetBranch("eventno");
	TBranch *b_sc1 = tree->GetBranch("Sc1RawADC");
	b_sc1->SetAddress(&sc1);

	//Set and fill sc1 histogram.
	TH1D *sc1_hist = new TH1D(config::sc1_beam_hist_name.data(), config::sc1_beam_hist_name.data(), config::bins, config::low_bin, config::high_bin);

	int event_index = 0;
	int noise_event_index = 0;
	int next_noise_event = noise_events[noise_event_index];
	while(b_event_num->GetEntry(event_index)) {
		if(event_index == next_noise_event) {
			next_noise_event = noise_events[++noise_event_index];
		} else {
			b_sc1->GetEntry(event_index);
			sc1_hist->Fill(sc1);
		}
		event_index++;
	}

	//Fit sc1 low peak.
	double mean, rms;
	int num_peaks;
	Double_t *x_peaks, *y_peaks;
	double x_peak, y_peak;
	TF1 *gauss = new TF1(config::sc1_beam_fit_name.data(), "gaus");
	TSpectrum *spec = new TSpectrum(2); //Hard-coded, fix
	num_peaks = spec->Search(sc1_hist, config::sc1_bkg_spec_sigma, "nobackground", config::sc1_bkg_spec_thresh);
	x_peaks = spec->GetPositionX();
	y_peaks = spec->GetPositionY();
	x_peak = x_peaks[0];
	y_peak = y_peaks[0];
	for(int i=1; i<num_peaks; i++) {
		if(x_peaks[i] < x_peak) {
			x_peak = x_peaks[i];
			y_peak = y_peaks[i];
		}
	}
	gauss->SetParameters(y_peak, x_peak, config::sc1_spec_sigma);
	gauss->SetRange(x_peak-config::sc1_spec_sigma*config::sc1_fit_range, x_peak+config::sc1_spec_sigma*config::sc1_fit_range);
	sc1_hist->Fit(gauss, "RQ");
	mean = gauss->GetParameter(1);
	rms = gauss->GetParameter(2);

	sc1_hist->Write();

	//Determine beam cut based on fit.
	double cut = mean + rms*config::sc1_beam_sigmas;

	//Determine beam events
	vector<int> beam_events;
	event_index = 0;
	noise_event_index = 0;
	next_noise_event = noise_events[noise_event_index];
	while(b_event_num->GetEntry(event_index)) {
		if(event_index == next_noise_event) {
			next_noise_event = noise_events[++noise_event_index];
		} else {
			b_sc1->GetEntry(event_index);
			if(sc1 > cut) {
				beam_events.push_back(event_index);
			}
		}
		event_index++;
	}


	//Set and save canvas showing fits and cuts.
	TCanvas *sc1_canvas = new TCanvas(config::sc1_beam_hist_name.data(), config::sc1_beam_hist_name.data(), config::canvas_x, config::canvas_y);

	TLine *cut_line = new TLine(cut, 0, cut, gauss->GetParameter(0));
	cut_line->SetLineColor(kRed);
	sc1_hist->Draw();
	cut_line->Draw();
	sc1_canvas->Write();


	//Delete objects from heap
	delete sc1_canvas;
	delete cut_line;
	delete sc1_hist;
	delete gauss;
	delete spec;


	return(beam_events);
}
