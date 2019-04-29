/*
 * filter_noise.cpp
 *
 * Get bad event numbers from 12 junk channels.
 *
 *  Created on: Apr 6, 2019
 *      Author: dylan
 */

#include <string>
#include <vector>

#include <TCanvas.h>
#include <TTree.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TH1.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TLine.h>

#include "../FNAL_Beam_Test_Source/config.h"

using namespace std;

vector<int> filter_noise(TTree *tree) {
	//Set branches from tree.
	float junk[config::junk_channels];
	TBranch *b_event_num = tree->GetBranch("eventno");
	TBranch *b_junk = tree->GetBranch("junk");
	b_junk->SetAddress(junk);

	//Set junk canvas.
	TCanvas *junk_canvas = new TCanvas("junk_Canvas" , "Noise Fits", 2200, 1400);
	junk_canvas->cd();
	junk_canvas->Divide(4,3);


	//Set junk histograms and lines
	TH1D *junk_hists[config::junk_channels];
//	TH1D *junk_noise_hists[config::junk_channels];
//	TH1D *junk_normal_hists[config::junk_channels];
	TLine *low_cuts[config::junk_channels];
	TLine *high_cuts[config::junk_channels];
	string name;
	for(int i=0; i<config::junk_channels; i++) {
		name = config::filter_hist_name + to_string(i);
		junk_hists[i] = new TH1D(name.data(), name.data(), config::bins, config::low_bin, config::high_bin);
//		name = config::filter_hist_name + " Noise " + to_string(i);
//		junk_noise_hists[i] = new TH1D(name.data(), name.data(), config::bins, config::low_bin, config::high_bin);
//		name = config::filter_hist_name + " Normal " + to_string(i);
//		junk_normal_hists[i] = new TH1D(name.data(), name.data(), config::bins, config::low_bin, config::high_bin);
	}

	//Fill junk histograms.
	int event_index = 0;
	while(b_event_num->GetEntry(event_index)) {
		b_junk->GetEntry(event_index);
		for(int i=0; i<config::junk_channels; i++) {
			junk_hists[i]->Fill(junk[i]);
		}
		event_index++;
	}

	//Get ADC channels outside of which events are considered noise.
	vector<vector<double>> junk_cuts;
	vector<double> cut = {0.0, 0.0};
	double mean, rms;
	Double_t *x_peak, *y_peak;
	TF1 *gauss = new TF1(config::filter_fit_name.data(), "gaus");
	for(int i=0; i<config::junk_channels; i++) {
		junk_canvas->cd(i+1);
		gPad->SetLogy();
		TSpectrum *spec = new TSpectrum(2); //Hard-coded, fix
		spec->Search(junk_hists[i], config::filter_spec_sigma, "nobackground", config::filter_spec_thresh);
		x_peak = spec->GetPositionX();
		y_peak = spec->GetPositionY();
		gauss->SetParameters(y_peak[0], x_peak[0], config::filter_spec_sigma);
		gauss->SetRange(x_peak[0]-config::filter_spec_sigma*config::filter_fit_range, x_peak[0]+config::filter_spec_sigma*config::filter_fit_range);
		junk_hists[i]->Fit(config::filter_fit_name.data(), "RQ");
		mean = gauss->GetParameter(1);
		rms = gauss->GetParameter(2);
		cut[0] = mean - rms * config::filter_rms_range;
		cut[1] = mean + rms * config::filter_rms_range;
		low_cuts[i] = new TLine(cut[0], 0, cut[0], gauss->GetParameter(0));
		low_cuts[i]->SetLineColor(0);
		high_cuts[i] = new TLine(cut[1], 0, cut[1], gauss->GetParameter(0));
		high_cuts[i]->SetLineColor(0);
		junk_hists[i]->GetXaxis()->SetRange(0,1000);
		junk_hists[i]->Draw();
		low_cuts[i]->Draw("SAME");
		high_cuts[i]->Draw("SAME");
		junk_cuts.push_back(cut);
		if(config::write_mode) {
			junk_hists[i]->Write();
		}
		delete spec;
	}

	//Write canvas
	if(config::write_mode) {
		junk_canvas->Write();
	}

	//Get number of noise coincidence for each event and determine noise events.
	vector<int> noise_events;
	int coincidence = 0;
	event_index = 0;
	while(b_event_num->GetEntry(event_index)) {
		b_junk->GetEntry(event_index);
		coincidence = 0;
		for(int i=0; i<config::junk_channels; i++) {
			if(junk[i] < junk_cuts[i][0] || junk[i] > junk_cuts[i][1]) {
				coincidence++;
			}
		}
		if(coincidence >= config::filter_coinc_thresh) {
			noise_events.push_back(event_index);
		}
		event_index++;
	}

	//Delete objects off of heap.
	for(int i=0; i<config::junk_channels; i++) {
		delete junk_hists[i];
		delete low_cuts[i];
		delete high_cuts[i];
	}
	delete gauss;
	delete junk_canvas;

	return(noise_events);
}


void plot_noise_vs_run(vector<double> noise_percent, vector<double> run_num) {
	TGraph *noise_vs_run = new TGraph((int)run_num.size(), run_num.data(), noise_percent.data());
	noise_vs_run->SetMarkerStyle(8);
	noise_vs_run->SetName("noise_vs_run");
	noise_vs_run->SetNameTitle("noise_vs_run");
	noise_vs_run->SetTitle("Noise vs Run Number");
	noise_vs_run->GetXaxis()->SetTitle("Run Number");
	noise_vs_run->GetYaxis()->SetTitle("Percentage of Events Flagged as Noise (%)");
	noise_vs_run->Write();
	delete noise_vs_run;
}


void plot_noise_corr(vector<vector<double>> noise_corr, vector<double> run_num) {
	TMultiGraph *multi_graph = new TMultiGraph();
	TMultiGraph *multi_graph11 = new TMultiGraph();
	vector<TGraph*> graphs;
	vector<int> colors = {1, 2, 3, 4, 6, 7, 11, 14, 28, 29, 39};
	reverse(colors.begin(), colors.end()); //Flip color order so that good colors are on 12, 11 etc.
	for(int i = 2; i<=12; i++) {
		graphs.push_back(new TGraph((int)run_num.size(), run_num.data(), noise_corr[i-2].data()));
		graphs.back()->SetMarkerStyle(8);
		graphs.back()->SetMarkerColor(colors[i-2]);
		graphs.back()->SetLineColor(colors[i-2]);
		graphs.back()->SetName((to_string(i) + " Coincidence Graph").data());
		graphs.back()->SetNameTitle((to_string(i) + " Coincidence Graph").data());
		graphs.back()->SetTitle((to_string(i) + " Coincidence Graph").data());
		graphs.back()->GetXaxis()->SetTitle("Run Number");
		graphs.back()->GetYaxis()->SetTitle("Percent More Noise Events for 1 Coincidence");
		multi_graph->Add(graphs.back(), "lp");
		if(i!=12) {multi_graph11->Add(graphs.back(), "lp");}
		graphs.back()->Write();
	}
	multi_graph->SetNameTitle("Coincidence Graphs", "Degree of Noise Correlation vs Run Number");
	multi_graph->GetXaxis()->SetTitle("Run Number");
	multi_graph->GetYaxis()->SetTitle("Percent More Noise Events for 1 Coincidence (%)");
	multi_graph->Write();

	multi_graph11->SetNameTitle("Coincidence Graphs11", "Degree of Noise Correlation (Excluding 12 Coincidence) vs Run Number");
	multi_graph11->GetXaxis()->SetTitle("Run Number");
	multi_graph11->GetYaxis()->SetTitle("Percent More Noise Events for 1 Coincidence (%)");
	multi_graph11->Write();

	TCanvas *multi_can = new TCanvas("multi_can", "Coincidence Graphs", config::canvas_x, config::canvas_y);
	multi_graph->Draw("alp");
	multi_can->BuildLegend();
	multi_can->Write();

	TCanvas *multi_can11 = new TCanvas("multi_can11", "Coincidence Graphs (Excluding 12 Coincidence)", config::canvas_x, config::canvas_y);
	multi_graph11->Draw("alp");
	multi_can11->BuildLegend();
	multi_can11->Write();

	//Delete Objects from heap
//	for(auto graph:graphs) {
//		delete graph;
//	}
	delete multi_graph; //Deletes graphs it owns too?
	delete multi_can;
	delete multi_can11;
}
