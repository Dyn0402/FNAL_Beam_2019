/*
 * LED_Stability.cpp
 *
 *  Created on: Apr 9, 2019
 *      Author: dylan
 */

#include <string>
#include <vector>

#include <TCanvas.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH1.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TSpectrum.h>

#include "../FNAL_Beam_Test_Source/config.h"
#include "../FNAL_Beam_Test_Source/Macros.h"

using namespace std;

vector<double> get_percent_list(vector<double> list) {
	vector<double> percent_list;
	double mean = get_mean(list);
	for(auto element:list) {
		percent_list.push_back( (element/mean - 1) * 100);
	}

	return(percent_list);
}

//Convert vector of events to vector of times, multiplying each event number by event_to_time conversion factor.
vector<double> event_to_time(vector<double> events, double event_to_time) {
	vector<double> time;

	for(auto event:events) {
		time.push_back(event*event_to_time);
	}

	return(time);
}


void get_led_stability(TTree *tree, vector<int> events, vector<vector<float>> ped_fits, int n_mv_avg) {
	int label_font = 62;
	//Set branches from tree.
	int event_num;
	float hCal[config::hCal_channels];
	float eCal[config::eCal_channels];
	TBranch *b_event_num = tree->GetBranch("eventno");
	TBranch *b_hCal = tree->GetBranch("HCalRawADC");
	TBranch *b_eCal = tree->GetBranch("ECalRawADC");
	b_event_num->SetAddress(&event_num);
	b_hCal->SetAddress(hCal);
	b_eCal->SetAddress(eCal);

	//Set histograms.
	string name;

	//Set hCal histograms.
	TH1D *hCal_hists[config::hCal_channels];
	for(int i=0; i<config::hCal_channels; i++) {
		name = "hCal LED Channel " + to_string(i);
		hCal_hists[i] = new TH1D(name.data(), name.data(), config::bins, config::low_bin, config::high_bin);
	}
	//Set eCal histograms.
	TH1D *eCal_hists[config::eCal_channels];
	for(int i=0; i<config::eCal_channels; i++) {
		name = "eCal LED Channel " + to_string(i);
		eCal_hists[i] = new TH1D(name.data(), name.data(), config::bins, config::low_bin, config::high_bin);
	}

	//Set vectors.
	vector<double> event_num_leds;
	vector<vector<double>> hCal_leds(config::hCal_channels, vector<double>{});
	vector<vector<double>> eCal_leds(config::eCal_channels, vector<double>{});
	vector<double> time_mv_avg;
	vector<vector<double>> hCal_mv_avg(config::hCal_channels, vector<double>{});
	vector<vector<double>> eCal_mv_avg(config::eCal_channels, vector<double>{});

	//Set first moving average points.
	for(int i=0; i<n_mv_avg; i++) {
		b_hCal->GetEntry(events[i]);
		b_eCal->GetEntry(events[i]);
	}

	//Fill channel histograms.
	for(int event_index:events) {
		b_hCal->GetEntry(event_index);
		b_eCal->GetEntry(event_index);

		//Fill channel histograms and vectors
		for(int i = 0; i<config::hCal_channels; i++) { //hCal
			hCal_hists[i]->Fill(hCal[i] - ped_fits[0][i+config::hCal_start_channel]);
			hCal_leds[i].push_back(hCal[i] - ped_fits[0][i+config::hCal_start_channel]);
		}
		for(int i = 0; i<config::eCal_channels; i++) { //eCal
			eCal_hists[i]->Fill(eCal[i] - ped_fits[0][i+config::eCal_start_channel]);
			eCal_leds[i].push_back(eCal[i] - ped_fits[0][i+config::eCal_start_channel]);
		}
		event_num_leds.push_back(event_index);
	}

	//Convert event numbers to time from start in hours assuming rate of 1Hz.
	vector<double> time_leds = event_to_time(event_num_leds, 1.0/3600.0);

	//Get n point moving averages for event number and hCal/eCal channels.
	time_mv_avg = get_moving_average(time_leds, n_mv_avg);
	for(int i = 0; i<config::hCal_channels; i++) { //hCal
		hCal_mv_avg[i] = get_moving_average(hCal_leds[i], n_mv_avg);
	}
	for(int i = 0; i<config::eCal_channels; i++) { //eCal
		eCal_mv_avg[i] = get_moving_average(eCal_leds[i], n_mv_avg);
	}

	//Make TGraphs
	TGraph *hCal_led_graphs[config::hCal_channels];
	TGraph *hCal_mv_avg_graph[config::hCal_channels];
	for(int i = 0; i<config::hCal_channels; i++) { //hCal
		hCal_led_graphs[i] = new TGraph((int)time_leds.size(), time_leds.data(), hCal_leds[i].data());
		hCal_led_graphs[i]->SetNameTitle(("hCal_led_vs_event channel " + to_string(i)).data());
		hCal_led_graphs[i]->SetTitle(("hCal LED Position vs Time Channel " + to_string(i)).data());
		hCal_led_graphs[i]->GetYaxis()->SetTitle("LED Peak Position (ADC Channel)");
		hCal_led_graphs[i]->GetYaxis()->SetTitleFont(label_font);
		hCal_led_graphs[i]->GetYaxis()->SetLabelFont(label_font);
		hCal_led_graphs[i]->GetXaxis()->SetTitle("Time from Start (hours)");
		hCal_led_graphs[i]->GetXaxis()->SetTitleFont(label_font);
		hCal_led_graphs[i]->GetXaxis()->SetLabelFont(label_font);
		hCal_mv_avg_graph[i] = new TGraph((int)time_mv_avg.size(), time_mv_avg.data(), hCal_mv_avg[i].data());
		hCal_mv_avg_graph[i]->SetNameTitle(("hCal_mv_avg channel " + to_string(i)).data());
		hCal_mv_avg_graph[i]->SetTitle(("hCal LED Position vs Time " + to_string(n_mv_avg) + " Point Moving Average Channel " + to_string(i)).data());
		hCal_mv_avg_graph[i]->SetLineColor(kRed);
		hCal_mv_avg_graph[i]->GetYaxis()->SetTitle("LED Peak Position (ADC Channel)");
		hCal_mv_avg_graph[i]->GetYaxis()->SetTitleFont(label_font);
		hCal_mv_avg_graph[i]->GetYaxis()->SetLabelFont(label_font);
		hCal_mv_avg_graph[i]->GetXaxis()->SetTitle("Time from Start (hours)");
		hCal_mv_avg_graph[i]->GetXaxis()->SetTitleFont(label_font);
		hCal_mv_avg_graph[i]->GetXaxis()->SetLabelFont(label_font);
	}
	TGraph *eCal_led_graphs[config::eCal_channels];
	TGraph *eCal_mv_avg_graph[config::eCal_channels];
	for(int i = 0; i<config::eCal_channels; i++) { //eCal
		eCal_led_graphs[i] = new TGraph((int)time_leds.size(), time_leds.data(), eCal_leds[i].data());
		eCal_led_graphs[i]->SetNameTitle(("eCal_led_vs_event channel " + to_string(i)).data());
		eCal_led_graphs[i]->SetTitle(("eCal LED Position vs Time Channel " + to_string(i)).data());
		eCal_led_graphs[i]->GetYaxis()->SetTitle("LED Peak Position (ADC Channel)");
		eCal_led_graphs[i]->GetYaxis()->SetTitleFont(label_font);
		eCal_led_graphs[i]->GetYaxis()->SetLabelFont(label_font);
		eCal_led_graphs[i]->GetXaxis()->SetTitle("Time from Start (hours)");
		eCal_led_graphs[i]->GetXaxis()->SetTitleFont(label_font);
		eCal_led_graphs[i]->GetXaxis()->SetLabelFont(label_font);
		eCal_mv_avg_graph[i] = new TGraph((int)time_mv_avg.size(), time_mv_avg.data(), eCal_mv_avg[i].data());
		eCal_mv_avg_graph[i]->SetNameTitle(("eCal_mv_avg channel " + to_string(i)).data());
		eCal_mv_avg_graph[i]->SetTitle(("eCal LED Position vs Time " + to_string(n_mv_avg) + " Point Moving Average Channel " + to_string(i)).data());
		eCal_mv_avg_graph[i]->SetLineColor(kRed);
		eCal_mv_avg_graph[i]->GetYaxis()->SetTitle("LED Peak Position (ADC Channel)");
		eCal_mv_avg_graph[i]->GetYaxis()->SetTitleFont(label_font);
		eCal_mv_avg_graph[i]->GetYaxis()->SetLabelFont(label_font);
		eCal_mv_avg_graph[i]->GetXaxis()->SetTitle("Time from Start (hours)");
		eCal_mv_avg_graph[i]->GetXaxis()->SetTitleFont(label_font);
		eCal_mv_avg_graph[i]->GetXaxis()->SetLabelFont(label_font);
	}

	//Make zero TF1
	TF1 *zero = new TF1("zero_eCal", "[0]");
	zero->SetParameter(0, 0);
	zero->SetLineColor(kBlue);
	zero->SetRange(hCal_led_graphs[0]->GetXaxis()->GetXmin(), hCal_led_graphs[0]->GetXaxis()->GetXmax()); //Sketchy definition.

	//Plot and write data to file.
	//HCal
	//Both Scattered data and moving average.
	TCanvas *hCal_can = new TCanvas("hCal_led_canvas", "HCal LED Drift", config::canvas_x, config::canvas_y);
	hCal_can->cd();
	hCal_can->Divide(4,4);

	//Just moving average.
	TCanvas *hCal_mv_avg_can = new TCanvas("hCal_mv_avg_canvas", ("HCal LED Drift " + to_string(n_mv_avg) + " Point Moving Average").data(), config::canvas_x, config::canvas_y);
	hCal_mv_avg_can->cd();
	hCal_mv_avg_can->Divide(4,4);


	TF1 *mean_hCal[config::hCal_channels];
	TGaxis *p_axis_hCal_raw[config::hCal_channels];
	TGaxis *p_axis_hCal[config::hCal_channels];
	double y_min, y_max, x_pos, percent_min, percent_max;
	for(int i=0; i<config::hCal_channels; i++) {
		//Draw raw graphs.
		hCal_can->cd(i+1); //cd to 4x4 pad.
		hCal_led_graphs[i]->Draw("AP"); //Draw raw graph.

		//Make and draw mean TF1
		mean_hCal[i] = new TF1(("mean_hCal" + to_string(i)).data(), "[0]", hCal_led_graphs[i]->GetXaxis()->GetXmin(), hCal_led_graphs[i]->GetXaxis()->GetXmax());
		mean_hCal[i]->SetParameter(0, get_mean(hCal_leds[i]));
		mean_hCal[i]->SetLineColor(kBlue);
		mean_hCal[i]->Draw("same");
		hCal_mv_avg_graph[i]->Draw("lsame");

		//Make and draw percent deviation axis.
		x_pos = hCal_led_graphs[i]->GetXaxis()->GetXmax();
		y_min = hCal_led_graphs[i]->GetYaxis()->GetXmin();
		y_max = hCal_led_graphs[i]->GetYaxis()->GetXmax();
		percent_min = (y_min / mean_hCal[i]->GetParameter(0) - 1) * 100;
		percent_max = (y_max / mean_hCal[i]->GetParameter(0) - 1) * 100;
		p_axis_hCal_raw[i] = new TGaxis(x_pos, y_min, x_pos, y_max, percent_min, percent_max, 510, "W+L", -0.8);
		p_axis_hCal_raw[i]->SetTitle("Percent Deviation from Mean");
		p_axis_hCal_raw[i]->SetTitleFont(label_font);
		p_axis_hCal_raw[i]->SetLabelFont(label_font);
		p_axis_hCal_raw[i]->Draw();

		//Make canvas for surpressed individual points.
		hCal_mv_avg_can->cd(i+1);
		gPad->SetGridx();
		hCal_mv_avg_graph[i]->Draw("AL");
		mean_hCal[i]->Draw("lsame");

		//Make and draw percent deviation axis.
		x_pos = hCal_mv_avg_graph[i]->GetXaxis()->GetXmax();
		y_min = hCal_mv_avg_graph[i]->GetYaxis()->GetXmin();
		y_max = hCal_mv_avg_graph[i]->GetYaxis()->GetXmax();
		percent_min = (y_min / mean_hCal[i]->GetParameter(0) - 1) * 100;
		percent_max = (y_max / mean_hCal[i]->GetParameter(0) - 1) * 100;
		p_axis_hCal[i] = new TGaxis(x_pos, y_min, x_pos, y_max, percent_min, percent_max, 510, "W+L", -0.8);
		p_axis_hCal[i]->SetTitle("Percent Deviation from Mean");
		p_axis_hCal[i]->SetTitleFont(label_font);
		p_axis_hCal[i]->SetLabelFont(label_font);
		p_axis_hCal[i]->Draw();
	}

	//ECal
	TCanvas *eCal_can = new TCanvas("eCal_led_canvas", "ECal LED Drift", config::canvas_x, config::canvas_y);
	eCal_can->cd();
	eCal_can->Divide(4,4);

	//Just moving average.
	TCanvas *eCal_mv_avg_can = new TCanvas("eCal_mv_avg_canvas", ("ECal LED Drift " + to_string(n_mv_avg) + " Point Moving Average").data(), config::canvas_x, config::canvas_y);
	eCal_mv_avg_can->cd();
	eCal_mv_avg_can->Divide(4,4);

	TF1 *mean_eCal[config::eCal_channels];
	TGaxis *p_axis_eCal_raw[config::eCal_channels];
	TGaxis *p_axis_eCal[config::eCal_channels];
	for(int i=0; i<config::eCal_channels; i++) {
		//Draw raw graphs.
		eCal_can->cd(i+1); //cd to 4x4 pad.
		eCal_led_graphs[i]->Draw("AP"); //Draw raw graph.

		//Make and draw mean TF1
		mean_eCal[i] = new TF1(("mean_eCal" + to_string(i)).data(), "[0]", eCal_led_graphs[i]->GetXaxis()->GetXmin(), eCal_led_graphs[i]->GetXaxis()->GetXmax());
		mean_eCal[i]->SetParameter(0, get_mean(eCal_leds[i]));
		mean_eCal[i]->SetLineColor(kBlue);
		mean_eCal[i]->Draw("same");
		eCal_mv_avg_graph[i]->Draw("lsame");

		x_pos = eCal_led_graphs[i]->GetXaxis()->GetXmax();
		y_min = eCal_led_graphs[i]->GetYaxis()->GetXmin();
		y_max = eCal_led_graphs[i]->GetYaxis()->GetXmax();
		percent_min = (y_min / mean_eCal[i]->GetParameter(0) - 1) * 100;
		percent_max = (y_max / mean_eCal[i]->GetParameter(0) - 1) * 100;
		p_axis_eCal_raw[i] = new TGaxis(x_pos, y_min, x_pos, y_max, percent_min, percent_max, 510, "W+L", -0.8);
		p_axis_eCal_raw[i]->SetTitle("Percent Deviation from Mean");
		p_axis_eCal_raw[i]->SetTitleFont(label_font);
		p_axis_eCal_raw[i]->SetLabelFont(label_font);
		p_axis_eCal_raw[i]->Draw();

		//Draw on canvas with raw data point suppressed.
		eCal_mv_avg_can->cd(i+1);
		gPad->SetGridx();
		eCal_mv_avg_graph[i]->Draw("AL");
		mean_eCal[i]->Draw("lsame");

		x_pos = eCal_mv_avg_graph[i]->GetXaxis()->GetXmax();
		y_min = eCal_mv_avg_graph[i]->GetYaxis()->GetXmin();
		y_max = eCal_mv_avg_graph[i]->GetYaxis()->GetXmax();
		percent_min = (y_min / mean_eCal[i]->GetParameter(0) - 1) * 100;
		percent_max = (y_max / mean_eCal[i]->GetParameter(0) - 1) * 100;
		p_axis_eCal[i] = new TGaxis(x_pos, y_min, x_pos, y_max, percent_min, percent_max, 510, "W+L", -0.8);
		p_axis_eCal[i]->SetTitle("Percent Deviation from Mean");
		p_axis_eCal[i]->SetTitleFont(label_font);
		p_axis_eCal[i]->SetLabelFont(label_font);
		p_axis_eCal[i]->Draw();
	}

	//Write objects to file
	for(int i = 0; i<config::hCal_channels; i++) { //hCal
		hCal_hists[i]->Write();
		hCal_led_graphs[i]->Write();
		hCal_mv_avg_graph[i]->Write();
	}
	for(int i = 0; i<config::eCal_channels; i++) { //eCal
		eCal_hists[i]->Write();
		eCal_led_graphs[i]->Write();
		eCal_mv_avg_graph[i]->Write();
	}
	hCal_can->Write();
	hCal_mv_avg_can->Write();
	eCal_can->Write();
	eCal_mv_avg_can->Write();

	//Save canvases as pdfs.
	hCal_can->SaveAs((config::led_stbl_path + "/hCal_raw_" + to_string(n_mv_avg) + ".pdf").data());
	hCal_mv_avg_can->SaveAs((config::led_stbl_path + "/hCal_mvavg_" + to_string(n_mv_avg) + ".pdf").data());
	eCal_can->SaveAs((config::led_stbl_path + "/eCal_raw_" + to_string(n_mv_avg) + ".pdf").data());
	eCal_mv_avg_can->SaveAs((config::led_stbl_path + "/eCal_mvavg_" + to_string(n_mv_avg) + ".pdf").data());

	//Delete heap objects
	for(int i = 0; i<config::hCal_channels; i++) { //hCal
		delete hCal_hists[i];
		delete hCal_led_graphs[i];
		delete hCal_mv_avg_graph[i];
		delete mean_hCal[i];
		delete p_axis_hCal[i];
	}
	for(int i = 0; i<config::eCal_channels; i++) { //eCal
		delete eCal_hists[i];
		delete eCal_led_graphs[i];
		delete eCal_mv_avg_graph[i];
		delete mean_eCal[i];
	}
	delete hCal_can;
	delete hCal_mv_avg_can;
	delete eCal_can;
	delete eCal_mv_avg_can;
}
