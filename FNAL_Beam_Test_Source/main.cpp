/*
 * main.cpp
 *
 *  Created on: Mar 11, 2019
 *      Author: dylan
 */


#include <iostream>

#include <TROOT.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TObjectTable.h>
#include <TGraphErrors.h>

#include "../FNAL_Beam_Test_Source/config.h"
#include "../FNAL_Beam_Test_Source/Macros.h"
#include "../FNAL_Beam_Test_Source/Plots.h"

using namespace std;

void Make_Root();
void Analysis();
void Plots();
void Noise_Events();
void LED_Stability(string run, int n_mv_avg);
void sc1_events();
void HCalibration();
void HCalibration20();
void HCalibration62();
void HCalibration130();
void ECalAnalysis();
void Correlate_Noise();
void Get_HECal_Scale136(string run);
void Get_HECal_Scale76();
void Get_Hadron_Res(string run);

int main() {
//	ECalAnalysis();
//	vector<string> runs = {"run136", "run137", "run138", "run139", "run140", "run141", "run142", "run143", "run144", "run145", "run146"};
//	vector<string> runs = {"run146"};
//	for(auto run:runs) {
//		Get_HECal_Scale136(run);
//		Get_Hadron_Res(run);
//	}
//	Get_HECal_Scale136();
//	Get_HECal_Scale76();
//	HCalibration130();
//	Correlate_Noise();
//	LED_Stability("run25", 600);
//	Make_Root();
	Noise_Events();

	cout << "donzo" << endl;

	return(0);
}

void Make_Root() {
	for(auto run:config::readout_runs) {
		string path_in = config::beam_path_raw + run + ".txt";
		string path_out = config::beam_path_root + run + ".root";
		cout << path_in << endl;
		ReadoutHistograms(path_in, path_out);
	}
}

//Need to rework
void Analysis() {
	TFile *f = new TFile((config::beam_path_root + config::file_name + ".root").data());
	TTree *tree = (TTree *)f->Get("T");
	string tfile_name = config::led_stbl_path + config::file_name + "_out.root";
	TFile *file = new TFile(tfile_name.data(), "RECREATE");
	vector<int> noise_events = filter_noise(tree);
	cout << noise_events.size() << endl;
	vector<double> sc1_cuts = get_sc1_cuts(tree, noise_events);
	cout << sc1_cuts[0] << " | " << sc1_cuts[1] << endl;
	vector<int> led_events = get_led_events(tree, noise_events);
	vector<int> bkg_events = get_bkg_events_sc1(tree, noise_events);
	cout << led_events.size() << endl;
	vector<vector<float>> ped_fits = get_pedestals(tree, led_events);
	get_hCal_iso(tree, led_events, ped_fits);
//	get_eCal_iso(tree, noise_events, led_events, ped_fits);
//	fit_eCal_iso(file);


	file->Close();
}


void HCalibration() {
	TFile *f = new TFile((config::beam_path_root + config::file_name + ".root").data());
	TTree *tree = (TTree *)f->Get("T");
	string tfile_name = config::hCalibration_path + config::file_name + "_hCalibration_out.root";
	TFile *file = new TFile(tfile_name.data(), "RECREATE");
	vector<int> noise_events = filter_noise(tree);
	cout << "Noise Events: " << noise_events.size() << endl;
	vector<int> bkg_events = get_bkg_events_sc1(tree, noise_events);
	cout << "Background Events: " << bkg_events.size() << endl;
	vector<int> beam_events = get_beam_events_sc1(tree, noise_events);
	cout << "Beam Events: " << beam_events.size() << endl;
	vector<vector<float>> ped_fits = get_pedestals(tree, bkg_events);
	vector<int> hadron_events = get_hadron_events(tree, beam_events, ped_fits);
	cout << "Hadron Events: " << hadron_events.size() << endl;
	get_hCal_iso(tree, hadron_events, ped_fits);
//	get_eCal_iso(tree, noise_events, led_events, ped_fits);
//	fit_eCal_iso(file);


	file->Close();
}

void HCalibration20() {
	TFile *f = new TFile((config::beam_path_root + "run20_21" + ".root").data());
	TTree *tree = (TTree *)f->Get("T");
	string tfile_name = config::hCalibration_path + "run20_21" + "_plot_test_out.root";
	TFile *file = new TFile(tfile_name.data(), "RECREATE");
	int num_events = get_num_events(tree);
	cout << "Total Events: " << num_events << endl;
	vector<int> noise_events = filter_noise(tree);
	cout << "Noise Events: " << noise_events.size() << endl;
	vector<int> clean_events = get_events_complement(noise_events, num_events);
	cout << "Clean Events: " << clean_events.size() << endl;
	vector<int> led_events = get_led_events(tree, clean_events);
	cout << "LED Events: " << led_events.size() << endl;
	vector<int> non_led_events = get_events_complement(led_events, num_events);
	cout << "Non-LED Events: " << non_led_events.size() << endl;
	vector<int> non_beam_events = get_cut(tree, clean_events, "Sc1RawADC", 0, 360);
	cout << "Non-Beam Events: " << non_beam_events.size() << endl;
	vector<int> beam_events = get_cut(tree, clean_events, "Sc1RawADC", 367, 4096);
	cout << "Beam Events: " << beam_events.size() << endl;
	vector<vector<float>> ped_fits = get_pedestals(tree, non_led_events);
	get_hCal_iso(tree, non_led_events, ped_fits);
//	get_eCal_iso(tree, noise_events, led_events, ped_fits);
//	fit_eCal_iso(file);


	file->Close();
}


void HCalibration130() {
	TFile *f = new TFile((config::beam_path_root + "run134_135" + ".root").data());
	TTree *tree = (TTree *)f->Get("T");
	string tfile_name = config::hCalibration_path + "run134_135" + "_hCalibration.root";
	TFile *file = new TFile(tfile_name.data(), "RECREATE");
	int num_events = get_num_events(tree);
	cout << "Total Events: " << num_events << endl;
	vector<int> noise_events = filter_noise(tree);
	cout << "Noise Events: " << noise_events.size() << endl;
	vector<int> clean_events = get_events_complement(noise_events, num_events);
	cout << "Clean Events: " << clean_events.size() << endl;
	vector<int> non_beam_events = get_cut(tree, clean_events, "Sc2RawADC", 0, 435);
	cout << "Non-beam Events: " << non_beam_events.size() << endl;
	vector<int> beam_events = get_cut(tree, clean_events, "Sc2RawADC", 435, 4096);
	cout << "Beam Events: " << beam_events.size() << endl;
	vector<int> led_events = get_cut_cal(tree, non_beam_events, "ECalRawADC", 2, 1700, 4096);
	cout << "LED Events: " << led_events.size() << endl;
	vector<int> non_led_events = get_events_complement(led_events, num_events);
	cout << "Non-LED Events: " << non_led_events.size() << endl;
	vector<int> ped_events = get_overlap_events(non_led_events, non_beam_events, num_events);
	cout << "Pedestal Events: " << ped_events.size() << endl;
	vector<vector<float>> ped_fits = get_pedestals(tree, ped_events);
	vector<int> clean_beam_events = get_overlap_events(non_led_events, beam_events, num_events);
	cout << "Clean Beam Events: " << clean_beam_events.size() << endl;
	get_hCal_iso(tree, clean_beam_events, ped_fits);
//	check_hcalibration(tree, clean_beam_events, ped_fits);

	file->Close();
}


void HCalibration62() {
	TFile *f = new TFile((config::beam_path_root + config::file_name + ".root").data());
	TTree *tree = (TTree *)f->Get("T");
	string tfile_name = config::hCalibration_path + config::file_name + "_hCalibration_out.root";
	TFile *file = new TFile(tfile_name.data(), "RECREATE");
	int num_events = get_num_events(tree);
	cout << "Total Events: " << num_events << endl;
	vector<int> noise_events = filter_noise(tree);
	cout << "Noise Events: " << noise_events.size() << endl;
	vector<int> clean_events = get_events_complement(noise_events, num_events);
	cout << "Clean Events: " << clean_events.size() << endl;
	vector<int> bkg_events = get_bkg_events_sc1(tree, noise_events);
	cout << "Background Events: " << bkg_events.size() << endl;
	vector<int> led_events = get_led_events(tree, bkg_events);
	vector<int> non_led_events = get_events_complement(led_events, num_events);
	vector<vector<float>> ped_fits = get_pedestals(tree, bkg_events);
	get_hCal_iso(tree, non_led_events, ped_fits);
//	get_eCal_iso(tree, noise_events, led_events, ped_fits);
//	fit_eCal_iso(file);


	file->Close();
}

void Plots() {
	hCal_Current_Calibration();
}

void Noise_Events() {
	vector<double> noise_percent;
	vector<int> noise_hold;
	vector<vector<double>> correlation_ratio(11, vector<double> {});
	vector<double> run_num;
	for(auto run:config::noise_runs) {
		string name = config::beam_path_root + run + ".root";
		TFile *file_in = new TFile(name.data(), "READ");
		TTree *tree = (TTree *)file_in->Get("T");
		int num_events = get_num_events(tree);
		string tfile_name = config::noise_path + run + "_noise.root";
		TFile *file_out = new TFile(tfile_name.data(), "RECREATE");
		config::filter_coinc_thresh = 1;
		vector<int> noise_events1 = filter_noise(tree);
		for(int i = 2; i<=12; i++) {
			config::filter_coinc_thresh = i;
			noise_hold = filter_noise(tree);
			if(noise_hold.size() > 0) {
				correlation_ratio[i-2].push_back(((float)noise_events1.size() / noise_hold.size() - 1.0) * 100);
			} else if(noise_events1.size() > 0) {
				correlation_ratio[i-2].push_back(-1.0);
			} else {
				correlation_ratio[i-2].push_back(0.0);
			}
		}
		if(noise_events1.size() > 0) {
			noise_percent.push_back((float)noise_events1.size() / num_events * 100);
		} else {
			noise_percent.push_back(0.0);
		}
		run_num.push_back(stod(run.substr(3,run.size()-7)));
		cout << run << ": \t" << noise_events1.size() << "\t|\t" << num_events << "\t|\t" << noise_percent.back() << "%" << endl;
		file_in->Close();
		file_out->Close();
	}

	string tfile_name = config::noise_path + "noise_vs_event.root";
	TFile *file_out = new TFile(tfile_name.data(), "RECREATE");
	plot_noise_vs_run(noise_percent, run_num);
	plot_noise_corr(correlation_ratio, run_num);
	file_out->Close();
}


void LED_Stability(string run, int n_mv_avg) {
	TFile *file_in = new TFile((config::beam_path_root + run + ".root").data());
	TTree *tree = (TTree *)file_in->Get("T");
	string tfile_name = config::led_stbl_path + run + "_LED_Stability_" + to_string(n_mv_avg) +".root";
	TFile *file_out = new TFile(tfile_name.data(), "RECREATE");
	int num_events = get_num_events(tree);
	cout << "Total Events: " << num_events << endl;
	vector<int> noise_events = filter_noise(tree);
	vector<int> clean_events = get_events_complement(noise_events, num_events);
	vector<int> led_events = get_led_events(tree, clean_events);
	cout << "LED Events: " << led_events.size() << endl;
	vector<int> non_led_events = get_events_complement(led_events, num_events);
	vector<vector<float>> ped_fits = get_pedestals(tree, non_led_events);
	get_led_stability(tree, led_events, ped_fits, n_mv_avg); //Add n_mv_avg to declaration/definition

	file_in->Close();
	file_out->Close();
}


void sc1_events() {
	TFile *file_in = new TFile((config::beam_path_root + config::file_name + ".root").data());
	TTree *tree = (TTree *)file_in->Get("T");
	string tfile_name = config::led_stbl_path + config::file_name + "_pedestal_plots.root";
	TFile *file_out = new TFile(tfile_name.data(), "RECREATE");
	vector<int> noise_events = filter_noise(tree);
	cout << noise_events.size() << endl;
	vector<double> sc1_cuts = get_sc1_cuts(tree, noise_events);
	cout << sc1_cuts[0] << " | " << sc1_cuts[1] << endl;
	vector<int> led_events = get_led_events(tree, noise_events);
	cout << led_events.size() << endl;
	file_in->Close();
	file_out->Close();
}


void ECalAnalysis() {
	TFile *file_in = new TFile((config::beam_path_root + config::file_name + ".root").data());
	TTree *tree = (TTree *)file_in->Get("T");
	string tfile_name = config::eCal_ele_scan_path + config::file_name + "_electron_scan.root";
	TFile *file_out = new TFile(tfile_name.data(), "RECREATE");
	int num_events = get_num_events(tree);
	cout << "Total events: " << num_events << endl;
	vector<int> noise_events = filter_noise(tree);
	cout << "Noise events: " << noise_events.size() << endl;
	vector<int> clean_events = get_events_complement(noise_events, num_events);
	cout << "Clean events: " << clean_events.size() << endl;
	vector<int> non_beam_events = get_cut(tree, clean_events, "Sc1RawADC", 0, 370);
	cout << "Non-Beam events: " << non_beam_events.size() << endl;
	vector<int> beam_events = get_cut(tree, clean_events, "Sc1RawADC", 370, 4096);
	cout << "Beam events: " << beam_events.size() << endl;
	vector<int> electron_events = get_cut(tree, beam_events, "Ce1RawADC", 420, 950);
	cout << "Electron events: " << electron_events.size() << endl;
	vector<vector<float>> ped_fits = get_pedestals(tree, non_beam_events);
	vector<int> hCal_zero_events = get_hCal_zero(tree, clean_events, ped_fits);
	cout << "hCal zero events: " << hCal_zero_events.size() << endl;
	vector<int> eCal_sum_events = get_overlap_events(electron_events, hCal_zero_events, num_events);
	cout << "Electron sum events: " << eCal_sum_events.size() << endl;
	get_eCal_sum(tree, eCal_sum_events, ped_fits);

//	for(auto i:electron_events) {
//		cout << i << "\t" << flush;
//	}

	file_in->Close();
	file_out->Close();
}


void ECalAnalysis2() {
	TFile *file_in = new TFile((config::beam_path_root + config::file_name + ".root").data());
	TTree *tree = (TTree *)file_in->Get("T");
	string tfile_name = config::eCal_ele_scan_path + config::file_name + "_electron_scan.root";
	TFile *file_out = new TFile(tfile_name.data(), "RECREATE");
	int num_events = get_num_events(tree);
	cout << "Total events: " << num_events << endl;
	vector<int> noise_events = filter_noise(tree);
	cout << "Noise events: " << noise_events.size() << endl;
	vector<int> clean_events = get_events_complement(noise_events, num_events);
	cout << "Clean events: " << clean_events.size() << endl;
	vector<int> non_beam_events = get_cut(tree, clean_events, "Sc1RawADC", 0, 354);
	cout << "Non-Beam events: " << non_beam_events.size() << endl;
	vector<int> beam_events = get_cut(tree, clean_events, "Sc1RawADC", 354, 4096);
	cout << "Beam events: " << beam_events.size() << endl;
	vector<int> electron_events = get_cut(tree, beam_events, "Ce1RawADC", 420, 950);
	cout << "Electron events: " << electron_events.size() << endl;
	vector<vector<float>> ped_fits = get_pedestals(tree, clean_events);
	vector<int> hCal_zero_events = get_hCal_zero(tree, clean_events, ped_fits);
	cout << "hCal zero events: " << hCal_zero_events.size() << endl;
	vector<int> eCal_sum_events = get_overlap_events(electron_events, hCal_zero_events, num_events);
	cout << "Electron sum events: " << eCal_sum_events.size() << endl;
	get_eCal_sum(tree, eCal_sum_events, ped_fits);

	for(auto i:electron_events) {
		cout << i << "\t" << flush;
	}

	file_in->Close();
	file_out->Close();
}


void Correlate_Noise() {
	TFile *f = new TFile((config::beam_path_root + "run131" + ".root").data());
	TTree *tree = (TTree *)f->Get("T");
	string tfile_name = config::hCalibration_path + "run131" + "_corr_noise.root";
	TFile *file = new TFile(tfile_name.data(), "RECREATE");
	int num_events = get_num_events(tree);
	cout << "Total Events: " << num_events << endl;
	vector<int> noise_events = filter_noise(tree);
	cout << "Noise Events: " << noise_events.size() << endl;
	vector<vector<float>> ped_fits = get_pedestals(tree, noise_events);

	file->Close();
}


void Get_HECal_Scale136(string run) {
	cout << run << endl;
	double low_rms = 1.5;
	double high_rms = 2.5;
	TFile *file_in = new TFile((config::beam_path_root + run + ".root").data());
	TTree *tree = (TTree *)file_in->Get("T");
	string tfile_name = config::heCal_scale_path + run + "_HECal_Expanded_Range" + to_string(low_rms) + "+" + to_string(high_rms) + ".root";
	TFile *file_out = new TFile(tfile_name.data(), "RECREATE");
	int num_events = get_num_events(tree);
//	cout << "Total events: " << num_events << endl;
	vector<int> noise_events = filter_noise(tree);
//	cout << "Noise events: " << noise_events.size() << endl;
	vector<int> clean_events = get_events_complement(noise_events, num_events);
//	cout << "Clean events: " << clean_events.size() << endl;
	vector<int> non_beam_events = get_cut(tree, clean_events, "Sc1RawADC", 0, 350);
//	cout << "Non-Beam events: " << non_beam_events.size() << endl;
	vector<int> beam_events = get_cut(tree, clean_events, "Sc1RawADC", 360, 500);
//	cout << "Beam events: " << beam_events.size() << endl;
//	vector<int> electron_events = get_cut(tree, beam_events, "Ce1RawADC", 420, 950);
//	cout << "Electron events: " << electron_events.size() << endl;
	vector<int> hadron_events = get_cut(tree, beam_events, "Ce1RawADC", 0, 350);
//	cout << "Hadron events: " << hadron_events.size() << endl;
	vector<vector<float>> ped_fits = get_pedestals(tree, non_beam_events);

	double eCal_scale_start = 0.25;
	double eCal_scale_end = 0.55;
	int num_scales = 200;
	double eCal_scale;
	double res_min = 10000.0; //Bad Hardcode
	double scale_err = (eCal_scale_end - eCal_scale_start) / num_scales;
	double scale_min = 0.0;
	vector<double> scales, scale_errs, reses, res_errs;
	double mean_min, mean_min_err, sigma_min, sigma_min_err, res_min_err, res, res_err;
	for(int i=0; i<=num_scales; i++) {
		eCal_scale = (eCal_scale_end - eCal_scale_start) * ((float)i)/num_scales  + eCal_scale_start;
		TF1 *gauss_bkg = get_HECal_scale(tree, hadron_events, ped_fits, eCal_scale, low_rms, high_rms);
		res = gauss_bkg->GetParameter(2) / gauss_bkg->GetParameter(1);
		res_err = res * TMath::Sqrt(TMath::Power(gauss_bkg->GetParError(2) / gauss_bkg->GetParameter(2), 2) + TMath::Power(gauss_bkg->GetParError(1) / gauss_bkg->GetParameter(1), 2));
		scales.push_back(eCal_scale);
		scale_errs.push_back(0.0);
		reses.push_back(res);
		res_errs.push_back(res_err);
		if(reses.back() < res_min) {
			res_min = reses.back();
			scale_min = scales.back();
			res_min_err = res_errs.back();
			mean_min = gauss_bkg->GetParameter(1);
			mean_min_err = gauss_bkg->GetParError(1);
			sigma_min = gauss_bkg->GetParameter(2);
			sigma_min_err = gauss_bkg->GetParError(2);
		}
		delete gauss_bkg;
	}
	TCanvas *graph_can = new TCanvas("Optimize_Canvas", "Resolution vs. ECal Weight", config::canvas_x, config::canvas_y);
	TGraphErrors *scale_graph = new TGraphErrors((int)scales.size(), scales.data(), reses.data(), scale_errs.data(), res_errs.data());
	scale_graph->SetNameTitle(("eCal_scale_opt_" + run).data());
	scale_graph->SetTitle(("ECal Scale Factor Optimization Run " + run + " Res_Min: " + to_string(res_min) + " Scale_Min: " + to_string(scale_min)).data());
	scale_graph->Draw();
	graph_can->Write();
	scale_graph->Write();
	delete graph_can;
	delete scale_graph;

	cout << scale_min << "\t" << scale_err << "\t" << res_min << "\t" << res_min_err << "\t" << mean_min << "\t" << mean_min_err << "\t" << sigma_min << "\t" << sigma_min_err << endl;

	file_in->Close();
	file_out->Close();
}


void Get_HECal_Scale76() {
	string run = "run85";
	double low_rms = 1.5;
	double high_rms = 2.5;
	TFile *file_in = new TFile((config::beam_path_root + run + ".root").data());
	TTree *tree = (TTree *)file_in->Get("T");
	string tfile_name = config::heCal_scale_path + run + "_HECal_Ting_cuts-" + to_string(low_rms) + "+" + to_string(high_rms) + ".root";
	TFile *file_out = new TFile(tfile_name.data(), "RECREATE");
	int num_events = get_num_events(tree);
	cout << "Total events: " << num_events << endl;
	vector<int> noise_events = filter_noise(tree);
	cout << "Noise events: " << noise_events.size() << endl;
	vector<int> clean_events = get_events_complement(noise_events, num_events);
	cout << "Clean events: " << clean_events.size() << endl;
	vector<int> non_beam_events = get_cut(tree, clean_events, "Sc1RawADC", 0, 350);
	cout << "Non-Beam events: " << non_beam_events.size() << endl;
	vector<int> beam_events = get_cut(tree, clean_events, "Sc1RawADC", 360, 500);
	cout << "Beam events: " << beam_events.size() << endl;
	vector<int> electron_events = get_cut(tree, beam_events, "Ce1RawADC", 420, 950);
	cout << "Electron events: " << electron_events.size() << endl;
	vector<int> hadron_events = get_cut(tree, beam_events, "Ce1RawADC", 0, 350);
	cout << "Hadron events: " << hadron_events.size() << endl;
	vector<vector<float>> ped_fits = get_pedestals(tree, non_beam_events);
	double eCal_scale_start = 0.25;
	double eCal_scale_end = 0.55;
	int num_scales = 200;
	double eCal_scale;
	double res_min = 10000.0; //Bad Hardcode
	double scale_min = 0.0;
	double scale_err = (eCal_scale_end - eCal_scale_start) / num_scales;
	double mean_min, mean_min_err, sigma_min, sigma_min_err, res_min_err, res, res_err;
	vector<double> scales, scale_errs, reses, res_errs;
	for(int i=0; i<=num_scales; i++) {
		eCal_scale = (eCal_scale_end - eCal_scale_start) * ((float)i)/num_scales  + eCal_scale_start;
		TF1 *gauss_bkg = get_HECal_scale(tree, hadron_events, ped_fits, eCal_scale, low_rms, high_rms);
		res = gauss_bkg->GetParameter(2) / gauss_bkg->GetParameter(1);
		res_err = res * TMath::Sqrt(TMath::Power(gauss_bkg->GetParError(2) / gauss_bkg->GetParameter(2), 2) + TMath::Power(gauss_bkg->GetParError(1) / gauss_bkg->GetParameter(1), 2));
		scales.push_back(eCal_scale);
		scale_errs.push_back(scale_err);
		reses.push_back(res);
		res_errs.push_back(res_err);
		if(reses.back() < res_min) {
			res_min = reses.back();
			scale_min = scales.back();
			res_min_err = res_errs.back();
			mean_min = gauss_bkg->GetParameter(1);
			mean_min_err = gauss_bkg->GetParError(1);
			sigma_min = gauss_bkg->GetParameter(2);
			sigma_min_err = gauss_bkg->GetParError(2);
		}
		delete gauss_bkg;
	}
	TCanvas *graph_can = new TCanvas("Optimize_Canvas", "Resolution vs. ECal Weight", config::canvas_x, config::canvas_y);
	TGraphErrors *scale_graph = new TGraphErrors((int)scales.size(), scales.data(), reses.data(), scale_errs.data(), res_errs.data());
	scale_graph->SetNameTitle(("eCal_scale_opt_" + run).data());
	scale_graph->SetTitle(("ECal Scale Factor Optimization Run " + run + " Res_Min: " + to_string(res_min) + " Scale_Min: " + to_string(scale_min)).data());
	scale_graph->Draw();
	graph_can->Write();
	scale_graph->Write();
	delete graph_can;
	delete scale_graph;

	cout << "For minimum resolution weight: " << endl;
	cout << "Weight: " << scale_min << "+-" << scale_err << endl;
	cout << "Res: " << res_min << "+-" << res_min_err << endl;
	cout << "Mean: " << mean_min << "+-" << mean_min_err << endl;
	cout << "Sigma: " << sigma_min << "+-" << sigma_min_err << endl;
	cout << scale_min << "\t" << scale_err << "\t" << res_min << "\t" << res_min_err << "\t" << mean_min << "\t" << mean_min_err << "\t" << sigma_min << "\t" << sigma_min_err << endl;

	cout << run << endl;

	file_in->Close();
	file_out->Close();
}


void Get_Hadron_Res(string run) {
	cout << run << endl;
	double low_rms = 1.5;
	double high_rms = 2.5;
	double eCal_weight = 0.35875;
	TFile *file_in = new TFile((config::beam_path_root + run + ".root").data());
	TTree *tree = (TTree *)file_in->Get("T");
	string tfile_name = config::heCal_scale_path + run + "_HECal_Expanded_Range_Ting_cuts-" + to_string(low_rms) + "+" + to_string(high_rms) + ".root";
	TFile *file_out = new TFile(tfile_name.data(), "RECREATE");
	int num_events = get_num_events(tree);
//	cout << "Total events: " << num_events << endl;
	vector<int> noise_events = filter_noise(tree);
//	cout << "Noise events: " << noise_events.size() << endl;
	vector<int> clean_events = get_events_complement(noise_events, num_events);
//	cout << "Clean events: " << clean_events.size() << endl;
	vector<int> non_beam_events = get_cut(tree, clean_events, "Sc1RawADC", 0, 350);
//	cout << "Non-Beam events: " << non_beam_events.size() << endl;
	vector<int> beam_events = get_cut(tree, clean_events, "Sc1RawADC", 360, 500);
//	cout << "Beam events: " << beam_events.size() << endl;
	vector<int> electron_events = get_cut(tree, beam_events, "Ce1RawADC", 420, 950);
//	cout << "Electron events: " << electron_events.size() << endl;
	vector<int> hadron_events = get_cut(tree, beam_events, "Ce1RawADC", 0, 350);
//	cout << "Hadron events: " << hadron_events.size() << endl;
	vector<vector<float>> ped_fits = get_pedestals(tree, non_beam_events);
	double res, res_err, mean, mean_err, sigma, sigma_err;
	TF1 *gauss_bkg = get_HECal_scale(tree, hadron_events, ped_fits, eCal_weight, low_rms, high_rms);
	res = gauss_bkg->GetParameter(2) / gauss_bkg->GetParameter(1);
	res_err = res * TMath::Sqrt(TMath::Power(gauss_bkg->GetParError(2) / gauss_bkg->GetParameter(2), 2) + TMath::Power(gauss_bkg->GetParError(1) / gauss_bkg->GetParameter(1), 2));
	mean = gauss_bkg->GetParameter(1);
	mean_err = gauss_bkg->GetParError(1);
	sigma = gauss_bkg->GetParameter(2);
	sigma_err = gauss_bkg->GetParError(2);
	delete gauss_bkg;
	cout << res << "\t" << res_err << "\t" << mean << "\t" << mean_err << "\t" << sigma << "\t" << sigma_err << endl;

	file_in->Close();
	file_out->Close();
}
