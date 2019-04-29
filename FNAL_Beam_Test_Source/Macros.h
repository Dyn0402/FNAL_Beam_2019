/*
 * Macros.h
 *
 *  Created on: Apr 6, 2019
 *      Author: dylan
 */

#ifndef MACROS_H_
#define MACROS_H_

#include <string>
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

//ReadoutHistograms.cpp
void ReadoutHistograms(string fname, string outfilename);

//filter_noise.cpp
vector<int> filter_noise(TTree *tree);
void plot_noise_vs_run(vector<double> noise_percent, vector<double> run_num);
void plot_noise_corr(vector<vector<double>> noise_corr, vector<double> run_num);

//functions.cpp
vector<int> get_events_complement(vector<int> events, int num_events);
vector<int> get_overlap_events(vector<int> events1, vector<int> events2, int num_events);
int get_num_events(TTree *tree);
vector<double> get_moving_average(vector<double> list, unsigned int n);
double get_mean(vector<double> list);

//sc1_cuts.cpp
vector<double> get_sc1_cuts(TTree *tree, vector<int> noise_events);
vector<int> get_bkg_events_sc1(TTree *tree, vector<int> noise_events);
vector<int> get_beam_events_sc1(TTree *tree, vector<int> noise_events);

//ce1_cuts.cpp
vector<int> get_electron_events(TTree *tree, vector<int> beam_events, vector<vector<float>> ped_fits);
vector<int> get_hadron_events(TTree *tree, vector<int> beam_events, vector<vector<float>> ped_fits);

//eCal_cuts.cpp
vector<int> get_led_events(TTree *tree, vector<int> events);
vector<int> get_led_events2(TTree *tree, vector<int> events);

//pedestals.cpp
vector<vector<float>> get_pedestals(TTree *tree, vector<int> ped_events);

//eCal_electron_scan.cpp
void get_eCal_sum(TTree *tree, vector<int> events, vector<vector<float>> ped_fits);
vector<int> get_hCal_zero(TTree *tree, vector<int> events, vector<vector<float>> ped_fits);
vector<int> get_cut(TTree* tree, vector<int> events, string branch, double lower, double upper);
vector<int> get_cut_cal(TTree* tree, vector<int> events, string branch, int channel, double lower, double upper);

//isolation.cpp
void get_hCal_iso(TTree *tree, vector<int> hadron_events, vector<vector<float>> ped_fits);
void get_hCal_iso2(TTree *tree, vector<int> hadron_events, vector<vector<float>> ped_fits);
void get_hCal_iso_dynamic(TTree *tree, vector<int> hadron_events, vector<vector<float>> ped_fits);
void get_eCal_iso(TTree *tree, vector<int> hadron_events, vector<vector<float>> ped_fits);
void fit_eCal_iso(TFile *file);
void check_hcalibration(TTree *tree, vector<int> events, vector<vector<float>> ped_fits);

//LED_Stability.cpp
void get_led_stability(TTree *tree, vector<int> events, vector<vector<float>> ped_fits, int n_mv_avg);

//HECal_scale.cpp
TF1 *get_HECal_scale(TTree *tree, vector<int> events, vector<vector<float>> ped_fits, double eCal_scale_factor, double low_rms, double high_rms);

//tree_pass_test.cpp
TH1D* pass_tree(TTree *tree);

void example();


#endif /* MACROS_H_ */
