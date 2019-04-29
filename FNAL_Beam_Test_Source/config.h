/*
 * config.h
 *
 *  Created on: Apr 6, 2019
 *      Author: dylan
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#include <string>
#include <vector>

using namespace std;

namespace config
{
	extern string beam_path_raw;
	extern string beam_path_root;
	extern string noise_path;
	extern string led_stbl_path;
	extern string hCalibration_path;
	extern string hCalibration_debug_path;
	extern string eCal_ele_scan_path;
	extern string heCal_scale_path;
	extern string file_name;

	extern string filter_hist_name;
	extern string filter_fit_name;
	extern string tfile_name;
	extern string sc1_hist_name;
	extern string sc1_fit_name;
	extern string sc1_bkg_hist_name;
	extern string sc1_bkg_fit_name;
	extern string sc1_beam_hist_name;
	extern string sc1_beam_fit_name;
	extern string ce1_electron_hist_name;
	extern string ce1_hadron_hist_name;
	extern string eCal_hist_name;
	extern string eCal_fit_name;
	extern string eCal_iso_hist_name;
	extern string hCal_iso_hist_name;

	extern bool write_mode;

	extern vector<int> eCal_cut_channels;
	extern vector<vector<int>> adj_channels;
	extern vector<float> hCal_weights;
	extern vector<float> eCal_weights;
	extern vector<string> noise_runs;
	extern vector<string> readout_runs;

	extern const int junk_channels;
	extern const int eCal_channels;
	extern const int hCal_channels;
	extern const int hod_channels;

	extern int hCal_start_channel;
	extern int eCal_start_channel;
	extern int hod_start_channel;
	extern int junk_start_channel;
	extern int sc1_channel;
	extern int sc2_channel;
	extern int ce1_channel;
	extern int ce2_channel;

	extern int bins;
	extern int canvas_x;
	extern int canvas_y;
	extern int filter_coinc_thresh;
	extern int eCal_coinc_thresh;
	extern int mv_avg_n;

	extern double low_bin;
	extern double high_bin;

	extern double filter_rms_range;
	extern double filter_fit_range;
	extern double filter_spec_sigma;
	extern double filter_spec_thresh;

	extern double ped_spec_thresh;

	extern double sc1_fit_range;
	extern double sc1_spec_sigma;
	extern double sc1_spec_thresh;
	extern double sc1_low_sigmas;
	extern double sc1_high_sigmas;
	extern double sc1_bkg_spec_sigma;
	extern double sc1_bkg_spec_thresh;
	extern double sc1_bkg_sigmas;
	extern double sc1_beam_sigmas;

	extern double ce1_electron_sigmas;
	extern double ce1_hadron_sigmas;

	extern int eCal_spec_peaks;
	extern double eCal_spec_sigma;
	extern double eCal_spec_thresh;
	extern double eCal_fit_range;
	extern double eCal_led_low_ADC;
	extern double eCal_rms_range;

	extern double iso_hCal_signal_sigma;
	extern double iso_hCal_adjsignal_sigma;
	extern int iso_hCal_rebinning;

	extern double iso_eCal_signal_sigma;
	extern double iso_eCal_adjsignal_sigma;
	extern double iso_eCal_adj_const;

	extern double eh_scale_hCal_signal_sigma;
	extern double eh_scale_eCal_signal_sigma;
}



#endif /* CONFIG_H_ */
