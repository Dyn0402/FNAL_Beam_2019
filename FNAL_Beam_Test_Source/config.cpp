/*
 * config.cpp
 *
 *  Created on: Apr 6, 2019
 *      Author: dylan
 */

#include <string>
#include <vector>

using namespace std;

namespace config
{
	string beam_path_raw = "/home/dylan/Dropbox/Research/Workspace2/FNAL_Beam_2019/Data/Beam_Runs/Raw/";
	string beam_path_root = "/home/dylan/Dropbox/Research/Workspace2/FNAL_Beam_2019/Data/Beam_Runs/Root/";
	string noise_path = "/home/dylan/Dropbox/Research/Workspace2/FNAL_Beam_2019/Analysis/Noise_Fits/";
	string led_stbl_path = "/home/dylan/Dropbox/Research/Workspace2/FNAL_Beam_2019/Analysis/LED_Stability/";
	string hCalibration_path = "/home/dylan/Dropbox/Research/Workspace2/FNAL_Beam_2019/Analysis/HCalibration/";
	string hCalibration_debug_path = "/home/dylan/Dropbox/Research/Workspace2/FNAL_Beam_2019/Analysis/HCalibration_Debug/";
	string eCal_ele_scan_path = "/home/dylan/Dropbox/Research/Workspace2/FNAL_Beam_2019/Analysis/ECal_Electron_Scan/";
//	string heCal_scale_path = "/home/dylan/Dropbox/Research/Workspace2/FNAL_Beam_2019/Analysis/HECal_scale/HEScale_04-16_hCalibration_clean/";
	string heCal_scale_path = "/home/dylan/Dropbox/Research/Workspace2/FNAL_Beam_2019/Analysis/HECal_scale/Hadron_Resolution_ECalWeight_0.35875/";
	string file_name = "run25";

	string filter_hist_name = "Junk";
	string filter_fit_name = "Gauss_Filter";
	string sc1_hist_name = "SC1";
	string sc1_fit_name = "Gauss_Sc1";
	string sc1_bkg_hist_name = "Sc1_bkg";
	string sc1_bkg_fit_name = "Gauss_Sc1_bkg";
	string sc1_beam_hist_name = "Sc1_beam";
	string sc1_beam_fit_name = "Gauss_Sc1_beam";
	string ce1_electron_hist_name = "Ce1_beam_events_electron_selection";
	string ce1_hadron_hist_name = "Ce1_beam_events_hadron_selection";
	string eCal_hist_name = "ECal Channel ";
	string eCal_fit_name = "Gauss_ECal";
	string hCal_iso_hist_name = "HCal Isolated Channel ";
	string eCal_iso_hist_name = "ECal Isolated Channel ";

	bool write_mode = true;

	vector<int> eCal_cut_channels {5, 6, 9, 10};
//	vector<int> eCal_cut_channels = {2, 4, 8, 11, 12, 13, 14, 15};
	vector<vector<int>> adj_channels {{1,4,5}, {0,2,4,5,6}, {1,3,5,6,7}, {2,6,7},
			{0,1,5,8,9}, {0,1,2,4,6,8,9,10}, {1,2,3,5,7,9,10,11}, {2,3,6,10,11},
			{4,5,9,12,13}, {4,5,6,8,10,12,13,14}, {5,6,7,9,11,13,14,15}, {6,7,10,14,15},
			{8,9,13}, {8,9,10,12,14}, {9,10,11,13,15}, {10,11,14}};
	vector<float> hCal_weights = {0.89348, 0.956762, 0.828132, 0.952171, 0.439308, 0.978776, 1.01891, 1.17377, 1.2098, 0.876184, 1.28611, 0.988631, 1.05044, 1.04991, 1.12791, 1.1697}; //Muon run 134 rough.
//	vector<float> hCal_weights = {0.825337496159864, 0.695614600375798, 0.963687673699176, 0.732770256271031, 0.736142888782516, 1.37158904, 1.1455043, 1.32926347977052, 1.09001507476656, 0.9962423, 0.48666436, 1.22153676885596, 0.945955747344055, 0.8490071944502, 0.872301952575213, 1.19172353878359}; //Muon runs 20 and 21.
	vector<float> eCal_weights = {.974462, .954643, .946785, 1.01851, 1.00689, .920467, .90186, 1.0301, 1.11047, 1.03755, 1.03422, 1.10631, 1.07239, .898038, 1.0314, .955892};
//	vector<string> noise_runs {"run18", "run19", "run20"};
	vector<string> noise_runs = {"run1","run2","run3","run4","run5","run6","run7","run8","run9","run10","run11","run12","run13","run14","run15","run16","run17","run18","run19","run20","run21","run22","run23","run24","run25","run26","run27","run28","run29","run30","run31","run32","run33","run34","run35","run36","run37","run38","run39","run40","run41","run42","run43","run44","run45","run46","run47","run48","run49","run50","run51","run52","run53","run54","run55","run56","run57","run58","run59","run60","run61","run62","run63","run64","run65","run66","run67","run68","run69","run70","run71","run72","run73","run74","run75","run76","run77","run78","run79","run80","run81","run82","run83","run84","run85","run86","run87","run88","run89","run90","run91","run92","run93","run94","run95","run96","run97","run98","run99","run100","run101","run102","run103","run104","run105","run106","run107","run108","run109","run110","run111","run112","run113","run114","run115","run116","run117","run118","run119","run120","run121","run122","run123","run124","run125","run126","run127","run128","run129","run130","run131","run132","run133","run134","run135","run136","run137","run138","run139","run140","run141","run142","run143","run144","run145","run146"};
	vector<string> readout_runs = {"run1","run2","run3","run4","run5","run6","run7","run8","run9","run10","run11","run12","run13","run14","run15","run16","run17","run18","run19","run20","run21","run22","run23","run24","run25","run26","run27","run28","run29","run30","run31","run32","run33","run34","run35","run36","run37","run38","run39","run40","run41","run42","run43","run44","run45","run46","run47","run48","run49","run50","run51","run52","run53","run54","run55","run56","run57","run58","run59","run60","run61","run62","run63","run64","run65","run66","run67","run68","run69","run70","run71","run72","run73","run74","run75","run76","run77","run78","run79","run80","run81","run82","run83","run84","run85","run86","run87","run88","run89","run90","run91","run92","run93","run94","run95","run96","run97","run98","run99","run100","run101","run102","run103","run104","run105","run106","run107","run108","run109","run110","run111","run112","run113","run114","run115","run116","run117","run118","run119","run120","run121","run122","run123","run124","run125","run126","run127","run128","run129","run130","run131","run132","run133","run134","run135","run136","run137","run138","run139","run140","run141","run142","run143","run144","run145","run146"};

	int junk_channels = 12; //Number of junk channels.
	int eCal_channels = 16; //Number of eCal channels.
	int hCal_channels = 16; //Number of hCal channels.
	int hod_channels = 16; //Number of hod channels.


	int hCal_start_channel = 0;
	int eCal_start_channel = 16;
	int hod_start_channel = 32;
	int junk_start_channel = 52;
	int sc1_channel = 48;
	int sc2_channel = 49;
	int ce1_channel = 50;
	int ce2_channel = 51;

	int bins = 4096;
	int canvas_x = 2200;
	int canvas_y = 1400;
	int filter_coinc_thresh = 1; //Number of junk channels flagged as noise required to flag event as noise event.
	int eCal_coinc_thresh = 1; //Number of eCal central channels flagged as LED required to flag event as LED event.
	int mv_avg_n = 5000; //Number of points to average for LED stability moving average.

	double low_bin = 0.5;
	double high_bin = 4096.5;

	double filter_rms_range = 5.0; //Number of sigma outside of noise pedestal fit required to be considred noise.
	double filter_fit_range = 3.0; //Number of filter_spec_sigma to fit junk pedestal.
	double filter_spec_sigma = 2.0; //TSpectrum initial guess for junk pedestal.
	double filter_spec_thresh = 0.3; //TSpectrum threshold for peak.

	double ped_spec_thresh = 0.3;

	double sc1_fit_range = 1.0; //Number of rms width
	double sc1_spec_sigma = 5.0;
	double sc1_spec_thresh = 0.01;
	double sc1_low_sigmas = 1.0; //Number of sigmas above fit low peak to set upper pulser cut.
	double sc1_high_sigmas = 3.0; //Number of sigmas above fit low peak to set lower beam cut.
	double sc1_bkg_spec_sigma = 5.0;
	double sc1_bkg_spec_thresh = 0.01;
	double sc1_bkg_sigmas = 3.5; //Number of sigmas above fit low peak to set background cut.
	double sc1_beam_sigmas = 4.5; //Number of sigmas above fit low peak to set beam cut.

	double ce1_electron_sigmas = 1.0; //Number of pedestal sigmas higher than pedestal mean to set lower cut for electron selection from ce1.
	double ce1_hadron_sigmas = 4.0; //Number of pedestal sigmas higher than pedestal mean to set upper cut for hadron selection from ce1.

	int eCal_spec_peaks = 4;
	double eCal_spec_sigma = 40.0;
	double eCal_spec_thresh = 0.3;
	double eCal_fit_range = 3.0; //Number of rms width
	double eCal_led_low_ADC = 250;//900;//1200;250
	double eCal_rms_range = 3.0;

	double iso_hCal_signal_sigma = 3.0;
	double iso_hCal_adjsignal_sigma = 3.0;
	int iso_hCal_rebinning = 16;

	double iso_eCal_signal_sigma = 1.0;
	double iso_eCal_adjsignal_sigma = 1.0;
	double iso_eCal_adj_const = 0;

	double eh_scale_hCal_signal_sigma = 3.0;
	double eh_scale_eCal_signal_sigma = 3.0;
}
