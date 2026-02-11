# -*- coding: utf-8 -*-
"""
Created on Wed May 7 12:11 2025

@author: Cristiana Pinheiro
"""
import sys
source_path = 'C:/Users/neuro/Desktop/Cristiana/paper_insula/githubINS/'
sys.path.append(source_path+'ft-insulaOS/ENIGMA')
sys.path.insert(0,source_path+"ft-atlas/src/scripts")
print('You need to change these paths during first use')
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mne
import numpy as np
import gzip
#from compute_atlasCP import main_CA
from py_tools.stats import get_binom_ci
import re
from enigmatoolbox.plotting import plot_subcortical
from PIL import Image
import statsmodels.stats.multitest as smm
from itertools import combinations
from matplotlib.colors import LinearSegmentedColormap
from statsmodels.stats.proportion import proportions_ztest
from dictROIs_Insula import idx, yticks
#from dictROIs_Amyg import idx, yticks

filter_default = {
                'age_range':(0, np.inf),
                'spike_count_max': np.inf,
                'spike_count_min':0,
                'spike_count_frac':0,
                'include_with_no_spike_count':True,
                'segmentations_stim':(0,1,2), #segmentation: nan / GM / WM (0,1,2)
                'segmentations_rec':(0,1,2),
                'dist_min':None,
                'only_symmetric_flag':False,
                'crf_exclusive_list':None,
                'allowed_centres_list':None,
                'fibre_condition':(),
                'intensity_range':(0, np.inf),
                'frequency_range':(0, 1.1), #(0, np.inf)
                'pulse_width_range':(0, np.inf),
                'only_different_electrodes_flag':False,
                'min_x_coordinate':0.,
                'allowed_contact_ids':None,
    }

#bas -400ms
"""surro_dict = { #probabilities
    100: 0.03,
    400: 0.08
}"""
#bas -200ms
surro_dict = { #probabilities
    50: 0.08,
    100: 0.15,
    150: 0.20,
    200: 0.20,
    250: 0.24,
    300: 0.24,
    350: 0.27,
    400: 0.29
}

# Ensure that you have dictROIs.py file if you need to merge parcels
def main(filter_default, resolution_stim = 1, resolution_rec = np.nan,
         stim_parc_name = 'MNI-insula', rec_parc_name = "Lausanne2018-scale2", plot_parc_name = "Lausanne2018-scale2",
         merged = [True, False], hemis_symmetrize = [True, "H"],
         ntotal_min = 50, time_window = ["max", 100],
         efferent = True, save = False):
    #variables
    # 'MNI-JulichBrain-3.0' 'MNI-insula'
    # 'Lausanne2008-125' 'Lausanne2008-60' 'Lausanne2008-33'
    # merged = [stim, rec]
    # hemis_symmetrize = [True, "H"] True/False, H/V (Horizontal/Vertical)
    # time_window = ["max", 100] max_peak_delay , 100ms
    trial = stim_parc_name + "/" + str(resolution_stim) + "/" + rec_parc_name + "/" + str(
        resolution_rec) + ""
    surro = surro_dict[time_window[1]]
    z = "5"
    min_value_impl = 2
    se_dir = source_path+'ft-insulaOS/'

    # select parcel groups
    parcels_idxStim, ytickStim, parcels_idxRec, ytickRec = select_parcel_groups(merged,
                                                                                stim_parc_name, rec_parc_name,
                                                                                resolution_stim, resolution_rec)

    output_directory_name = source_path+"ft-insulaOS/results/" + trial + "/"
    print(output_directory_name)
    if os.path.isdir(output_directory_name):
        print("Available data will be used")
        analyze(se_dir, time_window, z, ntotal_min, min_value_impl, surro, resolution_stim, resolution_rec,
                plot_parc_name, stim_parc_name, rec_parc_name, merged, hemis_symmetrize, efferent,
                parcels_idxStim, parcels_idxRec, ytickStim, ytickRec,
                output_directory_name, save)
    else:
        print("No data available")
        #access to raw data is not allowed
        """
        main_CA(filter_default, output_directory_name, se_dir, stim_parc_name, rec_parc_name,
                z, ntotal_min, time_window[1],
                merge_flag=merged, ROIsToMergeStim=parcels_idxStim, ROIsToMergeRec=parcels_idxRec,
                sym_flag=hemis_symmetrize,
                save_flag=True)
        main(filter_default, efferent=efferent, merged=merged, hemis_symmetrize=hemis_symmetrize, time_window =time_window, resolution_stim=resolution_stim,
             resolution_rec=resolution_rec, stim_parc_name=stim_parc_name, rec_parc_name=rec_parc_name, plot_parc_name=plot_parc_name, ntotal_min=ntotal_min)"""
    return

# FUNCTIONS
#______________________________________________________________________________________________________________________

"""
Compute the indices (`parcels_idx`) and labels (`ytick`) for stimulation and recording parcels
based on the provided stimulation and recording parcellation schemes.
"""
def select_parcel_groups(merged, stim_parc_name, rec_parc_name, resolution_stim, resolution_rec):
    if merged[0]==True:
        parcels_idxStim = idx[stim_parc_name][resolution_stim]
        ytickStim = yticks[stim_parc_name][resolution_stim]
    else:
        parcels_idxStim = []
        ytickStim = []
    if merged[1]==True:
        parcels_idxRec = idx[rec_parc_name][resolution_rec]
        ytickRec = yticks[rec_parc_name][resolution_rec]
    else:
        parcels_idxRec = []
        ytickRec = []
    return parcels_idxStim, ytickStim, parcels_idxRec, ytickRec

def analyze(se_dir, time_window, z, ntotal_min, min_value_impl, surro, resolution_stim, resolution_rec,
            plot_parc_name, stim_parc_name, rec_parc_name, merged, hemis_symmetrize, efferent,
            parcels_idxStim, parcels_idxRec, ytickStim, ytickRec,
            output_directory_name, save):
    # paths
    # mesh path
    if plot_parc_name == "MNI-JulichBrain-3.0":
        meshdirname = source_path+"ft-insulaOS/mne_data/MNE-sample-data/subjects/fsaverage"
    elif plot_parc_name[:8] == "Lausanne":
        meshdirname = source_path+"ft-insulaOS/mne_data/MNE-sample-data/subjects/cvs_avg35_inMNI152"
    else:
        print("meshdirname is not defined for this plot parcellation")
    # data path
    paths_discVAR, paths_contVAR = generate_paths(output_directory_name, time_window, z, ntotal_min)
    # order of measures_disc and measures_cont needs to match order of paths_discVAR and paths_contVAR, respectively
    measures_disc = ["probability", "NTotal", "NValue", "patients", "implantations"]
    measures_cont = ["peak_delay", "ampl"] #["peak_delay", "ampl", "speed"]
    features_cont = ["quantile_0.5", "mean", "std", "quantile_0.25", "quantile_0.75"]

    # read parcels names
    roi_stim = parcels_labels(se_dir, stim_parc_name, merged[0], hemis_symmetrize, False, idx=parcels_idxStim)
    roi_rec = parcels_labels(se_dir, rec_parc_name, merged[1], hemis_symmetrize, True, idx=parcels_idxRec)

    # read atlas data
    all_data_disc = []
    for i, j in enumerate(measures_disc):
        data = read_atlas_data(paths_discVAR[i], hemis_symmetrize)
        all_data_disc.append(data)
    all_data_cont = []
    for i, j in enumerate(measures_cont):
        for k in features_cont:
            data = read_atlas_data(paths_contVAR[i] + "nan" + k + ".txt", hemis_symmetrize)
            all_data_cont.append(data)

    # mask data
    all_data_disc_masked, all_data_cont_masked, mask, CI_masked = masking(all_data_disc, all_data_cont, stim_parc_name, rec_parc_name, hemis_symmetrize, ntotal_min, min_value_impl, surro)

    # plots
    overall_plot("Confidence Interval", CI_masked, 0, efferent,
                 plot_parc_name, stim_parc_name, rec_parc_name, meshdirname,
                 roi_stim, roi_rec, parcels_idxStim, parcels_idxRec,
                 output_directory_name, save)
    for i, j in enumerate(measures_disc):
        #if (i == 0):  # "probability","NTotal","NValue","patients","implantations"
        data = all_data_disc_masked[i]
        #"""
        overall_plot(j, data, surro, efferent,
                     plot_parc_name, stim_parc_name, rec_parc_name, meshdirname,
                     roi_stim, roi_rec, parcels_idxStim, parcels_idxRec,
                     output_directory_name, save)#"""
        print("overall_plot done: " + j)
    w = 0
    for i, j in enumerate(measures_cont):  # "peak_delay","ampl"
        for k2, k in enumerate(features_cont):  # "quantile_0.5", "mean", "std", "quantile_0.25", "quantile_0.75"
            if (i == 0) or (i == 1): #and (k2 == 0):
                data = all_data_cont_masked[w]
                #"""
                overall_plot(j + "_" + k, data, surro, efferent,
                             plot_parc_name, stim_parc_name, rec_parc_name, meshdirname,
                             roi_stim, roi_rec, parcels_idxStim, parcels_idxRec,
                             output_directory_name, save)#"""
            w = w + 1
            print("overall_plot done: " + j + "_" + k)

    # statistical differences between stimulation parcels
    if efferent == True:
        resolution = resolution_stim
    else:
        resolution = resolution_rec
    if resolution > 1:
        overall_stats(all_data_disc_masked, resolution, roi_stim, roi_rec, parcels_idxStim, parcels_idxRec,
                      plot_parc_name, stim_parc_name, rec_parc_name, hemis_symmetrize, efferent, surro, mask,
                      meshdirname, ytickStim, ytickRec, output_directory_name, save)
        print("overall_stats done")
    return

"""
Generate the file paths corresponding to all available data inside output_directory_name.
"""
def generate_paths(output_directory_name, time_window, z, ntotal_min):
    interm_path = (output_directory_name + "/" +
                   str(time_window[0]) + "_peak_delay_" + str(time_window[1]) +
                   ".0__zth" + z + "/feature_ampl_zth" + z +
                   "__N_with_value__min_value__" + str(ntotal_min))
    p_path = output_directory_name + "/probability.txt"
    ntotal_path = interm_path + "/feature_ampl_zth" + z + "/N_total.txt"
    nvalue_path = interm_path + "/feature_ampl_zth" + z + "/N_with_value.txt"
    npat_path =  interm_path +"/patient_name/count_unique_str.txt"
    no_impl_path = interm_path +"/implantation_name/count_unique_str.txt"
    pdelay_path = interm_path + "/feature_peak_delay_zth" + z + "/"
    ampl_path = interm_path + "/feature_ampl_zth" + z + "/"
    speed_path = interm_path + "/speed_euclidean_distance_pat_feature_peak_delay_zth" + z + "/"
    paths_disc = [p_path, ntotal_path, nvalue_path, npat_path, no_impl_path]
    paths_cont = [pdelay_path, ampl_path, speed_path]
    return paths_disc, paths_cont

"""
Load parcel names from supported parcellation schemes (Lausanne, Montreal,
and Julich) and modify them to reflect optional parcel merging and left–right hemispheric symmetry.
"""
def parcels_labels(se_dir, parc_name, merge_flag, sym_flag, rec, idx=[]):
    #read available rois
    roi = []
    labels_p = se_dir + "/parcellation_definitions/" + parc_name + ".txt"
    with open(labels_p, 'r') as file:
        for line in file:
            roi.extend(line.strip().split())
    roi = np.array(roi)
    """
    if 'Lausanne' in parc_name:
        roi = [roi_.replace("lh.", "_L") for roi_ in roi]
        roi = [roi_.replace("rh.", "_R") for roi_ in roi]"""

    #adjust rois to the analysis
    #roi0 = labels_uniformized(roi, parc_name, False)
    if merge_flag==True:
        roiM=labels_fitMerge(idx, roi)
    else:
        roiM=roi.copy()
    roiU = labels_uniformized(roiM, parc_name, merge_flag)
    if sym_flag[0] == True:
        if sym_flag[1] == "V":
            rec = not rec
        roiS = labels_fitHemiSym(rec, roiU)
    else:
        roiS = roiU.copy()

    return [roiS, roi]

"""
Generate a label for a group of merged parcels by selecting the label of the first
parcel in the group and adding an asterisk (*) to indicate merging.
"""
def labels_fitMerge(idx, roi):
    merged_list = []
    for i in idx:
        start = i[0]
        start_label = roi[start]
        if len(i) > 1:
            final_label = "*"
        else:
            final_label = " "
        merged_list.append(start_label+final_label)
    return merged_list

"""
Uniformize parcel naming by reorganizing parcel names so that all left-hemisphere
parcels precede right-hemisphere parcels, and by adding an 'L' (left) or 'R' (right) suffix to identify
the parcel’s hemisphere. Supported parcellation schemes: Lausanne, Montreal, and Julich.
"""
def labels_uniformized(roi, parc_name, merged):
    #Left - _L, Right - _R
    left = []
    right = []
    if 'Lausanne' in parc_name:  # remove brainstem and both unknowns
        if merged == True:
            roiU = roi.copy()
        else:
            if 'Lausanne2008' in parc_name:
                roiU = roi[:-2]
            elif 'Lausanne2018-scale1' in parc_name:
                roiU = roi[:-70]
            elif 'Lausanne2018-scale2' in parc_name:
                roiU = roi[:-116]
            elif 'Lausanne2018-scale3' in parc_name:
                roiU = roi[:-218]
            else:
                print("this parc_name is not defined in labels_uniformized function")
        left = ["lh.", "ctx-lh-", "Left-", "thal-lh-", "subc-lh-"]
        right = ["rh.", "ctx-rh-", "Right-", "thal-rh-", "subc-rh-"]
    elif 'JulichBrain' in parc_name:
        if merged == True:
            roiU = roi.copy()
        else:
            roiU = np.concatenate((roi[0::2], roi[1::2])) # first all left then all right
        roiU = [roi.replace("_", " ") for roi in roiU]
        left = [" L"]
        right = [" R"]
    elif "MNI-insula" in parc_name:
        roiU = roi.copy()
        left = ["L"]
        right = ["R"]
    else:  # just remove unknown
        print("this parc_name is not defined in labels_uniformized function")
        #roiU = roi.copy()

    list = []
    for i, j in enumerate(roiU):
        new = j
        for l in left:
            if l in j:
                new = j.replace(l, "") + "_L"
                break
        list.append(new)
    roiU = []
    for i, j in enumerate(list):
        new = j
        for r in right:
            if r in j:
                new = j.replace(r, "") + "_R"
                break
        roiU.append(new)

    return roiU

"""
Adapt parcel names for hemisphere symmetrization by substituting hemisphere labels
('L'/'R') with interaction labels ('intra', 'inter').

Symmetric left and right parcels are merged (e.g., L–L and R–R → intra; L–R and R–L → inter),
resulting in a twofold reduction in the number of parcels.
"""
def labels_fitHemiSym(rec, roi):
    new_roi = []
    if rec == True:
        for i, j in enumerate(roi):
            if j[-2:] == "_L":
                new_roi.append("intra." + j[:-2])
            else: #"_R"
                new_roi.append("inter." + j[:-2])
    else:
        roi = roi[:len(roi) // 2]
        for i, j in enumerate(roi):
            if j[-1] == "*":
                new_roi.append(j[:-3]+"*")
            else:
                new_roi.append(j[:-2])
    return new_roi

"""
Load data from files stored in the directory specified by `p_path`.
"""
def read_atlas_data(p_path, sym_flag):
    #read saved data
    atlas_data = []
    if os.path.exists(p_path):
        with open(p_path, 'r') as file:
            for line in file:
                if line[0:1]=="#":
                    continue
                else:
                    row = list(map(float, line.split()))  # Convert each element to an integer
                    atlas_data.append(row)
    else:
        atlas_data.append(None)
    atlas_data = np.array(atlas_data)

    #organize properly for the analysis
    if sym_flag[0] == True and sym_flag[1]=="V":
        atlas_data = np.concatenate((atlas_data[:,:atlas_data.shape[1]//2], atlas_data[:,atlas_data.shape[1]//2:]), axis=0)

    return atlas_data

"""
Create a parcel inclusion mask by selecting parcels that meet minimum response and
patient count criteria. The resulting brain coverage is computed as a percentage.

The mask is applied to probability, ntotal, nvalue, and npatients data, and it is
plotted.
"""
def masking(all_data_disc, all_data_cont, stim_parc_name, rec_parc_name, sym_flag, N_cond=50, min_value_impl=2, surro=0.15):
    prob = all_data_disc[0]
    NTotal = all_data_disc[1]
    NValue = all_data_disc[2]
    NPatients = all_data_disc[3]
    NImplantations = all_data_disc[4]

    #compute mask using number of minimum responses
    CI = get_binom_ci(prob, NTotal, 0.05, True)
    mask = (NTotal > N_cond) & (NPatients >= min_value_impl)
    mask2 = (NTotal > N_cond) & (NPatients >= min_value_impl) & (prob > surro)

    # print coverage as %
    coverage(sym_flag, mask)

    # apply mask
    all_data_disc[0] = np.where(mask, prob, np.nan)
    CI_masked = np.where(mask, CI, np.nan)
    all_data_disc[1] = np.where(mask, NTotal, np.nan)
    all_data_disc[2] = np.where(mask, NValue, np.nan)
    all_data_disc[3] = np.where(mask, NPatients, np.nan)
    all_data_disc[4] = np.where(mask, NImplantations, np.nan)
    for i in range(len(all_data_cont)):
        data_cont = all_data_cont[i].copy()
        all_data_cont[i] = np.where(mask2, data_cont, np.nan)

    # plot mask
    plot_mask(mask, CI_masked, stim_parc_name, rec_parc_name)
    print("Confidence interval: Max, Min, Mean")
    CImax = np.nanmax(CI_masked)
    CImin = np.nanmin(CI_masked)
    CImean = np.nanmean(CI_masked)
    print(CImax)
    print(CImin)
    print(CImean)

    return all_data_disc, all_data_cont, mask, CI_masked

"""
Plot the parcel mask and the confidence interval lengths of the probability estimates.
"""
def plot_mask(mask, CImasked, stim_parc_name, rec_parc_name):
    fig1, ax1 = plt.subplots(figsize=(10, 10))
    cax = ax1.matshow(mask, aspect='auto', cmap='Greys', interpolation='none', vmin=0, vmax=1)
    cbar = fig1.colorbar(cax)
    cbar.ax.tick_params(labelsize=14)
    ax1.set_title("Stim: " + stim_parc_name + "\n Rec: " + rec_parc_name, fontsize=16)
    ax1.set_xlabel('Recording Parcels', fontsize=16)
    ax1.set_ylabel('Stimulation Parcels', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.show()

    fig2, ax2 = plt.subplots(figsize=(10, 10))
    cax = ax2.matshow(CImasked, aspect='auto', vmin=0, vmax=0.15)
    cbar = fig2.colorbar(cax)
    cbar.set_label('Confidence Interval', fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    ax2.set_title("Stim: " + stim_parc_name + "\n Rec: " + rec_parc_name, fontsize=16)
    ax2.set_xlabel('Recording Parcels', fontsize=16)
    ax2.set_ylabel('Stimulation Parcels', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.show()

"""
Determine the percentage of brain coverage for intra- and inter-hemispheric parcels
based on the provided mask.
"""
def coverage(sym_flag, mask):
    rows, cols = mask.shape
    mid_row, mid_col = rows // 2, cols // 2
    q1 = mask[:mid_row, :mid_col]  # Top-left
    q2 = mask[:mid_row, mid_col:]  # Top-right
    q3 = mask[mid_row:, :mid_col]  # Bottom-left
    q4 = mask[mid_row:, mid_col:]  # Bottom-right
    if sym_flag[0] == False:
        intra = np.round(((np.sum(q1) + np.sum(q4)) / (q1.size + q4.size)), 3)
        inter = np.round(((np.sum(q2) + np.sum(q3)) / (q2.size + q3.size)), 3)
    else:
        if sym_flag[1] == "H":
            intra = np.round(((np.sum(q1) + np.sum(q3)) / (q1.size + q3.size)), 3)
            inter = np.round(((np.sum(q2) + np.sum(q4)) / (q2.size + q4.size)), 3)
        else:
            intra = np.round(((np.sum(q1) + np.sum(q2)) / (q1.size + q2.size)), 3)
            inter = np.round(((np.sum(q3) + np.sum(q4)) / (q3.size + q4.size)), 3)
    coverage = {
        "Intra Coverage": np.float64(intra * 100),
        "Inter Coverage": np.float64(inter * 100),
    }
    for quadrant, count in coverage.items():
        print(f"{quadrant}: {count} %")

"""
Plot a matrix and its associated brain maps, mapping matrix values onto brain parcels.
For efferent connectivity (efferent=True), one map is plotted per stimulation parcel; 
for afferent connectivity (efferent=False), one map is plotted per recording parcel.
"""
def overall_plot(measure, data_masked, surro, efferent,
                 plot_parc_name, stim_parc_name, rec_parc_name, meshdirname,
                 roi_stim, roi_rec, parcels_idxStim, parcels_idxRec,
                 output_directory_name, save, ytick=[]):
    if measure == "probability":
        vmin = surro
        vmax = np.nanmax(data_masked)
        vmax_l = vmax_list(data_masked,roi_stim[0], roi_rec[0], efferent)
    elif measure == "parcel_highestP":
        vmin = 0
        vmax = len(ytick) - 1
        vmax_l = [vmax]
    else:
        vmin = np.nanmin(data_masked)
        vmax = np.nanmax(data_masked)
        vmax_l = vmax_list(data_masked,roi_stim[0], roi_rec[0], efferent)
    #vmax=np.nanmax(vmax_l)
    vmax=vmax_l[0]
    if vmax > 0:
        # plot matrix
        plot_atlas(measure + "\n efferent: " + str(efferent), data_masked,
                   stim_parc_name, rec_parc_name, roi_stim[0], roi_rec[0],
                   vmax, vmin, output_directory_name, save, ytick=ytick)

        # plot brain maps from matrix
        if efferent == True:
            data = data_masked
            labels = roi_rec
            idxPlot = parcels_idxRec
            roi = roi_stim[0]
        else:
            data = np.array(data_masked).transpose()
            labels = roi_stim
            idxPlot = parcels_idxStim
            roi = roi_rec[0]
        for i in range(len(roi)):
            title = measure + "\n ROI: " + str(roi[i]) + "\n efferent: " + str(
                efferent) + " \n Stim: " + stim_parc_name + " \n Rec: " + rec_parc_name
            plot_brain_map(data[i, :], labels, roi[i], plot_parc_name,
                           meshdirname, output_directory_name, title, save,
                           idxPlot=idxPlot, measure=measure, vmax=vmax, vmin=vmin, ytick=ytick)
    return None

"""
Determine the maximum value among cortical parcels and a separate representative value for subcortical parcels ("SUBC" in label).
"""
def vmax_list(data_masked, roi_stim, roi_rec, efferent):
    if efferent == True:
        labels_cortex = roi_rec
        p = data_masked.copy()
    else:
        labels_cortex = roi_stim
        p = data_masked.transpose()
    indices = [i for i, val in enumerate(labels_cortex) if "SUBC" in val]
    subc_values = [p[0,i] for i in indices]
    row = p[0]
    other_values = [val for i, val in enumerate(row) if i not in indices]
    if subc_values == []:
        max_subc = 0
    else:
        max_subc = np.nanmax(subc_values)
    max_cort = np.nanmax(other_values)
    return [max_cort, max_subc]

"""
Plot a stimulation-to-recording parcel matrix and save the figure to file.
Rows correspond to stimulation parcels and columns to recording parcels.
"""
def plot_atlas(measure, prob_masked,
               stim_parc_name, rec_parc_name, roi_stim, roi_rec,
               vmax, vmin, path, save_flag, ytick=[]):
    fig2, ax2 = plt.subplots(figsize=(20, 20))
    ft = 14
    if "peak_delay" in measure:
        cmap = plt.cm.plasma
        #cmap = plt.cm.turbo_r
    else:
        cmap = plt.cm.plasma
        #cmap=plt.cm.turbo
    cmap = plt.cm.plasma.copy()
    cmap.set_under((0.8, 0.8, 0.8, 1.0))
    cax = ax2.matshow(prob_masked, aspect='auto', vmax=vmax, vmin=vmin, cmap=cmap)
    cbar = fig2.colorbar(cax, label=measure)
    cbar.ax.tick_params(labelsize=ft)
    cbar.set_label(measure, fontsize=ft)

    if ytick == []:
        ticks = np.linspace(vmin, vmax, 5)
        if "probability" in measure or "highest" in measure or "Confidence" in measure:
            ytick = [f"{t:.2f}" for t in ticks]
        else:
            ytick = [f"{round(t)}" for t in ticks]
    else:
        ticks = np.arange(vmin, vmax + 1)
    cbar.set_ticks(ticks)
    cbar.ax.set_yticklabels(ytick)

    ax2.set_title("Stim: " + stim_parc_name + " \n Rec: " + rec_parc_name, fontsize=ft)
    ax2.set_xlabel('Recording Parcels', fontsize=ft)
    ax2.set_ylabel('Stimulation Parcels', fontsize=ft)

    # Custom labels for rows and columns
    plt.xticks(fontsize=ft)
    plt.yticks(fontsize=ft)
    no = max(len(roi_stim), len(roi_rec))
    fontsize = max(2, min(8, int((8 - 0.02 * no) * 20/10)))
    yticks = np.arange(0, prob_masked.shape[0], 1)
    if isinstance(roi_stim, str):
        roi_stim = [roi_stim]
    plt.yticks(ticks=yticks, labels=roi_stim, fontsize=fontsize)
    xticks = np.arange(0, prob_masked.shape[1], 1)
    if isinstance(roi_rec, str):
        roi_rec = [roi_rec]
    plt.xticks(ticks=xticks, labels=roi_rec, rotation=90, fontsize=15)

    plt.tight_layout()
    if save_flag == True:
        safe_filename = re.sub(r'[<>:"/\\|?*\n\r\t]', '_', measure)
        plt.savefig(path + "Matrix_" + safe_filename + "_" + str(np.round(vmax,2)) + ".svg")
    plt.show()
    return None

"""
Plot brain maps and save the figures.
Cortical parcels visualized with MNE FreeSurfer meshes.
- Lausanne parcellation: subcortical parcels with ENIGMA templates.
- Other parcellations: standard 2D plot for subcortical parcels.
Connectivity type (efferent/afferent) determines whether recording or stimulation parcels are mapped.
"""
def plot_brain_map(mat_p, labels, roi, plot_parc_name,
                   meshdirname, path, title, save_flag,
                   idxPlot=[], measure="probability", vmax=1, vmin=0, ytick=[]):
    figL = mne.viz.create_3d_figure(size=(400, 400), bgcolor=(255, 255, 255), smooth_shading=None, handle=None,
                                    scene=True, show=False)
    figR = mne.viz.create_3d_figure(size=(400, 400), bgcolor=(255, 255, 255), smooth_shading=None, handle=None,
                                    scene=True, show=False)
    Brain = mne.viz.get_brain_class()
    brainL = Brain("", hemi='lh', surf='inflated', subjects_dir=meshdirname, figure=figL, cortex='#64646400')
    brainR = Brain("", hemi='rh', surf='inflated', subjects_dir=meshdirname, figure=figR, cortex='#64646400')
    labels_path = meshdirname + "/label/" + plot_parc_name

    fig1, axs = plt.subplots(4, 2, figsize=(23, 30))
    fig1.suptitle(title, fontsize=40)
    gridspec.GridSpec(4, 2, height_ratios=[1, 1, 1, 1])

    # costum plot parcellation
    mat_p_cost, labels_cost = costum_plot_parcellation(mat_p, labels, idxPlot=idxPlot)
    #print(labels)
    #print(labels_cost)

    # plot cortical
    plot_cortical(measure, mat_p_cost, labels_cost, labels_path,
                  brainL, brainR, ytick, vmin, vmax, fig1, axs)

    # plot subcortical
    axs[3, 0].axis('off')
    axs[3, 1].axis('off')
    axs[2, 0].axis('off')
    axs[2, 1].axis('off')
    if 'Lausanne' in plot_parc_name and len(labels[0])>2:
        plot_subcortical_enigma(measure, mat_p_cost, vmin, vmax, fig1, plot_parc_name)
    else:
        indexes = [i for i, item in enumerate(labels[0]) if "SUB" in item]
        if indexes == []:
            print("No subcortical parcels found")
        else:
            plot_subcortical_standard(measure, mat_p_cost, vmin, vmax, fig1, axs, labels_cost, ytick)
    fig1.tight_layout()

    if save_flag == True:
        safe_filename = re.sub(r'[<>:"/\\|?*\n\r\t]', '_', measure)
        safe_roi = re.sub(r'[*]', '_', roi)
        fig1.savefig(path + "BrainMap_" + safe_roi + "_" + safe_filename + "_" + str(np.round(vmax,2)) + ".svg")
    fig1.show()
    plt.close(fig1)
    brainL.close()
    brainR.close()
    return None

"""
Adapt the data array to the available meshes for visualization.
When parcels are merged, each mesh in the merged group is assigned the same value.
"""
def costum_plot_parcellation(mat_p, labels, idxPlot=[]):
    if idxPlot==[]:
        mat_p = mat_p
        labels = labels[0]
    else:
        labels = labels[1]
        old_p = mat_p
        new_p = []
        for l, ln in enumerate(labels):
            group_index = next((l_ for l_, group in enumerate(idxPlot) if l in group), None)
            if group_index is not None:
                new_p.append(old_p[group_index])
            else:
                new_p.append(np.nan)
        mat_p = np.array(new_p)
        labels = [roi.replace("_", " ") for roi in labels]
        labels = [roi.replace(" L", "_L") for roi in labels]
        labels = [roi.replace(" R", "_R") for roi in labels]
        #print(labels)
    return mat_p, labels

"""
Generate cortical surface plots showing lateral and medial views for the left and right brain hemispheres.
"""
def plot_cortical(measure, mat_p, labels, labels_path,
                  brainL, brainR, ytick, vmin, vmax, fig1, axs, colorbar='on', colormap='plasma'):
    # associate probability values with colors
    if "peak_delay" in measure:
        cmap = plt.get_cmap('plasma')
        #cmap = plt.get_cmap('turbo_r')
    else:
        cmap = plt.get_cmap('plasma')
        #cmap = plt.get_cmap('turbo')
    print(measure)
    colors = measureToColor(mat_p, vmin, vmax, measure, cmap)

    # associate colors with parcels
    colorToParcel(mat_p, labels, colors, brainL, brainR, labels_path)

    #lateral view
    brainR.show_view(view="lat")
    img_R_lat = brainR.screenshot(time_viewer=False)
    img_R_lat = img_R_lat[50:-50, :, :]
    brainL.show_view(view="lat")
    img_L_lat = brainL.screenshot(time_viewer=False)
    img_L_lat = img_L_lat[50:-50, :, :]
    axs[0, 0].imshow(img_L_lat)
    axs[0, 0].axis('off')
    axs[0, 1].imshow(img_R_lat)
    axs[0, 1].axis('off')
    # medial view
    brainR.show_view(view="med")
    img_R_med = brainR.screenshot(time_viewer=False)
    img_R_med = img_R_med[50:-50, :, :]
    brainL.show_view(view="med")
    img_L_med = brainL.screenshot(time_viewer=False)
    img_L_med = img_L_med[50:-50, :, :]
    axs[1, 0].imshow(img_L_med)
    axs[1, 0].axis('off')
    axs[1, 1].imshow(img_R_med)
    axs[1, 1].axis('off')

    # colorbar
    if colorbar=="on":
        #cm = plt.colormaps[colormap]
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        cb_ax = fig1.add_axes([0.45, 0.38, 0.025, 0.4])  # [left, bottom, width, height]
        cbar1 = fig1.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cb_ax, orientation='vertical')
        cbar1.ax.tick_params(labelsize=40)
        if ytick == []:
            ticks = np.linspace(vmin, vmax, 5)
            if "probability" in measure or "highest" in measure or "Confidence" in measure:
                ytick = [f"{t:.2f}" for t in ticks]
            else:
                ytick = [f"{round(t)}" for t in ticks]
        else:
            ticks = np.arange(vmin, vmax + 1)
        cbar1.set_ticks(ticks)
        cbar1.ax.set_yticklabels(ytick)
    return None

"""
Map numerical data values to colors for plotting.
Probabilities below surrogates appear in white.
Non-significant measures appear in black.
"""
def measureToColor(mat_p, vmin, vmax, measure, cmap):
    cmaplist = [cmap(i) for i in range(cmap.N)]
    #print(cmaplist)
    values = np.linspace(vmin, vmax, 255)
    colors = []
    for i, x in enumerate(mat_p):
        if "probability" in measure and x < vmin:
            colors.append((0.8, 0.8, 0.8, 1.0)) #white below surrogate
        elif (measure=="parcel_highestP" and x==-1) or (measure=="pvalue" and x>0.05):
            colors.append((0.2, 0.2, 0.2, 1.0)) #black if no statistical differences were found
        elif measure=="pvalue" and np.isnan(x):
            colors.append((0.392, 0.392, 0.392, 0.0)) #medium gray
        else:
            j = np.argmin(np.abs(values - x))
            color0=cmaplist[j]
            colors.append(color0)
    return colors

"""
Add color to meshes for visualization.
Each mesh is assigned a color based on the data value associated with it.
Mesh labels are saved inside mne_data folder with filenames following the convention:
- Left hemisphere: 'lh.<parcel_name>.label'
- Right hemisphere: 'rh.<parcel_name>.label'
where <parcel_name> corresponds to the name of each parcel.
"""
def colorToParcel(mat_p, labels, colors, brainL, brainR, labels_path):
    os.chdir(labels_path)
    for c,r in enumerate(labels):
        #print(r)
        if ('intra.' in r):
            lab=r.replace('intra.', "")
            lab= "lh." + lab + ".label"
            if os.path.exists(labels_path+"/"+lab):
                if not np.isnan(mat_p[c]):
                    brainL.add_label(lab, color=colors[c], alpha=colors[c][-1])
        elif ('_L' in r):
            lab = r.replace('_L', "")
            #lab = lab.replace(' ', "_")
            lab = "lh." + lab + ".label"
            if os.path.exists(labels_path+"/"+lab):
                if not np.isnan(mat_p[c]):
                    brainL.add_label(lab, color=colors[c], alpha=colors[c][-1])
        elif ('inter.' in r):
            lab = r.replace('inter.', "")
            lab = "rh." + lab + ".label"
            if os.path.exists(labels_path+"/"+lab):
                if not np.isnan(mat_p[c]):
                    brainR.add_label(lab, color=colors[c], alpha=colors[c][-1])
        elif ('_R' in r):
            lab = r.replace('_R', "")
            #lab = lab.replace(' ', "_")
            lab = "rh." + lab + ".label"
            if os.path.exists(labels_path+"/"+lab):
                if not np.isnan(mat_p[c]):
                    brainR.add_label(lab, color=colors[c], alpha=colors[c][-1])

"""
Visualize subcortical parcels by plotting ENIGMA subcortical meshes.
Parcels with values below surrogates or masked or non-significant are not displayed. 
"""
def plot_subcortical_enigma(measure, mat_p, vmin, vmax, fig1, plot_parc_name):
    if (plot_parc_name[:12] == "Lausanne2008"):
        if (plot_parc_name == "Lausanne2008-60"):
            add=0
            add2=0
        elif (plot_parc_name == "Lausanne2008-125"):
            add=54
            add2=105
        elif (plot_parc_name == "Lausanne2008-33"):
            add=-23
            add2=-46
        else:
            print("plot_subcortical function is not defined for this parcellation")
        # L_accumbens, L_amygdala, L_caudate, L_hippocampus,L_pallidun, L_putamen, L_thalamus, 0-6
        sub_cort = np.empty(14)
        sub_cort[:] = np.nan
        # sub_cort[0] = mat_p[62+add]
        sub_cort[1] = mat_p[63 + add]
        # sub_cort[2] = mat_p[61+add]
        sub_cort[3] = mat_p[60 + add]
        # sub_cort[4] = mat_p[58+add]
        # sub_cort[5] = mat_p[59+add]
        # sub_cort[6] = mat_p[57+add]
        # R_accumbens, R_amygdala, R_caudate, R_hippocampus, R_pallidun, R_putamen, R_thalamus, 7-13
        # sub_cort[7] = mat_p[126+add2]
        sub_cort[8] = mat_p[127 + add2]
        # sub_cort[9] = mat_p[125+add2]
        sub_cort[10] = mat_p[124 + add2]
        # sub_cort[11] = mat_p[122+add2]
        # sub_cort[12] = mat_p[123+add2]
        # sub_cort[13] = mat_p[121+add2]
    elif (plot_parc_name[:12] == "Lausanne2018"):
        if (plot_parc_name == "Lausanne2018-scale2"):
            add2=0
            add=0
        elif (plot_parc_name == "Lausanne2018-scale3"):
            add2=51
            add=102
        elif (plot_parc_name == "Lausanne2018-scale1"):
            add2=-23
            add=-46
        else:
            print("plot_subcortical function is not defined for this parcellation")
        # L_accumbens, L_amygdala, L_caudate, L_hippocampus,L_pallidun, L_putamen, L_thalamus, 0-6
        sub_cort = np.empty(14)
        sub_cort[:] = np.nan
        # sub_cort[0] = mat_p[137+add]
        sub_cort[1] = mat_p[138 + add]
        # sub_cort[2] = mat_p[134+add]
        sub_cort[3] = mat_p[139 + add]
        # sub_cort[4] = mat_p[136+add]
        # sub_cort[5] = mat_p[135+add]
        # R_accumbens, R_amygdala, R_caudate, R_hippocampus, R_pallidun, R_putamen, R_thalamus, 7-13
        # sub_cort[7] = mat_p[67+add2]
        sub_cort[8] = mat_p[68 + add2]
        # sub_cort[9] = mat_p[64+add2]
        sub_cort[10] = mat_p[69 + add2]
        # sub_cort[11] = mat_p[66+add2]
        # sub_cort[12] = mat_p[65+add2]
    else:
        print("plot_subcortical function is not defined for this parcellation")

    if "probability" in measure:
        sub_cort[sub_cort < vmin] = np.nan
    if measure == "parcel_highestP":
        sub_cort[sub_cort == -1] = np.nan
    if "peak_delay" in measure:
        cmap = 'plasma'
        #cmap = 'turbo_r'
    else:
        cmap='plasma'
        #cmap = 'turbo'
    fig_sc = plot_subcortical(array_name=sub_cort, size=(1000, 500), cmap=cmap, color_bar=False, ventricles=False,
                              color_range=(vmin, vmax), screenshot=True, filename= 'subcort.png')
    ax_custom = fig1.add_axes([0.12, 0.2, 0.8, 0.3])
    img_subcortical = Image.open(fig_sc)
    ax_custom.imshow(np.array(img_subcortical))
    ax_custom.axis('off')

"""
Create a 2D visualization of subcortical data.
"""
def plot_subcortical_standard(measure, mat_p, vmin, vmax, fig1, axs, labels, ytick=[]):
    axs[3, 0].remove()
    axs[3, 1].remove()
    axs2 = fig1.add_subplot(4, 2, (7, 8))
    indexes = [i for i, item in enumerate(labels) if "SUB" in item]
    mat_p=mat_p[indexes]
    labels = [labels[i] for i in indexes]
    axs2.plot(mat_p, 'o', markersize=16)
    axs2.set_facecolor('darkgray')
    axs2.grid(True, linewidth=3)

    axs2.set_xticks(np.arange(len(labels)))
    axs2.set_xticklabels(labels, rotation=90, ha='right', fontsize=14)
    axs2.set_xlabel("Parcels", fontsize=14)
    axs2.tick_params(axis='y', labelsize=14)
    if measure=="parcel_highestP":
        axs2.set_ylim(-1, vmax)
    elif measure=="directionality":
        axs2.set_ylim(-1, 1)
    else:
        axs2.set_ylim(0, vmax)
    axs2.axhline(y=vmin, color='lightgray', linestyle='--', linewidth=6)

    if ytick == []:
        ticks = np.linspace(vmin, vmax, 5)
        ytick = [f"{t:.2f}" for t in ticks]
    else:
        ticks = np.arange(vmin, vmax + 1)
    axs2.set_yticks(ticks)
    axs2.set_yticklabels(ytick)

"""
Analyze insula parcel dominance based on probability data.
This function performs pairwise statistical analysis between insula parcels
and constructs a dominance matrix (idx_all) to identify dominant insula parcels
connected to each brain region. It also generates figures to summarize the results.
Figures:
- Plot the highest probability insula parcels.
- Display the corresponding effect sizes.
- Show the dominance matrix (idx_all).
- Generate a bar plot of the percentage of dominant parcels per brain system
 (occipital, temporal, parietal, frontal, limbic, central). Only supported for Lausanne parcellation.
"""
def overall_stats(all_data_disc_masked, resolution, roi_stim, roi_rec, parcels_idxStim, parcels_idxRec,
                  plot_parc_name, stim_parc_name, rec_parc_name, hemis_symmetrize, efferent, surro, mask,
                  meshdirname, ytickStim, ytickRec, output_directory_name, save):
    if efferent == True:
        prob_masked = all_data_disc_masked[0].transpose()
        NTotal_masked = all_data_disc_masked[1].transpose()
        NValue_masked = all_data_disc_masked[2].transpose()
        ytick = ytickStim
    else:
        prob_masked = all_data_disc_masked[0]
        NTotal_masked = all_data_disc_masked[1]
        NValue_masked = all_data_disc_masked[2]
        ytick = ytickRec

    if hemis_symmetrize[0] == False:
        slcs = [slice(0, resolution),slice(resolution, None)]
    else:
        slcs = [slice(0, resolution)]

    for slc in slcs:
        pval_l, effect_l, idx_maxp_l, idx_all_l = stats(NTotal_masked[:, slc],NValue_masked[:, slc],
                                                        prob_masked[:, slc],resolution,surro)
        if efferent == True:
            idx_maxp_l = idx_maxp_l.transpose()
            effect_l = effect_l.transpose()
            idx_all_l = idx_all_l
            prob = prob_masked.transpose()
            prob = prob[slc]
            mask_l = mask[slc,:]
            roi_stim_lv1 = roi_stim[0][0] #[roi_stim[0][0], roi_stim[1][0]]
            roi_rec_lv1 = roi_rec
            roi_stim_l = roi_stim[0][slc]
            roi_rec_l = roi_rec[0]
            if hemis_symmetrize[0] == False:
                roi = roi_stim[0]
            else:
                roi_intra = ["intra." + item for item in roi_stim[0]]
                roi_inter = ["inter." + item for item in roi_stim[0]]
                roi = roi_intra + roi_inter
            rois_systems = roi_rec[0]
        else:
            idx_maxp_l = idx_maxp_l
            effect_l = effect_l
            idx_all_l = idx_all_l.transpose()
            prob = prob_masked[:, slc]
            mask_l = mask[:, slc]
            roi_rec_lv1 = roi_rec[0][0] #[roi_rec[0][0], roi_rec[1][0]]
            roi_stim_lv1 = roi_stim
            roi_rec_l = roi_rec[0][slc]
            roi_stim_l = roi_stim[0]
            if hemis_symmetrize[0] == False:
                roi = roi_rec[0]
            else:
                roi_intra = ["intra." + item for item in roi_rec[0]]
                roi_inter = ["inter." + item for item in roi_rec[0]]
                roi = roi_intra + roi_inter
            rois_systems = roi_stim[0]
        #"""
        overall_plot("parcel_highestP", idx_maxp_l, surro, efferent,
                     plot_parc_name, stim_parc_name, rec_parc_name, meshdirname,
                     roi_stim_lv1, roi_rec_lv1, parcels_idxStim, parcels_idxRec,
                     output_directory_name, save, ytick=ytick)
        overall_plot("effect_highestP", effect_l, surro, efferent,
                     plot_parc_name, stim_parc_name, rec_parc_name, meshdirname,
                     roi_stim_lv1, roi_rec_lv1, parcels_idxStim, parcels_idxRec,
                     output_directory_name, save, ytick=[])#"""
        plot_atlas("similar_highestP" + "\n efferent: " + str(efferent), np.where(mask_l, idx_all_l, np.nan),
                   stim_parc_name, rec_parc_name, roi_stim_l, roi_rec_l, 2, 0,
                   output_directory_name, save, ytick=["Not similar/No sig diff", "Similar to Highest", "Highest P"])
        if plot_parc_name == 'Lausanne2018-scale2':
            if efferent == True:
                idx_all_l = idx_all_l.copy()
                prob = prob_masked.copy()
            else:
                idx_all_l = idx_all_l.transpose()
                prob = prob_masked.transpose()
            ROI_dominance_by_brainSystems(idx_all_l, roi, save, "ROI_Dominance", output_directory_name)
            #results = HighestDominantP_within_brainSystem(idx_all_l, prob, roi, rois_systems)
        else:
            print("ROI_dominance_by_brainSystems function is not prepared for this plot_parc_name")

"""
Analyze insula parcel dominance based on probability data.
This function performs pairwise statistical analysis between insula parcels
and constructs a dominance matrix (idx_all) to identify dominant insula parcels
connected to each brain region.
"""
def stats(NTotal, NValue, prob, resolution, surro):
    pval_posthoc=[]
    effect_posthoc=[]
    pval_l1 = []
    effect_l1 = []
    idx_maxp = []
    idx_all = np.zeros((resolution, prob.shape[0]))
    pairs = list(combinations(range(resolution), 2))
    #runs for each cortex parcel
    for j in range(prob.shape[0]):
        #finds the p-value and effect after comparing each pair of insula parcels
        #considers binomial distribution and independent proportions
        pval_pairs = []
        eff_pairs = []
        vals = prob[j, 0:len(prob[j, :])]
        not_nan_count = np.sum(~np.isnan(vals))
        count_above_surro = np.sum(vals > surro)
        if not_nan_count < 2 or count_above_surro < 2:
            pval_pairs.extend([np.nan] * len(pairs))
            eff_pairs.extend([np.nan] * len(pairs))
            pval_l1.append(np.nan)
            effect_l1.append(np.nan)
            idx_maxp.append(np.nan)
            idx_all[:, j] = np.nan
        else:
            pval_pairs, eff_pairs=sigtest_pairs(prob, NValue, NTotal, pairs, j, surro)
            pval_l1, effect_l1, idx_maxp, idx_all=sigtest_DominantParcel(prob, pval_pairs, eff_pairs, pairs, j,
                                                                         pval_l1, effect_l1, idx_maxp, idx_all)
        pval_posthoc.append(np.array(pval_pairs))
        effect_posthoc.append(np.array(eff_pairs))
    pval_l1 = np.array(pval_l1)
    pval_l1 = pval_l1[:, np.newaxis]
    effect_l1 = np.array(effect_l1)
    effect_l1 = effect_l1[:, np.newaxis]
    pval_posthoc = np.array(pval_posthoc)
    effect_posthoc = np.array(effect_posthoc)
    idx_maxp = np.array(idx_maxp)
    idx_maxp = idx_maxp[:, np.newaxis]
    idx_all = np.array(idx_all)
    return pval_l1, effect_l1, idx_maxp, idx_all

"""
Analyze pairwise differences between insula parcels (only supports probability data).
- Compute proportion z-tests for each pair of parcels (NaN is given if both are masked or below surrogates).
- Calculate effect size according to Cohen's h.
- Apply Benjamini-Hochberg correction to control false discovery rate of p-values.
"""
def sigtest_pairs(prob,NValue,NTotal,pairs,j,surro):
    pval_pairs = []
    eff_pairs = []
    for i, k in pairs:
        if ~np.isnan(prob[j, i]) and ~np.isnan(prob[j, k]):
            if prob[j, i] < surro and prob[j, k] < surro:
                pval_pairs.append(np.nan)
                eff_pairs.append(np.nan)
            else:
                count = np.array([NValue[j, i], NValue[j, k]]).T
                nobs = np.array([NTotal[j, i], NTotal[j, k]]).T
                z_stat, p_val_pair = proportions_ztest(count, nobs)
                #eff = abs((prob[j, i] - prob[j, k])/np.nanmax([prob[j, i],prob[j, k]]))
                eff = abs(2 * (np.arcsin(np.sqrt(prob[j, i])) - np.arcsin(np.sqrt(prob[j, k]))))
                pval_pairs.append(p_val_pair)
                eff_pairs.append(eff)
        else:
            pval_pairs.append(np.nan)
            eff_pairs.append(np.nan)
    pval_pairs = np.array(pval_pairs)
    eff_pairs = np.array(eff_pairs)
    valid_mask = ~np.isnan(pval_pairs)
    pvals_clean = pval_pairs[valid_mask]
    fdr_corrected_clean = smm.fdrcorrection(pvals_clean)[1]  # Benjamini/Hochberg correction
    pval_pairs = np.full_like(pval_pairs, np.nan, dtype=float)
    pval_pairs[valid_mask] = fdr_corrected_clean
    return pval_pairs, eff_pairs

"""
Create a dominance matrix for insula parcels based on statistical comparisons: idx_all.
For each brain region:
1. Identify the insula parcel with the highest probability.
2. Determine which insula parcels are statistically similar to the highest.
3. Encode a matrix (idx_all) where:
   - 2: highest probability insula parcel
   - 1: insula parcels statistically similar to the highest
   - 0: remaining insula parcels
4. If no differences are found between any insula parcel pairs, all entries are 0.
"""
def sigtest_DominantParcel(prob, pval_pairs, eff_pairs, pairs, j, pval_l1, effect_l1, idx_maxp, idx_all):
    if np.any(pval_pairs < 0.05):
        idx = np.nanargmax(prob[j, :])
        n = -1
        flagdif = False
        pvalue_temp = []
        effect_temp = []
        for i, k in pairs:
            n = n + 1
            if ((i == idx) or (k == idx)) and (pval_pairs[n] < 0.05):
                flagdif = True
                pvalue_temp.append(pval_pairs[n])
                effect_temp.append(eff_pairs[n])
        if flagdif==True:
            idx_maxp.append(idx)
            idx_all[idx, j] = 2
            effect_l1.append(np.nanmax(effect_temp))
            pval_l1.append(pvalue_temp[np.nanargmin(effect_temp)])
            n = -1
            for i, k in pairs:
                n = n + 1
                if (i == idx) and (pval_pairs[n] > 0.05):
                    idx_all[k, j] = 1
                elif (k == idx) and (pval_pairs[n] > 0.05):
                    idx_all[i, j] = 1
        else:
            pval_l1.append(np.nanmin(pval_pairs))
            effect_l1.append(np.nan)
            idx_maxp.append(-1) #all similar to the highest prob
            idx_all[:, j] = 0
    else:
        pval_l1.append(np.nanmin(pval_pairs))
        effect_l1.append(np.nan)
        idx_maxp.append(-1)  # all similar to the highest prob
        idx_all[:, j] = 0
    return pval_l1, effect_l1, idx_maxp, idx_all

"""
Compute and visualize dominant insula parcels proportions per brain system.
- Supported brain systems: occipital, temporal, parietal, frontal, limbic, central.
- Dominant parcels are identified from the dominance matrix as:
    - Highest-probability insula parcel
    - Insula parcels statistically similar to the highest
- Calculates the percentage of dominant insular connections in each system.
- Generates a bar plot showing percentages for each insula parcel and brain system.
- Only supported for Lausanne parcellation.
"""
def ROI_dominance_by_brainSystems(idx_all, roi, save_flag, title, path):
    idx_all[idx_all == 2] = 1
    # Define region order
    regions = ["occipital", "temporal", "parietal", "frontal", "limbic", "central"]
    atlas = 'Lausanne2018-scale2' #'Lausanne2008-60'

    # Left hemisphere
    brainSystemsL = idx[atlas][6][:6]
    left_vals = {
        r: np.nansum(idx_all[:, np.array(brainSystemsL[i])], axis=1) / len(brainSystemsL[i]) * 100
        for i, r in enumerate(regions)
    }
    # Right hemisphere
    brainSystemsR = idx[atlas][6][6:]
    right_vals = {
        r: np.nansum(idx_all[:, np.array(brainSystemsR[i])], axis=1) / len(brainSystemsR[i]) * 100
        for i, r in enumerate(regions)
    }
    # Combine hemispheres
    combined = {r: np.concatenate([left_vals[r], right_vals[r]]) for r in regions}

    ftsize = 20
    fig, axs = plt.subplots(1,1,figsize=(20, 14))
    cmap = plt.cm.plasma
    colors = [cmap(i) for i in np.linspace(0, 1, len(regions))]
    bottom = np.zeros_like(combined[regions[0]], dtype=float)
    for label, color in zip(regions, reversed(colors)):
        vals = combined[label]
        axs.bar(roi, vals, bottom=bottom, label=label.capitalize(), color=color)
        bottom += vals
    axs.set_ylabel("Cumulative % Dominant Regions", fontsize=ftsize)
    axs.tick_params(axis='y', labelsize=ftsize)
    axs.set_xlabel("ROIs", labelpad=40, fontsize=ftsize)
    axs.set_xticklabels(roi, fontsize=ftsize, rotation=90)
    plt.legend(fontsize=ftsize)

    text_colors = ["black","black","black","white","white","white"]
    for i in range(len(roi)):
        bottom = 0
        c=-1
        for region in regions:
            val = combined[region][i]
            c = c + 1
            if int(val) != 0:
                axs.text(roi[i], bottom + val / 2, str(int(val)),
                        ha='center', va='center', color=text_colors[c], fontsize=ftsize)
            bottom += val

    plt.tight_layout()
    if save_flag == True:
        safe_title = re.sub(r'[*]', '_', title)
        plt.savefig(path + safe_title + ".svg")
    plt.show()

def HighestDominantP_within_brainSystem(idx_all, prob, roi, rois_systems):
    idx_all[idx_all == 2] = 1
    prob[idx_all != 1] = np.nan
    brainSystemsL = idx['Lausanne2018-scale2'][6][0:6] #'Lausanne2008-60'

    system_order = ["occipital", "temporal", "parietal", "frontal", "limbic", "central"]
    results = {name: [] for name in system_order}

    for i in range(len(roi)//2):
        for sys_name, sys_idx in zip(system_order, brainSystemsL):
            cols = np.array(sys_idx)
            if np.isnan(prob[i, cols]).all():
                results[sys_name].append(0)
            else:
                best_idx = np.nanargmax(prob[i, cols], axis=0)
                results[sys_name].append(rois_systems[cols[best_idx]])
    return results

#%%
if __name__ == "__main__":
    main(filter_default)