# -*- coding: utf-8 -*-
"""
Created on Wed May 14 10:20 2025

@author: Cristiana Pinheiro
"""
import sys
source_path = 'C:/Users/neuro/Desktop/Cristiana/'
sys.path.append(source_path+'ft-insula/data/ENIGMA')
print('You need to change these paths during first use')
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mne
import numpy as np
import re
from enigmatoolbox.plotting import plot_subcortical
from PIL import Image
import statsmodels.stats.multitest as smm
from statsmodels.stats.proportion import proportions_ztest
from scipy.stats import spearmanr
from analyze_atlas import parcels_labels, read_atlas_data, masking, colorToParcel, surro_dict, select_parcel_groups, generate_paths, costum_plot_parcellation, plot_subcortical_standard, overall_plot
import matplotlib.colors as mcolors
import matplotlib as mpl

def main(title="EffvsAff_", merged = [[True, False], [False, True]], hemis_symmetrize = [[True, "H"], [True, "V"]],
         resolution_stim = [1, np.nan], resolution_rec = [np.nan, 1],
         stim_parc_name = ['MNI-insula',"Lausanne2018-scale2"], rec_parc_name = ["Lausanne2018-scale2",'MNI-insula'],
         plot_parc_name = "Lausanne2018-scale2", time_window = [["max", 100], ["max", 100]], ntotal_min = [50,50], save = False):
    # Variables
    # 'MNI-JulichBrain-3.0' 'MNI-insula' 'Lausanne2008-125' 'Lausanne2008-60' 'Lausanne2008-33'
    # title = LeftvsRight_     EffvsAff_   AllMatter_vs_Gray_
    # merged = stim, rec ; stim, rec
    trial1 = "INSULA/" + stim_parc_name[0] + "/" + str(resolution_stim[0]) + "/" + rec_parc_name[0] + "/" + str(
        resolution_rec[0]) + "/efferent_True"
    trial2 = "INSULA/" + stim_parc_name[1] + "/" + str(resolution_stim[1]) + "/" + rec_parc_name[1] + "/" + str(
        resolution_rec[1]) + "/efferent_False"
    surro = surro_dict[time_window[0][1]]
    z = "5"
    se_dir = source_path+'/ft-insula/data/real_data/se_v3'
    # __________________________________________________________________________________________________________________

    # select parcels groups
    parcels_idxStim0, ytickStim0, parcels_idxRec0, ytickRec0 = select_parcel_groups(merged[0],
                                                                                stim_parc_name[0], rec_parc_name[0],
                                                                                resolution_stim[0], resolution_rec[0])
    parcels_idxStim1, ytickStim1, parcels_idxRec1, ytickRec1 = select_parcel_groups(merged[1],
                                                                                    stim_parc_name[1], rec_parc_name[1],
                                                                                    resolution_stim[1], resolution_rec[1])
    parcels_idxStim = [parcels_idxStim0, parcels_idxStim1]
    ytickStim = [ytickStim0, ytickStim1]
    parcels_idxRec = [parcels_idxRec0, parcels_idxRec1]
    ytickRec = [ytickRec0, ytickRec1]

    output_directory_name1 = source_path + "ft-insula/data/results/" + trial1 + "/"
    print("output_directory_name1", output_directory_name1)
    output_directory_name2 = source_path + "ft-insula/data/results/" + trial2 + "/"
    print("output_directory_name2", output_directory_name2)
    output_directory_name = [output_directory_name1, output_directory_name2]
    if os.path.isdir(output_directory_name1) & os.path.isdir(output_directory_name2):
        print("Available data will be used")
        compare(se_dir, time_window, z, ntotal_min, surro, resolution_stim, resolution_rec,
                plot_parc_name, stim_parc_name, rec_parc_name, merged, hemis_symmetrize, title,
                parcels_idxStim, parcels_idxRec, ytickStim, ytickRec,
                output_directory_name, save)
    else:
        print("No data available")

# FUNCTIONS
#______________________________________________________________________________________________________________________

def compare(se_dir, time_window, z, ntotal_min, surro, resolution_stim, resolution_rec,
            plot_parc_name, stim_parc_name, rec_parc_name, merged, hemis_symmetrize, title,
            parcels_idxStim, parcels_idxRec, ytickStim, ytickRec,
            output_directory_name, save):
    # paths
    # mesh path
    if plot_parc_name == "MNI-JulichBrain-3.0":
        meshdirname =  source_path + "/ft-insula/data/mne_data/MNE-sample-data/subjects/fsaverage"
    elif plot_parc_name[:8] == "Lausanne":
        meshdirname = source_path + "/ft-insula/data/mne_data/MNE-sample-data/subjects/cvs_avg35_inMNI152"
    else:
        print("meshdirname is not defined for this plot parcellation")
    # data path
    paths_discVAR0, paths_contVAR0 = generate_paths(output_directory_name[0], time_window[0], z, ntotal_min[0])
    paths_discVAR1, paths_contVAR1 = generate_paths(output_directory_name[1], time_window[1], z, ntotal_min[1])
    # order of measures_disc and measures_cont needs to match order of paths_discVAR and paths_contVAR, respectively
    measures_disc = ["probability", "NTotal", "NValue", "patients", "implantations"]
    measures_cont = ["peak_delay", "ampl"]
    features_cont = ["quantile_0.5", "mean", "std", "quantile_0.25", "quantile_0.75"]

    # read parcels names
    roi_stim0 = parcels_labels(se_dir, stim_parc_name[0], merged[0][0], hemis_symmetrize[0], False, idx=parcels_idxStim[0])
    roi_rec0 = parcels_labels(se_dir, rec_parc_name[0], merged[0][1], hemis_symmetrize[0], True, idx=parcels_idxRec[0])
    roi_stim1 = parcels_labels(se_dir, stim_parc_name[1], merged[1][0], hemis_symmetrize[1], False, idx=parcels_idxStim[1])
    roi_rec1 = parcels_labels(se_dir, rec_parc_name[1], merged[1][1], hemis_symmetrize[1], True, idx=parcels_idxRec[1])

    # read atlas data
    all_data_disc0 = []
    all_data_disc1 = []
    for i, j in enumerate(measures_disc):
        print(paths_discVAR1[i])
        data = read_atlas_data(paths_discVAR0[i], hemis_symmetrize[0])
        all_data_disc0.append(data)
        data = read_atlas_data(paths_discVAR1[i], hemis_symmetrize[1])
        all_data_disc1.append(data)
    all_data_cont0 = []
    all_data_cont1 = []
    #"""
    for i, j in enumerate(measures_cont):
        for k in features_cont:
            data = read_atlas_data(paths_contVAR0[i] + "nan" + k + ".txt", hemis_symmetrize[0])
            all_data_cont0.append(data)
            data = read_atlas_data(paths_contVAR1[i] + "nan" + k + ".txt", hemis_symmetrize[1])
            all_data_cont1.append(data)#"""

    # mask data
    all_data_disc_masked0, all_data_cont_masked0, mask0, CI0 = masking(all_data_disc0, all_data_cont0, stim_parc_name[0], rec_parc_name[0], hemis_symmetrize[0], ntotal_min[0])
    all_data_disc_masked1, all_data_cont_masked1, mask1, CI1 = masking(all_data_disc1, all_data_cont1, stim_parc_name[1], rec_parc_name[1], hemis_symmetrize[1], ntotal_min[1])

    # plots
    if title == "EffvsAff_":
        efferent = [True, False]
    else:
        efferent = [True, True]
    for i, j in enumerate(measures_disc):
        if (i == 0):  # "probability","NTotal","NValue","patients","implantations"
            #"""
            data = all_data_disc_masked0[i]
            overall_plot(j, data, surro, efferent[0],
                         plot_parc_name, stim_parc_name[0], rec_parc_name[0], meshdirname,
                         roi_stim0, roi_rec0, parcels_idxStim[0], parcels_idxRec[0],
                         output_directory_name[0], save)
            data = all_data_disc_masked1[i]
            overall_plot(j, data, surro, efferent[1],
                         plot_parc_name, stim_parc_name[1], rec_parc_name[1], meshdirname,
                         roi_stim1, roi_rec1, parcels_idxStim[1], parcels_idxRec[1],
                         output_directory_name[1], save)#"""
    #"""
    w = 0
    for i, j in enumerate(measures_cont):  # "peak_delay","ampl","speed"
        for k2, k in enumerate(features_cont):  # "quantile_0.5", "mean", "std", "quantile_0.25", "quantile_0.75"
            if (i == 1) and (k2 == 0):
                data = all_data_cont_masked0[w]
                overall_plot(j + "_" + k, data, surro, efferent[0],
                             plot_parc_name, stim_parc_name[0], rec_parc_name[0], meshdirname,
                             roi_stim0, roi_rec0, parcels_idxStim[0], parcels_idxRec[0],
                             output_directory_name[0], save)
                data = all_data_cont_masked1[w]
                overall_plot(j + "_" + k, data, surro, efferent[1],
                             plot_parc_name, stim_parc_name[1], rec_parc_name[1], meshdirname,
                             roi_stim1, roi_rec1, parcels_idxStim[1], parcels_idxRec[1],
                             output_directory_name[1], save)
            w = w + 1#"""

    #correlations
    """
    sig1 = all_data_disc_masked0[0][0,:] #all_data_disc0[0][0,:]    all_data_disc_masked0[0][0,:]  all_data_cont0[5][0,:]  all_data_cont0[0][0,:]
    sig2 = all_data_disc_masked1[0][0,:] #all_data_cont1[5][0,:]   all_data_cont1[0][0,:]   all_data_disc_masked1[0][0,:]
    corrv2, p_value_corrv2 = spearmanr(sig1, sig2, nan_policy='omit')
    print("Correlation: r and p-value")
    print(corrv2)
    print(p_value_corrv2)"""

    # statistical differences between stimulation parcels
    #"""
    all_data_disc_masked = [all_data_disc_masked0, all_data_disc_masked1]
    roi_stim = [roi_stim0, roi_stim1]
    roi_rec = [roi_rec0, roi_rec1]
    mask = [mask0, mask1]
    overall_stats_compare(all_data_disc_masked, resolution_stim[0], resolution_rec[0], roi_stim, roi_rec, parcels_idxStim, parcels_idxRec,
                  plot_parc_name, title, surro, meshdirname, output_directory_name, save)#"""
    return

"""
Analyze pairwise differences between two atlases and visualize directionality of probability.
This function performs statistical comparisons between two parcellated atlases 
(probability data only) and generates brain maps to visualize directional measure between atlas pairs.
Types of comparisons:
- Efferent vs afferent atlases.
- Right vs left hemispheres (comparing the same two atlases).
- Atlas generated with grey and white matter contacts vs atlas with only grey matter contacts.
"""
def overall_stats_compare(all_data_disc_masked, resolution_stim, resolution_rec, roi_stim, roi_rec, parcels_idxStim, parcels_idxRec,
                  plot_parc_name, title, surro, meshdirname, output_directory_name, save):
    corr_ = []
    red_blue = mcolors.LinearSegmentedColormap.from_list("red_blue", ["red", "blue"])
    mpl.colormaps.register(cmap=red_blue, name="red_blue")
    if title == "LeftvsRight_":
        range_ = int(len(roi_stim[0][0])/2)
    else:
        range_ = len(roi_stim[0][0])
    for i in range(range_):
        if title == "EffvsAff_":
            p0 = all_data_disc_masked[0][0][i, :]
            p1 = all_data_disc_masked[1][0][:, i]
            nt0 = all_data_disc_masked[0][1][i, :]
            nt1 = all_data_disc_masked[1][1][:, i]
            nv0 = all_data_disc_masked[0][2][i, :]
            nv1 = all_data_disc_masked[1][2][:, i]
            roi = roi_stim[0][0]
            labels = roi_rec[0]
            idxPlot = parcels_idxRec[0]
        elif title == "LeftvsRight_":
            resolution_rec = int(len(roi_rec[0][0])/2)
            p0 = all_data_disc_masked[0][0][i, :]
            p1 = all_data_disc_masked[1][0][resolution_stim + i, :]
            p1 =  np.concatenate([p1[resolution_rec:], p1[:resolution_rec]])
            nt0 = all_data_disc_masked[0][1][i, :]
            nt1 = all_data_disc_masked[1][1][resolution_stim + i, :]
            nt1 = np.concatenate([nt1[resolution_rec:], nt1[:resolution_rec]])
            nv0 = all_data_disc_masked[0][2][i, :]
            nv1 = all_data_disc_masked[1][2][resolution_stim + i, :]
            nv1 = np.concatenate([nv1[resolution_rec:], nv1[:resolution_rec]])
            roi = roi_stim[0][0]
            labels = roi_rec[0]
            idxPlot = parcels_idxRec[0]
        elif title == "AllMatter_vs_Gray_":
            p0 = all_data_disc_masked[0][0][i, :]
            p1 = all_data_disc_masked[1][0][i, :]
            nt0 = all_data_disc_masked[0][1][i, :]
            nt1 = all_data_disc_masked[1][1][i, :]
            nv0 = all_data_disc_masked[0][2][i, :]
            nv1 = all_data_disc_masked[1][2][i, :]
            roi = roi_stim[0][0]
            labels = roi_rec[0]
            idxPlot = parcels_idxRec[0]
        pval_pairs, eff_pairs, corr = stats_between_atlas(p0, p1, nt0, nt1, nv0, nv1, surro)
        corr_.append(corr)
        title2 = "Directionality \n ROI: " + str(roi[i]) + " \n Comparison: " + title
        plot_brain_map(pval_pairs, eff_pairs, labels, roi[i], plot_parc_name,
                       meshdirname, output_directory_name[0], title2, save,
                       idxPlot=idxPlot, measure="directionality", vmax=1, vmin=-1, ytick=[])
    print("Correlation: Mean and Std")
    print(np.nanmean(corr_))
    print(np.nanstd(corr_))

"""
Analyze pairwise differences between two atlases (only supports probability data).
- Compute proportion z-tests for each pair of parcels (NaN is given if both are masked or below surrogates).
- Calculate effect size as a measure of directionality. A measure of directionality was computed by 
subtracting the probabilities being compared corresponding to the same parcel pair and 
then normalized (-1:+1) this subtraction by the maximum probability.   
- Apply Benjamini-Hochberg correction to control false discovery rate of p-values.
"""
def stats_between_atlas(prob1, prob2, NTotal1, NTotal2, NValue1, NValue2, surro):
    pval_pairs=[]
    eff_pairs=[]
    for j in range(NTotal1.shape[0]):
        if ~np.isnan(prob1[j]) & ~np.isnan(prob2[j]):
            if prob1[j] < surro and prob2[j] < surro:
                pval_pairs.append(np.nan)
                eff_pairs.append(np.nan)
            else:
                count = np.array([NValue1[j], NValue2[j]]).T
                nobs = np.array([NTotal1[j], NTotal2[j]]).T
                z_stat, p_val_pair = proportions_ztest(count, nobs)
                #eff=2 * np.arcsin(np.sqrt(prob1[j])) - 2 * np.arcsin(np.sqrt(prob2[j])) #ref: Cohen (1988), doi: 10.4300/JGME-D-12-00156.1
                eff=(prob1[j]-prob2[j])/np.nanmax([prob1[j],prob2[j]])
                pval_pairs.append(p_val_pair)
                eff_pairs.append(eff)
        else:
            pval_pairs.append(np.nan)
            eff_pairs.append(np.nan)
    pval_pairs = np.array(pval_pairs, dtype=float)
    eff_pairs = np.array(eff_pairs, dtype=float)
    valid_mask = ~np.isnan(pval_pairs)
    pvals_clean = pval_pairs[valid_mask]
    fdr_corrected_clean = smm.fdrcorrection(pvals_clean)[1] #Benjamini/Hochberg correction
    pval_pairs = np.full_like(pval_pairs, np.nan, dtype=float)
    pval_pairs[valid_mask] = fdr_corrected_clean
    reciprocal_perc=(np.sum(pvals_clean>0.05)/(len(pvals_clean)))*100
    print(f"% reciprocal parcel connectivity: {reciprocal_perc:.2f}")
    corr, p_value_corr = spearmanr(prob1, prob2, nan_policy='omit')
    print(f"Spearman correlation: {corr:.4f}")
    print(f"Spearman p-value: {p_value_corr}")
    return pval_pairs, eff_pairs, corr

"""
Plot directionality brain maps and save the figures.
Cortical parcels visualized with MNE FreeSurfer meshes.
- Lausanne parcellation: subcortical parcels with ENIGMA templates.
- Other parcellations: standard 2D plot for subcortical parcels.
"""
def plot_brain_map(pvalue, effect, labels, roi, plot_parc_name,
                   meshdirname, path, title, save_flag,
                   idxPlot=[], measure="directionality", vmax=1, vmin=-1, ytick=[]):
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
    pval_cost, labels_cost = costum_plot_parcellation(pvalue, labels, idxPlot=idxPlot)
    eff_cost, labels_cost = costum_plot_parcellation(effect, labels, idxPlot=idxPlot)

    # plot cortical
    plot_cortical(measure, eff_cost, pval_cost, labels_cost, labels_path,
                  brainL, brainR, ytick, vmin, vmax, fig1, axs, colorbar='on')

    # plot subcortical
    axs[3, 0].axis('off')
    axs[3, 1].axis('off')
    axs[2, 0].axis('off')
    axs[2, 1].axis('off')
    if 'Lausanne' in plot_parc_name and len(labels[0])>2:
        plot_subcortical_enigma(eff_cost, pval_cost, vmin, vmax, fig1, plot_parc_name)
    else:
        indexes = [i for i, item in enumerate(labels[0]) if "SUB" in item]
        if indexes == []:
            print("No subcortical parcels found")
        else:
            eff_cost[pval_cost > 0.5] = np.nan
            plot_subcortical_standard(measure, eff_cost, vmin, vmax, fig1, axs, labels_cost, ytick)
    fig1.tight_layout()

    if save_flag == True:
        safe_filename = re.sub(r'[<>:"/\\|?*\n\r\t]', '_', measure)
        safe_roi = re.sub(r'[*]', '_', roi)
        fig1.savefig(path + safe_roi + "_" + safe_filename + ".svg")
    fig1.show()
    plt.close(fig1)
    brainL.close()
    brainR.close()
    return None

"""
Visualize subcortical parcels by plotting ENIGMA subcortical meshes.
Parcels with values masked or non-significant are not displayed. 
"""
def plot_subcortical_enigma(eff, pval, vmin, vmax, fig1, plot_parc_name):
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
        sub_cort = np.empty((14, 2));
        sub_cort[:] = np.nan
        #sub_cort[0] = [pval[62 + add], eff[62 + add]]
        sub_cort[1] = [pval[63 + add], eff[63 + add]]
        #sub_cort[2] = [pval[61 + add], eff[61 + add]]
        sub_cort[3] = [pval[60 + add], eff[60 + add]]
        #sub_cort[4] = [pval[58 + add], eff[58 + add]]
        #sub_cort[5] = [pval[59 + add], eff[59 + add]]
        #sub_cort[6] = [pval[57 + add], eff[57 + add]]
        # R_accumbens, R_amygdala, R_caudate, R_hippocampus, R_pallidun, R_putamen, R_thalamus, 7-13
        #sub_cort[7] = [pval[126 + add2], eff[126 + add2]]
        sub_cort[8] = [pval[127 + add2], eff[127 + add2]]
        #sub_cort[9] = [pval[125 + add2], eff[125 + add2]]
        sub_cort[10] = [pval[124 + add2], eff[124 + add2]]
        #sub_cort[11] = [pval[122 + add2], eff[122 + add2]]
        #sub_cort[12] = [pval[123 + add2], eff[123 + add2]]
        #sub_cort[13] = [pval[121 + add2], eff[121 + add2]]
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
        sub_cort = np.empty((14, 2));
        sub_cort[:] = np.nan
        # sub_cort[0] = [pval[137 + add], eff[137 + add]]
        sub_cort[1] = [pval[138 + add], eff[138 + add]]
        # sub_cort[2] = [pval[134 + add], eff[134 + add]]
        sub_cort[3] = [pval[139 + add], eff[139 + add]]
        # sub_cort[4] = [pval[136 + add], eff[136 + add]]
        # sub_cort[5] = [pval[135 + add], eff[135 + add]]
        # R_accumbens, R_amygdala, R_caudate, R_hippocampus, R_pallidun, R_putamen, R_thalamus, 7-13
        # sub_cort[7] = [pval[67 + add2], eff[67 + add2]]
        sub_cort[8] = [pval[68 + add2], eff[68 + add2]]
        # sub_cort[9] = [pval[64 + add2], eff[64 + add2]]
        sub_cort[10] = [pval[69 + add2], eff[69 + add2]]
        # sub_cort[11] = [pval[66 + add2], eff[66 + add2]]
        # sub_cort[12] = [pval[65 + add2], eff[65 + add2]]
    else:
        print("plot_subcortical function is not defined for this parcellation")

    scort = np.empty(14)
    scort[:] = np.nan
    values = np.linspace(vmin, vmax, 255)
    for i in range(len(scort)):
        if sub_cort[i][0] > 0.05:
            scort[i] = np.nan
        elif np.isnan(sub_cort[i][0]):
            scort[i] = np.nan
        else:
            scort[i] = values[np.argmin(np.abs(values - sub_cort[i][1]))]
    fig_sc = plot_subcortical(array_name=scort, size=(1000, 500), cmap="seismic_r", color_bar=False, ventricles=False,
                              color_range=(vmin, vmax), screenshot=True, filename= 'subcort.png',
                              nan_color=(0.2, 0.2, 0.2, 0.0), background=(0.392, 0.392, 0.392))
    ax_custom = fig1.add_axes([0.12, 0.2, 0.8, 0.3])
    img_subcortical = Image.open(fig_sc)
    ax_custom.imshow(np.array(img_subcortical))
    ax_custom.axis('off')
    return None

"""
Generate cortical surface plots showing lateral and medial views for the left and right brain hemispheres.
"""
def plot_cortical(measure, eff, pval, labels, labels_path,
                  brainL, brainR, ytick, vmin, vmax, fig1, axs, colorbar='on'):
    # associate probability values with colors
    cmap = mcolors.LinearSegmentedColormap.from_list("red_blue", ["red", "blue"])
    cmap = plt.cm.seismic_r
    colors = measureToColor(eff, pval, vmin, vmax, measure, cmap)

    # associate colors with parcels
    colorToParcel(pval, labels, colors, brainL, brainR, labels_path)

    #lateral view
    brainR.show_view(view="lat")
    img_R_lat = brainR.screenshot(time_viewer=True)
    img_R_lat = img_R_lat[50:-50, :, :]
    brainL.show_view(view="lat")
    img_L_lat = brainL.screenshot(time_viewer=True)
    img_L_lat = img_L_lat[50:-50, :, :]
    axs[0, 0].imshow(img_L_lat)
    axs[0, 0].axis('off')
    axs[0, 1].imshow(img_R_lat)
    axs[0, 1].axis('off')
    # medial view
    brainR.show_view(view="med")
    img_R_med = brainR.screenshot(time_viewer=True)
    img_R_med = img_R_med[50:-50, :, :]
    brainL.show_view(view="med")
    img_L_med = brainL.screenshot(time_viewer=True)
    img_L_med = img_L_med[50:-50, :, :]
    axs[1, 0].imshow(img_L_med)
    axs[1, 0].axis('off')
    axs[1, 1].imshow(img_R_med)
    axs[1, 1].axis('off')

    # colorbar
    if colorbar == "on":
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        cb_ax = fig1.add_axes([0.45, 0.48, 0.025, 0.4])  # [left, bottom, width, height]
        cbar1 = fig1.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cb_ax, orientation='vertical')
        cbar1.ax.tick_params(labelsize=40)
        if ytick == []:
            ticks = np.linspace(vmin, vmax, 5)
            ytick = [f"{t:.2f}" for t in ticks]
        else:
            ticks = np.arange(vmin, vmax + 1)
        cbar1.set_ticks(ticks)
        cbar1.ax.set_yticklabels(ytick)
    return None

"""
Map directionality values to colors for plotting.
Non-significant measures appear in black.
"""
def measureToColor(eff, pval, vmin, vmax, measure, cmap):
    cmaplist = [cmap(i) for i in range(cmap.N)]
    values = np.linspace(vmin, vmax, 255)
    colors = []
    for i, x in zip(eff, pval):
        if measure=="directionality" and x>0.05:
            colors.append((0.2, 0.2, 0.2, 1.0)) #black if no statistical differences were found
        elif measure=="directionality" and np.isnan(x):
            colors.append((0.392, 0.392, 0.392, 0.0)) #medium gray
        else:
            j = np.argmin(np.abs(values - i))
            color0=cmaplist[j]
            colors.append(color0)
    return colors

#%%
if __name__ == "__main__":
    main()