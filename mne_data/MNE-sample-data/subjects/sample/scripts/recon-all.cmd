
 mri_convert /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/T1_002.mgz /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/sample/mri/orig/001.mgz 

#--------------------------------------------
#@# MotionCor Thu Sep 15 14:10:39 EDT 2011

 cp /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/sample/mri/orig/001.mgz /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/sample/mri/rawavg.mgz 


 mri_convert /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/sample/mri/rawavg.mgz /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/sample/mri/orig.mgz --conform 


 mri_add_xform_to_header -c /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/sample/mri/transforms/talairach.xfm /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/sample/mri/orig.mgz /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/sample/mri/orig.mgz 

#--------------------------------------------
#@# Talairach Thu Sep 15 14:10:56 EDT 2011

 talairach_avi --i orig.mgz --xfm transforms/talairach.auto.xfm 


 cp transforms/talairach.auto.xfm transforms/talairach.xfm 

#--------------------------------------------
#@# Talairach Failure Detection Thu Sep 15 14:11:47 EDT 2011

 talairach_afd -T 0.005 -xfm transforms/talairach.xfm 


 awk -f /usr/local/freesurfer/stable5_1_0/bin/extract_talairach_avi_QA.awk /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/sample/mri/transforms/talairach_avi.log 

#--------------------------------------------
#@# Nu Intensity Correction Thu Sep 15 14:11:47 EDT 2011

 mri_nu_correct.mni --i orig.mgz --o nu.mgz --uchar transforms/talairach.xfm --n 2 

#--------------------------------------------
#@# Intensity Normalization Thu Sep 15 14:16:26 EDT 2011

 mri_normalize -g 1 nu.mgz T1.mgz 

#--------------------------------------------
#@# Skull Stripping Thu Sep 15 14:18:36 EDT 2011

 mri_em_register -skull nu.mgz /usr/local/freesurfer/stable5_1_0/average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta 


 mri_watershed -T1 -brain_atlas /usr/local/freesurfer/stable5_1_0/average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta T1.mgz brainmask.auto.mgz 


 cp brainmask.auto.mgz brainmask.mgz 

#-------------------------------------
#@# EM Registration Thu Sep 15 14:32:05 EDT 2011

 mri_em_register -uns 3 -mask brainmask.mgz nu.mgz /usr/local/freesurfer/stable5_1_0/average/RB_all_2008-03-26.gca transforms/talairach.lta 

#--------------------------------------
#@# CA Normalize Thu Sep 15 14:54:30 EDT 2011

 mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz /usr/local/freesurfer/stable5_1_0/average/RB_all_2008-03-26.gca transforms/talairach.lta norm.mgz 

#--------------------------------------
#@# CA Reg Thu Sep 15 14:56:19 EDT 2011

 mri_ca_register -nobigventricles -T transforms/talairach.lta -align-after -mask brainmask.mgz norm.mgz /usr/local/freesurfer/stable5_1_0/average/RB_all_2008-03-26.gca transforms/talairach.m3z 

#--------------------------------------
#@# CA Reg Inv Thu Sep 15 18:17:28 EDT 2011

 mri_ca_register -invert-and-save transforms/talairach.m3z 

#--------------------------------------
#@# Remove Neck Thu Sep 15 18:18:49 EDT 2011

 mri_remove_neck -radius 25 nu.mgz transforms/talairach.m3z /usr/local/freesurfer/stable5_1_0/average/RB_all_2008-03-26.gca nu_noneck.mgz 

#--------------------------------------
#@# SkullLTA Thu Sep 15 18:20:10 EDT 2011

 mri_em_register -skull -t transforms/talairach.lta nu_noneck.mgz /usr/local/freesurfer/stable5_1_0/average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta 

#--------------------------------------
#@# SubCort Seg Thu Sep 15 18:41:01 EDT 2011

 mri_ca_label -align -nobigventricles norm.mgz transforms/talairach.m3z /usr/local/freesurfer/stable5_1_0/average/RB_all_2008-03-26.gca aseg.auto_noCCseg.mgz 


 mri_cc -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz -lta /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/sample/mri/transforms/cc_up.lta sample 

#--------------------------------------
#@# Merge ASeg Thu Sep 15 18:57:47 EDT 2011

 cp aseg.auto.mgz aseg.mgz 

#--------------------------------------------
#@# Intensity Normalization2 Thu Sep 15 18:57:47 EDT 2011

 mri_normalize -aseg aseg.mgz -mask brainmask.mgz norm.mgz brain.mgz 

#--------------------------------------------
#@# Mask BFS Thu Sep 15 19:01:22 EDT 2011

 mri_mask -T 5 brain.mgz brainmask.mgz brain.finalsurfs.mgz 

#--------------------------------------------
#@# WM Segmentation Thu Sep 15 19:01:25 EDT 2011

 mri_segment brain.mgz wm.seg.mgz 


 mri_edit_wm_with_aseg -keep-in wm.seg.mgz brain.mgz aseg.mgz wm.asegedit.mgz 


 mri_pretess wm.asegedit.mgz wm norm.mgz wm.mgz 

#--------------------------------------------
#@# Fill Thu Sep 15 19:03:56 EDT 2011

 mri_fill -a ../scripts/ponscc.cut.log -xform transforms/talairach.lta -segmentation aseg.auto_noCCseg.mgz wm.mgz filled.mgz 

#--------------------------------------------
#@# Tessellate lh Thu Sep 15 19:04:45 EDT 2011

 mri_pretess ../mri/filled.mgz 255 ../mri/norm.mgz ../mri/filled-pretess255.mgz 


 mri_tessellate ../mri/filled-pretess255.mgz 255 ../surf/lh.orig.nofix 


 rm -f ../mri/filled-pretess255.mgz 


 mris_extract_main_component ../surf/lh.orig.nofix ../surf/lh.orig.nofix 

#--------------------------------------------
#@# Smooth1 lh Thu Sep 15 19:04:54 EDT 2011

 mris_smooth -nw -seed 1234 ../surf/lh.orig.nofix ../surf/lh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 lh Thu Sep 15 19:04:59 EDT 2011

 mris_inflate -no-save-sulc ../surf/lh.smoothwm.nofix ../surf/lh.inflated.nofix 

#--------------------------------------------
#@# QSphere lh Thu Sep 15 19:05:38 EDT 2011

 mris_sphere -q -seed 1234 ../surf/lh.inflated.nofix ../surf/lh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology lh Thu Sep 15 19:10:43 EDT 2011

 cp ../surf/lh.orig.nofix ../surf/lh.orig 


 cp ../surf/lh.inflated.nofix ../surf/lh.inflated 


 mris_fix_topology -mgz -sphere qsphere.nofix -ga -seed 1234 sample lh 


 mris_euler_number ../surf/lh.orig 


 mris_remove_intersection ../surf/lh.orig ../surf/lh.orig 


 rm ../surf/lh.inflated 

#--------------------------------------------
#@# Make White Surf lh Thu Sep 15 19:39:07 EDT 2011

 mris_make_surfaces -noaparc -whiteonly -mgz -T1 brain.finalsurfs sample lh 

#--------------------------------------------
#@# Smooth2 lh Thu Sep 15 19:45:07 EDT 2011

 mris_smooth -n 3 -nw -seed 1234 ../surf/lh.white ../surf/lh.smoothwm 

#--------------------------------------------
#@# Inflation2 lh Thu Sep 15 19:45:13 EDT 2011

 mris_inflate ../surf/lh.smoothwm ../surf/lh.inflated 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 ../surf/lh.inflated 


#-----------------------------------------
#@# Curvature Stats lh Thu Sep 15 19:47:25 EDT 2011

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/lh.curv.stats -F smoothwm sample lh curv sulc 

#--------------------------------------------
#@# Sphere lh Thu Sep 15 19:47:33 EDT 2011

 mris_sphere -seed 1234 ../surf/lh.inflated ../surf/lh.sphere 

#--------------------------------------------
#@# Surf Reg lh Thu Sep 15 20:42:23 EDT 2011

 mris_register -curv ../surf/lh.sphere /usr/local/freesurfer/stable5_1_0/average/lh.average.curvature.filled.buckner40.tif ../surf/lh.sphere.reg 

#--------------------------------------------
#@# Jacobian white lh Thu Sep 15 21:13:47 EDT 2011

 mris_jacobian ../surf/lh.white ../surf/lh.sphere.reg ../surf/lh.jacobian_white 

#--------------------------------------------
#@# AvgCurv lh Thu Sep 15 21:13:50 EDT 2011

 mrisp_paint -a 5 /usr/local/freesurfer/stable5_1_0/average/lh.average.curvature.filled.buckner40.tif#6 ../surf/lh.sphere.reg ../surf/lh.avg_curv 

#-----------------------------------------
#@# Cortical Parc lh Thu Sep 15 21:13:53 EDT 2011

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 sample lh ../surf/lh.sphere.reg /usr/local/freesurfer/stable5_1_0/average/lh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs ../label/lh.aparc.annot 

#--------------------------------------------
#@# Make Pial Surf lh Thu Sep 15 21:14:40 EDT 2011

 mris_make_surfaces -white NOWRITE -mgz -T1 brain.finalsurfs sample lh 

#--------------------------------------------
#@# Surf Volume lh Thu Sep 15 21:26:18 EDT 2011

 mris_calc -o lh.area.mid lh.area add lh.area.pial 


 mris_calc -o lh.area.mid lh.area.mid div 2 


 mris_calc -o lh.volume lh.area.mid mul lh.thickness 

#-----------------------------------------
#@# Parcellation Stats lh Thu Sep 15 21:26:19 EDT 2011

 mris_anatomical_stats -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab sample lh white 

#-----------------------------------------
#@# Cortical Parc 2 lh Thu Sep 15 21:26:37 EDT 2011

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 sample lh ../surf/lh.sphere.reg /usr/local/freesurfer/stable5_1_0/average/lh.destrieux.simple.2009-07-29.gcs ../label/lh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 lh Thu Sep 15 21:27:31 EDT 2011

 mris_anatomical_stats -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.a2009s.stats -b -a ../label/lh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab sample lh white 

#--------------------------------------------
#@# Tessellate rh Thu Sep 15 21:27:50 EDT 2011

 mri_pretess ../mri/filled.mgz 127 ../mri/norm.mgz ../mri/filled-pretess127.mgz 


 mri_tessellate ../mri/filled-pretess127.mgz 127 ../surf/rh.orig.nofix 


 rm -f ../mri/filled-pretess127.mgz 


 mris_extract_main_component ../surf/rh.orig.nofix ../surf/rh.orig.nofix 

#--------------------------------------------
#@# Smooth1 rh Thu Sep 15 21:28:00 EDT 2011

 mris_smooth -nw -seed 1234 ../surf/rh.orig.nofix ../surf/rh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 rh Thu Sep 15 21:28:05 EDT 2011

 mris_inflate -no-save-sulc ../surf/rh.smoothwm.nofix ../surf/rh.inflated.nofix 

#--------------------------------------------
#@# QSphere rh Thu Sep 15 21:28:44 EDT 2011

 mris_sphere -q -seed 1234 ../surf/rh.inflated.nofix ../surf/rh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology rh Thu Sep 15 21:33:36 EDT 2011

 cp ../surf/rh.orig.nofix ../surf/rh.orig 


 cp ../surf/rh.inflated.nofix ../surf/rh.inflated 


 mris_fix_topology -mgz -sphere qsphere.nofix -ga -seed 1234 sample rh 


 mris_euler_number ../surf/rh.orig 


 mris_remove_intersection ../surf/rh.orig ../surf/rh.orig 


 rm ../surf/rh.inflated 

#--------------------------------------------
#@# Make White Surf rh Thu Sep 15 21:58:48 EDT 2011

 mris_make_surfaces -noaparc -whiteonly -mgz -T1 brain.finalsurfs sample rh 

#--------------------------------------------
#@# Smooth2 rh Thu Sep 15 22:04:50 EDT 2011

 mris_smooth -n 3 -nw -seed 1234 ../surf/rh.white ../surf/rh.smoothwm 

#--------------------------------------------
#@# Inflation2 rh Thu Sep 15 22:04:56 EDT 2011

 mris_inflate ../surf/rh.smoothwm ../surf/rh.inflated 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 ../surf/rh.inflated 


#-----------------------------------------
#@# Curvature Stats rh Thu Sep 15 22:07:10 EDT 2011

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/rh.curv.stats -F smoothwm sample rh curv sulc 

#--------------------------------------------
#@# Sphere rh Thu Sep 15 22:07:18 EDT 2011

 mris_sphere -seed 1234 ../surf/rh.inflated ../surf/rh.sphere 

#--------------------------------------------
#@# Surf Reg rh Thu Sep 15 23:55:45 EDT 2011

 mris_register -curv ../surf/rh.sphere /usr/local/freesurfer/stable5_1_0/average/rh.average.curvature.filled.buckner40.tif ../surf/rh.sphere.reg 

#--------------------------------------------
#@# Jacobian white rh Fri Sep 16 00:40:13 EDT 2011

 mris_jacobian ../surf/rh.white ../surf/rh.sphere.reg ../surf/rh.jacobian_white 

#--------------------------------------------
#@# AvgCurv rh Fri Sep 16 00:40:16 EDT 2011

 mrisp_paint -a 5 /usr/local/freesurfer/stable5_1_0/average/rh.average.curvature.filled.buckner40.tif#6 ../surf/rh.sphere.reg ../surf/rh.avg_curv 

#-----------------------------------------
#@# Cortical Parc rh Fri Sep 16 00:40:19 EDT 2011

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 sample rh ../surf/rh.sphere.reg /usr/local/freesurfer/stable5_1_0/average/rh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs ../label/rh.aparc.annot 

#--------------------------------------------
#@# Make Pial Surf rh Fri Sep 16 00:41:07 EDT 2011

 mris_make_surfaces -white NOWRITE -mgz -T1 brain.finalsurfs sample rh 

#--------------------------------------------
#@# Surf Volume rh Fri Sep 16 00:53:07 EDT 2011

 mris_calc -o rh.area.mid rh.area add rh.area.pial 


 mris_calc -o rh.area.mid rh.area.mid div 2 


 mris_calc -o rh.volume rh.area.mid mul rh.thickness 

#-----------------------------------------
#@# Parcellation Stats rh Fri Sep 16 00:53:07 EDT 2011

 mris_anatomical_stats -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab sample rh white 

#-----------------------------------------
#@# Cortical Parc 2 rh Fri Sep 16 00:53:26 EDT 2011

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 sample rh ../surf/rh.sphere.reg /usr/local/freesurfer/stable5_1_0/average/rh.destrieux.simple.2009-07-29.gcs ../label/rh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 rh Fri Sep 16 00:54:21 EDT 2011

 mris_anatomical_stats -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.a2009s.stats -b -a ../label/rh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab sample rh white 

#--------------------------------------------
#@# Cortical ribbon mask Fri Sep 16 00:54:41 EDT 2011

 mris_volmask --label_left_white 2 --label_left_ribbon 3 --label_right_white 41 --label_right_ribbon 42 --save_ribbon sample 

#--------------------------------------------
#@# ASeg Stats Fri Sep 16 01:08:16 EDT 2011

 mri_segstats --seg mri/aseg.mgz --sum stats/aseg.stats --pv mri/norm.mgz --empty --excludeid 0 --excl-ctxgmwm --supratent --subcortgray --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --etiv --surf-wm-vol --surf-ctx-vol --totalgray --ctab /usr/local/freesurfer/stable5_1_0/ASegStatsLUT.txt --subject sample 

#-----------------------------------------
#@# AParc-to-ASeg Fri Sep 16 01:15:21 EDT 2011

 mri_aparc2aseg --s sample --volmask 


 mri_aparc2aseg --s sample --volmask --a2009s 

#-----------------------------------------
#@# WMParc Fri Sep 16 01:18:44 EDT 2011

 mri_aparc2aseg --s sample --labelwm --hypo-as-wm --rip-unknown --volmask --o mri/wmparc.mgz --ctxseg aparc+aseg.mgz 


 mri_segstats --seg mri/wmparc.mgz --sum stats/wmparc.stats --pv mri/norm.mgz --excludeid 0 --brain-vol-from-seg --brainmask mri/brainmask.mgz --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --subject sample --surf-wm-vol --ctab /usr/local/freesurfer/stable5_1_0/WMParcStatsLUT.txt --etiv 

#--------------------------------------------
#@# BA Labels lh Fri Sep 16 01:34:07 EDT 2011
INFO: fsaverage subject does not exist in SUBJECTS_DIR
INFO: Creating symlink to fsaverage subject...

 cd /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects; ln -s /usr/local/freesurfer/stable5_1_0/subjects/fsaverage; cd - 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/lh.BA1.label --trgsubject sample --trglabel ./lh.BA1.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/lh.BA2.label --trgsubject sample --trglabel ./lh.BA2.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/lh.BA3a.label --trgsubject sample --trglabel ./lh.BA3a.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/lh.BA3b.label --trgsubject sample --trglabel ./lh.BA3b.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/lh.BA4a.label --trgsubject sample --trglabel ./lh.BA4a.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/lh.BA4p.label --trgsubject sample --trglabel ./lh.BA4p.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/lh.BA6.label --trgsubject sample --trglabel ./lh.BA6.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/lh.BA44.label --trgsubject sample --trglabel ./lh.BA44.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/lh.BA45.label --trgsubject sample --trglabel ./lh.BA45.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/lh.V1.label --trgsubject sample --trglabel ./lh.V1.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/lh.V2.label --trgsubject sample --trglabel ./lh.V2.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/lh.MT.label --trgsubject sample --trglabel ./lh.MT.label --hemi lh --regmethod surface 


 mris_label2annot --s sample --hemi lh --ctab /usr/local/freesurfer/stable5_1_0/average/colortable_BA.txt --l lh.BA1.label --l lh.BA2.label --l lh.BA3a.label --l lh.BA3b.label --l lh.BA4a.label --l lh.BA4p.label --l lh.BA6.label --l lh.BA44.label --l lh.BA45.label --l lh.V1.label --l lh.V2.label --l lh.MT.label --a BA --maxstatwinner --noverbose 


 mris_anatomical_stats -mgz -f ../stats/lh.BA.stats -b -a ./lh.BA.annot -c ./BA.ctab sample lh white 

#--------------------------------------------
#@# BA Labels rh Fri Sep 16 01:36:39 EDT 2011

 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/rh.BA1.label --trgsubject sample --trglabel ./rh.BA1.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/rh.BA2.label --trgsubject sample --trglabel ./rh.BA2.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/rh.BA3a.label --trgsubject sample --trglabel ./rh.BA3a.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/rh.BA3b.label --trgsubject sample --trglabel ./rh.BA3b.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/rh.BA4a.label --trgsubject sample --trglabel ./rh.BA4a.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/rh.BA4p.label --trgsubject sample --trglabel ./rh.BA4p.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/rh.BA6.label --trgsubject sample --trglabel ./rh.BA6.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/rh.BA44.label --trgsubject sample --trglabel ./rh.BA44.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/rh.BA45.label --trgsubject sample --trglabel ./rh.BA45.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/rh.V1.label --trgsubject sample --trglabel ./rh.V1.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/rh.V2.label --trgsubject sample --trglabel ./rh.V2.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects/fsaverage/label/rh.MT.label --trgsubject sample --trglabel ./rh.MT.label --hemi rh --regmethod surface 


 mris_label2annot --s sample --hemi rh --ctab /usr/local/freesurfer/stable5_1_0/average/colortable_BA.txt --l rh.BA1.label --l rh.BA2.label --l rh.BA3a.label --l rh.BA3b.label --l rh.BA4a.label --l rh.BA4p.label --l rh.BA6.label --l rh.BA44.label --l rh.BA45.label --l rh.V1.label --l rh.V2.label --l rh.MT.label --a BA --maxstatwinner --noverbose 


 mris_anatomical_stats -mgz -f ../stats/rh.BA.stats -b -a ./rh.BA.annot -c ./BA.ctab sample rh white 

#--------------------------------------------
#@# Ex-vivo Entorhinal Cortex Label lh Fri Sep 16 01:39:13 EDT 2011
INFO: lh.EC_average subject does not exist in SUBJECTS_DIR
INFO: Creating symlink to lh.EC_average subject...

 cd /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects; ln -s /usr/local/freesurfer/stable5_1_0/subjects/lh.EC_average; cd - 


 mris_spherical_average -erode 1 -orig white -t 0.4 -o sample label lh.entorhinal lh sphere.reg lh.EC_average lh.entorhinal_exvivo.label 


 mris_anatomical_stats -mgz -f ../stats/lh.entorhinal_exvivo.stats -b -l ./lh.entorhinal_exvivo.label sample lh white 

#--------------------------------------------
#@# Ex-vivo Entorhinal Cortex Label rh Fri Sep 16 01:39:30 EDT 2011
INFO: rh.EC_average subject does not exist in SUBJECTS_DIR
INFO: Creating symlink to rh.EC_average subject...

 cd /autofs/cluster/fusion/gramfort/work/data/bem_giorgio/subjects; ln -s /usr/local/freesurfer/stable5_1_0/subjects/rh.EC_average; cd - 


 mris_spherical_average -erode 1 -orig white -t 0.4 -o sample label rh.entorhinal rh sphere.reg rh.EC_average rh.entorhinal_exvivo.label 


 mris_anatomical_stats -mgz -f ../stats/rh.entorhinal_exvivo.stats -b -l ./rh.entorhinal_exvivo.label sample rh white 

