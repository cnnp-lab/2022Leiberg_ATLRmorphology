# 2023Leiberg_ATLRmorphology

Analysis code for the paper "Effects of anterior temporal lobe resection on cortical morphology".

## Dependencies:

The script voxel_featuers_excludeTL.m requires the MATLAB libraries [CorticalFoldingAnalysisTools](https://github.com/cnnp-lab/CorticalFoldingAnalysisTools) and [iso2mesh](http://iso2mesh.sourceforge.net).
The script run_analysis.m requires the MATLAB library [SurfStat](https://www.math.mcgill.ca/keith/surfstat/).
The scripts run on MATLAB versions starting with R2016b.

## Details on how to run:

Data used for the analysis can be found on zenodo, [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7957831.svg)](https://doi.org/10.5281/zenodo.7957831).

The script **fs_voxel_features_excludeTL.m** is a version of the [script](https://github.com/KarolineLeiberg/folding_pointwise) published along with the paper [Local Morphological Measures Confirm that Folding within Small Partitions of the Human Cortex Follows Universal Scaling Law](https://doi.org/10.1007/978-3-030-87234-2_65), adapted to allow the exclusion of the temporal lobe region from the morphological analysis. We ran the script on individuals with temporal lobe epilepsy pre- and postoperatively excluding the ipsilateral temporal lobe, and we used it to process all controls twice, with and without temporal lobes. The script requires the FreeSurfer annotation files ```?.aparc.annot```, in addition to the files used in the previous version ```?h.pial```, ```?h.thickness```, ```?h.pial-outer-smoothed``` (see the [FreeSurfer lGI website](https://surfer.nmr.mgh.harvard.edu/fswiki/LGI)), and an expansion of the ```?h.white``` surface by 0.5 saved as ```?h.mid``` (using FreeSurfer ```mris_expand -thickness ?h.white 0.5 ?h.mid```).

After running **fs_voxel_features_excludeTL.m**, we registered each subject's data to the fsaverage using the FreeSurfer function ```mri_surf2surf```. We then collated the subjects' data into one file for each measure, hemisphere, and onset side (to use healthy controls with their left temporal lobe excluded when correcting left onset patients etc). We also produced corresponding meta data files containing age, sex, site/scanning protocol, and whether the subject is a control, in the same order as the data.

These data and meta data files can be found in the folder *not_corrected* on [zenodo](https://doi.org/10.5281/zenodo.7957831). They are required to run the script **pointwiseGAM.R**. Line 7 of the script can be edited to run for different measures, onset side, and hemisphere. (Paths on line 9,10,38 can be edited for different input/ output folders.)

After running the GAM correction, we used the freesurfer function [mris_apply_reg](http://ftp.nmr.mgh.harvard.edu/pub/dist/freesurfer/tutorial_packages/OSX/freesurfer/docs/html/mris_apply_reg.help.xml.html) with the [Xhemi](https://surfer.nmr.mgh.harvard.edu/fswiki/Xhemi) interhemispheric registration to convert left hemisphere data to right hemisphere and vice versa for individuals with right onset epilepsy. We then collated (using [SurfStat](https://www.math.mcgill.ca/keith/surfstat/) functions ```SurfStatListDir``` & ```SurfStatReadData```) data of all subjects with TLE into pre- and postoperative tables for each measure to use in the **run_analysis.m** script. The data can be found in the folder *age_sex_corrected* on [zenodo](https://doi.org/10.5281/zenodo.7957831). The script requires the [surfstat](https://www.math.mcgill.ca/keith/surfstat/) folder in its path, as well as the fsaverage FreeSurfer folder. (Paths on line 11 and 14 can be edited for different input/ output folders.)
