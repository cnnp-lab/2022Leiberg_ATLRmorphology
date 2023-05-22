# 2023Leiberg_ATLRmorphology

Analysis code for the paper "Effects of anterior temporal lobe resection on cortical morphology".

## Dependencies:

The script voxel_featuers_excludeTL.m requires the MATLAB libraries [CorticalFoldingAnalysisTools](https://github.com/cnnp-lab/CorticalFoldingAnalysisTools) and [iso2mesh](http://iso2mesh.sourceforge.net).
The script run_analysis.m requires the MATLAB library [SurfStat](https://www.math.mcgill.ca/keith/surfstat/).
The scripts run on MATLAB versions starting with R2016b.

## Details on how to run:

The script voxel_features_excludeTL.m is an adaptation of the [script](https://github.com/KarolineLeiberg/folding_pointwise) published along with the paper [Local Morphological Measures Confirm that Folding within Small Partitions of the Human Cortex Follows Universal Scaling Law](https://doi.org/10.1007/978-3-030-87234-2_65) to allow the exclusion of the temporal lobe region from the morphological analysis. We ran the script on individuals with temporal lobe epilepsy excluding the ipsilateral temporal lobe, and we used it to process all controls twice, with and without temporal lobes. The script requires the FreeSurfer annotation files ```?.aparc.annot```, in addition to the files used in the previous version ```?h.pial```, ```?h.thickness```, ```?h.pial-outer-smoothed``` (see the [FreeSurfer lGI website](https://surfer.nmr.mgh.harvard.edu/fswiki/LGI)), and an expansion of the ```?h.white``` surface by 0.5 saved as ```?h.mid``` (using FreeSurfer ```mris_expand -thickness ?h.white 0.5 ?h.mid```).

The script pointwiseGAM.R requires an input of all subjects' voxel-wise data in a table, and a meta data file containing age, sex, site/scanning protocol, and whether the subject is a control, in the same order as the data.

After running the GAM correction, we used the freesurfer function [mris_apply_reg](http://ftp.nmr.mgh.harvard.edu/pub/dist/freesurfer/tutorial_packages/OSX/freesurfer/docs/html/mris_apply_reg.help.xml.html) with the [Xhemi](https://surfer.nmr.mgh.harvard.edu/fswiki/Xhemi) interhemispheric registration to convert left hemisphere data to right hemisphere and vice versa for individuals with right onset epilepsy. We collated all data of subjects with right and left onset into a single table for each measure to use in the run_analysis.m script.
