READ ME

Enclosed in this project are the MATLAB based scripts and data for the Boylan et al. study titled, “Patterns of neural sharpening and generalization emerge in the infant brain following novel object associative learning” 

Scripts/ Code
	bootstrappedBayesianTopographicModels.m
		- This script contains the analysis code for the bootstrapped Bayesian topographic models
		- Input variables: allData, ages_alphabetical
		- This script also includes code for:
			- 2D topographies with multi-colored electrode markers denoting various levels of BF > 3 
			- Regionalized topography figure
			- All data tables 
			- Summary model comparison by region figure
	This script requires 2 other functions:
		- bootstrap2BF_z.m
			- This function computes BFs following the bootstrapped distribution creation
		- topoplotMB.m
			- This function is a slightly modified topoplot from EEGLab. I edited it to allow for many colors to be used as electrode markers


Data Variables in Boylan_etal_DCN_bayesianModels.mat
	- ages_alphabetical: vector with each subj age group, in the order of the allData mat
	- allData: contains topographic SNR data for all subjs
		- 4 dimensions: 129 (electrodes) x 51 (infants) x 4 (TS,GS1-3) x 3 (label)
	- BF10_models_avgLab: BF values across the topography, 129 (elec) x 2 (model) x 3 (age group)
	- BF10_models_byLabel: BF values across the topography, 129 (elec) x3 (label) x 2 (model) x 3 (age group)
	- locsEEGLAB129HCL: EEGLAB structure of the 129 montage
	- regions3, regions5: cell array with the elecs included in the regions across the topography
		- regions5: 1 = frontal, 2 = central, 3 = temporal, 4 = pariental, 5 = occipital
		- regions3: 1 = frontal/central, 2 = temporal/parietal, 3 = occipital

NOTE: specific distributions used in this project are too large to be shared, but can be re-created. 
The shared “BF…” variables are the ones used in the results of the project 

