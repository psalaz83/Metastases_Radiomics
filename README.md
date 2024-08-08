# Metastases_Radiomics - Author: P salazar 2024. Tested with RStudio 2023.09.0 Build 463 and R version ‘4.3.0’
R-code, data and example results - metastases radiomics FPCA

1. Figure1_FPCA_Metastases_3FPC.R
Input file: metastases_io_files
Output file: Figure1_MetastasesFPCA.pdf
Compute Functional Principal Components and plots shown in Figure 1 of metastases_radiomic paper:
Predefined and data-driven CT radiomics predict recurrence-free and overall survival in patients with pulmonary metastases treated with stereotactic body radiotherapy - PLOS One to appear 2024

2. Figure2A2B_KLPlots_RFS_ssf4_Entropy_Age.R
Computation of Kaplan-Meier (KM) curves corresponding to Recurrence free survival (RFS) for the variables: ssf4_Entropy and Age, for the same paper.
Input file: Metas111_v3_all_variables.csv
Output files: KMPlot_RFS_age.pdf and KMPlot_RFS_ssf4_entropy.pdf

3. Figure2B-2C_KM_Plots_OS_LungPrimLoc_Size.R
Computation of KM curves corresponding to Overall Survival (OS) for the variables: lung primary location and Size (tumor size), for the same paper.
Input file: Metas_v4_all_variables.csv
Output files: KM_Plot_OS_lungOrigin.pdf and KM_Plot_OS_size.pdf

4. Figure3A_3B_SegmentRegression_ssf4Entropy_size_ssf0_ROI_Area_Log10.R
   Input file: Metas_v4_all_variables.csv
   Output test files: Fig3A_SegmentedRegression_ssf4entropy_vs_size.pdf Fig3B_SegmentedRegression_ssf4entropy_vs_log10ROI_Area.pdf
   Figure 3A. Find breakpoint with segmented linear regression for ssf4_entropy vs. size (cm)
   Figure 3B. Find breakpoint with segmented linear regression for ssf0_entropy vs. log10ssf0Total (that is the number of voxels in the segmented metastasis)
   Davis tests for breakpoint significance are included

5. Figure_Suppl_1_FPCA_MetasPeritumor_3FPC.R
   Input file: All_metasPeritumor_fda_out_v3.csv
   Output file: Suppl_Figure1_PerimetastasesFPCA.pdf and Suppl_Figure1B_Perimetastases_RowCurvesandMean.pdf
   Compute Functional Principal Components and plots shown in Figure 1 of peri-metastases_radiomic paper. It also visualize the row curves and computes two 
   type of mean curves: the classic Euclidean L2 cross sectional mean curve and the Wasserstein-Frechet mean curve.

6. SmoothDensityCurve_Metas_fda.R
   Input files: Alldata_111tumor_filenames_v3.csv and MET_0526_1A_output2.csv, MET_0527_1A_output2.csv, MET_0530_1A_output2.csv, MET_0541_1A_output2.csv
   Output file: MET_0526_1A_output_output_smoothC.csv, MET_0527_1A_output_output_smoothC.csv, MET_0530_1A_output_output_smoothC.csv, MET_0541_1A_output_smoothC.csv
   See: folder: ThreeSmoothCurves for input curves and output curves.
   Computes a smooth curve from a raw curve using a smoothing method adapted to density distribution (which have constraint of area sum up to one).

8. Compute_MedianFollowUpTime.R
   Input file: Metas111_v3_all_variables.csv
   Short script to compute median follow-up time in months and days taking into account the survival effect. We use the 'reverse' Kaplan-Meier method: from D.Moore's 
   book: 'Applied survival analysis using R'. Chap. 3.3.
   
