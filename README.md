AG3C
====


	Included files and description


	run_analysis.py
		python script that executes the analysis pipeline.
	
	GO_analysis/ChartReport.py
		Runs a chart report request on the David web service
	GO_analysis/DAVIDWebService_Client.py
		Client script for David web service query
	GO_analysis/DAVIDWEBREADME.txt
		Readme file for david web service scripts
	GO_analysis/FuncClust.py
		Runs a functional cluster request on the David web service
	GO_analysis/go_rise.py
		Calculates the average time it takes for all proteins in a GO term
		to begin changing (t1 in main txt).
	GO_analysis/id_convert.p
		Converts Accession to Entrez ID
	GO_analysis/mean_of_genes_in_go_term_mRNA.py
		Averages the gene expression profiles that are contained in a given GO term. 
	GO_analysis/readcsv.py
		C simple script to read in txt delimeted files



	data/Oi_name_map.p
		Pickle file, gives the apex OI values indexed by their common names.
	data/Test_data.csv 
		Gives some test RNA raw count data. Unit test for RNA normalization and data formatting
	data/Test_data_answers.xlsx
		Gives the test RNA raw count data answers to normalization and data formatting
	data/Test_name_map.p
		Maps the accession to entrez IDs for querying DAVID webservice
	data/emili_orphan_predictions_ecoli.csv
		Emili lab computational predictions of gene functions downloaded from publication <cite>.
	data/prot_data.csv
		Raw protein counts
	data/rna_data.csv
		Raw rna counts
	data/rna_name_map.csv
		Maps the RNA (ECB number) to Entrez accession.


	data_processing/cluster_kmeans.py
		Preform k-means clustering and plot the results
	data_processing/format_data.py
		Processes count data as described in main text. normalizes counts using read depth.
	data_processing/format_data_deseq.py
		Process count data. normalize counts using DEseq results
	data_processing/format_data_norm_to_length.py
		Process count data.normalize counts by transcript length (for absolute comparison)
	data_processing/format_data_wapex.py
		Process count data. Normalize counts using APEX OI values (for absolute comparison)
	run_analysis.py
	
	
	Figures/
		Contains figures of publication

	fit_time_course/time_course_fit.py
		Script that fit's piecewise continuous curve to expression profiles. Uses custom implimination of Differential Evolution algorthim.
	fit_time_course/classify_fits.py
		Sort fits based upon parameters estimated from fitting the piecewise continous curve.


	operon_analysis/operon_analysis_and_cluster_prot.ipynb
		ipython notebook script for creating correlation of expression profiles within operons

	operon_analysis/OpersonSet2.txt
		List of operons downloaded from RegulonDB.



