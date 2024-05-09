README: 

CODEBOX folder:
	USE pd_paper.yml FOR RECREATING THE CONDA ENVIRONMENT FOR RUNNING THE CODE
	01_Gen_Data.sh runs code/01_Gen_Data.py; generates the data used for the simulations from the raw Gaia data; takes about 2 hours to run
	run_02_batch_pas_even.sh runs code/02_PDM_AOV_SL_multi.py; calculates the PDM, AOV, and SL statistics for evenly spaced observations with homoscedastic error under H0; takes about 7 hours to run
	run_02_batch_pas_uneven.sh runs code/02_PDM_AOV_SL_multi.py; calculates the PDM, AOV, and SL statistics for unevenly spaced observations with homoscedastic error under H0; takes about 7 hours to run
	run_02_batch_pas_hetero_trans.sh runs code/02_PDM_AOV_SL_multi.py; calculates the PDM, AOV, and SL statistics for unevenly spaced observations with heteroscedastic error (from not scaling) under H0; takes about 7 hours to run
	run_02_batch_pas_hetero_sim.sh runs code/02_PDM_AOV_SL_multi.py; calculates the PDM, AOV, and SL statistics for unevenly spaced observations with heteroscedastic error (growing over time) error under H0; takes about 7 hours to run
	run_02_batch_pas_h1.sh runs runs code/02_PDM_AOV_SL_multi.py; calculates the PDM, AOV, and SL statistics for unevenly spaced observations with homoscedastic error under H1; takes about 7 hours to run
	run_03_ls.sh runs code/03_LSPeriodogramSims.py; calculates LS power statistics for all scenarios; takes about 3 hours to run
	run_04_paslls.sh runs code/04_plot_PDM_AOV_SL_LS.py; makes the majority of the plots in the paper; takes about 15 minutes to run

RAW folder:

	Gaia documentation: https://gea.esac.esa.int/archive/

	raw/gdr3_vari_eclipsing_binary: contains the list files created by Mowlavi 2022 which fits Gaussians to the observed light curves
		Files downloaded from http://cdn.gea.esac.esa.int/Gaia/gdr3/Variability/vari_eclipsing_binary/
		To download these files using command line: 
			1) Used terminal to get a list of files: curl http://cdn.gea.esac.esa.int/Gaia/gdr3/Variability/vari_eclipsing_binary/ --output list.txt
			2) Used Excel to clean and prepend "http://cdn.gea.esac.esa.int/Gaia/gdr3/Variability/vari_eclipsing_binary/" to each name
			3) Saved new urls in list.txt
			4) Used terminal in raw/gdr3_vari_eclipsing_binary folder: wget -i list.txt
			5) Unzipped using gunzip *.csv.gz in terminal (WARNING: about 1 GB unzipped)
			6) Rezipped after running 01_Gen_Data.py to save space


	raw/gdr3_epoch_photometry: contains the Gaia DR3 photometry data; 
		this will be merged with the ones used in Mowlavi to determine the average number of observations per star

		Files downloaded from http://cdn.gea.esac.esa.int/Gaia/gdr3/Photometry/epoch_photometry/
		To download these files using command line: 
			1) Used terminal to get a list of files: curl http://cdn.gea.esac.esa.int/Gaia/gdr3/Photometry/epoch_photometry/ --output list.txt
			2) Used Excel to clean and prepend "http://cdn.gea.esac.esa.int/Gaia/gdr3/Photometry/epoch_photometry/" to each name
			3) Saved new urls in list.txt
			4) Used terminal in raw/gdr3_epoch_photometry folder: wget -i list.txt
			5) Unzipped using gunzip *.csv.gz in terminal (WARNING: about 300 GB unzipped)
			6) Rezipped after running 01_Gen_Data.py to save space


CODE folder:
	
	01_Gen_Data.py: Takes about 2 hours to run; generates the data for the simulations


	02_PDM_AOV_SL_multi.py: Calculates the PDM, AOV, and SL statistics for a single scenario considered, e.g. evenly sampled observations under H0


	03_LSPeriodogramSims.py: Calculates the LS power for all scenarios considered


	04_plot_PDM_AOV_SL_LS.py: Plots the histograms and periodograms for the AOV, PDM, SL, and LS power statistics


	05_LombScarglePeriodogramFAPGraph.py: Creates Figure 2 in the manuscript


	06_OtherFigs.py: Creates Figure 1 in the manuscript


INTER folder contains the intermediate output from all the code described above

OUTPUT folder contains all the finalized output from the code above