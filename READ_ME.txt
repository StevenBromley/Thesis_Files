This directory contains two sub-folders, each related to the 'Gas Jet' and 'Gold Spectroscopy' projects.
In most files, the file format is listed on the top line after a python comment character (#)
The files stored within each are described. For ANY questions about the contents or use of the listed scripts, please contact
sjbroml@g.clemson.edu.

Directory 1: 'gas_jet':
--
-folder 'high_p_designs':
	-contains solidworks parts + assembly for the components of the 'high pressure gas jet' design
-folder 'low_p_designs':
	-contains solidworks parts + assembly for the components of the 'low pressure gas jet' design
-script 'jet calcs_script' (python script):
	-estimates jet operating conditions (on-axis density, temp, etc.)
		-(optional, does not need to be run, e.g. commented out): estimates # of particles
		leaving the skimmer based on data in the file 'jet_inp.csv'
-file 'jet_inp.csv':
	-contains 3 pressures (see Ch 3 of dissertation: P_N (gas reservoir), P_cc (nozzle chamber), and P_PKR (jet dump))
--
	
Directory 2: 'gold_analysis'
--
-folder 'raw_data': all gold raw (x-ray contaminated) spectra files
-folder 'raw_data_noxrays_nosum': all gold x-ray-removed spectra files)
-files:
		-'master shot log': .csv file summarizing the CTH gold and nickel discharges
			-each numerical listing corresponds to a probe and depth with a unique shot number
			-a 'clean' frame in the x-ray-contaminated spectra is listed in parantheses (outdated inclusion but left in)
		-'regions_with_wavelngths': file describing the wavelength ranges used for this work. Format:
			- Col 0 (Placeholder); index of each line is used in 'countr' variable in multiple scripts (see below)
			- Col 1: Approximate lower wavelength of each wavelength window
			- Col 2: Approximate upper wavelength of each wavelength window
		-'shots.csv': .csv file containing all shot numbers used for the final gold analysis
			-each row corresponds to a single wavelength range. Each entry in each row is a shot number (see 'master shot log')
		-'table_auI_level_info.txt' and 'table_auII_level_info.txt':
			-contains level information for Au I and Au II w/ LaTeX formatting characters
			-Format: Col 0 (Numerical Label), Col 1 (configuration), Col 2 (J Value), Col 3 (Energy in cm^-1), Col 4 (Reference)
-scripts:
	-'steve_lib.py': master library contains many functions used in the codes below. Heavily suggested to import before running any of these scripts
	-'Calculating New Lines from Levels.py':
		-calculates all allowed E1 transitions (NOTE: Does not apply Delta L = +/-1, Delta S = 0) and
		-compares to a line list.
		-May be modified to *not* compare to a line list; may be modified to include Delta L and Delta S selection rules
	
	-'Generate Peak List from raw spectra files.py'
		-carries out peak finding on a user-specified probe depth. Must be modified to work on your system (adjust directory variables)
			-requires 'Peak Utils' codes. See Appendix D2 of dissertation.
		-generates line list. NOTE: 'Counts' output likely to be inaccurate. Use integrated intensity function / code for this purpose.
	-'Intensity Scatter Plot Code.py'
		-takes in a user-defined peak list and plots the intensity of each peak as a function of time for every probe and depth
		-Must be modified to run on your PC (change file paths & line list variables)
		-saves each plot in a user-specified directory
	-'Master Plot Script.py'
		-simple plotting script to plot the spectra in a single, user-specified wavelength range (controlled by the variable 'countr' from 0 to 16')
	-'peak_find.py'
		-sample script showing the use of the peak finding code
	-'peak_intensity.py'
		-sample script showing the use of the integrated intensity code. Constructed from intensity function within steve_lib.py
	-'plotting 187 - 800 nm.py'
		-plots one probe & depth across the entire wavelength range at a single user-specified frame.
	-'x-ray_filter.py'
		-sample script showing the use of the x-ray filter. Used to process all raw spectra files.
	-'LOPT_levels_script.py'
		-script for porting output of 'LOPT' code to a LateX-compilable text file. All LaTeX code must be user-generated, and the 
		table is imported with e.g. '\input{AuII_level_table.txt}'
		-requires a level file, e.g. 'table_auII_level_info.txt' above
	-'LOPT_lines_script.py'
		-script for porting output of 'LOPT' code to a LaTeX-compilable text file. All LaTeX code must be generated, and the table is 
		imported with e.g. '\input{AuII_lines_table.txt}'
		-requires a level file, e.g. 'table_auII_level_info.txt' above
		-requires a 'references' file (see comments within script)
		-requires a 'comments' file (see comments within script)
		-intensities are calculated with a function within steve_lib.py 
		
		