**Modeling Breeding Distribution of the Red-backed Shrike (Lanius collurio) in Iberia**
-
This repository represents the bulk of the content used to model the breeding distribution of the Red-backed Shrike in this master thesis.
The following repository contains:

  1. R script "_biomod2_github.R_"
  2. Occurrence data in "_Final qgis products_"
       .Iberian Peninsula independent census data
     
       .EBBA2, GBIF, and combined (FULL) occurrence data sets using the migration route pseudo-absence sampling design
     
       .remaining pseudo-absence sampling designs considered using the combined occurrence data set, named poly10, poly5 and squared
  3. Environmental data in folders for "current" conditions and global change future scenarios (A1B & A2 for carbon emissions; MR & CS for global circulation models)

**How to use**
-

Users are only required to download the R script "_biomod2_github.R_". The script automatically downloads the necessary files to run the analysis into a temporary folder in the user's drive.

This analysis also requires the user to download the necessary R packages from R-CRAN repository listed in the first few lines of the script. Packages between lines 7-24 are required to run the full modeling procedure, while packages between lines 25-29 are only used for more in-depth analysis and graphical outputs.

The user should also run line 52 to obtain the same results as the ones provided in the thesis.

Between lines 60-94 of the script the user should select **ONLY ONE** of the provided sampling area designs to run the analysis. A total of 4 sampling area designs were considered:
  1. "Original" area based on the migration route (lines 60-67)
  2. "Polygon 5º" which adds a 5º coordinate buffer to the known breeding range (lines 69-76)
  3. "Polygon 10º" which adds a 10º coordinate buffer to the known breeding range (lines 78-85)
  4. "Rectangular" baed on the limits of the lat/long values of recorded presences (lines 87-94)

Afterwards the code between lines 97-109 is required to save the Iberian Peninsula census data to an object.

An additional **optional** step appears between lines 111-122 where the user can remove the Iberian Peninsula from the calibration dataset. This is not necessary to run the analysis and should be skipped unless the user wants to check the results from removing the Iberian Peninsula occurrence data from calibration.

The code between lines 124-989 can be run afterwards without any further inputs from the user. This runs the whole analysis performed in this thesis.

All code after line 990 is used to provide a more in-depth analysis of the models, as well as create more detailed graphical outputs.

**Outputs**
-

The script provided creates a temporary "Results" folder in the user's drive, found in "User\AppData\Local\Temp\". To access this folder it may be required to open the "Run" ("Executar" in portuguese) program and put in "appdata" as the folder you wish to open.

The R script produces some files inside the temporary folder, including:
  1. ".txt" files containing model accuracy metrics, variable importance for fitted models and predicted changes in the distribution caused by global change
  2. ".asc" files that contain probability estimates for different models, binary presence-absence predictions and predicted range shifts in future scenarios

Additionally, graphical outputs are also produced within the R Studio software "Plots" window. Some of these are automatically produced when running the main analysis, while others require users to further run code after line 990.
