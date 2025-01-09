**Modeling Breeding Distribution of the Red-backed Shrike (Lanius collurio) in Iberia**
-

The following repository contains:

  1. R script "_biomod2_github.R_"
  2. Occurrence data in "_Final qgis products_"
  3. Environmental data in folders for "current" conditions and global change future scenarios (A1B & A2 for carbon emissions; MR & CS for global circulation models)

**How to use**
-

Only the R script "_biomod2_github.R_" is required to download. The script automatically downloads the necessary files to run the analysis into a temporary folder in the user's drive.

This analysis also requires the user to download the necessary R packages from R-CRAN repository listed in the first few lines of the script. Packages between lines 7-24 are required to run the full modeling procedure, while packages between lines 25-29 are only used for more in-depth analysis and graphical outputs.

Between lines 60-94 of the script the user should select **only one** of the provided sampling area designs to run the analysis. A total of 4 sampling area designs were considered:
  1. "Original" area based on the migration route (lines 60-67)
  2. "Polygon 5º" which adds a 5º coordinate buffer to the known breeding range (lines 69-76)
  3. "Polygon 10º" which adds a 10º coordinate buffer to the known breeding range (lines 78-85)
  4. "Rectangular" baed on the limits of the lat/long values of recorded presences (lines 87-94)

