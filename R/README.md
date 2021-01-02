THis folder contains all files needed to run the code.
_____
**Initial Setup**:<br/>
One needs **R code** and **RStudio** as requirements. <br/> Zhao uses Windows 10, so please forgive me for lack of knowledge on other platforms. <br/>
After installation, one can run **package_verification.R** to check for needed R packages. If a package is not found, it will try to install it. <br/>
**Energy grid data** for structures that will be used for machine learning. 
_____
**Before Machine Learning**<br/>
Energy Grid machine learning uses tabulated energy points. In its first step, it converts energy values stored in energy grid files for each structure into a fine histogram with bin size equals to 0.1 kJ/mol. <br/>
This process is enabled by **save_h2_hists.R**<br/>
Using Vanilla (base) R, in command prompt (Windows), type the following command:<br/>
```
Rscript --vanilla R\save_h2_hists.R whatever.rds Energies\ <keyword>
```
**'Energies\'** is the directory where the energy grid files are stored. <br/>
Keywords specifies the maximum and minimum bounds for the histogram, using **'autotune'** will automatically detect the minimum energy from the energy grid files. <br/>
The fine histograms are stored in **whatever.rds**, ready for machine learning. <br/>
_____
The wrapper for single component prediction is **Propane_ML_Histogram.R**, a similar wrapper for mixture fitting is **Read_XeKr.R**<br/>
It reads the data from "**SI.xlsx**" from "**All_data**" folder. <br/>
The inputs one needs to define in this wrapper is: <br/>
1. **'molecule_name'**: "Ethane", "Propane", "Xe", "Kr" <br/>
2. **'Temperature'**: "273K", "298K" <br/>
3. **'Pressure'**: "1Bar", "5Bar", "10Bar" (For more, see different sheets in SI.xlsx) <br/>
4. **.rds** file that stores fine histograms from energy grids. (Look for the previous section for details.)<br/>
_____
**How Energy Histograms Are Used as Features**<br/>
The fine histograms are aggregated in when performing machine learning. <br/>
There are mainly two ways of aggregation:<br/>
1. Equally sized bins, for example, bin size equals 1kJ/mol. <br/>
2. Automatic aggregation. The code decides the bin bounds depending on the weight of energy regions: if a region has more counts, then the code use more bins to describe it. <br/>
_____
**Machine Learning**<br/>
The code first uses LASSO and then RandomForest for the fitting. <br/>
Additional features are also tested, for example, R score between energy histograms, textural properties, and energy statistics from energy grid files. <br/>
_____
**Where does output go???**<br/>
During the machine learning run, the outputs (training/testing data, parity plots, histograms) are stored in **Results/** folder. <br/>
Subfolders are named and created according to the molecule, condition, and grid probe. <br/>
_____
