THis folder contains all files needed to run the code.
_____
**Initial Setup**:<br/>
One needs **R code** and **RStudio** as requirements. Zhao uses Windows 10, so please forgive me for lack of knowledge on other platforms. 
After installation, one can run **package_verification.R** to check for needed R packages. If a package is not found, it will try to install it. 
_____
The wrapper is **Propane_ML_Histogram.R**
It reads the data from "**SI.xlsx**" from "**All_data**" folder. <br/>
The inputs one needs to define in this file is: <br/>
1. **'molecule_name'**: "Ethane", "Propane", "Xe", "Kr" <br/>
2. **'Temperature'**: "273K", "298K" <br/>
3. **'Pressure'**: "1Bar", "5Bar", "10Bar" (For more, see different sheets in SI.xlsx) <br/>
_____

