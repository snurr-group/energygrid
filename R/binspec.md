# Histogram bin specifications of how to transform an energy grid

## Planning the structure of a "binspec" class

Class: list


## Data

df: SAME FORMAT AS BEFORE


## Methods

* Convert grid(s)

## Constructors

One big constructor?
* to/from/step/width/above/below
	* bounds_from_params in regression_and_features.R
* also a default color scheme based on the other params: see color.R: color_from_binloc for a map
	* though maybe that should remain as its own function, with binspec as an argument


# Other notes
## Deprecated files/functions to delete
(COMPLETED)* copy_training_data.R: irrelevant artifact from the old NN hMOF file structure
(COMPLETED)* csd_formulas.Rds was generated in the Makefile by R/filter_ccdc_chemistry.R.  Remove the R and Makefile lines
(COMPLETED)* many_plots.R: delete it
(DELETED)* plot_diagnostics.R: plot_bin_z_vs_y has hardcoded information and no longer used
(DELETED)* scotty_figs.R
* What is the purpose of save_test_train_data.R? (It is used to save the test and train data to a csv file
Easier to reuse for other purpose, like test a new model like random forest
It is newly created)


## Files that are probably less useful
(MOVED)* hyper_tuned.R shoud be moved to an example MSDE folder, not as part of the general project.  It's highly specific to that one figure.
(MOVED to NOTEBOOK FOLDER)* load_data.R should be moved into Notebooks/setup_data.R

## Files that need to be reorganized
* get_energy_stats.R is still more convoluted than it should be. 
  (DID)Refactor by lumping tidy_energy_hists into energy_stats, testing, 
  then refactoring bin_width,min_max as a new_binspec which might just call metric_from_hists
* plot_diagnostics.R: binspec should be part of the model, not something explicitly feed with the model.
	* Honestly the model should be a full object with a predict.glmnet_mod (or predict.grid_model) function so that the nuances of pred_grid are no longer necessary
* regression_and_features.R has lots of good content
	* bin_loc_from_spec should be part of the binspec object and pre-computed
	

# Overall ideas on the types of necessary functions
* fit/predict a model from input?  But then is it the responsibility of the caller or function itself to transform grid + other data into the right format?
* The binspec of course
* Probably a griddata class type
* NOTE: is there a way to easily check input types in R, that it's the correct grid/binspec/other?




