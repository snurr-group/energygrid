.PHONY: all backup test

RUN_R := export R_LIBS_USER="C:/Users/Benjamin/Documents/R/win-library/3.4" && Rscript --vanilla

# Can also consider restructuring the repo, e.g. http://drivendata.github.io/cookiecutter-data-science/
# Variable tutorial: https://ftp.gnu.org/old-gnu/Manuals/make-3.79.1/html_chapter/make_6.html

all:
	@echo "Sample make file for experimentation.  Still needs work.  Only backup implemented"

backup:
	rsync -av --exclude=".*" --exclude="BigData/" --delete . ../../Box\ Sync/Projects/GitBackups/EnergyGrid

BigData/init_h2_hist_vals.Rds:
	$(RUN_R) R/save_h2_hists.R $@ BigData/10k-hMOFs/part1/CIF_FILES BigData/10k-hMOFs/part2/CIF_FILES

BigData/Robj/ccdc_hist_vals.Rds:
	$(RUN_R) R/save_h2_hists.R $@ BigData/CSD/csddata/CIF_FILES

BigData/ccdc_500_hist_vals.Rds:
	$(RUN_R) R/save_h2_hists.R $@ BigData/500cifs-CCDC
	# Process test subset of CCDC MOFs from Scotty 2017-11-07

BigData/Robj/tobacco_h2.Rds:
	${RUN_R} R/save_h2_hists.R $@ BigData/Convergence/tob-1.0-complete

BigData/Robj/tobacco_ch4.Rds:
	${RUN_R} R/save_h2_hists.R $@ BigData/tobacco-20171114/ch4grids use_ch4

BigData/Robj/hmof_h2.Rds:
	${RUN_R} R/save_h2_hists.R $@ BigData/2500hmof-data/h2grids

BigData/Robj/hmof_ch4.Rds:
	${RUN_R} R/save_h2_hists.R $@ BigData/2500hmof-data/ch4grids use_ch4

BigData/Robj/ad_hoc.Rds:
	${RUN_R} R/save_h2_hists.R $@ BigData/AdHoc
	
BigData/Convergence/top_20171122.Rds:
	rmdir --ignore-fail-on-non-empty BigData/Convergence/top_h05_20171122/*  # remove empty directories
	${RUN_R} R/save_h2_hists.R $@ BigData/Convergence/top_h05_20171122

