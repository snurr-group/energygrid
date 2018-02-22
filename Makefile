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
	$(RUN_R) R/save_h2_hists.R $@ BigData/Mateo/EnergyGrid/csd_h2

BigData/Robj/ccdc_hist_160k.Rds:
	$(RUN_R) R/save_h2_hists.R $@ BigData/Mateo/EnergyGrid/CSD-160K

BigData/ccdc_500_hist_vals.Rds:
	$(RUN_R) R/save_h2_hists.R $@ BigData/500cifs-CCDC
	# Process test subset of CCDC MOFs from Scotty 2017-11-07

BigData/Robj/tobacco_h2.Rds:
	${RUN_R} R/save_h2_hists.R $@ BigData/Mateo/EnergyGrid/TobH2

BigData/Robj/tobacco_ch4.Rds:
	${RUN_R} R/save_h2_hists.R $@ BigData/tobacco-20171114/ch4grids use_ch4

BigData/Robj/hmof_h2.Rds:
	${RUN_R} R/save_h2_hists.R $@ BigData/Mateo/EnergyGrid/2500hMOF

BigData/Robj/hmof_ch4.Rds:
	${RUN_R} R/save_h2_hists.R $@ BigData/Mateo/EnergyGrid/Methane2500 use_ch4

BigData/Robj/hmof_160k.Rds:
	${RUN_R} R/save_h2_hists.R $@ BigData/Mateo/EnergyGrid/Grids160K

BigData/Robj/adhoc_h2.Rds:
	${RUN_R} R/save_h2_hists.R $@ BigData/Mateo/EnergyGrid/AdHocGrids

BigData/Robj/adhoc_ch4.Rds:
	${RUN_R} R/save_h2_hists.R $@ BigData/Mateo/EnergyGrid/AdHocCH4 use_ch4

BigData/Robj/ad_hoc.Rds:
	${RUN_R} R/save_h2_hists.R $@ BigData/AdHoc
	
BigData/Convergence/top_20171122.Rds:
	rmdir --ignore-fail-on-non-empty BigData/Convergence/top_h05_20171122/*  # remove empty directories
	${RUN_R} R/save_h2_hists.R $@ BigData/Convergence/top_h05_20171122

BigData/Robj/hyper_tuned.Rds: R/hyper_tuned.R Notebooks/setup_data.R
	${RUN_R} $<

BigData/csd_formula.txt:
	#export BABEL_DATADIR=/cygdrive/c/Users/Benjamin/AppData/Roaming/OpenBabel-2.4.1/data; \
	# Need to make sure the data files have the linux line endings:
	export BABEL_DATADIR=/cygdrive/c/Users/Benjamin/Git/mofid/src/ob_datadir/; \
	rm $@; \
	cd /cygdrive/c/Users/Benjamin/Git/EnergyGrid/BigData/cifs_from_csddata; \
	obabel *.cif -ab -otxt --append "\tFORMULA" >> ../csd_formula.txt

BigData/Robj/csd_formulas.Rds: BigData/csd_formula.txt
	${RUN_R} R/filter_ccdc_chemistry.R

BigData/TopCCDC/top.smi: Figs/ccdc_screening.tex
	top3d=$$(tail -n +3 $< | cut -f1 -d' '); \
	rm -r BigData/TopCCDC; \
	mkdir BigData/TopCCDC; \
	mkdir -p Test/; \
	for structure in $$top3d; \
		do cif="$${structure}_stripped.cif"; \
		cp "BigData/cifs_from_csddata/$${cif}" BigData/TopCCDC; \
		python ../mofid/Python/extract_moffles.py "BigData/TopCCDC/$${cif}" >> $@; \
		done ;\
	obabel $@ -O BigData/TopCCDC/top.svg

