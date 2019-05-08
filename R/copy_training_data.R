# Prepare hMOFs for LJ grid calculations and NN training data

HMOFS_PATH <- "C:/Users/Benjamin/Desktop/Jiayi/Files/Dataset Comparison/hMOF"
TRAINING_COUNT <- 10000  # Let's start the proof-of-concept with 10000 MOFs. 
OUT_DIR <- "Python/gridcalc/CIF_FILES"


extract_hmof_id <- function(file_name) {
  # Extract an hMOF ID number from its filename
  hmof_id <- sub("hypotheticalMOF_(\\d+)_.*", "\\1", file_name)
}



if (length(list.files(OUT_DIR)) != 0) {
  stop(paste("Output directory", OUT_DIR, "is not empty!"))
}

# Randomly select hMOFs for analysis
set.seed(20160912)
hMOFs <- list.files(HMOFS_PATH)
train <- sample(1:length(hMOFs), TRAINING_COUNT)

# Copy files and set up the Python directory structure
for (i in train) {
  mof_name <- hMOFs[i]
  mof_id <- extract_hmof_id(mof_name)
  file_name <- paste0("h", mof_id)
  out_path <- paste0(OUT_DIR, "/", file_name)
  dir.create(out_path)
  file.copy(paste0(HMOFS_PATH, "/", mof_name), paste0(out_path, "/", file_name, ".cif"))
  
  # Also save the i, j, k, m, cat ID information
  cif_basename <- sub(".cif", "", mof_name, fixed=TRUE)
  cat("", file=paste0(out_path, "/", cif_basename, ".txt"))
}

