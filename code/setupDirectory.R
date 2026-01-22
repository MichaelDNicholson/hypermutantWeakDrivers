# set up directory, once pulled from github


projectRoot <- paste0(normalizePath("."),"/") 

# Image sub directories
subdirs <- c("krasPublicData", "krasSimPlots", "krasMutationalBias", 
             "pcawg_trinuc_plots", "approxSimOutput", "paramRegimesAnalytic")

# Path to the main images directory
images_dir <- file.path(projectRoot, "images")

# Create the main images directory if it doesn't exist
if (!dir.exists(images_dir)) {
  dir.create(images_dir)
}

# Loop over subdirectories and create each one
for (sub in subdirs) {
  dir_path <- file.path(images_dir, sub)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path)
  }
}

cat("Created images directory and subdirectories successfully.\n")
