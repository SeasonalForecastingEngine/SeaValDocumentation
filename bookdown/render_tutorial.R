### This script renders the SeaVal documentation

# It essentially knits all .rmd files in this folder, starting with index.rmd and then following the naming order (so the file names should start with a number)

### adjust to your setup: ###
setwd('~/pkg/SeaValDocumentation/bookdown/')
data_dir = '/nr/project/stat/CONFER/Data/SeaVal/example_data/202102/'
#chirps_fn = "/nr/project/stat/CONFER/Data/CHIRPS_prec_upscaled.csv"
output_dir = '/nr/common/www/virtual/files.nr.no/htdocs/samba/CONFER/SeaVal/'
# this directory is shared at http://files.nr.no/samba/CONFER/SeaVal/

#save(data_dir,chirps_fn,file = '_data_dir.RData') # just for loading it in the bookdown
save(data_dir,file = '_data_dir.RData') # just for loading it in the bookdown


bookdown::render_book('index.rmd',
                      'bookdown::gitbook',
                      new_session = T)


# We need to copy paste, the output_dir option of render_book does not work, see https://github.com/rstudio/bookdown/issues/804


file.copy(list.files('./_book/',full.names = TRUE), to = output_dir, recursive = TRUE)

