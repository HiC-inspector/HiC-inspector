install.packages('gplots',dependencies=TRUE, repos='http://cran.rstudio.com/');
install.packages('MCMCpack',dependencies=TRUE, repos='http://cran.rstudio.com/');
install.packages('RColorBrewer',dependencies=TRUE, repos='http://cran.rstudio.com/'); 

install.packages("BiocManager",dependencies=TRUE, repos='http://cran.rstudio.com/');
BiocManager::install(c("Heatplus", "HiTC"))
