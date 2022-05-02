doc:
	Rscript -e "library(roxygen2); roxygenize('.'); devtools::build_manual()"