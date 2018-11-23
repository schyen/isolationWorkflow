isolationWorkflow


install instructions:
```
# if devtools package not yet installed, install devtools first
install.packages('devtools')

# loading devtools library
library(devtools)

# download and install isolationWorkflow
install_github("schyen/isolationWorkflow", build_vignette = TRUE)

# you may need to install bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DECIPHER", version = "3.8")

# load isolationWorkflow package
library('isolationWorkflow')
```
To view the vignette:
`browseVignettes(package = "isolationWorkflow")`
And click `html`

Or you can see the instructions in the Wiki tab
