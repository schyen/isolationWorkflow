isolationWorkflow


install instructions:
```
# if devtools package not yet installed, install devtools first
install_packages('devtools')

# loading devtools library
library(devtools)

# download and install isolationWorkflow
install_github("schyen/isolationWorkflow", build_vignette = TRUE)

# you may need to install bioconductor dependencies
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

source("https://bioconductor.org/biocLite.R")
biocLite("DECIPHER")

# load isolationWorkflow package
library('isolationWorkflow')
```
To view the vignette:
`browseVignettes(package = "isolationWorkflow")`
And click `html`

Or you can see the instructions in the Wiki tab
