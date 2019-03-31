# isolationWorkflow

R package to help with workflow of building a strain library during isolation experiments. Allows you to preview, trim Sanger sequences before organizing into a fasta file that is ready to BLAST. Reads Blast results, and aligns sequences in order to build and curate a strain library

Install Instructions:
```
# if devtools package not yet installed, install devtools first
install.packages('devtools')

# you may need to install bioconductor dependencies if you have not done so previously
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DECIPHER", version = "3.8")

# loading devtools library
library(devtools)

# download and install isolationWorkflow
install_github("schyen/isolationWorkflow", build_vignette = TRUE)

# load isolationWorkflow package
library('isolationWorkflow')
```
To view the vignette:
`browseVignettes(package = "isolationWorkflow")`
And click `html`

Or you can see the instructions in the Wiki tab
