![FuSviz_overview](logo_white_back.png)

### Overview

A shiny app to visualise, interpret and prioritise genomic/transcriptomic structural variations (SVs) of multiple samples. It provides multiple solutions and an interactive user interface to investigate the prevalence and recurrence of SVs and their relevant genes in a cohort of cases. The tool is designed for combining SVs called from DNA-seq (e.g. whole genome or target exome) and RNA-seq to illustrate a biological implication of SVs to the host genes and interested genomic regions in context of various annotations, and it can also integrate a mutation profile (SNVs/Indels) to reveal a connection between small variants and complex genomic aberrations. 

<br>
<img align="center" width = "100%" height = "100%" border = 2, src="FuSViz_overview.gif"/>

### Getting a quick start

#### Installation & launch

    if (! require('remotes')) install.packages('remotes')
    remotes::install_github('senzhaocode/FuSViz')
    
    source(file.path(system.file("app", package = "FuSViz"), "global.R"), local = TRUE, chdir = TRUE)
    FuSViz_app()

##### NOTE: some libraries need to be installed in Linux properly before FuSViz setup.

1. Install `openssl` a dependency for R package 'httpr'

    * For Debian or Ubuntu OS, `sudo apt-get install -y libssl-dev`
    * For Fedora, CentOS or RHEL OS, `sudo yum install openssl-devel`
    
#### Usage & manual

A full description of FuSViz documentation is available: &nbsp;&nbsp; [![Documentation Status](https://readthedocs.org/projects/fusviz-docs/badge/?version=latest)](https://fusviz-docs.readthedocs.io/en/latest/index.html)

### Contact

t.cytotoxic AT gmail.com

 
