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

1. Install a software library [OpenSSL](https://www.openssl.org) - a depedency of R package ['openssl'](https://github.com/jeroen/openssl)

    * For Debian or Ubuntu OS, `sudo apt-get install -y libssl-dev`
    * For Fedora, CentOS or RHEL OS, `sudo yum install openssl-devel`
    * If root privillege is not available, users have to download source code and install at their HOME directory. For example,
    
            ./Configure --prefix=/tsd/p1532/home/p1532-senz/openssl --openssldir=/tsd/p1532/home/p1532-senz/openssl/ssl
            make && install
            C_INCLUDE_PATH=/tsd/p1532/home/p1532-senz/openssl/include
            export C_INCLUDE_PATH
            LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/tsd/p1532/home/p1532-senz/openssl/lib
            export LD_LIBRARY_PATH
            
    * Install R package 'openssl': `install.packages("openssl")`
2. install libjpeg
    
#### Usage & manual

A full description of FuSViz documentation is available: &nbsp;&nbsp; [![Documentation Status](https://readthedocs.org/projects/fusviz-docs/badge/?version=latest)](https://fusviz-docs.readthedocs.io/en/latest/index.html)

### Contact

t.cytotoxic AT gmail.com

 
