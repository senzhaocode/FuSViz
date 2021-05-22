![FuSviz_overview](logo_white_back.png)

### Overview

A shiny app to visualise, interpret and prioritise genomic/transcriptomic structural variations (SVs) of multiple samples. It provides multiple solutions and an interactive user interface to investigate the prevalence and recurrence of SVs and their relevant genes in a cohort of cases. The tool is designed for combining SVs called from DNA-seq (e.g. whole genome or target exome) and RNA-seq to illustrate a biological implication of SVs to the host genes and interested genomic regions in context of various annotations, and it can also integrate a mutation profile (SNVs/Indels) to reveal a connection between small variants and complex genomic aberrations. 

<br>
<img align="center" width = "100%" height = "100%" border = 2, src="FuSViz_overview.gif"/>

### Getting a quick start

#### Prerequisite

* R (>=4.0.0): https://www.r-project.org/; [RStudio](https://rstudio.com/products/rstudio/download/#download) is recommended but not mandatory. 
* For windows users, if an earlier version of R (< 4.0) is present in the system, please uninstall it firstly and make sure only R >=4.0 is available.

#### Installation

    if (! require('remotes')) install.packages('remotes')
    remotes::install_github('senzhaocode/FuSViz')

#### IMPORTANT NOTE for Linux OS, some libraries need to be installed properly before setup FuSViz.

1. Install software library [OpenSSL](https://www.openssl.org) - a dependency of R package [openssl](https://cran.r-project.org/web/packages/openssl/index.html)

    * For **Debian or Ubuntu**: `sudo apt-get install -y libssl-dev`; For **Fedora, CentOS or RHEL**: `sudo yum install openssl-devel`
    * If root privillege is not available, users have to download [source code](https://github.com/openssl/openssl) and install at $HOME directory. For example,
    
            ./Configure --prefix=/OpenSSL_path --openssldir=/OpenSSL_path/ssl
            make && make install
            C_INCLUDE_PATH=/OpenSSL_path/include
            export C_INCLUDE_PATH
            LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/OpenSSL_path/lib
            export LD_LIBRARY_PATH
            
    * Install R package [openssl](https://cran.r-project.org/web/packages/openssl/index.html): `install.packages("openssl")`

2. Install software library [libxml2](http://xmlsoft.org) - a dependency of R package [xml2](https://cran.r-project.org/web/packages/XML/index.html)

    * For **Debian or Ubuntu**: `sudo apt-get install libxml2-dev`; For **Fedora, CentOS or RHEL**: `sudo yum install libxml2-devel`
    * If root privillege is not available, users have to download [source code](http://xmlsoft.org/downloads.html) and install at $HOME directory. For example,
    
            ./configure --prefix=/libxml2_path # if ./configure file does not exist, please run ./autogen.sh --prefix=/libxml2_path instead.
            make && make install
            C_INCLUDE_PATH=/libxml2_path/include
            export C_INCLUDE_PATH
            CPLUS_INCLUDE_PATH=/libxml2_path/include
            export CPLUS_INCLUDE_PATH
            LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/libxml2_path/lib
            export LD_LIBRARY_PATH
    
    * Install R package [xml2](https://cran.r-project.org/web/packages/XML/index.html): `install.packages("xml2")`

3. Install software library [libjpeg](https://ijg.org) - a dependency of R package [jpeg](https://cran.r-project.org/web/packages/jpeg/index.html)

    * For **Debian or Ubuntu**: `sudo apt-get install libjpeg-dev`; For **Fedora, CentOS or RHEL**: `sudo yum install libjpeg-turbo-devel`
    * If root privillege is not available, users have to download [source code](https://ijg.org) and install at $HOME directory. For example,
    
            ./Configure --prefix=/libjpeg_path --libdir=/libjpeg_path/lib --includedir=/libjpeg_path/include
            make && make install
            C_INCLUDE_PATH=/libjpeg_path/include
            export C_INCLUDE_PATH
            LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/libjpeg_path/lib
            export LD_LIBRARY_PATH
            
     * Install R package [jpeg](https://cran.r-project.org/web/packages/jpeg/index.html): `install.packages("jpeg")`

4. Install software library [libpng](https://libpng.sourceforge.io) - a dependency of R package [png](https://cran.r-project.org/web/packages/png/index.html)

    * For **Debian or Ubuntu**: `sudo apt-get install libpng-dev`; For **Fedora, CentOS or RHEL**: `sudo yum install libpng-devel`
    * If root privillege is not available, users have to download [source code](https://libpng.sourceforge.io) and install at $HOME directory. For example,

            ./Configure --prefix=/libpng_path
            make && make install
            C_INCLUDE_PATH=/libpng_path/include
            export C_INCLUDE_PATH
            LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/libpng_path/lib
            export LD_LIBRARY_PATH
            
     * Install R package [png](https://cran.r-project.org/web/packages/png/index.html): `install.packages("png")`
 
5. Install software library [libcurl](https://curl.se/libcurl/) - a dependency of R package [RCurl](https://cran.r-project.org/web/packages/RCurl/index.html)

    * Need for **Debian or Ubuntu**: `sudo apt install libcurl4-openssl-dev`

6. If users have a problem to install the latest version of one dependency R package **stringi** (v1.6.2) for Linux OS, please install the earlier one (e.g. v1.5.3) before setup FuSViz: `install.packages("https://cran.r-project.org/src/contrib/Archive/stringi/stringi_1.5.3.tar.gz", repos=NULL, type="source")`

#### Launch FuSViz App

    source(file.path(system.file("app", package = "FuSViz"), "global.R"), local = TRUE, chdir = TRUE)
    FuSViz_app()

#### Usage & manual

A full description of FuSViz documentation is available: &nbsp;&nbsp; [![Documentation Status](https://readthedocs.org/projects/fusviz-docs/badge/?version=latest)](https://fusviz-docs.readthedocs.io/en/latest/index.html)

### Contact

t.cytotoxic AT gmail.com

 
