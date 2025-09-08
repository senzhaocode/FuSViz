FROM rocker/shiny:4.4.1

LABEL maintainer="Sen ZHAO <t.cytotoxic@gmail.com>"

RUN dpkg --configure -a

RUN apt update && apt install -y -f \
	build-essential \
	sudo \
	unzip \
	wget \
	curl \
	git \
	libbz2-dev \
	libgsl0-dev \
	liblzma-dev \
	libglpk-dev \
	libncurses5-dev \
	libperl-dev \
	zlib1g-dev \
	libcurl4-openssl-dev \
	libxt-dev \
	libcairo2-dev \
	libsqlite3-dev \
	libpng-dev \
	libjpeg-dev \
	libxml2-dev \
	libssl-dev \
	libssh2-1-dev \
    && rm -rf /var/lib/apt/lists/*

ENV LIBRARY_PATH="$LIBRARY_PATH:/usr/local/lib/R/lib/"
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/R/lib/"

WORKDIR /tmp
ARG htsversion=1.19
RUN curl -L https://github.com/samtools/htslib/releases/download/${htsversion}/htslib-${htsversion}.tar.bz2 | tar xj && \
    (cd htslib-${htsversion} && ./configure --enable-plugins --with-plugin-path='$(libexecdir)/htslib:/usr/libexec/htslib' && make install) && \
    ldconfig && \
    curl -L https://github.com/samtools/samtools/releases/download/${htsversion}/samtools-${htsversion}.tar.bz2 | tar xj && \
    (cd samtools-${htsversion} && ./configure --with-htslib=system && make install) && \
    curl -L https://github.com/samtools/bcftools/releases/download/${htsversion}/bcftools-${htsversion}.tar.bz2 | tar xj && \
    (cd bcftools-${htsversion} && ./configure --enable-libgsl --enable-perl-filters --with-htslib=system && make install)

RUN install2.r -e remotes
RUN install2.r -e devtools
RUN R -e "install.packages('BiocManager'); BiocManager::install('Rhtslib'); BiocManager::install('GenomicFeatures'); BiocManager::install('AnnotationDbi'); BiocManager::install('GenomicRanges'); BiocManager::install('IRanges'); BiocManager::install('S4Vectors'); BiocManager::install('GenomicAlignments'); BiocManager::install('Rsamtools'); BiocManager::install('Biostrings'); BiocManager::install('Gviz');"
# RUN installGithub.r "lchiffon/wordcloud2"
# RUN installGithub.r "senzhaocode/FuSViz"
RUN wget -t 0 -c "https://github.com/senzhaocode/FuSViz/archive/refs/tags/v2.0.0.tar.gz" && R -e "remotes::install_local('v2.0.0.tar.gz', dependencies=T)"
RUN rm -rf /tmp/bcftools* && rm -rf /tmp/htslib-* && rm -rf /tmp/samtools-* && rm -rf /tmp/file* && rm -rf /tmp/*.tar.gz

RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" >> /usr/local/lib/R/etc/Rprofile.site

# Set Volume
RUN mkdir /data && chmod 777 /data
VOLUME /data

# Final clean
RUN apt-get clean autoclean
RUN rm -rf /var/tmp/*
RUN rm -rf /tmp/downloaded_packages

# set portal
EXPOSE 3838

CMD ["R", "-e", "shiny::runApp(file.path(system.file('app', package='FuSViz')))"]

