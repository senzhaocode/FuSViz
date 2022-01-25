FROM rocker/shiny:4.1

LABEL maintainer="Sen ZHAO <t.cytotoxic@gmail.com>"

RUN apt-get update && apt-get install -y \
    sudo \
	unzip \
	wget \
	libbz2-dev \
	zlib1g-dev \
	git \
    libcurl4-openssl-dev \
	libcurl4-gnutls-dev \
	libxt-dev \
	libcairo2-dev \
	libsqlite3-dev \
    libpng-dev \
	libjpeg-dev \
	libxml2-dev \
    libssl-dev \
    libssh2-1-dev \
    && rm -rf /var/lib/apt/lists/*

RUN RUN installGithub.r "senzhaocode/FuSViz" && rm -rf /tmp/downloaded_packages/

RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" > /usr/lib/R/etc/Rprofile.site

RUN addgroup --system senzhao && adduser --system --ingroup senzhao senzhao

WORKDIR /home/senzhao

COPY inst/app/*.R .

RUN chown senzhao:senzhao -R /home/senzhao

USER senzhao

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/home/senzhao')"]