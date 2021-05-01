FROM ubuntu:18.04
LABEL software="Auspice from PHO BCC dev"
LABEL software.version="2.25.1"
LABEL description="A web-based viewer for phylodynamics based on Nextstrain"
LABEL website="https://repo.coreontario.ca/COVID-19/nextstrain"
LABEL maintainer="Public Health Ontario"

ENTRYPOINT ["Building the Auspice phylogenetic viewer"]
RUN apt-get update && brew instll node && \
    apt-get upgrade -y && \
    apt-get install -y git && \
    git clone https://github.com/matt-sd-watson/auspice && cd auspice && npm install argparse && \
	npm install express && npm install && npm install --global . 
