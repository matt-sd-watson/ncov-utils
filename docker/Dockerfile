FROM node:latest
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

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y apt-utils && apt-get install -y git && \
    apt-get install -y software-properties-common && \
    apt-get install -y curl && apt-get install -y sudo && \
    sudo apt-get update -y && \
    sudo apt install -y aptitude && \
    sudo aptitude install -y libnode-dev && \                        
    sudo aptitude install -y libnode64 && \
    sudo aptitude install -y node-gyp && \
    sudo apt-get update -y && \
    sudo apt-get upgrade -y && \
    sudo aptitude install -y npm && \
    git clone https://github.com/matt-sd-watson/auspice.git && cd auspice && npm install && npm install --global . && \
    which auspice

ENV PATH="/usr/local/bin/auspice:$PATH"

WORKDIR /app
