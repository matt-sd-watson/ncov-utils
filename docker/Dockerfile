FROM ubuntu
LABEL software="Auspice from PHO BCC dev"
LABEL software.version="2.25.1"
LABEL description="A web-based viewer for phylodynamics based on Nextstrain"
LABEL website="https://repo.coreontario.ca/COVID-19/nextstrain"
LABEL maintainer="Public Health Ontario"

ENTRYPOINT ["Building the Auspice phylogenetic viewer"]
RUN apt-get update && apt-get install -y git && \
    apt-get install -y curl && apt-get install -y sudo && \
    sudo apt-get install -y apt-utils && \
    curl -sL https://deb.nodesource.com/setup_15.x | sudo -E bash - && \
    sudo apt-get install -y nodejs && \
    sudo apt-get update && \
    sudo apt-get upgrade && \
    apt-get install -y npm && \
    git clone https://github.com/nextstrain/auspice.git && cd auspice && npm install argparse && \
	npm install express && npm install && npm install --global . 
