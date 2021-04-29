LABEL software="Auspice from PHO BCC dev"
LABEL software.version="2.25.1"
LABEL description="A web-based viewer for phylodynamics based on Nextstrain"
LABEL website="https://repo.coreontario.ca/COVID-19/nextstrain"
LABEL maintainer-"Public Health Ontario"

RUN git clone https://repo.coreontario.ca/COVID-19/nextstrain.git && cd nextstrain && npm install argparse && \
	npm install express && npm install && npm install --global . 

ENV PATH /usr/local/auspice/:$PATH
