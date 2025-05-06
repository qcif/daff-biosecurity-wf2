FROM python:3.12

ENV COX1_HMM_PROFILE=/hmm/pf00115.hmm

RUN apt update && apt install -y nano less

# Install taxonkit
RUN wget https://github.com/shenwei356/taxonkit/releases/latest/download/taxonkit_linux_amd64.tar.gz \
    && tar -zxvf taxonkit_linux_amd64.tar.gz \
    && mv taxonkit /usr/local/bin/ \
    && rm taxonkit_linux_amd64.tar.gz

# Install HMMSearch
RUN wget http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz \
    && tar -xzf hmmer-3.4.tar.gz \
    && rm hmmer-3.4.tar.gz \
    && cd hmmer-3.4 \
    && ./configure --prefix=/opt/hmmer \
    && make \
    && make install \
    && strip /opt/hmmer/bin/hmmsearch \
    && find /opt/hmmer/bin ! -name 'hmmsearch' -type f -delete \
    && cp /opt/hmmer/bin/hmmsearch /usr/local/bin \
    && strip /usr/local/bin/hmmsearch \
    && cd .. \
    && rm -rf hmmer-3.4/

# Install HMM profile for COX1
RUN wget https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00115?annotation=hmm -O pf00115.hmm.gz \
    && gunzip pf00115.hmm.gz \
    && mkdir /hmm \
    && mv pf00115.hmm $COX1_HMM_PROFILE

COPY . /app

RUN pip install -r /app/requirements.txt
