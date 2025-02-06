FROM python:3.12

# Install taxonkit
RUN wget https://github.com/shenwei356/taxonkit/releases/latest/download/taxonkit_linux_amd64.tar.gz \
    && tar -zxvf taxonkit_linux_amd64.tar.gz \
    && mv taxonkit /usr/local/bin/ \
    && rm taxonkit_linux_amd64.tar.gz

COPY . /app

RUN pip install -r /app/requirements.txt
