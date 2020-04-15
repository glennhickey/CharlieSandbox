FROM python:3
WORKDIR /usr/src/app
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install cython
RUN python3 -m pip install IPython==5.0
RUN python3 -m pip install requests pyvcf hgvs pyhgvs pyfaidx biocommons.seqrepo pyensembl
RUN apt-get -y update
RUN apt-get install -y rsync tabix bcftools vcftools
COPY . .
CMD [ "/bin/bash" ]

