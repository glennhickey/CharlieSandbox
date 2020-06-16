FROM ubuntu:20.10

WORKDIR /usr/src/app

RUN wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz .
RUN tar xvzf Eagle_v2.4.1.tar.gz
COPY ./Eagle_v2.4.1/eagle .
COPY ./Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz .
COPY ./Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz .
RUN rm -fr Eagle_v2.4.1

CMD [ "/bin/bash" ]

