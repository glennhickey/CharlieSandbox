FROM ubuntu:20.10

WORKDIR /usr/src/app

COPY Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz .
COPY Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz .
COPY Eagle_v2.4.1/eagle .
COPY Eagle_v2.4.1/src/ .

CMD [ "/bin/bash" ]

