FROM openjdk:8-slim
COPY --from=python:3.7 / /
WORKDIR /usr/src/app
RUN mkdir -p /usr/share/man/man1 && \
    apt-get update && apt-get install -y rsync tabix bcftools vcftools \
    && rm -rf /var/lib/apt/lists/* && \
    pip3 --no-cache-dir install cython IPython==5.0 requests pyvcf hgvs pyhgvs pyfaidx biocommons.seqrepo pyensembl numpy scipy matplotlib jupyter pandas sympy nose
COPY . .
CMD [ "/bin/bash" ]

