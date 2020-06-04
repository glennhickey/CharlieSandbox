FROM conda/miniconda3
WORKDIR /usr/src/app
RUN conda create -n hail python'>=3.6,<3.8'
SHELL ["conda", "run", "-n", "hail", "/bin/bash", "-c"]
RUN pip install --upgrade pip
RUN pip install cython
RUN pip install IPython==5.0
RUN conda install -c bioconda pyvcf
RUN pip install "urllib3<1.25,>=1.21.1" "requests<2.21.1,>=2.21.0"
RUN pip install hail
RUN pip install hgvs pyhgvs pyfaidx biocommons.seqrepo pyensembl numpy "scipy<1.4,>1.2" "tqdm==4.42.1" matplotlib ipython jupyter "pandas<0.26,>0.24" sympy nose
RUN apt-get -y update
RUN apt-get install -y rsync tabix bcftools vcftools openjdk-8-jre-headless
COPY . .
CMD [ "/bin/bash" ]

