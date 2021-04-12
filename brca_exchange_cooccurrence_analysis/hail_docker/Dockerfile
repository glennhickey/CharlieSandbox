FROM python:3.6.9-slim-stretch
RUN mkdir -p /usr/share/man/man1 && \
    apt-get update && apt-get install -y \
    openjdk-8-jre-headless bcftools tabix \
    && rm -rf /var/lib/apt/lists/* && \
    pip3 --no-cache-dir install hail==0.2.64

ENTRYPOINT ["/bin/bash"]

