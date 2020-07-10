#!/bin/sh
wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz
tar xvzf Eagle_v2.4.1.tar.gz
docker build -t quay.io/cmarkello/eagle .
#rm -fr Eagle_v2.4.1 Eagle_v2.4.1.tar.gz
