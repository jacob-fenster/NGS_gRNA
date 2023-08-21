#this is to build a container than can run the BBMap suite

#load ubuntu linux 20.04 on an amd64 platform
FROM --platform=linux/amd64 ubuntu:20.04
WORKDIR /app
#command to not have interactive mode to solve java install
ARG DEBIAN_FRONTEND=noninteractive
# Install base utilities
RUN apt-get update \
    && apt-get install -y build-essential \
    && apt-get install -y wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* 
#Install java
RUN apt-get update \
    && apt-get install -y openjdk-11-jre-headless 
#copy bbmap script suite to the /opt linux folder
COPY software/bbmap /opt
#add the /opt directory to the PATH variable
ENV PATH=/opt:$PATH