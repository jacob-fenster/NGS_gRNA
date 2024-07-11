# Use an ARM64 architecture version of Ubuntu 22.04 for apple silicon chips
FROM --platform=linux/arm64/v8 ubuntu:22.04
#make the working directory /app
WORKDIR /app
# Command to disable interactive mode during package installation
ARG DEBIAN_FRONTEND=noninteractive
# Install base utilities
RUN apt-get update \
    && apt-get install -y build-essential \
    && apt-get install -y wget \
    && apt-get install -y zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev \
    && apt-get clean
# Download and install Python 3.11.4 from source code
RUN wget https://www.python.org/ftp/python/3.11.4/Python-3.11.4.tgz \
    && tar -xvf Python-3.11.4.tgz \
    && cd Python-3.11.4 \
    && ./configure --enable-optimizations \
    && make altinstall \
    && cd .. \
    && rm -rf Python-3.11.4 \
    && rm Python-3.11.4.tgz
# Update alternatives and Install pip
RUN update-alternatives --install /usr/bin/python3 python3 /usr/local/bin/python3.11 1 \
    && wget https://bootstrap.pypa.io/get-pip.py \
    && python3 get-pip.py \
    && rm get-pip.py
#install python packages
RUN pip install pandas==2.1.1 \
    && pip install biopython==1.81 && pip install matplotlib==3.8.0
#Install java
RUN apt-get update \
    && apt-get install -y openjdk-11-jre-headless 
#copy bbmap software suite folder to the /opt linux folder
COPY software/bbmap /opt
#add the /opt directory to the PATH variable
ENV PATH=/opt:$PATH
RUN chmod +x /opt/bbmerge.sh
#remove chache of software lists
RUN rm -rf /var/lib/apt/lists/*