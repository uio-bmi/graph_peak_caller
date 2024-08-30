FROM ubuntu:20.04
ENV PIP_BREAK_SYSTEM_PACKAGES 1
ENV PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata parallel
RUN apt install --assume-yes git python3 python3-pip pkg-config build-essential wget
RUN pip3 install matplotlib
RUN git clone https://github.com/cgroza/graph_peak_caller.git \
RUN cd graph_peak_caller  && pip3 install .

RUN wget https://github.com/vgteam/vg/releases/download/v1.58.0/vg -O /usr/bin/vg
RUN chmod +x /usr/bin/vg

ENTRYPOINT [ "/bin/bash", "-l", "-c" ]
