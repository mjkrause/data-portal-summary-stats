FROM python:3.6.9-slim-stretch
SHELL ["/bin/bash", "-c"]

RUN mkdir /build
RUN mkdir /build/src

WORKDIR /build

COPY . /build
RUN sudo apt install build-essential libxml2-dev libglpk-dev libgmp3-dev libblas-dev liblapack-dev libarpack2-dev python-dev
RUN pip install --upgrade pip \
    && pip install -r requirements.txt

ENTRYPOINT [ "python", "./data_portal_summary_stats.py"]
