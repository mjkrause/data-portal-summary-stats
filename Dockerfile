FROM python:3.6.9-slim-stretch
SHELL ["/bin/bash", "-c"]

RUN mkdir /build
RUN mkdir /build/src

WORKDIR /build

COPY requirements.txt \
     data_portal_summary_stats.py /build/

COPY src/* /build/src/

RUN pip install --upgrade pip \
    && pip install -r requirements.txt

ENTRYPOINT [ "python", "./data_portal_summary_stats.py"]
