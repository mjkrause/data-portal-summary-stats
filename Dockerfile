FROM python:3.6.9-slim-stretch
SHELL ["/bin/bash", "-c"]
RUN mkdir /build
WORKDIR /build
COPY . /build
RUN pip install --upgrade pip \
    && pip install -r requirements.txt
ENTRYPOINT [ "python", "./data_portal_summary_stats.py"]
