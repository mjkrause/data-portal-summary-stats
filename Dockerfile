FROM python:3.6.9-slim-stretch
SHELL ["/bin/bash", "-c"]
RUN mkdir /build
WORKDIR /build
COPY . /build
RUN pip install -r requirements.txt
# Run the following script when container launches.
CMD [ "python", "./data_portal_summary_stats.py" ]