FROM openfaas/of-watchdog:0.5.3 as watchdog
FROM disarm/geopandas-r-base:1.0.0

COPY --from=watchdog /fwatchdog /usr/bin/fwatchdog
RUN chmod +x /usr/bin/fwatchdog

RUN mkdir -p /app

WORKDIR /app

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY index.py .
COPY config.py .
COPY preprocess_helpers.py .
COPY function function

# Populate example here - i.e. "cat", "sha512sum" or "node index.js"
ENV fprocess="python index.py"
# Set to true to see request in function logs
# ENV combine_output='false'
# ENV write_debug="true"
ENV write_timeout=6000
ENV read_timeout=6000
ENV exec_timeout=6000

EXPOSE 8080

HEALTHCHECK --interval=3s CMD [ -e /tmp/.lock ] || exit 1
CMD [ "fwatchdog" ]
