#!/bin/bash

# Run docker
docker ps -q | xargs -r docker stop
docker ps -aq | xargs -r docker rm
docker images -q | xargs -r docker rmi
docker build -t rendrick27/snakemake .
docker push rendrick27/snakemake:latest
docker run -it --name app-container rendrick27/snakemake /bin/bash
