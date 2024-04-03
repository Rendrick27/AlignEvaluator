#!/bin/bash

# Run docker without publishing
docker ps -q | xargs -r docker stop
docker ps -aq | xargs -r docker rm
docker images -q | xargs -r docker rmi
docker build -t rendrick27/snakemake .
docker run -it --name app-container rendrick27/snakemake /bin/bash


////

# Stop all running containers
docker ps -q | ForEach-Object { docker stop $_ }

# Remove all containers
docker ps -aq | ForEach-Object { docker rm $_ }

# Remove all images
docker images -q | ForEach-Object { docker rmi $_ }

# Build the Docker image
docker build -t rendrick27/snakemake .

# Run the Docker container
docker run -it --name app-container rendrick27/snakemake /bin/bash
