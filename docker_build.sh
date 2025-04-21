#!/usr/bin/env bash

set -e

IMAGE=neoformit/daff-taxonomic-assignment

read -p "Enter the tag for the image (default: latest): " TAG
TAG=${TAG:-latest}

# Build the Docker image
docker build -t $IMAGE:$TAG .
docker tag $IMAGE:$TAG $IMAGE:latest

# if -p in args, push the image to the registry
if [[ $* == *-p* ]]; then
  docker push $IMAGE:$TAG
  docker push $IMAGE:latest
fi
