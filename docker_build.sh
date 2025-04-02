#!/usr/bin/env bash

IMAGE=neoformit/daff-taxonomic-assignment

read -p "Enter the tag for the image (default: latest): " TAG
TAG=${TAG:-latest}

# Build the Docker image
docker build -t $IMAGE:$TAG .

# if -p in args, push the image to the registry
if [[ $* == *-p* ]]; then
  docker push $IMAGE:$TAG
fi
