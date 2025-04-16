#!/usr/bin/env bash

set -e

IMAGE=neoformit/daff-taxonomic-assignment

read -p "Enter the tag for the image (default: latest): " TAG
TAG=${TAG:-latest}

read -p "Quick build? (update code only) [Y/n]: " QUICK_BUILD
if [[ $QUICK_BUILD =~ ^[Yy]$ ]]; then
  DOCKERFILE=Dockerfile.code_only
else
  DOCKERFILE=Dockerfile
fi

# Build the Docker image
docker build -t $IMAGE:$TAG -f $DOCKERFILE .

# if -p in args, push the image to the registry
if [[ $* == *-p* ]]; then
  docker tag $IMAGE:$TAG $IMAGE:latest
  docker push $IMAGE:$TAG
  docker push $IMAGE:latest
fi
