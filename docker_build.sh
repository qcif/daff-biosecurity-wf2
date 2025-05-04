#!/usr/bin/env bash

set -e

IMAGE=neoformit/daff-taxonomic-assignment

# Check for -t argument and set TAG if provided
while getopts "t:" opt; do
  case $opt in
    t)
      TAG=$OPTARG
      ;;
    *)
      ;;
  esac
done

if [[ -z $TAG ]]; then
  # Prompt for the tag if not provided
  read -p "Enter the tag for the image (default: latest): " TAG
  TAG=${TAG:-latest}
fi

# Build the Docker image
docker build -t $IMAGE:$TAG .
docker tag $IMAGE:$TAG $IMAGE:latest

# if -p in args, push the image to the registry
if [[ $* == *-p* ]]; then
  docker push $IMAGE:$TAG
  docker push $IMAGE:latest
fi
