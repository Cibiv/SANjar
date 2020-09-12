#!/bin/bash

# Create shiny user & group with the uid/gid of the user running the container
export LD_PRELOAD='/usr/lib/libnss_wrapper.so'
export NSS_WRAPPER_PASSWD="$(mktemp)"
export NSS_WRAPPER_GROUP="$(mktemp)"
echo "shiny:x:$(id -u):$(id -g):Shiny Server:/state:/bin/false" > "$NSS_WRAPPER_PASSWD"
echo "shiny:x:$(id -g):" > "$NSS_WRAPPER_GROUP"

# Create necessary directories
mkdir -p /state/log /state/bookmarks

# Run server
exec $*
