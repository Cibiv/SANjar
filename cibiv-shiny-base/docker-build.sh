#!/bin/bash
function die() {
	echo $*
	exit 1
}

DOCKER_REGISTRY="minsky:5000"
CONTAINER_NAME="cibiv-shiny-base-r40"

img_tag="$CONTAINER_NAME:$(date +'%Y-%m-%d')"
echo "*** Building image $"
docker build -t "$img_tag" "$(dirname $BASH_SOURCE)" || die "image build failed"
img_hash="$(docker images -q "$img_tag")"

publish_tag="$DOCKER_REGISTRY/$img_tag"
echo "*** Tagging image $img_tag (hash $img_hash) for publication as $publish_tag"
docker tag "$img_tag" "$publish_tag" || die "tagging for publication failed"

echo "*** Publishing $publish_tag"
docker push "$publish_tag" || die "image push failed"
