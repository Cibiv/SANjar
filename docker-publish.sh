#!/bin/bash
function die() {
	echo $*
	exit 1
}


DOCKER_REGISTRY="minsky:5000"
CONTAINER_NAME="sanmodelexplorer"
APPLICATION_SERVER="minsky"

if [ "$(git symbolic-ref --short HEAD)" != "master" ]; then
	echo "Currently checked out branch must be 'master'" >&2
        exit 1
fi

if [ $(git status --porcelain | wc -l) != 0 ]; then
	echo "Working copy contains uncommitted changes" >&2
	exit 1
fi

echo "*** Fetching tags"
git fetch --tags || die "Unable to fetch tags"

rev=$(git tag --points-at HEAD | sed -n 's/^rev\([0-9]*\)$/\1/p')
if [ "$rev" == "" ]; then
	rev=$[ $(git tag | sed -n 's/^rev\([0-9]*\)$/\1/p' | sort -n | tail -n1) + 1 ]
	echo "*** Publishing new revision $rev"
	echo "*** Tagging HEAD as rev$rev and pushing tag"
	git tag rev$rev || die "Unable to tag current version as rev$rev"
	git push origin master rev$rev || die "Unable to push rev$rev to origin"
else
	echo "*** Publishing existing revision $rev"
fi

img_tag="$CONTAINER_NAME:r$rev"
echo "*** Building image $img_tag"
docker build --build-arg CONTAINER_REVISION="rev$rev" -t "$img_tag" "$(dirname $BASH_SOURCE)" || exit 1
img_hash="$(docker images -q "$img_tag")"

publish_rev_tag="$DOCKER_REGISTRY/$img_tag"
publish_latest_tag="$DOCKER_REGISTRY/$CONTAINER_NAME:latest"
echo "*** Tagging image $img_tag (hash $img_hash) for publication as $publish_rev_tag and as $publish_latest_tag"
docker tag "$img_tag" "$publish_rev_tag" || exit 1
docker tag "$img_tag" "$publish_latest_tag" || exit 1

echo "*** Publishing $publish_rev_tag"
docker push "$publish_rev_tag" || exit 1
docker push "$publish_latest_tag" || exit 1

echo "*** Pulling new image on $APPLICATION_SERVER, shinyproxy should use it automatically"
ssh "$APPLICATION_SERVER" docker pull --quiet "$publish_latest_tag" || exit 1
