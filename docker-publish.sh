#!/bin/bash
function die() {
	echo $*
	exit 1
}

if [ "$(git symbolic-ref --short HEAD)" != "master" ]; then
	echo "Currently checked out branch must be 'master'" >&2
        exit 1
fi

if [ $(git status --porcelain | wc -l) != 0 ]; then
	echo "Working copy contains uncommitted changes" >&2
	exit 1
fi

git fetch --tags || die "Unable to fetch tags"
rev=$[ $(git tag | sed -n 's/^rev\([0-9]*\)$/\1/p' | sort -n | tail -n1) + 1 ]
echo "*** Publishing revision $rev"

echo "*** Tagging revision"
git tag rev$rev || die "Unable to tag current version as rev$rev"
git push origin master rev$rev

img_tag="sanmodelexplorer:r$rev"
echo "*** Building image $img_tag"
docker build --build-arg CONTAINER_REVISION="rev$rev" -t "$img_tag" "$(dirname $BASH_SOURCE)" || exit 1
img_hash="$(docker images -q "$img_tag")"

publish_rev_tag="minsky:5000/$img_tag"
publish_latest_tag="minsky:5000/sanmodelexplorer:latest"
echo "*** Tagging image $img_tag (hash $img_hash) for publication as $publish_rev_tag and as $publish_latest_tag"
docker tag "$img_tag" "$publish_rev_tag" || exit 1
docker tag "$img_tag" "$publish_latest_tag" || exit 1

echo "*** Publishing $publish_rev_tag"
docker push "$publish_rev_tag" || exit 1
docker push "$publish_latest_tag" || exit 1

echo "*** Pulling new image on minsky, shinyproxy should use it automatically"
ssh minsky docker pull --quiet "$publish_latest_tag" || exit 1
