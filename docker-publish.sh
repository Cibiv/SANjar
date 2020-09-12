#!/bin/bash
function die() {
	echo $*
	exit 1
}

if [ $(git status --porcelain | wc -l) != 0 ]; then
	echo "Working copy contains uncommitted changes" >&2
	exit 1
fi

git fetch --tags || die "Unable to fetch tags"
rev=$[ $(git tag | sed -n 's/^rev\([0-9]*\)$/\1/p' | sort -n | tail -n1) + 1 ]
echo "*** Publishing revision $rev"

echo "*** Tagging revision"
git tag rev$rev || die "Unable to tag current version as rev$rev"
git push origin rev$rev

img_tag="sanmodelexplorer:r$rev"
echo "*** Building image $img_tag"
docker build -t "$img_tag" "$(dirname $BASH_SOURCE)" || exit 1
img_hash="$(docker images -q "$img_tag")"

publish_tag="minsky:5000/$img_tag"
echo "*** Tagging image $img_tag (hash $img_hash) for publication as $publish_tag"
docker tag "$img_tag" "$publish_tag" || exit 1

echo "*** Publishing $publish_tag"
docker push "$publish_tag" || exit 1

echo "*** Disabling and stopping all sanmodelexplorer instances on minsky"
ssh minsky sudo systemctl stop "sanmodelexplorer@*" || exit 1
ssh minsky sudo systemctl disable "sanmodelexplorer@*" || exit 1

echo "*** Enabling and starting sanmodelexplorer@r$rev on minsky"
ssh minsky sudo systemctl enable "sanmodelexplorer@r$rev" || exit 1
ssh minsky sudo systemctl start "sanmodelexplorer@r$rev" || exit 1
