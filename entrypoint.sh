#!/bin/bash

# Create /state/parametersets which is symlinked to /sanmodelexplorer/parametersets/ (see Dockerfile)
mkdir -p /state/parametersets

(cd /sanmodelexplorer/parametersets.dist; find . -type f -name "*.rd" -print0) | \
while read -d $'\0' -r ps ; do
	if ! test -e "/sanmodelexplorer/parametersets/$ps"; then
		cp "/sanmodelexplorer/parametersets.dist/$ps" "/sanmodelexplorer/parametersets/"
	fi
done

# Run server
exec $*