#!/bin/bash

# Create /state/parametersets which is symlinked to /sanmodelexplorer/parametersets/ (see Dockerfile)
mkdir -p /state/parametersets

(cd /sanmodelexplorer/parametersets.dist; find . -type f -name "*.rd" -print0) | \
while read -d $'\0' -r ps ; do
	if ! test -e "/sanmodelexplorer/parametersets/$ps"; then
		cp "/sanmodelexplorer/parametersets.dist/$ps" "/sanmodelexplorer/parametersets/"
	fi
done

# Copy default dataset specification unless load_by_default.overwrite indicates we shouldn't
if ! test -e "/sanmodelexplorer/parametersets/load_by_default.overwrite" && \
     test -e "/sanmodelexplorer/parametersets/load_by_default.txt"; then
	cp "/sanmodelexplorer/parametersets.dist/load_by_default.txt" "/sanmodelexplorer/parametersets/"
fi

# Run server
exec $*
