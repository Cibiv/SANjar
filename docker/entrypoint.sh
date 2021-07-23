#!/bin/bash

# Create /state/parametersets which is symlinked to /sanmodelexplorer/parametersets/ (see Dockerfile)
mkdir -p /state/parametersets

(cd /sanmodelexplorer/parametersets.dist; find . -type f -name "*.rd" -print0) | \
while read -d $'\0' -r ps ; do
	if ! test -e "/sanmodelexplorer/parametersets/$ps"; then
		cp "/sanmodelexplorer/parametersets.dist/$ps" "/sanmodelexplorer/parametersets/"
	fi
done

# Copy default dataset specification unless load_by_default.override indicates we shouldn't
if ! test -e "/sanmodelexplorer/parametersets/load_by_default.override"; then
	cp "/sanmodelexplorer/parametersets.dist/load_by_default.txt" "/sanmodelexplorer/parametersets/"
fi

# Run serve
# NOTE: We use script here to run the server under a tty. This is a dirty hack that ensures
# that progress message are output by progress_bar, which we then intercept and turn into
# shiny progress updates.
cmd=""
for arg in "$@"; do
	cmd="$cmd '$(printf "%s\n" "$arg" | sed "s/'/'\\\\''/g")'"
done
exec script -qefc "$cmd" /dev/null
