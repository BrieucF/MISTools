#! /bin/bash

scripts=`find $1/ -name "*run_limits.sh"`

for script in $scripts; do
    dir=$(dirname $script)
    script=$(basename $script)
    echo "Computing limits with ${script}"
    pushd $dir &> /dev/null
    . $script
    popd &> /dev/null
done

