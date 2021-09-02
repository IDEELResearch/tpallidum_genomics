 for f in *mm.fail.txt; do cat $f| awk '{print $1}' | sort -T /tmp --parallel=24|uniq > ${f}.tmp; done
