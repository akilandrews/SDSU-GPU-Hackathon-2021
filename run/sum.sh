#!/bin/bash
# Summarize results of a run

run() {
for n in *.log; do
  i=`basename $n .log`
  sed -n -e "s/Time: /$i /p;" $n
done
}

run | sort -n -k 1,1
