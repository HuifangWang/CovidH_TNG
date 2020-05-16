#!/bin/bash

csv=$1
grep -v "^#" $csv | cut -f1-18 -d, | tr , " " | column -t > .tmp
nl .tmp | head -n1
nl .tmp | tail -n20



