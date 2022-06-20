#!/usr/bin/bash
curl -X get "https://localhost:7020/optimise/general?doOpt=false&datafile=$3" --output $1
#echo $pp
sed -i "s/^\[//;s/\]$//" $1
curl -X POST -H "Content-type: application/json" -d @$1 https://localhost:7020/optimise/$2