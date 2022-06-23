#!/usr/bin/bash
res=${4:-result}
echo $res
curl -X get "https://localhost:7020/optimise/general?doOpt=false&datafile=$3&round=0" --output $1
sed  "s/^\[//;s/\]$//" $1  | tr "," "\n" | sed "/null$/d" | sed "/doOpt/s/false/true/" | tr "\n" "," > kk
mv kk $1

curl -X POST -H "Content-type: application/json" -d @$1 https://localhost:7020/optimise/$2 --output $res

