sed -i "s/^\[//;s/\]$//" $1
tr "," "\n" < $1 |sed -i "/null$/d" |tr "\n" ","  > kk
mv kk $1
echo curl -X POST -H "\"Content-type: application/json\"" -d > t2  #need to get rid of \ here
cat $1 >> t2
tr "\n" " " < t2 > $1
rm t2
echo " https://localhost:7020/optimise/general ">> $1
sed -i "s/{/\"{/;s/}/}\"/" $1
sed -i "s/\"/\\\\\"/g" $1
sed -i "s/\\\\\"{/\"{/" $1
sed -i "s/}\\\\\"/}\"/" $1
