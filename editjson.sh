sed -i "s/^\[//;s/\]$//" $1
sed "s/,\"/,\n\"/g" $1 > t1
sed -i "/null,/d" t1
tr "\n" " " < t1 > $1
rm t1
sed -i "s/, /,/g" $1
echo curl -X POST -H "\"Content-type: application/json\"" -d > t2  #need to get rid of \ here
cat $1 >> t2
tr "\n" " " < t2 > $1
rm t2
echo " https://localhost:7020/optimise/general ">> $1
sed -i "s/{/\"{/;s/}/}\"/" $1
sed -i "s/\"/\\\\\"/g" $1
sed -i "s/\\\\\"{/\"{/" $1
sed -i "s/}\\\\\"/}\"/" $1
