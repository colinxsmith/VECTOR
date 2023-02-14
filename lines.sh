#!/usr/bin/bash
cat $1 | tr "," "\n" | sed "s/{/{\n/" | sed "s/null}/null\n}/"  |   sed "/null$/d" | tr "\n" ","|sed "s/,\"/,\n\"/g" | sed "s/{,/{/" | sed "s/,}/}/" > kk
mv kk $1
