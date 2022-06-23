#!/usr/bin/bash
cat $1 | sed  "s/^\[//;s/\]$//"| tr "," "\n" | sed "/null$/d" | tr "\n" ","