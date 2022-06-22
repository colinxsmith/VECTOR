@echo off
type %1 | \cygwin\bin\sed  "s/^\[//;s/\]$//"| \cygwin\bin\tr "," "\n" | \cygwin\bin\sed "/null$/d" | \cygwin\bin\tr "\n" ","