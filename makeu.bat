@echo off
set res=%4
if "%4"=="" set res=result
curl -X get "http://localhost:7777/optimise/general?doOpt=false&datafile=%3&round=0" --output %1
\cygwin\bin\sed -i "s/^\[//;s/\]$//" %1
type %1 | \cygwin\bin\tr "," "\n" | \cygwin\bin\sed "/null$/d" | \cygwin\bin\sed "/doOpt/s/false/true/" | \cygwin\bin\tr "\n" "," > kk
move kk %1

curl -X POST -H "Content-type: application/json" -d @%1 http://localhost:7777/optimise/%2 --output %res%

