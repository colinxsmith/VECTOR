echo off
type %1 | \cygwin\bin\tr "," "\n" | \cygwin\bin\sed "s/{/{\n/"  |   \cygwin\bin\sed "/null$/d" | \cygwin\bin\sed "/doOpt/s/false/true/" | \cygwin\bin\tr "\n" ","|\cygwin\bin\sed "s/,\"/,\n\"/g" | \cygwin\bin\sed "s/{,/{/" > kk
move kk %1