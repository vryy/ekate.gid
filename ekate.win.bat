REM    OutputFile: $2/$1.info
REM    ErrorFile: $2/$1.err
REM    delete previous result file
DEL "%2\*.post.res"
DEL "%2\*.post.msh"
DEL "%2\*.post.bin"
REM COPY "%3\sed.exe" "%2"
REM COPY "%3\regex2.dll" "%2"
REM COPY "%3\libintl3.dll" "%2"
REM COPY "%3\libiconv2.dll" "%2"
ECHO "Please install Python-2 and add to the path"
PYTHON "%3\clean_mdpa.py" "%2\%1.dat" "%2\%1.bak"
DEL "%2\%1.dat"
REN "%2\%1.bak" "%2\%1.dat"
REM renaming Kratos input files
REN "%2\%1.dat" "%2\%1.mdpa"
REN "%2\%1-1.dat" "%2\%1.inp"
REN "%2\%1-2.dat" "%2\%1.py"
REN "%2\%1-3.dat" "%2\%1_include.py"
REN "%2\%1-4.dat" "%2\%1_layers.py"
REN "%2\%1-5.dat" "%2\%1_parallel_include.py"
ECHO >> "%2\%1.ess"
TYPE "%2\%1.ess" >> "%2\%1.py"
REM "%2\SED.EXE" "s\rEpLaCeMeNtStRiNg\%1\g" < "%2\%1.py" > "%2\%1.py_changed"
"%3\SED.EXE" -e "s\rEpLaCeMeNtStRiNg\%1\g" -e "s/__DaTe__/%(date)/g" -e "s/__TiMe__/%(time)/g" < "%2\%1.py" > "%2\%1.py_changed"
REN "%2\%1.py_changed" "%2\%1.py"
CD "%2"
CD ..

