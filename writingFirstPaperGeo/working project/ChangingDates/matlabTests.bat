@echo on
echo Run MATLAB tests
matlab -automation -wait -r "cd ('%cd%'); matlabTests" -logfile "matlab.log"
set EL=%ERRORLEVEL%
type "matlab.log"

IF %EL% NEQ 0 GOTO Error
cd %WORKING_DIRECTORY%
echo Build succeeded
EXIT /B 0

:Error
cd %WORKING_DIRECTORY%
echo Build failed
EXIT /B 1