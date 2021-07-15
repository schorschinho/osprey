matlab -nojvm -nosplash -nodisplay -wait -logfile matlab-output.txt -r pushTest
type matlab-output.txt

set content=
for /f "delims=" %%i in (ExitCode.txt) do set content=%content% %%i

exit %content%
