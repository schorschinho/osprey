% assume it always fails
exit_code = 1;

% MATLAB runs the test
[result,rt] = unitTest(0,'PRESS',0)

% check the result
if sum(result{1}.Failed) == 0 && sum(result{1}.Incomplete) == 0
    exit_code = 0;
end

% write the ExitCode
fid = fopen('ExitCode.txt','w');
fprintf(fid,'%d',exit_code);
fclose(fid);

% Ensure that we ALWAYS call exit that is always a success so that CI doesn't stop
exit
