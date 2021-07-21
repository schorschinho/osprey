% assume it always fails
exit_code = 1;

% MATLAB runs the test
[result,rt] = UnitTest('SinglePRESS',2,1);

% check the result
if sum(rt{1}.Failed) == 0 && sum(rt{1}.Incomplete) == 0 &&sum(rt{2}.Failed) == 0 && sum(rt{2}.Incomplete) == 0
    exit_code = 0;
end
