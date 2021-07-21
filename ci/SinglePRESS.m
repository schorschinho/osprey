function [exit_code] = SinglePRESS
    % assume it always fails
    exit_code = 1;

    % MATLAB runs the test
    [result,rt] = UnitTest('SinglePRESS',0,0);

    % check the result
    if sum(rt{1}.Failed) == 0 && sum(rt{1}.Incomplete) == 0
        exit_code = 0;
    end
end