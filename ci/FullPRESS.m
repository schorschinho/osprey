% MATLAB runs the test
[result,rt] = UnitTest('PRESS',1,1);

assert(sum(rt{1}.Failed) == 0 && sum(rt{1}.Incomplete) == 0 &&sum(rt{2}.Failed) == 0 && sum(rt{2}.Incomplete) == 0)
