
    [result,rt] = UnitTest('SinglePRESS',0,0);

    assert(sum(rt{1}.Failed) == 0 && sum(rt{1}.Incomplete) == 0)
