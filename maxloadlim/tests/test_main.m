addpath('../');
% Run specific tests as well?
run_specific = 0;
%% All tests
testcase = matlab.unittest.TestSuite.fromClass(?Test_maxloadlim);
res = run(testcase);

if run_specific
    %% Specific tests only
    testcase = matlab.unittest.TestSuite.fromMethod(?Test_maxloadlim,'testVarGen_case39');%'testAgainstCPF_case39');%'testAgainstMatpowerCPF_case39');
    res_2 = run(testcase);
    %%
end
