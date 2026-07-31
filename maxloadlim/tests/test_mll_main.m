function success = test_mll_main(verbose, exit_on_fail)
%TEST_MLL_MAIN Run all MaxLoadLim tests.
%   TEST_MLL_MAIN
%   TEST_MLL_MAIN(VERBOSE)
%   TEST_MLL_MAIN(VERBOSE, EXIT_ON_FAIL)
%   SUCCESS = TEST_MLL_MAIN(...)
%
%   Runs all of the MaxLoadLim tests. If VERBOSE is true (false by default),
%   it prints the details of the individual tests. If EXIT_ON_FAIL is true
%   (false by default), it will exit MATLAB or Octave with a status of 1
%   unless T_RUN_TESTS returns ALL_OK.
%
%   See also T_RUN_TESTS.

%   MATPOWER
%   Copyright (c) 2015-2019, Power Systems Engineering Research Center (PSERC)
%   by Camille Hamon
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER/mx-maxloadlim.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mx-maxloadlim/ for more info.

if nargin < 2
    exit_on_fail = 0;
    if nargin < 1
        verbose = 0;
    end
end

tests = {};

tests{end+1} = 't_cpf_case9';
tests{end+1} = 't_cpf_case39';
tests{end+1} = 't_varGen_case9';
tests{end+1} = 't_varGen_case39';

%% warning for bug in Octave 11.x
is_octave11 = have_feature('octave') && floor(have_feature('octave', 'vnum')) == 11;
if is_octave11
    persistent show_octave11_warning_once;
    if isempty(show_octave11_warning_once)
        show_octave11_warning_once = 1;
        warning(sprintf('\n###############################################################################\n#  GNU Octave 11.x has a bug (https://savannah.gnu.org/bugs/index.php?68227)  #\n#  that results in lots of warnings when running MATPOWER. One workaround is  #\n#  to turn off all warnings using:  warning off                               #\n#  (applied automatically and temporarily for test_mll_main)                  #\n###############################################################################\n'));
    end
    w = warning();
    warning('off');
end

%% run the tests
all_ok = t_run_tests( tests, verbose );

%% handle success/failure
if exit_on_fail && ~all_ok
    exit(1);
end
if nargout
    success = all_ok;
end
