function success = test_sdp_pf(verbose, exit_on_fail)
%TEST_SDP_PF  Run all SDP_PF tests.
%   TEST_SDP_PF
%   TEST_SDP_PF(VERBOSE)
%   TEST_SDP_PF(VERBOSE, EXIT_ON_FAIL)
%   SUCCESS = TEST_SDP_PF(...)
%
%   Runs all of the SDP_PF tests. If VERBOSE is true (false by default),
%   it prints the details of the individual tests. If EXIT_ON_FAIL is true
%   (false by default), it will exit MATLAB or Octave with a status of 1
%   unless T_RUN_TESTS returns ALL_OK.
%
%   See also T_RUN_TESTS.

%   MATPOWER
%   Copyright (c) 2004-2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER/mx-sdp_pf.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mx-sdp_pf/ for more info.

if nargin < 2
    exit_on_fail = 0;
    if nargin < 1
        verbose = 0;
    end
end

tests = {};

if have_feature('mosek') || have_feature('sdpt3') || have_feature('sedumi')
    tests{end+1} = 't_opf_sdpopf';
    tests{end+1} = 't_insolvablepf';
    tests{end+1} = 't_insolvablepf_limitQ';
    tests{end+1} = 't_insolvablepfsos';
    tests{end+1} = 't_insolvablepfsos_limitQ';
end
tests{end+1} = 't_testglobalopt';

%% warning for bug in Octave 11.x
is_octave11 = have_feature('octave') && floor(have_feature('octave', 'vnum')) == 11;
if is_octave11
    persistent show_octave11_warning_once;
    if isempty(show_octave11_warning_once)
        show_octave11_warning_once = 1;
        warning(sprintf('\n###############################################################################\n#  GNU Octave 11.x has a bug (https://savannah.gnu.org/bugs/index.php?68227)  #\n#  that results in lots of warnings when running MATPOWER. One workaround is  #\n#  to turn off all warnings using:  warning off                               #\n#  (applied automatically and temporarily for test_sdp_pf)                    #\n###############################################################################\n'));
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
