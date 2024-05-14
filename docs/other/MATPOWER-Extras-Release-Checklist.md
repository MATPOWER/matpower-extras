MATPOWER Extras Release Checklist
=================================

Pre-release
-----------

- Check issue trackers for to do items for:
  - [`mx-maxloadlim`](https://github.com/CamilleH/Max-Load-Lim-matpower/issues)
  - [`mx-reduction`](https://github.com/MATPOWER/mx-reduction/issues)
  - [`mx-sdp_pf`](https://github.com/MATPOWER/mx-sdp_pf/issues)
  - [`mx-se`](https://github.com/MATPOWER/mx-se/issues)
  - [`mx-simulink_matpower`](https://github.com/MATPOWER/mx-simulink_matpower/issues)
  - [`mx-syngrid`](https://github.com/MATPOWER/mx-syngrid/issues)
  - [`matpower-extras`](https://github.com/MATPOWER/matpower-extras/issues)
- For each of the following, (1) update `master` branch, (2) push to GitHub,
  (3) merge into `release` branch (4) tag with version number, (5) push to
  GitHub, again:
  - `mx-maxloadlim` - submit upstream pull requests, if necessary
  - `mx-reduction`
  - `mx-sdp_pf` - follow `mx-syngrid/untracked/SDP_PF-Release-Checklist.md`
  - `mx-se`
  - `mx-simulink_matpower`
  - `mx-syngrid` - follow `mx-syngrid/untracked/SynGrid-Release-Checklist.md`
- Create & checkout new `prep-for-release` branch from latest `master`.
- Get released subrepos (use `--branch=release` if different from `master`):
  - `git subrepo pull --branch=master maxloadlim`
  - `git subrepo pull --branch=master reduction`
  - `git subrepo pull --branch=master sdp_pf`
  - `git subrepo pull --branch=master se`
  - `git subrepo pull --branch=main simulink_matpower`
  - `git subrepo pull --branch=master syngrid`
- Add release notice with date and version in `smartmarket/SM_CHANGES.md`.
- Commit all changes to `prep-for-release`.
- Push `prep-for-release` to GitHub.
- Make sure CI checks are ok.


Release
-------
- Merge latest `prep-for-release` into `master`.
- Tag with version number, e.g. `8.0`.
- Push `master` to GitHub.
- Publish new release on GitHub: https://github.com/MATPOWER/matpower-extras/releases/new


Post-release
------------
- Merge latest `master` into `release`.
- Push `release` to GitHub.
