SDP_PF_ Release Checklist
=========================


Pre-release
-----------
- Check [SDP_PF issue tracker](https://github.com/MATPOWER/mx-sdp_pf/issues)
  for to do items.
- Create & checkout new `prep-for-release` branch from latest `master`.
- Update version number and date in:
  - `sdp_pf_ver.m`
  - Copyright line in `LICENSE`.
- Add release notice with date and version in `CHANGES.md`.
- Commit all changes to `prep-for-release`.
- Push `prep-for-release` to GitHub.
- Make sure CI checks are ok.


Release
-------
- Merge latest `prep-for-release` into `master`.
- Tag with version number, e.g. `1.0.1`. (actually, use MATPOWER version for tag)
- Push `master` to GitHub.


Post-release
------------
- Merge latest `master` into `release`.
- Push `release` to GitHub.
