SynGrid Release Checklist
=========================


Pre-release
-----------
- Check [SynGrid issue tracker](https://github.com/MATPOWER/mx-syngrid/issues)
  and `untracked/SynGrid-To-Do-List.md` for to do items.
- Create & checkout new `prep-for-release` branch from latest `master`.
- Release notes:
  - Make sure Release History in Appendix A of `SynGrid-manual.tex` is
    up-to-date.
  - Create `docs/relnotes/SynGrid-Release-Notes-#.#.md` document from
    Appendix A of `SynGrid-manual.tex`.
- Update date in Copyright line in:
  - `LICENSE`
- Update version number and date in:
  - `sgver.m`
  - `docs/relnotes/SynGrid-Release-Notes-#.#.md`
  - `docs/src/SynGrid-manual/SynGrid-manual.tex`
    - title page
    - copyright (front page and LICENSE text)
    - Appendix A Release History
    - `\syngridver` (update `\mpver`, `\mipsver`, and `\mostver` too)
- In `README.md` and `docs/src/SynGrid-manual/SynGrid-manual.tex`
  - update output of `test_syngrid` in Section 2.2
- In `docs/src/SynGrid-manual/SynGrid-manual.tex`
  - check for any highlighting `\hl`
- Create new DOI for this version of the User's Manual
  - Go to https://doi.org/10.5281/zenodo.3238679
    - Click "New Version" to reserve new DOI for new version
  - Make updates for current version specific citations:
    - version number (3 places)
    - latest version DOI, current is: 10.5281/zenodo.3251099  10.5281/zenodo.11194332
      - (update here each time)
    ... in the following places ...
    - CITATION file
    - Citing ... section of README.md
    - Citing ... section of User's Manual
    - Citing ... section of website (not currently here)
  - Make updates for non-version specific citations:
    - search everywhere for 10.5281/zenodo.3238679 and update year
      - User's Manual
      - search citations in all other projects being updated simultaneously
        (does not exist as of 6/20/19)
    - search everywhere for 10.5281/zenodo.3236535 and update year (MATPOWER Software)
      - User's Manual
    - search everywhere for 10.5281/zenodo.3236519 and update year (MATPOWER User's Manual)
      - User's Manual
- Create `SynGrid-manual.pdf` from `SynGrid-manual.tex` and  move to `docs`.
- Add release notice with date and version in `CHANGES.md`.
- Commit all changes to `prep-for-release`.
- Push `prep-for-release` to GitHub.
- Make sure CI checks are ok.
- Make copy of `docs/SynGrid-manual.pdf` named `SynGrid-manual-x.x.pdf`
    - copy to `docs` directory of `matpower.org-static` git repo
      - update `SynGrid-manual.pdf` symlink on `https://matpower.org/docs/` to point
        to new `SynGrid-manual-x.x.pdf` (replaces existing current version)
        - `cd dev/projects/matpower.org-static/docs`
        - `rm SynGrid-manual.pdf`
        - `ln -s ./SynGrid-manual-x.x.pdf SynGrid-manual.pdf`
      - commit & push, then pull to matpower.org
    - upload `SynGrid-manual-x.x.pdf` to Zenodo and finish entry for "New Version"
      - update:
        - Publication date
        - Version
        - Identifiers:
          - version number in "identical to"
  - add link on `https://matpower.org/doc/manuals/` page


Release
-------
- Merge latest `prep-for-release` into `master`.
- Tag with version number, e.g. `1.3`.
- Push `master` to GitHub.
- Publish new release on GitHub: https://github.com/MATPOWER/mx-syngrid/releases/new
  - use (possibly shortened) contents of `docs/relnotes/SynGrid-Release-Notes-#.#.md`
- Update matpower.org website
  - update in matpower.org-static Git repo, then pull to matpower.org
    - in docs, update symlink for current manual to point to new versions
        - SynGrid-manual.pdf


Post-release
------------
- Merge latest `master` into `release`.
- Push `release` to GitHub.
- In manual
  - update version to dev version
  - comment out date line so it uses current date
