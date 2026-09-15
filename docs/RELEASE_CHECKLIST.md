# Release checklist

Before publishing the empirical ECM extension as `v1.3.0`:

- [ ] Run the full unit-test suite from a clean checkout.
- [ ] Verify the checksums in `data/empirical_ecm/SOURCE_NOTICE.md`.
- [ ] Confirm that a full ECM reproduction matches all committed result files.
- [ ] Review `git diff` and ensure no exploratory files or credentials are staged.
- [ ] Create a GitHub release tagged `v1.3.0`.
- [ ] Optionally archive the release on Zenodo to obtain a DOI.
- [ ] If a manuscript/preprint is posted, add its citation near the top of the README.
- [ ] Keep the claim boundary in `docs/CLAIMS_AND_LIMITATIONS.md` intact unless new evidence supports a change.
