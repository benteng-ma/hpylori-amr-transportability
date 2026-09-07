# Public release checklist

This is a reusable checklist, not a certification that a release has occurred.

- Confirm the author list, actual contributions, release permission and citation metadata.
- Confirm the intended version and set its actual release date in CITATION.cff.
- Run the regression tests and record the operating system, passing tests and explicit skips. The Snippy /tmp safety test requires Linux.
- Check the current presentation update against its frozen source tables and verify the file manifest.
- Exclude credentials, private correspondence, raw third-party sequence files and local runtime caches.
- Keep v1.0.0 and its DOI immutable. Create a new version only after the candidate is finalized.
- After Zenodo actually publishes the new version, verify the record and creators, then add its real version DOI to the citation metadata and manuscript as appropriate. Do not assign the old version DOI to new files.
