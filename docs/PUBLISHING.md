# Publishing PySlicer

**Build** tests PySlicer and produces a wheel, source distribution, and SHA-256 manifest. **Publish** promotes the artifact from the successful `master` push build for the tagged commit. PyPI and GitHub Releases receive the same distribution files.

## One-time setup

1. In GitHub repository settings, create the `pypi` environment and allow deployments only from the `master` branch.
2. In your PyPI account, add a [pending trusted publisher](https://docs.pypi.org/trusted-publishers/creating-a-project-through-oidc/) for the first release. For an existing PyPI project, use its **Publishing** settings instead.

   | Setting | Value |
   | --- | --- |
   | PyPI project | `pyslicer` |
   | Owner | `isque03` |
   | Repository | `pyslicer` |
   | Workflow | `publish.yml` |
   | Environment | `pypi` |

The publish job uses GitHub's OIDC identity. You do not need a PyPI API token in GitHub.

## Build and publish

Build runs on pull requests and pushes to `master`. It tests Python 3.10 and 3.14 on Linux and macOS, checks package metadata, and installs both built distributions in fresh environments on both operating systems. The installed CLI smoke test loads a YAML profile and generates G-code and an HTML preview from a cube. The **Build passed** check requires every job in this workflow to succeed before a pull request can merge. CI retains the `pyslicer-distributions` artifact for 30 days.

For a stable release, set `[project].version` in `pyproject.toml`, merge the release changes, and identify the successful **push** Build for that exact commit. The following examples use `0.2.0` and a clean checkout of `master`. The package verification command requires Python 3.11 or later.

```bash
git pull --ff-only
VERSION=0.2.0
SHA=$(git rev-parse HEAD)
TAG="v$VERSION"

gh run list --workflow build.yml --branch master --event push --commit "$SHA" --status success
# Set BUILD_RUN to the successful run ID shown above.
CANDIDATE=$(mktemp -d)
gh run download "$BUILD_RUN" --name pyslicer-distributions --dir "$CANDIDATE"
python3 scripts/release_artifacts.py verify "$CANDIDATE" "$SHA" "$VERSION"

git tag -a "$TAG" "$SHA" -m "PySlicer $VERSION"
git push origin "refs/tags/$TAG"
gh workflow run publish.yml --ref master -f tag="$TAG"
```

Publish requires an annotated `vMAJOR.MINOR.PATCH` tag that matches `pyproject.toml` and points into `master` history. It selects the exact successful push Build and verifies both distributions against the manifest and their embedded metadata. The PyPI job uploads only those files. After checking PyPI's published file hashes, the workflow creates a stable GitHub Release with the wheel, sdist, and manifest and verifies their asset hashes.

## Verify distribution

Watch the Publish run complete, then check both public endpoints. `$CANDIDATE`, `$VERSION`, and `$TAG` refer to the verified artifact and tag from above.

```bash
gh run list --workflow publish.yml --branch master --event workflow_dispatch
# Set PUBLISH_RUN to the new Publish run ID.
gh run watch "$PUBLISH_RUN" --exit-status
python3 scripts/release_artifacts.py pypi "$CANDIDATE" --wait
gh release view "$TAG" --json tagName,isDraft,isPrerelease,assets,url \
  | python3 scripts/release_artifacts.py github "$CANDIDATE"
```

To check the installed public package, create a fresh virtual environment, run `python -m pip install "pyslicer==$VERSION"`, and try `pyslicer --help`. Build tests verify installed software and generated files. Physical print acceptance depends on the printer and its profile.

## Recovery

- **Build failed or the artifact expired before publication.** Fix the failure before tagging, or rerun the Build for the same tagged SHA if only the artifact expired. Download and verify the new artifact. Keep release tags immutable.
- **PyPI has both matching distribution files and GitHub release creation failed.** Dispatch Publish again with the original Build artifact available. The workflow will compare PyPI's hashes and skip the upload. If GitHub retained a draft, use the exact original files to finish uploading its assets, publish the draft as stable, and rerun verification.
- **PyPI accepted only one distribution.** Dispatch Publish again with the original Build artifact available. Publish verifies the existing file against that artifact before resuming the upload.
- **PyPI has different hashes or extra files for the version.** Investigate the registry state and use a new version and tag for a corrected release. Keep published tags and files immutable.
- **PyPI rejects the publisher identity.** Check the project, owner, repository, workflow filename, branch, and environment in the GitHub run and PyPI settings. A pending publisher becomes a trusted publisher after the first accepted upload.
