import io
import json
import tarfile
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

from scripts import release_artifacts


class ReleaseArtifactTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.bundle = Path(temp.name)
        self.packages = self.bundle / "packages"
        self.packages.mkdir()
        self.source_sha = "a" * 40
        self.version = "0.2.0"
        self.wheel = self.packages / "pyslicer-0.2.0-py3-none-any.whl"
        self.sdist = self.packages / "pyslicer-0.2.0.tar.gz"
        self.metadata = b"Metadata-Version: 2.4\nName: pyslicer\nVersion: 0.2.0\n"
        self.write_wheel(self.metadata)
        member = tarfile.TarInfo("pyslicer-0.2.0/PKG-INFO")
        member.size = len(self.metadata)
        with tarfile.open(self.sdist, "w:gz") as archive:
            archive.addfile(member, io.BytesIO(self.metadata))
        self.manifest = release_artifacts.describe(self.bundle, self.source_sha, self.version)
        (self.bundle / "manifest.json").write_text(json.dumps(self.manifest), encoding="utf-8")

    def write_wheel(self, metadata):
        with zipfile.ZipFile(self.wheel, "w") as archive:
            archive.writestr("pyslicer-0.2.0.dist-info/METADATA", metadata)

    def test_verify_same_build(self):
        release_artifacts.verify(self.bundle, self.source_sha, self.version)

    def test_reject_manifest_from_another_commit(self):
        with self.assertRaisesRegex(ValueError, "tagged source"):
            release_artifacts.verify(self.bundle, "b" * 40, self.version)

    def test_reject_modified_wheel(self):
        with zipfile.ZipFile(self.wheel, "a") as archive:
            archive.writestr("changed.txt", "changed")
        with self.assertRaisesRegex(ValueError, "package hashes"):
            release_artifacts.verify(self.bundle, self.source_sha, self.version)

    def test_reject_wrong_package_metadata(self):
        self.write_wheel(self.metadata.replace(b"0.2.0", b"0.3.0"))
        with self.assertRaisesRegex(ValueError, "version does not match"):
            release_artifacts.verify(self.bundle, self.source_sha, self.version)

    def test_pypi_status_for_a_new_version_and_matching_existing_files(self):
        with mock.patch.object(release_artifacts, "pypi_files", return_value=None):
            self.assertEqual(release_artifacts.registry_status(self.manifest), "missing")
        with mock.patch.object(release_artifacts, "pypi_files", return_value=self.manifest["files"]):
            self.assertEqual(release_artifacts.registry_status(self.manifest), "published")

    def test_pypi_status_distinguishes_matching_partial_from_different_release(self):
        matching = {self.wheel.name: self.manifest["files"][self.wheel.name]}
        with mock.patch.object(release_artifacts, "pypi_files", return_value=matching):
            self.assertEqual(release_artifacts.registry_status(self.manifest), "partial")
        different = {self.wheel.name: "b" * 64}
        with mock.patch.object(release_artifacts, "pypi_files", return_value=different):
            with self.assertRaisesRegex(ValueError, "outside the selected Build"):
                release_artifacts.registry_status(self.manifest)

    def test_pypi_waits_for_the_complete_release(self):
        partial = {self.wheel.name: self.manifest["files"][self.wheel.name]}
        with (
            mock.patch.object(release_artifacts, "pypi_files", side_effect=[None, partial, self.manifest["files"]]),
            mock.patch.object(release_artifacts.time, "sleep"),
        ):
            self.assertEqual(release_artifacts.registry_status(self.manifest, wait_for_publish=True), "published")

    def test_github_requires_stable_release_with_the_same_files(self):
        files = {**self.manifest["files"], "manifest.json": release_artifacts.sha256(self.bundle / "manifest.json")}
        release = {
            "isDraft": False,
            "isPrerelease": False,
            "tagName": "v0.2.0",
            "assets": [{"name": name, "digest": f"sha256:{digest}"} for name, digest in files.items()],
        }
        release_artifacts.verify_github(self.bundle, release)
        with self.assertRaisesRegex(ValueError, "stable"):
            release_artifacts.verify_github(self.bundle, {**release, "isDraft": True})
        with self.assertRaisesRegex(ValueError, "assets"):
            release_artifacts.verify_github(self.bundle, {**release, "assets": release["assets"][:-1]})
