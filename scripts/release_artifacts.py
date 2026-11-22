"""Record and verify the exact package files promoted by GitHub Actions."""

import argparse
import hashlib
import json
import sys
import tarfile
import time
import tomllib
import zipfile
from email.parser import BytesParser
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import urlopen


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def package_metadata(path):
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as archive:
            names = [name for name in archive.namelist() if name.endswith(".dist-info/METADATA")]
            if len(names) == 1:
                return BytesParser().parsebytes(archive.read(names[0]), headersonly=True)
    else:
        with tarfile.open(path, "r:gz") as archive:
            members = [
                entry for entry in archive.getmembers()
                if entry.isfile() and entry.name.count("/") == 1 and entry.name.endswith("/PKG-INFO")
            ]
            if len(members) == 1:
                with archive.extractfile(members[0]) as stream:
                    return BytesParser().parse(stream, headersonly=True)
    raise ValueError(f"Expected a single package metadata file in {path.name}")


def describe(bundle, source_sha, version):
    packages = sorted((bundle / "packages").iterdir())
    expected = {f"pyslicer-{version}-py3-none-any.whl", f"pyslicer-{version}.tar.gz"}
    if {path.name for path in packages} != expected:
        raise ValueError(f"Expected exactly these distributions: {', '.join(sorted(expected))}")

    for path in packages:
        metadata = package_metadata(path)
        if metadata.get("Name") != "pyslicer" or metadata.get("Version") != version:
            raise ValueError(f"Package name or version does not match in {path.name}")

    return {"source_sha": source_sha, "version": version,
            "files": {path.name: sha256(path) for path in packages}}


def read_manifest(bundle):
    return json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))


def verify(bundle, source_sha, version):
    if read_manifest(bundle) != describe(bundle, source_sha, version):
        raise ValueError("Build manifest does not match the tagged source and package hashes")


def pypi_files(version):
    url = f"https://pypi.org/pypi/pyslicer/{version}/json?check={time.time_ns()}"
    try:
        with urlopen(url, timeout=15) as response:
            package = json.load(response)
    except HTTPError as error:
        if error.code == 404:
            return None
        raise
    if package["info"]["version"] != version:
        raise ValueError("PyPI returned an unexpected package version")
    if any(file.get("yanked") for file in package["urls"]):
        raise ValueError("This PyPI version has yanked files")
    return {file["filename"]: file["digests"]["sha256"] for file in package["urls"]}


def registry_status(manifest, wait_for_publish=False):
    delays = (0, 3, 6, 12, 24, 30) if wait_for_publish else (0,)
    for delay in delays:
        if delay:
            time.sleep(delay)
        files = pypi_files(manifest["version"])
        if files == manifest["files"]:
            return "published"
        if files is not None and not files.items() <= manifest["files"].items():
            raise ValueError("This PyPI version contains files or hashes outside the selected Build")
    if wait_for_publish:
        raise ValueError("PyPI does not have the complete set of verified distributions")
    return "missing" if files is None else "partial"


def verify_github(bundle, release):
    manifest = read_manifest(bundle)
    expected = {**manifest["files"], "manifest.json": sha256(bundle / "manifest.json")}
    actual = {asset["name"]: asset["digest"] for asset in release["assets"]}
    if release["tagName"] != f"v{manifest['version']}" or release["isDraft"] or release["isPrerelease"]:
        raise ValueError("GitHub release is not a stable published release for this version")
    if actual != {name: f"sha256:{digest}" for name, digest in expected.items()}:
        raise ValueError("GitHub release assets do not match the selected Build")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for command in ("record", "verify", "pypi", "github"):
        subparser = commands.add_parser(command)
        subparser.add_argument("bundle", type=Path)
        if command in ("record", "verify"):
            subparser.add_argument("source_sha")
        if command == "verify":
            subparser.add_argument("version")
        if command == "pypi":
            subparser.add_argument("--wait", action="store_true")
    args = parser.parse_args()
    try:
        if args.command == "record":
            version = tomllib.loads(Path("pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
            manifest = describe(args.bundle, args.source_sha, version)
            (args.bundle / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        elif args.command == "verify":
            verify(args.bundle, args.source_sha, args.version)
        elif args.command == "pypi":
            print(registry_status(read_manifest(args.bundle), args.wait))
        else:
            verify_github(args.bundle, json.load(sys.stdin))
    except (OSError, KeyError, ValueError, tarfile.TarError, zipfile.BadZipFile) as error:
        parser.exit(1, f"Release verification failed: {error}\n")


if __name__ == "__main__":
    main()
