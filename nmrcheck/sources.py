"""Public data sources and the code that fetches them."""
import hashlib
import shutil
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import requests


@dataclass(frozen=True)
class Source:
    name: str
    filename: str
    url: str
    md5: str
    license: str
    citation: str
    homepage: str
    required: bool = True


SOURCES = {
    "nmrshiftdb2": Source(
        name="nmrshiftdb2",
        filename="nmrshiftdb2withsignals.sd",
        url="https://sourceforge.net/projects/nmrshiftdb2/files/data/"
            "nmrshiftdb2withsignals.sd/download",
        md5="",
        license="nmrshiftdb2 Data License (open content)",
        citation="Kuhn S., Schloerer N.E., Magn. Reson. Chem. 2015, 53, 582",
        homepage="https://nmrshiftdb.nmr.uni-koeln.de",
    ),
    "nmrexp": Source(
        name="nmrexp",
        filename="NMRexp_10to24_1_1004.parquet",
        url="https://zenodo.org/records/17296666/files/"
            "NMRexp_10to24_1_1004.parquet?download=1",
        md5="2f9ed8bc533364dfb9dce96c83703dd4",
        license="CC BY 4.0",
        citation="Wang J.-J. et al., Sci. Data 2025; Zenodo 10.5281/zenodo.17296666",
        homepage="https://zenodo.org/records/17296666",
    ),
    # NMRBank has no stable direct download link and is CC BY-NC, so it is
    # never fetched automatically. Drop the CSV into the data directory to
    # have it indexed.
    "nmrbank": Source(
        name="nmrbank",
        filename="NMRBank.csv",
        url="",
        md5="",
        license="CC BY-NC 3.0",
        citation="Wang Q. et al., Chem. Sci. 2025, 16, 11548",
        homepage="https://doi.org/10.1039/D4SC08802F",
        required=False,
    ),
}


def source_path(data_dir, name):
    return Path(data_dir) / name / SOURCES[name].filename


def file_md5(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


class FetchError(RuntimeError):
    pass


def _download(url, dest, retries=3, timeout=60):
    part = dest.with_name(dest.name + ".part")
    last = None
    for attempt in range(retries):
        offset = part.stat().st_size if part.exists() else 0
        headers = {"Range": f"bytes={offset}-"} if offset else {}
        try:
            with requests.get(url, stream=True, headers=headers, timeout=timeout,
                              allow_redirects=True) as resp:
                if resp.status_code == 416:
                    break
                resp.raise_for_status()
                if offset and resp.status_code != 206:
                    offset = 0
                total = int(resp.headers.get("Content-Length", 0)) + offset
                done = offset
                tick = time.time()
                with open(part, "ab" if offset else "wb") as out:
                    for block in resp.iter_content(1 << 20):
                        out.write(block)
                        done += len(block)
                        if time.time() - tick > 5:
                            tick = time.time()
                            pct = f" ({100 * done / total:.0f}%)" if total else ""
                            print(f"  {dest.name}: {done >> 20} MB{pct}",
                                  file=sys.stderr)
            break
        except requests.RequestException as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    else:
        raise FetchError(f"download failed for {url}: {last}")
    part.replace(dest)


def fetch(name, data_dir, from_dir=None, offline=False, log=None):
    """Make sure the raw file for a source exists. Returns its path or None."""
    src = SOURCES[name]
    dest = source_path(data_dir, name)
    note = log or (lambda *a, **k: None)

    if not dest.exists() and from_dir:
        local = Path(from_dir) / src.filename
        if local.exists():
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(local, dest)

    if not dest.exists():
        if not src.url:
            note("source_unavailable", name,
                 f"{src.filename} not present and the source has no automatic "
                 "download; skipped")
            return None
        if offline:
            note("source_unavailable", name, "offline mode and file not present")
            return None
        dest.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading {name} from {src.url}", file=sys.stderr)
        try:
            _download(src.url, dest)
        except FetchError as exc:
            note("source_unavailable", name, str(exc))
            return None

    if src.md5:
        stamp = dest.with_name(dest.name + ".md5")
        digest = stamp.read_text().strip() if stamp.exists() else ""
        if not digest:
            digest = file_md5(dest)
            stamp.write_text(digest)
        if digest != src.md5:
            note("checksum_mismatch", name,
                 f"expected md5 {src.md5}, got {digest}; file used as is")
    return dest


def fetch_all(data_dir, from_dir=None, offline=False, log=None, only=None):
    found = {}
    for name in SOURCES:
        if only and name not in only:
            continue
        path = fetch(name, data_dir, from_dir=from_dir, offline=offline, log=log)
        if path is not None:
            found[name] = path
    return found
