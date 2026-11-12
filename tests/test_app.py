import functools
import hashlib
import http.server
import json
import os
import threading

import pytest

from conftest import write_sources
from nmrcheck import cli, sources
from nmrcheck.sources import SOURCES, Source, fetch, source_path
from nmrcheck.web import create_app


@pytest.fixture
def served(tmp_path):
    """Serve a folder over HTTP so downloads can be exercised offline."""
    root = tmp_path / "remote"
    root.mkdir()
    handler = functools.partial(http.server.SimpleHTTPRequestHandler,
                                directory=str(root))
    handler.log_message = lambda *a, **k: None
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    yield root, f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()


def fake_source(name, url, md5=""):
    old = SOURCES[name]
    return Source(name=name, filename=old.filename, url=url, md5=md5,
                  license=old.license, citation=old.citation,
                  homepage=old.homepage)


def test_download_and_checksum(served, tmp_path, monkeypatch):
    root, base = served
    payload = b"x" * 300000
    (root / "file.bin").write_bytes(payload)
    monkeypatch.setitem(SOURCES, "nmrexp", fake_source(
        "nmrexp", f"{base}/file.bin", hashlib.md5(payload).hexdigest()))
    notes = []
    path = fetch("nmrexp", tmp_path / "data", log=lambda *a: notes.append(a))
    assert path.read_bytes() == payload
    assert notes == []
    assert not path.with_name(path.name + ".part").exists()


def test_checksum_mismatch_is_logged(served, tmp_path, monkeypatch):
    root, base = served
    (root / "file.bin").write_bytes(b"abc")
    monkeypatch.setitem(SOURCES, "nmrexp",
                        fake_source("nmrexp", f"{base}/file.bin", "0" * 32))
    notes = []
    fetch("nmrexp", tmp_path / "data", log=lambda *a: notes.append(a))
    assert notes[0][0] == "checksum_mismatch"


def test_checksum_cache_is_dropped_when_file_changes(tmp_path, monkeypatch):
    good = b"published content"
    monkeypatch.setitem(SOURCES, "nmrexp", fake_source(
        "nmrexp", "", hashlib.md5(good).hexdigest()))
    path = source_path(tmp_path, "nmrexp")
    path.parent.mkdir(parents=True)
    path.write_bytes(good)
    stamp = path.with_name(path.name + ".md5")

    notes = []
    fetch("nmrexp", tmp_path, offline=True, log=lambda *a: notes.append(a))
    assert notes == []
    assert stamp.read_text() == hashlib.md5(good).hexdigest()

    # while the cache is newer than the file it is trusted and nothing is hashed
    calls = []
    real_md5 = sources.file_md5
    monkeypatch.setattr(sources, "file_md5",
                        lambda p, *a: calls.append(p) or real_md5(p, *a))
    fetch("nmrexp", tmp_path, offline=True, log=lambda *a: notes.append(a))
    assert calls == [] and notes == []

    # the file is replaced by a different download: the cache must not hide it
    bad = b"truncated or tampered"
    path.write_bytes(bad)
    later = stamp.stat().st_mtime + 60
    os.utime(path, (later, later))
    fetch("nmrexp", tmp_path, offline=True, log=lambda *a: notes.append(a))
    assert calls == [path]
    assert [n[0] for n in notes] == ["checksum_mismatch"]
    assert stamp.read_text() == hashlib.md5(bad).hexdigest()


def test_failed_download_is_logged(served, tmp_path, monkeypatch):
    _, base = served
    monkeypatch.setitem(SOURCES, "nmrexp",
                        fake_source("nmrexp", f"{base}/missing.bin"))
    monkeypatch.setattr("nmrcheck.sources.time.sleep", lambda s: None)
    notes = []
    assert fetch("nmrexp", tmp_path / "data",
                 log=lambda *a: notes.append(a)) is None
    assert notes[0][0] == "source_unavailable"


def test_offline_and_manual_sources(tmp_path):
    notes = []
    assert fetch("nmrexp", tmp_path, offline=True,
                 log=lambda *a: notes.append(a)) is None
    assert fetch("nmrbank", tmp_path, log=lambda *a: notes.append(a)) is None
    assert [n[0] for n in notes] == ["source_unavailable"] * 2


def test_cli_end_to_end(served, tmp_path, monkeypatch, capsys):
    """fetch -> build -> check with nothing prepared beforehand."""
    root, base = served
    files = write_sources(root)
    for name, path in files.items():
        flat = root / path.name
        flat.write_bytes(path.read_bytes())
        monkeypatch.setitem(SOURCES, name, fake_source(name, f"{base}/{path.name}"))
    data = tmp_path / "data"
    out = tmp_path / "reports"
    code = cli.main(["check", "c1ccccc1", "--data-dir", str(data),
                     "-o", str(out), "-s", "CDCl3"])
    assert code == 0
    assert source_path(data, "nmrexp").exists()
    assert (data / "index.sqlite").exists()
    assert "reported" in capsys.readouterr().out
    report = json.loads((out / "c1ccccc1.json").read_text())
    assert report["verdict"]["status"] == "reported"
    assert (out / "c1ccccc1.md").read_text().startswith("# NMR novelty check")

    # second run reuses the index and needs no network
    monkeypatch.setitem(SOURCES, "nmrexp", fake_source("nmrexp", "http://127.0.0.1:9/x"))
    assert cli.main(["check", "CCO", "--data-dir", str(data), "-o", str(out),
                     "--offline"]) == 0


def test_cli_invalid_structure(data_dir, tmp_path):
    assert cli.main(["check", "???", "--data-dir", str(data_dir),
                     "-o", str(tmp_path)]) == 2


def test_cli_from_dir(tmp_path):
    files = write_sources(tmp_path / "staging")
    flat = tmp_path / "downloads"
    flat.mkdir()
    for path in files.values():
        (flat / path.name).write_bytes(path.read_bytes())
    data = tmp_path / "data"
    assert cli.main(["build", "--data-dir", str(data), "--from-dir", str(flat),
                     "--offline"]) == 0
    assert (data / "index.sqlite").exists()


def test_web_form_and_api(index):
    client = create_app(index).test_client()
    assert b"<form" in client.get("/").data
    page = client.get("/?target=c1ccccc1&impurities=Cc1ccccc1,CCO&solvent=CDCl3")
    assert page.status_code == 200
    assert b"<strong>reported</strong>" in page.data
    assert b"Other species" in page.data
    api = client.get("/api/check?target=c1ccccc1")
    assert api.get_json()["verdict"]["status"] == "reported"
    assert client.get("/api/check").status_code == 400
    assert client.get("/api/check?target=C>C").status_code == 400
    assert b"class=\"error\"" in client.get("/?target=C>C").data
