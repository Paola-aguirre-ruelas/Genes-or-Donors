#!/usr/bin/env python3
"""Stream the public H5AD once and verify its size and SHA-256 without writing it to disk.

Sends exactly one HTTP GET (standard-library ``http.client``; redirects are not followed and a
failed transfer is not retried), reads the response body into one fixed 1 MiB buffer, updates a
SHA-256 digest and a byte counter for every chunk and then discards the chunk. The file is never
written to disk and never parsed. Only a small JSON verification record is written (``--record``).

    python -I scripts/verify_public_h5ad.py [--record PATH]

Exit status 0 means PUBLIC_SOURCE_SHA256_VERIFIED; 1 means PUBLIC_SOURCE_MISMATCH or
TRANSFER_INCOMPLETE (details in the record).
"""
import argparse
import hashlib
import http.client
import json
import platform
import sys
import time
from datetime import datetime, timezone
from urllib.parse import urlsplit

URL = "https://exampledata.scverse.org/pertpy/stephenson_2021_subsampled.h5ad"
EXPECTED_BYTES = 700_549_709
EXPECTED_SHA256 = "2345f3c512552892bb8327f4779dc69972faccac3b00dcf90c0e9d9ae5ab738e"
CHUNK_BYTES = 1 << 20
TIMEOUT_S = 120
KEEP_HEADERS = ("content-length", "content-type", "content-encoding", "last-modified", "etag",
                "location")
METHOD = ("one HTTP GET with Python standard-library http.client (no redirects followed, no retry, "
          "Accept-Encoding: identity); body read into one fixed 1 MiB buffer with readinto(); "
          "hashlib.sha256 updated per chunk; chunks discarded; nothing written to disk")


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def verify(url: str) -> dict:
    parts = urlsplit(url)
    rec = {"status": None, "url": url, "method": METHOD, "started_utc": now(),
           "expected_bytes": EXPECTED_BYTES, "expected_sha256": EXPECTED_SHA256,
           "observed_bytes": 0, "observed_sha256": None, "http_status": None,
           "response_headers": {}, "written_to_disk": False, "phase": "connect", "error": None,
           "python": platform.python_version()}
    digest = hashlib.sha256()
    buf = bytearray(CHUNK_BYTES)
    view = memoryview(buf)
    t0 = time.monotonic()
    conn = http.client.HTTPSConnection(parts.hostname, parts.port or 443, timeout=TIMEOUT_S)
    try:
        conn.request("GET", parts.path, headers={"User-Agent": "compact-panels-verify/1.0",
                                                 "Accept-Encoding": "identity"})
        rec["phase"] = "response"
        resp = conn.getresponse()
        rec["http_status"] = resp.status
        rec["response_headers"] = {k: v for k, v in ((h.lower(), v) for h, v in resp.getheaders())
                                   if k in KEEP_HEADERS}
        if resp.status != 200:
            rec["error"] = f"HTTP {resp.status} {resp.reason}; not followed"
            rec["status"] = "TRANSFER_INCOMPLETE"
            return rec
        rec["phase"] = "body"
        next_report = 100 << 20
        while True:
            n = resp.readinto(buf)
            if not n:
                break
            digest.update(view[:n])
            rec["observed_bytes"] += n
            if rec["observed_bytes"] >= next_report:
                el = time.monotonic() - t0
                print(f"{rec['observed_bytes'] / 1e6:8.1f} MB  {el:6.1f} s  "
                      f"{rec['observed_bytes'] / 1e6 / el:5.1f} MB/s", file=sys.stderr, flush=True)
                next_report += 100 << 20
        rec["phase"] = "complete"
        rec["observed_sha256"] = digest.hexdigest()
        ok = rec["observed_bytes"] == EXPECTED_BYTES and rec["observed_sha256"] == EXPECTED_SHA256
        rec["status"] = "PUBLIC_SOURCE_SHA256_VERIFIED" if ok else "PUBLIC_SOURCE_MISMATCH"
        return rec
    except Exception as exc:  # report, never retry
        rec["error"] = f"{type(exc).__name__}: {exc}"
        rec["status"] = "TRANSFER_INCOMPLETE"
        if rec["observed_bytes"]:
            rec["sha256_of_received_prefix"] = digest.hexdigest()
        return rec
    finally:
        conn.close()
        rec["finished_utc"] = now()
        rec["elapsed_s"] = round(time.monotonic() - t0, 1)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--record", help="write the JSON verification record to this path")
    a = ap.parse_args()
    rec = verify(URL)
    text = json.dumps(rec, indent=2) + "\n"
    sys.stdout.write(text)
    sys.stdout.flush()
    if a.record:
        with open(a.record, "w", encoding="utf-8") as fh:
            fh.write(text)
    return 0 if rec["status"] == "PUBLIC_SOURCE_SHA256_VERIFIED" else 1


if __name__ == "__main__":
    sys.exit(main())
