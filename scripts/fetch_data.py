#!/usr/bin/env python3
"""Download official PE samples, sky maps, and ZTF flare light curves.

Preferred entry points (via uv)::

    uv run python scripts/fetch_data.py --help
    uv run python scripts/fetch_data.py all --events GW190521 GW231123
    uv run python scripts/fetch_data.py pe --events GW190521
    uv run python scripts/fetch_data.py skymaps --events GW190521 --bulk gwtc2.1
    uv run python scripts/fetch_data.py zenodo-pe --catalogs gwtc5
    uv run python scripts/fetch_data.py everything          # ALL official PE + skymaps + ZTF
    uv run python scripts/fetch_data.py everything --dry-run # sizes only

Bulk Zenodo totals (approx., 2026-11)::

    gwtc2.1  ~23 GB   Zenodo 6513631
    gwtc3    ~14 GB   Zenodo 22685054
    gwtc4.0  ~14 GB   Zenodo 17602505
    gwtc4.1  ~15 GB   Zenodo 20275769
    gwtc5    ~55 GB   Zenodo 20348005 + 20348006
    skymaps  ~1 GB    tarballs inside the PE records above

Outputs land under ``data/`` (override with ``EAGW_DATA_ROOT`` or ``--data-root``).
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _default_data_root() -> Path:
    env = os.environ.get("EAGW_DATA_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    return ROOT / "data"


DATA = _default_data_root()
PE_DIR = DATA / "pe"
SKY_DIR = DATA / "skymaps"
ZTF_DIR = DATA / "ztf"

# Official Zenodo PE data releases (IGWN / LVK). Skymap tarballs live inside these.
ZENODO_PE: dict[str, dict] = {
    "gwtc2.1": {
        "records": ["6513631"],
        "label": "GWTC-2.1 PE (+ PESkyMaps tarball)",
        "skymap_tarball": "IGWN-GWTC2p1-v2-PESkyMaps.tar.gz",
        "approx_gb": 22.7,
    },
    "gwtc3": {
        "records": ["22685054"],
        "label": "GWTC-3 PE (+ PESkyLocalizations tarball)",
        "skymap_tarball": "IGWN-GWTC3p0-v3-PESkyLocalizations.tar.gz",
        "approx_gb": 14.1,
    },
    "gwtc4.0": {
        "records": ["17602505"],
        "label": "GWTC-4.0 PE (+ Archived_Skymaps tarball)",
        "skymap_tarball": "IGWN-GWTC4p0-38214bd95_724-Archived_Skymaps.tar.gz",
        "approx_gb": 14.4,
    },
    "gwtc4.1": {
        "records": ["20275769"],
        "label": "GWTC-4.1 PE (+ Archived_Skymaps tarball)",
        "skymap_tarball": "IGWN-GWTC4p1-18965dda8_5-Archived_Skymaps.tar.gz",
        "approx_gb": 14.5,
    },
    "gwtc5": {
        "records": ["20348005", "20348006"],
        "label": "GWTC-5.0 PE parts 1–2 (+ Archived_Skymaps in part 1)",
        "skymap_tarball": "IGWN-GWTC5p0-29ebe06b7_25-Archived_Skymaps.tar.gz",
        "skymap_record": "20348005",
        "approx_gb": 55.4,
    },
}

# Short names used in docs → GWOSC commonName / API slug when they differ.
EVENT_ALIASES = {
    "GW231123": "GW231123_135430",
    "GW190521": "GW190521",  # also GW190521_030229 in some releases
}

GWOSC_CATALOGS = [
    "GWTC-5.0",
    "GWTC-4.1",
    "GWTC-4.0",
    "GWTC-3-confident",
    "GWTC-2.1-confident",
    "GWTC-1-confident",
]

# AGN J124942.3+344929 / ZTF19abanrhr (Graham et al. 2020)
ZTF19ABANRHR = {
    "object": "ZTF19abanrhr",
    "ra": 192.426223,
    "dec": 34.824839,
    "radius_deg": 0.002,  # ~7 arcsec
}


def _set_data_root(root: Path) -> None:
    global DATA, PE_DIR, SKY_DIR, ZTF_DIR
    DATA = root.resolve()
    PE_DIR = DATA / "pe"
    SKY_DIR = DATA / "skymaps"
    ZTF_DIR = DATA / "ztf"


def _download(url: str, dest: Path, overwrite: bool = False) -> Path:
    """Download ``url`` to ``dest`` atomically, with retries.

    Writes to a ``.partial`` file and renames only on success, and checks the
    body against ``Content-Length``, so a truncated transfer can never be
    mistaken for a complete cached file.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not overwrite and dest.stat().st_size > 0:
        print(f"  exists, skip: {dest}")
        return dest
    print(f"  GET {url}")
    print(f"  ->  {dest}")
    tmp = dest.with_suffix(dest.suffix + ".partial")
    last: Exception | None = None
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(url, timeout=300) as resp, tmp.open("wb") as fh:
                expected = resp.headers.get("Content-Length")
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            got = tmp.stat().st_size
            if expected is not None and int(expected) != got:
                raise OSError(f"truncated download: {got} of {expected} bytes")
            tmp.replace(dest)
            return dest
        except Exception as exc:  # noqa: BLE001 - retry transient network errors
            last = exc
            tmp.unlink(missing_ok=True)
            print(f"  download failed ({exc}); retry {attempt}/3")
            time.sleep(2 * attempt)
    raise RuntimeError(f"could not download {url}") from last


def _validate_gzip(path: Path) -> None:
    """Raise unless ``path`` is a complete gzip stream (catches truncation)."""
    with gzip.open(path, "rb") as fh:
        while fh.read(1 << 20):
            pass


def _event_name_candidates(event: str) -> list[str]:
    names = [EVENT_ALIASES.get(event, event)]
    if event not in names:
        names.append(event)
    seen: set[str] = set()
    out: list[str] = []
    for name in names:
        if name not in seen:
            seen.add(name)
            out.append(name)
    return out


def _load_eventapi_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=60) as resp:
        payload = json.load(resp)
    events = payload.get("events") or payload
    if isinstance(events, dict) and events:
        return next(iter(events.values()))
    if isinstance(events, dict):
        raise RuntimeError(f"Empty events payload from {url}")
    return events


def _gwosc_event_json(event: str) -> dict:
    """Fetch GWOSC Event API JSON for a catalog event."""
    last_err: Exception | None = None

    for name in _event_name_candidates(event):
        v2 = f"https://gwosc.org/api/v2/events/{urllib.parse.quote(name)}?format=json"
        try:
            with urllib.request.urlopen(v2, timeout=60) as resp:
                meta = json.load(resp)
            versions = meta.get("versions") or []
            if not versions:
                continue
            detail_url = versions[-1]["detail_url"]
            with urllib.request.urlopen(detail_url, timeout=60) as resp:
                detail = json.load(resp)
            jsonurl = detail.get("jsonurl") or detail.get("eventapi_json_url")
            if jsonurl:
                return _load_eventapi_json(jsonurl)
            params_url = detail.get("parameters_url")
            if not params_url:
                params_url = detail_url.rstrip("/").split("?")[0] + "/parameters?format=json"
            with urllib.request.urlopen(params_url, timeout=60) as resp:
                params_payload = json.load(resp)
            parameters: dict = {}
            for row in params_payload.get("results") or []:
                parameters[row["name"]] = {
                    "data_url": row.get("data_url") or "",
                    "is_preferred": bool(row.get("is_preferred")),
                    "pipeline_type": row.get("pipeline_type"),
                    "waveform_family": row.get("waveform_family"),
                    "links": {
                        link.get("label"): link.get("url")
                        for link in (row.get("links") or [])
                        if isinstance(link, dict) and link.get("label") and link.get("url")
                    },
                }
            return {
                "commonName": meta.get("name") or name,
                "catalog.shortName": versions[-1].get("catalog"),
                "parameters": parameters,
            }
        except Exception as exc:  # noqa: BLE001
            last_err = exc

    for name in _event_name_candidates(event):
        for cat in GWOSC_CATALOGS:
            url = f"https://gwosc.org/eventapi/json/{cat}/{name}/"
            try:
                return _load_eventapi_json(url)
            except Exception as exc:  # noqa: BLE001
                last_err = exc
                continue

    raise RuntimeError(f"Could not resolve {event} on GWOSC Event API") from last_err


def list_catalog_events(catalog: str) -> list[str]:
    """Return commonNames for every event in a GWOSC catalog JSON dump."""
    url = f"https://gwosc.org/eventapi/json/{catalog}/"
    with urllib.request.urlopen(url, timeout=120) as resp:
        payload = json.load(resp)
    events = payload.get("events") or {}
    names: list[str] = []
    for _key, rec in events.items():
        if isinstance(rec, dict):
            names.append(str(rec.get("commonName") or _key.split("-")[0]))
        else:
            names.append(str(_key))
    # Unique, stable order
    seen: set[str] = set()
    out: list[str] = []
    for n in sorted(names):
        if n not in seen:
            seen.add(n)
            out.append(n)
    return out


def _is_pe_data_url(url: str) -> bool:
    low = url.lower()
    return url.startswith("http") and any(
        token in low for token in ("pedata", ".h5", ".hdf5", "posterior", "/content")
    )


def _pick_urls(event_record: dict) -> dict[str, str]:
    """Extract PE HDF5 and skymap URLs from a GWOSC event record."""
    out: dict[str, str] = {}

    for key, val in event_record.items():
        if not isinstance(val, str):
            continue
        low = key.lower()
        if "hdf5" in low and val.startswith("http"):
            out.setdefault("hdf5", val)
        if ("skymap" in low or "fits" in low) and val.endswith(".fits") and val.startswith("http"):
            out.setdefault("skymap", val)
        if low in {"hdf5_url", "pesummary_hdf5_url"} and val.startswith("http"):
            out["hdf5"] = val
        if low in {"skymap_fits_url", "skymap_url"} and val.startswith("http"):
            out["skymap"] = val

    parameters = event_record.get("parameters") or {}
    if isinstance(parameters, dict):
        candidates: list[tuple[int, int, str, dict]] = []
        for name, entry in parameters.items():
            if not isinstance(entry, dict):
                continue
            data_url = entry.get("data_url") or ""
            if _is_pe_data_url(data_url):
                candidates.append(
                    (
                        0 if entry.get("is_preferred") else 1,
                        0 if entry.get("pipeline_type") == "pe" else 1,
                        name,
                        entry,
                    )
                )
        candidates.sort()
        if candidates:
            _pref, _pipe, _name, entry = candidates[0]
            out["hdf5"] = entry["data_url"]
            links = entry.get("links") or {}
            skymap_url = None
            if isinstance(links, dict):
                skymap_url = links.get("skymap")
            elif isinstance(links, list):
                for link in links:
                    if isinstance(link, dict) and link.get("label") == "skymap":
                        skymap_url = link.get("url")
                        break
            if isinstance(skymap_url, str) and skymap_url.startswith("http"):
                out["skymap"] = skymap_url
                if any(tok in skymap_url.lower() for tok in (".tar.gz", "skymaps", "peskymaps")):
                    out["skymap_tarball"] = skymap_url

    return out


def _skymap_search_tokens(event: str, record: dict) -> list[str]:
    tokens = []
    for value in (
        record.get("commonName"),
        record.get("grace_id"),
        event,
        EVENT_ALIASES.get(event),
    ):
        if value and value not in tokens:
            tokens.append(str(value))
    for token in list(tokens):
        if "_" in token:
            short = token.split("_")[0]
            if short not in tokens:
                tokens.append(short)
    return tokens


def _extract_skymap_from_tarball(
    tarball: Path,
    event: str,
    record: dict,
    dest: Path,
) -> Path:
    """Extract one event FITS skymap from an official PESkyMaps tarball."""
    tokens = [t.lower() for t in _skymap_search_tokens(event, record)]
    with tarfile.open(tarball, "r:*") as tar:
        members = [
            m
            for m in tar.getmembers()
            if m.isfile()
            and (
                m.name.lower().endswith(".fits")
                or m.name.lower().endswith(".fits.gz")
            )
        ]
        matched = [
            m
            for m in members
            if any(token in m.name.lower().replace("-", "_") for token in tokens)
            or any(token in Path(m.name).name.lower().replace("-", "_") for token in tokens)
        ]
        if not matched:
            raise RuntimeError(
                f"No FITS member matching {tokens} in {tarball.name} "
                f"({len(members)} FITS files)"
            )

        def rank(member: tarfile.TarInfo) -> tuple[int, int, str]:
            name = member.name.lower()
            # Prefer Mixed maps; prefer uncompressed .fits when both exist.
            prefer_mixed = 0 if "mixed" in name and "mixed+xo4a" not in name else 1
            prefer_ext = 0 if name.endswith(".fits") else 1
            return (prefer_mixed, prefer_ext, name)

        member = sorted(matched, key=rank)[0]
        print(f"  extracting {member.name}")
        dest.parent.mkdir(parents=True, exist_ok=True)
        # Keep .fits.gz suffix if that is what the archive stores.
        out_path = dest
        if member.name.lower().endswith(".fits.gz") and not str(dest).endswith(".gz"):
            out_path = dest.with_suffix(dest.suffix + ".gz")
        with tar.extractfile(member) as src, out_path.open("wb") as out:
            assert src is not None
            out.write(src.read())
        if out_path != dest:
            print(f"  note: archive stores gzip FITS -> {out_path.name}")
        return out_path
    return dest


def _skymap_dest_exists(dest: Path) -> Path | None:
    """Return existing skymap path (.fits or .fits.gz) if present."""
    for path in (dest, dest.with_suffix(dest.suffix + ".gz")):
        if path.exists() and path.stat().st_size > 0:
            return path
    return None


def _dir_has_skymap_fits(path: Path) -> bool:
    """True if directory contains extracted .fits / .fits.gz skymaps."""
    if not path.exists():
        return False
    return any(path.rglob("*.fits")) or any(path.rglob("*.fits.gz"))


def _catalog_key_from_record(record: dict) -> str | None:
    """Map a GWOSC catalog.shortName to a ZENODO_PE key."""
    raw = str(record.get("catalog.shortName") or record.get("catalog") or "").lower()
    mapping = (
        ("gwtc-5", "gwtc5"),
        ("gwtc5", "gwtc5"),
        ("gwtc-4.1", "gwtc4.1"),
        ("gwtc4.1", "gwtc4.1"),
        ("gwtc-4.0", "gwtc4.0"),
        ("gwtc4.0", "gwtc4.0"),
        ("gwtc-3", "gwtc3"),
        ("gwtc3", "gwtc3"),
        ("gwtc-2.1", "gwtc2.1"),
        ("gwtc2.1", "gwtc2.1"),
        ("gwtc-2", "gwtc2.1"),
    )
    for needle, key in mapping:
        if needle in raw:
            return key
    return None


def fetch_pe_gwosc(events: list[str], overwrite: bool = False) -> list[Path]:
    """Download official PE HDF5 files via GWOSC Event API links."""
    saved: list[Path] = []
    for event in events:
        print(f"[pe] {event}")
        try:
            record = _gwosc_event_json(event)
            urls = _pick_urls(record)
            if "hdf5" not in urls:
                print(f"  WARNING: no hdf5 URL in GWOSC record; keys={sorted(record)}")
                print("  Use Zenodo bulk PE releases instead (`zenodo-pe` / `everything`).")
                continue
            label = record.get("commonName") or event
            dest = PE_DIR / f"{label}_PEDataRelease.hdf5"
            saved.append(_download(urls["hdf5"], dest, overwrite=overwrite))
        except Exception as exc:  # noqa: BLE001 — keep going for other events
            print(f"  ERROR fetching PE for {event}: {exc}")
    return saved


def fetch_skymaps_gwosc(
    events: list[str],
    overwrite: bool = False,
    bulk_fallback: str | None = None,
) -> list[Path]:
    """Download official PE sky-map FITS via GWOSC Event API links.

    If an event has no direct skymap URL (common for GWTC-4/5), falls back to the
    Zenodo archived skymap tarball for ``bulk_fallback`` or the event's catalog.
    """
    saved: list[Path] = []
    for event in events:
        print(f"[skymap] {event}")
        try:
            record = _gwosc_event_json(event)
            urls = _pick_urls(record)
            label = record.get("commonName") or event
            dest = SKY_DIR / f"{label}_skymap.fits"
            existing = _skymap_dest_exists(dest)
            if existing is not None and not overwrite:
                print(f"  exists, skip: {existing}")
                saved.append(existing)
                continue
            if "skymap_tarball" in urls or (
                "skymap" in urls
                and any(tok in urls["skymap"].lower() for tok in (".tar.gz", "skymaps"))
            ):
                tarball_url = urls.get("skymap_tarball") or urls["skymap"]
                url_path = Path(urllib.parse.urlparse(tarball_url).path)
                if url_path.name == "content" and url_path.parent.name:
                    tarball_name = url_path.parent.name
                else:
                    tarball_name = url_path.name or "skymaps.tar.gz"
                tarball = _download(tarball_url, SKY_DIR / tarball_name, overwrite=False)
                saved.append(_extract_skymap_from_tarball(tarball, event, record, dest))
                continue
            if "skymap" in urls:
                saved.append(_download(urls["skymap"], dest, overwrite=overwrite))
                continue

            # No direct URL — extract from Zenodo archived skymap tarball.
            catalog = bulk_fallback or _catalog_key_from_record(record)
            if catalog is None:
                print(f"  WARNING: no skymap URL in GWOSC record; keys={sorted(record)}")
                print("  Use Zenodo archived skymap tarballs instead (`--bulk` / `everything`).")
                continue
            print(f"  no direct skymap URL; falling back to Zenodo bulk '{catalog}'")
            tarball_path = fetch_skymaps_zenodo_bulk(catalog)
            # fetch_skymaps_zenodo_bulk returns extract dir; locate the .tar.gz beside it
            meta = ZENODO_PE[catalog]
            archive = SKY_DIR / meta["skymap_tarball"]
            if not archive.exists():
                print(f"  ERROR: expected archive missing: {archive}")
                continue
            saved.append(_extract_skymap_from_tarball(archive, event, record, dest))
        except Exception as exc:  # noqa: BLE001
            print(f"  ERROR fetching skymap for {event}: {exc}")
            if bulk_fallback:
                try:
                    print(f"  retry via Zenodo bulk '{bulk_fallback}'")
                    fetch_skymaps_zenodo_bulk(bulk_fallback)
                    meta = ZENODO_PE[bulk_fallback]
                    archive = SKY_DIR / meta["skymap_tarball"]
                    dest = SKY_DIR / f"{event}_skymap.fits"
                    record = {"commonName": EVENT_ALIASES.get(event, event)}
                    saved.append(
                        _extract_skymap_from_tarball(archive, event, record, dest)
                    )
                except Exception as exc2:  # noqa: BLE001
                    print(f"  ERROR bulk fallback for {event}: {exc2}")
    return saved


def zenodo_record_size_gb(record_id: str) -> tuple[float, int]:
    """Query Zenodo API for total file size (GB) and file count."""
    url = f"https://zenodo.org/api/records/{record_id}"
    with urllib.request.urlopen(url, timeout=60) as resp:
        payload = json.load(resp)
    files = payload.get("files") or []
    total = sum(f.get("size", 0) for f in files)
    return total / 1e9, len(files)


def fetch_zenodo_record(record_id: str, dest_dir: Path, files: list[str] | None = None) -> None:
    """Download a Zenodo record with zenodo_get (must be installed in the uv env)."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-m", "zenodo_get", record_id, "-o", str(dest_dir)]
    if files:
        for f in files:
            cmd.extend(["-g", f])
    print("  running:", " ".join(cmd))
    subprocess.check_call(cmd)


#: Backwards-compatible alias: notebooks import either spelling.
fetch_skymaps = fetch_skymaps_gwosc


def fetch_skymaps_zenodo_bulk(catalog: str) -> Path:
    """Download archived sky-map tarball from the official Zenodo PE data release."""
    if catalog not in ZENODO_PE:
        raise ValueError(f"Unknown catalog {catalog}; choose from {sorted(ZENODO_PE)}")
    meta = ZENODO_PE[catalog]
    name = meta["skymap_tarball"]
    record = meta.get("skymap_record") or meta["records"][0]
    SKY_DIR.mkdir(parents=True, exist_ok=True)
    url = (
        f"https://zenodo.org/api/records/{record}/files/"
        f"{urllib.parse.quote(name)}/content"
    )
    dest = SKY_DIR / name
    # A truncated archive must never be reused from the cache.
    if dest.exists():
        try:
            _validate_gzip(dest)
        except Exception as exc:  # noqa: BLE001 - corrupt cache
            print(f"  cached archive is incomplete ({exc}); re-downloading")
            dest.unlink()
    _download(url, dest)
    _validate_gzip(dest)
    extract_dir = SKY_DIR / catalog
    if not _dir_has_skymap_fits(extract_dir):
        print(f"  extracting -> {extract_dir}")
        extract_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(dest, "r:*") as tar:
            tar.extractall(extract_dir)
    return extract_dir


def fetch_zenodo_pe_catalogs(catalogs: list[str], dry_run: bool = False) -> None:
    """Download full official PE Zenodo records for the named catalogs."""
    total_gb = 0.0
    for cat in catalogs:
        if cat not in ZENODO_PE:
            raise SystemExit(f"Unknown catalog '{cat}'. Choose from: {', '.join(ZENODO_PE)}")
        meta = ZENODO_PE[cat]
        print(f"[zenodo-pe] {cat}: {meta['label']}")
        for rid in meta["records"]:
            try:
                gb, nfiles = zenodo_record_size_gb(rid)
            except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
                gb, nfiles = meta["approx_gb"] / max(len(meta["records"]), 1), -1
                print(f"  (size query failed for {rid}: {exc}; using approx {gb:.1f} GB)")
            total_gb += gb
            print(f"  record {rid}: ~{gb:.2f} GB ({nfiles} files)")
            if dry_run:
                continue
            fetch_zenodo_record(rid, PE_DIR / f"zenodo_{rid}")
    print(f"[zenodo-pe] total ~{total_gb:.1f} GB across {len(catalogs)} catalog(s)")


def fetch_ztf_irsa(
    ra: float | None = None,
    dec: float | None = None,
    radius_deg: float | None = None,
    object_name: str = "ZTF19abanrhr",
    overwrite: bool = False,
) -> Path:
    """Download ZTF light curves from the IRSA ZTF Light Curve API (official)."""
    ra = ZTF19ABANRHR["ra"] if ra is None else ra
    dec = ZTF19ABANRHR["dec"] if dec is None else dec
    radius_deg = ZTF19ABANRHR["radius_deg"] if radius_deg is None else radius_deg

    query = urllib.parse.urlencode(
        {
            "POS": f"CIRCLE {ra} {dec} {radius_deg}",
            "BAD_CATFLAGS_MASK": "32768",
            "FORMAT": "csv",
        }
    )
    url = f"https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?{query}"
    dest = ZTF_DIR / f"{object_name}_irsa.csv"
    _download(url, dest, overwrite=overwrite)

    lasair = ZTF_DIR / f"{object_name}_lasair.html"
    try:
        _download(
            f"https://lasair-ztf.lsst.ac.uk/objects/{object_name}/",
            lasair,
            overwrite=overwrite,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"  (optional Lasair page skipped: {exc})")

    meta = ZTF_DIR / f"{object_name}_meta.json"
    meta.write_text(
        json.dumps(
            {
                "object": object_name,
                "ra": ra,
                "dec": dec,
                "radius_deg": radius_deg,
                "host_agn": "J124942.3+344929",
                "source": "IRSA ZTF light-curve API",
                "api_docs": "https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html",
                "csv": str(dest.relative_to(DATA.parent)) if DATA.parent == ROOT else str(dest),
            },
            indent=2,
        )
        + "\n"
    )
    print(f"  wrote {meta}")
    return dest


def fetch_ztf_ztfquery(object_name: str = "ZTF19abanrhr") -> Path:
    """Optional: same light curve via ztfquery (IRSA under the hood)."""
    from ztfquery import lightcurve

    ra, dec, radius = (
        ZTF19ABANRHR["ra"],
        ZTF19ABANRHR["dec"],
        ZTF19ABANRHR["radius_deg"] * 3600.0,
    )
    print(f"[ztf/ztfquery] {object_name} @ ({ra}, {dec})")
    lcq = lightcurve.LCQuery.from_position(ra, dec, radius)
    dest = ZTF_DIR / f"{object_name}_ztfquery.csv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    lcq.data.to_csv(dest, index=False)
    print(f"  -> {dest}")
    return dest


def fetch_everything(
    catalogs: list[str] | None = None,
    include_skymaps: bool = True,
    include_ztf: bool = True,
    dry_run: bool = False,
) -> None:
    """Pull all official PE Zenodo releases (+ skymap tarballs + ZTF19abanrhr)."""
    catalogs = catalogs or list(ZENODO_PE.keys())
    print(f"[everything] catalogs={catalogs}")
    print(f"[everything] data root={DATA}")
    fetch_zenodo_pe_catalogs(catalogs, dry_run=dry_run)
    if include_skymaps:
        for cat in catalogs:
            if "skymap_tarball" not in ZENODO_PE[cat]:
                continue
            print(f"[everything/skymaps] {cat}")
            if dry_run:
                meta = ZENODO_PE[cat]
                print(f"  would fetch {meta['skymap_tarball']} from record "
                      f"{meta.get('skymap_record') or meta['records'][0]}")
                continue
            # Full PE download already includes the tarball; still extract a copy under skymaps/
            meta = ZENODO_PE[cat]
            name = meta["skymap_tarball"]
            # Prefer already-downloaded copy inside zenodo_* tree
            found = None
            for rid in meta["records"]:
                candidate = PE_DIR / f"zenodo_{rid}" / name
                if candidate.exists():
                    found = candidate
                    break
            if found is not None:
                dest = SKY_DIR / name
                SKY_DIR.mkdir(parents=True, exist_ok=True)
                if not dest.exists():
                    print(f"  linking/copying {found} -> {dest}")
                    try:
                        dest.hardlink_to(found)
                    except OSError:
                        import shutil

                        shutil.copy2(found, dest)
                extract_dir = SKY_DIR / cat
                if not _dir_has_skymap_fits(extract_dir):
                    print(f"  extracting -> {extract_dir}")
                    extract_dir.mkdir(parents=True, exist_ok=True)
                    with tarfile.open(dest, "r:*") as tar:
                        tar.extractall(extract_dir)
            else:
                fetch_skymaps_zenodo_bulk(cat)
    if include_ztf:
        print("[everything/ztf] ZTF19abanrhr")
        if not dry_run:
            fetch_ztf_irsa(object_name="ZTF19abanrhr")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--data-root",
        type=Path,
        default=None,
        help="Override data directory (default: $EAGW_DATA_ROOT or ./data)",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    pe = sub.add_parser("pe", help="Official PE HDF5 via GWOSC Event API")
    pe.add_argument("--events", nargs="+", default=["GW190521"])
    pe.add_argument(
        "--catalog",
        default=None,
        help="If set, fetch PE for every event in this GWOSC catalog "
        "(e.g. GWTC-5.0). Overrides --events.",
    )
    pe.add_argument("--overwrite", action="store_true")

    sky = sub.add_parser("skymaps", help="Official PE sky-map FITS via GWOSC Event API")
    sky.add_argument("--events", nargs="+", default=["GW190521"])
    sky.add_argument("--overwrite", action="store_true")
    sky.add_argument(
        "--bulk",
        choices=["none", *ZENODO_PE.keys()],
        default="none",
        help="Also download/extract Zenodo archived skymap tarball for a catalog",
    )

    ztf = sub.add_parser("ztf", help="ZTF flare light curves from IRSA")
    ztf.add_argument("--object", default="ZTF19abanrhr")
    ztf.add_argument("--ra", type=float, default=None)
    ztf.add_argument("--dec", type=float, default=None)
    ztf.add_argument("--radius-deg", type=float, default=None)
    ztf.add_argument("--ztfquery", action="store_true", help="Also fetch via ztfquery")
    ztf.add_argument("--overwrite", action="store_true")

    bulk = sub.add_parser("zenodo-pe", help="Bulk official PE from Zenodo (large!)")
    bulk.add_argument(
        "--catalogs",
        nargs="+",
        default=["gwtc5"],
        choices=list(ZENODO_PE.keys()),
        help="Which catalog releases to pull",
    )
    bulk.add_argument(
        "--records",
        nargs="+",
        default=None,
        help="Raw Zenodo record IDs (overrides --catalogs)",
    )
    bulk.add_argument("--dry-run", action="store_true", help="Print sizes only")

    all_p = sub.add_parser("all", help="PE + skymaps (GWOSC) + ZTF19abanrhr (IRSA)")
    all_p.add_argument(
        "--events",
        nargs="+",
        default=["GW190521"],
        help="Event names (GW231123 resolves to GW231123_135430)",
    )
    all_p.add_argument("--overwrite", action="store_true")

    everything = sub.add_parser(
        "everything",
        help="ALL official Zenodo PE releases + skymap tarballs + ZTF (~120 GB)",
    )
    everything.add_argument(
        "--catalogs",
        nargs="+",
        default=list(ZENODO_PE.keys()),
        choices=list(ZENODO_PE.keys()),
    )
    everything.add_argument("--no-skymaps", action="store_true")
    everything.add_argument("--no-ztf", action="store_true")
    everything.add_argument("--dry-run", action="store_true")

    lst = sub.add_parser("list-events", help="List GWOSC catalog event names")
    lst.add_argument(
        "--catalog",
        default="GWTC-5.0",
        help="GWOSC catalog short name (default GWTC-5.0)",
    )

    sizes = sub.add_parser("sizes", help="Show Zenodo PE release sizes")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.data_root is not None:
        _set_data_root(args.data_root)
    else:
        _set_data_root(_default_data_root())
    DATA.mkdir(parents=True, exist_ok=True)

    if args.cmd == "pe":
        events = list_catalog_events(args.catalog) if args.catalog else args.events
        if args.catalog:
            print(f"[pe] catalog {args.catalog}: {len(events)} events")
        fetch_pe_gwosc(events, overwrite=args.overwrite)
    elif args.cmd == "skymaps":
        bulk = None if args.bulk == "none" else args.bulk
        # Always attempt bulk first when requested so per-event GWOSC failures
        # (or events with no skymap URL, e.g. GW231123) can still extract FITS.
        if bulk is not None:
            try:
                fetch_skymaps_zenodo_bulk(bulk)
            except Exception as exc:  # noqa: BLE001
                print(f"  WARNING: bulk skymap download failed: {exc}")
        fetch_skymaps_gwosc(
            args.events,
            overwrite=args.overwrite,
            bulk_fallback=bulk,
        )
    elif args.cmd == "ztf":
        fetch_ztf_irsa(
            ra=args.ra,
            dec=args.dec,
            radius_deg=args.radius_deg,
            object_name=args.object,
            overwrite=args.overwrite,
        )
        if args.ztfquery:
            fetch_ztf_ztfquery(args.object)
    elif args.cmd == "zenodo-pe":
        if args.records:
            for rec in args.records:
                if args.dry_run:
                    gb, n = zenodo_record_size_gb(rec)
                    print(f"record {rec}: ~{gb:.2f} GB ({n} files)")
                else:
                    fetch_zenodo_record(rec, PE_DIR / f"zenodo_{rec}")
        else:
            fetch_zenodo_pe_catalogs(args.catalogs, dry_run=args.dry_run)
    elif args.cmd == "all":
        fetch_pe_gwosc(args.events, overwrite=args.overwrite)
        fetch_skymaps_gwosc(args.events, overwrite=args.overwrite)
        fetch_ztf_irsa(object_name="ZTF19abanrhr", overwrite=args.overwrite)
    elif args.cmd == "everything":
        fetch_everything(
            catalogs=args.catalogs,
            include_skymaps=not args.no_skymaps,
            include_ztf=not args.no_ztf,
            dry_run=args.dry_run,
        )
    elif args.cmd == "list-events":
        try:
            for name in list_catalog_events(args.catalog):
                print(name)
        except BrokenPipeError:
            return 0
    elif args.cmd == "sizes":
        grand = 0.0
        for cat, meta in ZENODO_PE.items():
            print(f"{cat}: {meta['label']}")
            for rid in meta["records"]:
                gb, n = zenodo_record_size_gb(rid)
                grand += gb
                print(f"  {rid}: {gb:.2f} GB ({n} files)")
        print(f"TOTAL ~{grand:.1f} GB")
    else:
        raise SystemExit(f"unknown command {args.cmd}")
    print("Done.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        raise SystemExit(0)
