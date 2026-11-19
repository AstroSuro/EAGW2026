import argparse
import gzip
import hashlib
import io
import re
import tarfile
from pathlib import Path

from astropy.io import fits
import astropy_healpix as ah
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
SOURCE = ('https://zenodo.org/records/20348005/files/'
          'IGWN-GWTC5p0-29ebe06b7_25-Archived_Skymaps.tar.gz')
MD5 = '6a4b5d6fc3f62d350be9b73dae1796c8'
PATTERN = re.compile(r'-(GW\d{6}_\d{6})-(.+)_Skymap_PEDataRelease\.fits\.gz$')
PREFERENCE = ['Mixed', 'IMRPhenomXPHM_SpinTaylor', 'SEOBNRv5PHM', 'IMRPhenomXPNR', 'NRSur7dq4']


def peak(raw):
    # The release includes both gzip-compressed and uncompressed FITS bytes,
    # even though both kinds use the .fits.gz filename suffix.
    content = gzip.decompress(raw) if raw.startswith(b'\x1f\x8b') else raw
    with fits.open(io.BytesIO(content), memmap=False) as hdus:
        table, header = hdus[1].data, hdus[1].header
        if header.get('COORDSYS') != 'C':
            raise ValueError('Sky map is not in equatorial coordinates')
        if 'UNIQ' in table.names:
            level, pixels = ah.uniq_to_level_ipix(table['UNIQ'])
            nside = 2**level
            density = np.asarray(table['PROBDENSITY'], dtype=float)
            mass = density * (4*np.pi/(12*nside*nside))
            i = np.argmax(density)
            lon, lat = ah.healpix_to_lonlat(pixels[i], int(nside[i]), order='nested')
        else:
            mass = np.asarray(table['PROB'], dtype=float).ravel()
            nside = ah.npix_to_nside(len(mass))
            ordering = header.get('ORDERING', '').strip()
            if ordering not in ('RING', 'NESTED'):
                raise ValueError(f'Unknown HEALPix ordering: {ordering}')
            lon, lat = ah.healpix_to_lonlat(int(np.argmax(mass)), nside, order=ordering.lower())
        if not np.isfinite(mass).all() or (mass < 0).any() or abs(mass.sum()-1) > 1e-4:
            raise ValueError(f'Invalid sky-map probabilities: sum={mass.sum()}')
        return float(lon.deg), float(lat.deg)


def save_peak_map(peaks, path):
    from matplotlib.figure import Figure

    figure = Figure(figsize=(10, 6), layout='constrained')
    ax = figure.add_subplot(projection='mollweide')
    ax.scatter(np.deg2rad(180 - peaks.RA), np.deg2rad(peaks.DEC), s=18, color='#275c91')
    ax.set_xticks(np.deg2rad(np.arange(-150, 180, 30)),
                  labels=[f'{ra}°' for ra in range(330, 0, -30)])
    ax.set(xlabel='Right ascension', ylabel='Declination',
           title=f'GWTC-5.0 sky-localization peaks ({len(peaks)} events)')
    ax.grid(alpha=0.3)
    figure.savefig(path, format='png', dpi=200)


def read_peaks(path, names):
    chosen = {}
    with tarfile.open(path, 'r:gz') as archive:
        for member in archive:
            match = PATTERN.search(member.name)
            if not member.isfile() or not match:
                continue
            event, waveform = match.groups()
            if event not in names or waveform not in PREFERENCE:
                continue
            rank = PREFERENCE.index(waveform)
            if event not in chosen or rank < chosen[event][0]:
                chosen[event] = (rank, waveform, member.name, archive.extractfile(member).read())
    if not chosen:
        raise ValueError('No matching GWTC-5.0 maps found; existing output was not changed')
    rows = []
    for event, (_, waveform, member, raw) in sorted(chosen.items()):
        ra, dec = peak(raw)
        rows.append(dict(event=event, RA=ra, DEC=dec, map_release='GWTC-5.0 stable release 9',
                         waveform=waveform, position_type='center of maximum sky-probability-density pixel',
                         skymap_member=member, skymap_sha256=hashlib.sha256(raw).hexdigest(),
                         archive_url=SOURCE, archive_md5=MD5))
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    default_archive = ROOT/'gw_data'/SOURCE.rsplit('/', 1)[-1]
    parser.add_argument('archive', type=Path, nargs='?', default=default_archive,
                        help=f'Official sky-map tar.gz (default: gw_data/{default_archive.name})')
    parser.add_argument('--catalogue', type=Path, default=ROOT/'gw_data/gwtc_catalogue.csv',
                        help='Input event catalogue (default: gw_data/gwtc_catalogue.csv)')
    parser.add_argument('--output', type=Path, default=ROOT/'gw_data/gwtc_skymap_peaks.csv',
                        help='Output peak coordinates (.csv) or whole-sky plot (.png)')
    args = parser.parse_args()
    if not args.archive.is_file():
        parser.error(f'Archive file not found: {args.archive}. Download {SOURCE} '
                     f'to {default_archive} or pass the path to the downloaded file.')
    if not args.catalogue.is_file():
        parser.error(f'Catalogue file not found: {args.catalogue}')
    catalogue = pd.read_csv(args.catalogue)
    if not {'event', 'release'}.issubset(catalogue.columns):
        parser.error('--catalogue needs event and release columns. Use gw_data/gwtc_catalogue.csv.')
    if args.output.is_dir():
        parser.error('--output needs a filename, e.g. skymaps/gwtc_skymap_peaks.csv')
    digest = hashlib.md5()
    with args.archive.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    if digest.hexdigest() != MD5:
        raise ValueError('Archive MD5 does not match the official stable release 9 archive')
    names = set(catalogue.loc[catalogue.release.eq('GWTC-5.0'), 'event'])
    peaks = read_peaks(args.archive, names)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + '.tmp')
    if args.output.suffix.lower() == '.png':
        save_peak_map(peaks, temporary)
    else:
        peaks.to_csv(temporary, index=False, float_format='%.9f')
    temporary.replace(args.output)
    print(f'Saved {len(peaks)} GWTC-5.0 sky-map peaks to {args.output}; '
          f'{len(names)-len(peaks)} catalogue events have no selected map in the archive.')


if __name__ == '__main__':
    main()
