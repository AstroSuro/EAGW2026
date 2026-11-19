"""Prepare GW events and their published candidate optical counterparts.

Additional GW catalogues from the Gravitational Wave Open Science Center:
    GWTC-4.0 (O4a): https://gwosc.org/GWTC-4.0/
    GWTC-5.0 (O4b): https://gwosc.org/GWTC-5.0/
"""

from pathlib import Path
import hashlib
import json

from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
import astropy_healpix as ah
import contourpy
import numpy as np
import pandas as pd

FOLDER = Path(__file__).resolve().parent
MAPS = FOLDER / 'gw_data' / 'skymaps'
SOURCES = {
    'GWTC-2.1 v2': 'https://zenodo.org/api/records/6513631/files/IGWN-GWTC2p1-v2-PESkyMaps.tar.gz/content',
    'GWTC-2 v7': 'https://dcc.ligo.org/public/0169/P2000223/007/all_skymaps.tar',
    'GWTC-3 v2': 'https://zenodo.org/api/records/8177023/files/IGWN-GWTC3p0-v2-PESkyLocalizations.tar.gz/content',
    'GWTC-4.0': 'https://zenodo.org/records/17602505/files/IGWN-GWTC4p0-38214bd95_724-Archived_Skymaps.tar.gz',
    'GWTC-5.0': 'https://zenodo.org/records/20348005/files/IGWN-GWTC5p0-29ebe06b7_25-Archived_Skymaps.tar.gz',
    'GraceDB Bilby v0': 'https://gracedb.ligo.org/api/v2/superevents/S230922g/files/Bilby.fits.gz,0',
}
DOWNGRADED = {'GW190424_180648', 'GW190909_114149'}


def load_map(path):
    """Return native density/mass/area and a vectorized sky-density lookup.

    Multi-order pixels are located using their contiguous NESTED pixel ranges.
    This preserves native probability masses without rasterizing large maps.
    """
    with fits.open(path) as hdus:
        table = hdus[1].data
        header = hdus[1].header
        if header.get('COORDSYS') != 'C':
            raise ValueError(f'{path}: expected equatorial coordinates')
        if 'UNIQ' in table.names:
            order, ipix = ah.uniq_to_level_ipix(table['UNIQ'])
            nside = 2**order
            area = 4*np.pi/(12*nside*nside)
            density = np.asarray(table['PROBDENSITY'], float).copy()
            mass = density*area
            max_order = int(order.max())
            starts = ipix << (2*(max_order-order))
            stops = (ipix+1) << (2*(max_order-order))
            sorting = np.argsort(starts)
            starts, stops = starts[sorting], stops[sorting]
            assert starts[0] == 0 and stops[-1] == 12*4**max_order
            assert np.array_equal(starts[1:], stops[:-1]), 'Incomplete or overlapping MOC'

            def lookup(ra, dec):
                pix = ah.lonlat_to_healpix(ra*u.deg, dec*u.deg, 2**max_order, order='nested')
                return sorting[np.searchsorted(starts, pix, side='right')-1]

            peak = int(np.argmax(density))
            lon, lat = ah.healpix_to_lonlat(int(ipix[peak]), int(nside[peak]), order='nested')
        else:
            mass = np.asarray(table['PROB'], float).ravel().copy()
            nside = ah.npix_to_nside(len(mass))
            area = np.full(len(mass), 4*np.pi/len(mass))
            density = mass/area
            ordering = header['ORDERING'].strip()
            if ordering not in ('NESTED', 'RING'):
                raise ValueError(f'Unknown ordering {ordering}')
            ordering = ordering.lower()

            def lookup(ra, dec):
                return ah.lonlat_to_healpix(ra*u.deg, dec*u.deg, nside,
                                           order=ordering)

            lon, lat = ah.healpix_to_lonlat(int(np.argmax(density)), nside,
                                           order=ordering)
    total = mass.sum()
    assert np.isfinite(density).all() and (density >= 0).all()
    assert abs(total-1) < 1e-4, f'Probability normalization failed: {total}'
    mass /= total
    density /= total
    ranking = np.argsort(-density, kind='stable')
    cumulative = np.cumsum(mass[ranking])
    thresholds, areas = {}, {}
    for level in (0.5, 0.9):
        count = np.searchsorted(cumulative, level) + 1
        thresholds[level] = float(density[ranking[count-1]])
        areas[level] = float(area[ranking[:count]].sum()*(180/np.pi)**2)
    return density, lookup, thresholds, areas, (lon.deg, lat.deg), total


def choose_map(event):
    if event == 'S230922g':
        return MAPS/'S230922g_Bilby.fits.gz', 'GraceDB Bilby v0', 'Bilby'
    if event in DOWNGRADED:
        waveform = 'NRSur7dq4' if event.startswith('GW190424') else 'SEOBNRv4PHM'
        return MAPS/f'{event}_C01:{waveform}.fits', 'GWTC-2 v7', waveform
    key = 'GW190521_030229' if event == 'GW190521' else event
    release = 'GWTC-3 v2' if event.startswith('GW20') else 'GWTC-2.1 v2'
    candidates = list(MAPS.glob(f'*{key}*IMRPhenomXPHM*.fits'))
    if len(candidates) != 1:
        raise ValueError(f'Expected one IMRPhenomXPHM map for {event}; found {candidates}')
    return candidates[0], release, 'IMRPhenomXPHM'


def main():
    pairs = pd.read_csv(FOLDER/'gw_data/published_candidate_pairs.csv')
    coordinates = SkyCoord(pairs.position_sexagesimal.tolist(), frame='icrs')
    pairs['RA'], pairs['DEC'] = coordinates.ra.deg, coordinates.dec.deg
    pairs['position_type'] = 'Proposed optical counterpart; not a GW localization coordinate'
    pairs['retrieved_date'] = '2026-11-18'
    # Regular angular grid is used only for displaying contours, not for areas.
    ra = np.linspace(0, 360, 1441)
    dec = np.linspace(-90, 90, 721)
    ra_grid, dec_grid = np.meshgrid(ra, dec)
    events, contours = [], {}
    for event, group in pairs.groupby('event', sort=False):
        path, release, waveform = choose_map(event)
        density, lookup, thresholds, areas, peak, total = load_map(path)
        grid = density[lookup(ra_grid, dec_grid)]
        generator = contourpy.contour_generator(x=ra, y=dec, z=grid, name='serial')
        contours[event] = {}
        for level, threshold in thresholds.items():
            lines = generator.lines(threshold)
            assert lines, f'Missing {level:.0%} contour for {event}'
            contours[event][str(round(100 * level))] = [line.round(6).tolist() for line in lines]

        pixels = lookup(group.RA.to_numpy(), group.DEC.to_numpy())
        pairs.loc[group.index, 'inside_selected_map_90pct_sky'] = density[pixels] >= thresholds[0.9]
        events.append(dict(event=event, RA=peak[0], DEC=peak[1],
            position_type='Center of maximum probability-density HEALPix cell; not counterpart position',
            sky_area_50_deg2=areas[0.5], sky_area_90_deg2=areas[0.9],
            gw_status=('Downgraded in GWTC-2.1; historical candidate' if event in DOWNGRADED
                       else 'Public BBH alert' if event == 'S230922g' else 'Confident GWTC event'),
            counterpart_status='Proposed; unconfirmed',
            counterparts='; '.join(group.counterpart),
            map_release=release, waveform=waveform,
            skymap_file=str(path.relative_to(FOLDER)), skymap_url=SOURCES[release],
            skymap_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            native_probability_sum=total,
            reference=group.reference.iloc[0],
            source_url=group.source_url.iloc[0], retrieved_date='2026-11-18'))
        print(f'{event}: 90% area = {areas[0.9]:.1f} deg², {path.name}', flush=True)
    pd.DataFrame(events).to_csv(FOLDER/'gw_events.csv', index=False, float_format='%.8g')
    pairs.to_csv(FOLDER/'gw_agn_candidate_pairs.csv', index=False, float_format='%.9f')
    (FOLDER/'gw_events.txt').write_text('\n'.join(e['event'] for e in events)+'\n')
    (FOLDER/'gw_data/gw_contours.json').write_text(json.dumps(dict(
        grid_spacing_degrees=0.25, method='Native HEALPix density thresholds; marching-squares display contours',
        contours=contours), separators=(',', ':')))
    print(f'Saved {len(events)} GW events, {pairs.counterpart.nunique()} counterparts, {len(pairs)} pairs.')


if __name__ == '__main__':
    main()
