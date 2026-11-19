import argparse
import json
from html import escape
from pathlib import Path

import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio

FOLDER = Path(__file__).resolve().parent


def mollweide(ra, dec):
    """RA increases leftward; RA=180 degrees lies at the map center."""
    ra, dec = np.broadcast_arrays(np.asarray(ra, float), np.asarray(dec, float))
    longitude = np.deg2rad(180 - ra)
    latitude = np.deg2rad(dec)
    # Bisection is stable at the poles, where Newton's derivative vanishes.
    lower = np.full(latitude.shape, -np.pi / 2)
    upper = np.full(latitude.shape, np.pi / 2)
    for _ in range(55):
        theta = (lower + upper) / 2
        below = 2 * theta + np.sin(2 * theta) < np.pi * np.sin(latitude)
        lower = np.where(below, theta, lower)
        upper = np.where(below, upper, theta)
    theta = np.where(np.abs(dec) == 90, np.sign(dec) * np.pi / 2, (lower + upper) / 2)
    return 2 * np.sqrt(2) / np.pi * longitude * np.cos(theta), np.sqrt(2) * np.sin(theta)


def read_positions(path, name_column):
    data = pd.read_csv(path)
    names = {col.lower(): col for col in ('RA', 'DEC', name_column)}
    data = data.rename(columns=lambda col: names.get(col.lower(), col))
    for col in names.values():
        if col not in data:
            raise ValueError(f'{path}: missing column {col!r}')
    for col in ('RA', 'DEC'):
        data[col] = pd.to_numeric(data[col])
    valid = data.RA.between(0, 360) & data.DEC.between(-90, 90)
    if not valid.all():
        raise ValueError(f'{path}: {int((~valid).sum())} invalid coordinate rows; expected degrees')
    if data.empty or data[name_column].isna().any():
        raise ValueError(f'{path}: empty table or missing object names')
    data['RA'] %= 360
    return data


def project_lines(lines):
    """Separate projected paths with gaps so Plotly draws each line independently."""
    xs, ys = [], []
    for ra, dec in lines:
        x, y = mollweide(ra, dec)
        xs.extend([*x, None])
        ys.extend([*y, None])
    return xs, ys


def add_grid(fig):
    meridians = [(ra, np.linspace(-90, 90, 181)) for ra in range(0, 361, 30)]
    parallels = [(np.linspace(0, 360, 361), dec) for dec in range(-60, 61, 30)]
    xs, ys = project_lines(meridians + parallels)
    fig.add_trace(go.Scatter(x=xs, y=ys, mode='lines', line=dict(color='#d7dfeb', width=0.7),
                             hoverinfo='skip', showlegend=False))
    angle = np.linspace(0, 2 * np.pi, 721)
    fig.add_trace(go.Scatter(x=2*np.sqrt(2)*np.cos(angle), y=np.sqrt(2)*np.sin(angle),
                             mode='lines', line=dict(color='#8091a8', width=1),
                             hoverinfo='skip', showlegend=False))
    for ra in range(30, 360, 30):
        x, y = mollweide(ra, 0)
        fig.add_annotation(x=float(x), y=float(y), text=f'{ra}°', showarrow=False,
                           yshift=-12, font=dict(size=10, color='#45556c'), bgcolor='rgba(255,255,255,0.8)')
    for dec in (-60, -30, 30, 60):
        x, y = mollweide(360, dec)
        fig.add_annotation(x=float(x), y=float(y), text=f'{dec:+d}°', showarrow=False,
                           xshift=-22, font=dict(size=11, color='#45556c'))


def add_catalogue(fig, data, label, color, size, name_column='AGN_name', gw=False):
    x, y = mollweide(data.RA, data.DEC)
    details = []
    for _, row in data.iterrows():
        lines = [f'<b>{escape(str(row[name_column]))}</b>',
                 f'RA: {row.RA:.6f}°', f'Dec: {row.DEC:+.6f}°']
        for col, title in [('z', 'Redshift'), ('t0', 'Flare peak (MJD)'), ('source', 'Source catalogue')]:
            if col in row and pd.notna(row[col]):
                lines.append(f'{title}: {escape(str(row[col]))}')
        if gw:
            lines.append('GW point estimate; uncertainty not shown')
        details.append('<br>'.join(lines))
    fig.add_trace(go.Scattergl(x=x, y=y, mode='markers', name=f'{label} ({len(data):,})',
        text=details, hovertemplate='%{text}<extra>%{fullData.name}</extra>',
        marker=dict(size=size, color=color, opacity=0.95 if gw else 0.65,
                    symbol='diamond' if gw else 'circle',
                    line=dict(width=1 if gw else 0, color='#111827'))))


def build_figure(coarse, refined, gw=None):
    fig = go.Figure()
    add_grid(fig)
    # The refined catalogue is a subset: draw each object only once.
    other = coarse.loc[~coarse.AGN_name.isin(refined.AGN_name)]
    add_catalogue(fig, other, 'AGN coarse only', '#4c89bd', 3)
    add_catalogue(fig, refined, 'AGN refined', '#d77816', 5)
    plane = SkyCoord(l=np.linspace(0, 360, 1441)*u.deg, b=np.zeros(1441)*u.deg,
                     frame='galactic').icrs
    order = np.argsort(plane.ra.deg)
    x, y = mollweide(plane.ra.deg[order], plane.dec.deg[order])
    fig.add_trace(go.Scatter(x=x, y=y, mode='lines', name='Galactic plane',
                            line=dict(color='#8b5e9e', width=1.5, dash='dot'), hoverinfo='skip'))
    if gw is not None:
        add_catalogue(fig, gw, 'GW candidates', '#d92e4b', 13, 'event', gw=True)
    fig.update_layout(template='plotly_white', height=760,
        margin=dict(l=55, r=40, t=35, b=80), paper_bgcolor='white', plot_bgcolor='white',
        font=dict(family='Arial, sans-serif', color='#22324a'),
        xaxis=dict(range=[-3.18, 3.12], visible=False, constrain='domain'),
        yaxis=dict(range=[-1.62, 1.62], visible=False, scaleanchor='x', scaleratio=1,
                   constrain='domain'),
        dragmode='pan', hovermode='closest',
        legend=dict(orientation='h', x=0.5, xanchor='center', y=-0.04, yanchor='top',
                    itemclick='toggle', itemdoubleclick='toggleothers'))
    return fig


def add_published_candidates(fig, events, pairs, contours):
    """Keep GW locations, localization contours, and optical positions distinct."""
    # Keep the general GWTC layer visible for "All", but hide it when the
    # user selects AGN-only or an individual curated event.
    membership = [({'__gwtc__'} if (trace.name or '').startswith('GWTC sky-map peaks') else None)
                  for trace in fig.data]
    colors = ['#c82836', '#7f7f7f', '#7b3fb2', '#008271', '#d85f00',
              '#245ad1', '#999999', '#c51b7d', '#726900', '#009cba']
    for index, row in events.iterrows():
        color = colors[index % len(colors)]
        name = row.event + (' [downgraded]' if 'Downgraded' in row.gw_status else '')
        for level, dash in [('50', 'solid'), ('90', 'dot')]:
            paths = contours['contours'][row.event][level]
            xs, ys = project_lines(np.asarray(path).T for path in paths)
            fig.add_trace(go.Scatter(x=xs, y=ys, mode='lines', name=f'{row.event}: {level}%',
                legendgroup=row.event, showlegend=False,
                line=dict(color=color, width=1.6 if level=='50' else 1.2, dash=dash),
                hovertemplate=f'{escape(row.event)} · {level}% sky credible contour<extra></extra>'))
            membership.append({row.event})
        x, y = mollweide(row.RA, row.DEC)
        text = (f'<b>{escape(row.event)}</b><br>GW probability-density peak (pixel center)'
                f'<br>RA {row.RA:.3f}° · Dec {row.DEC:+.3f}°'
                f'<br>50% area: {row.sky_area_50_deg2:.1f} deg²'
                f'<br>90% area: {row.sky_area_90_deg2:.1f} deg²'
                f'<br>{escape(row.gw_status)}<br>{escape(row.map_release)} · {escape(row.waveform)}'
                f'<br>Proposed counterparts: {escape(row.counterparts)}')
        fig.add_trace(go.Scatter(x=[float(x)], y=[float(y)], mode='markers', name=name,
            legendgroup=row.event, marker=dict(symbol='diamond', size=11, color=color,
            line=dict(width=1, color='white')), text=[text], hovertemplate='%{text}<extra></extra>'))
        membership.append({row.event})
    for host, group in pairs.groupby('counterpart', sort=False):
        row = group.iloc[0]
        x, y = mollweide(row.RA, row.DEC)
        text = (f'<b>{escape(host)}</b><br>Proposed optical counterpart; association unconfirmed'
                f'<br>RA {row.RA:.6f}° · Dec {row.DEC:+.6f}°<br>Redshift: {row.z}'
                f'<br>GW events: {escape("; ".join(group.event))}<br>{escape(row.reference)}')
        fig.add_trace(go.Scatter(x=[float(x)], y=[float(y)], mode='markers',
            name='Proposed optical counterparts', legendgroup='hosts',
            showlegend=host == pairs.counterpart.iloc[0], marker=dict(symbol='star', size=14,
            color='#f9d648', line=dict(width=1.4, color='#282130')),
            text=[text], hovertemplate='%{text}<extra></extra>'))
        membership.append(set(group.event))
    options = [('All candidates', None), ('AGN catalogues only', '')]
    options += [(event, event) for event in events.event]
    buttons = []
    for label, chosen in options:
        visible = [members is None or chosen is None or chosen in members
                   for members in membership]
        # The first host may be hidden when filtering: move its legend to a visible host.
        showlegend = [trace.showlegend for trace in fig.data]
        shown_host = False
        for i, trace in enumerate(fig.data):
            if trace.legendgroup == 'hosts':
                showlegend[i] = bool(visible[i] and not shown_host)
                shown_host |= visible[i]
        buttons.append(dict(label=label, method='update', args=[
            dict(visible=visible, showlegend=showlegend), {}]))
    fig.update_layout(height=880, margin=dict(l=55, r=40, t=75, b=170),
        updatemenus=[dict(buttons=buttons, direction='down', x=0.01, y=1.11,
                         xanchor='left', yanchor='top', active=0)],
        legend=dict(groupclick='togglegroup', y=-0.08, font=dict(size=11)))


def render_page(plot, status, gw, published, coarse_count, refined_count):
    html = '''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>AGN + GW · Interactive Mollweide sky map</title>
<style>body{margin:0;background:#f2f5f9;color:#22324a;font:15px/1.6 Arial,sans-serif}
main{max-width:1500px;margin:24px auto;padding:24px;background:white;border-radius:14px}
h1{font-size:26px;margin:0 0 5px}p{margin:6px 0}.muted{color:#58677b}
.status{background:#f3f6fb;padding:10px 14px;border-radius:7px;margin-top:12px}
footer{font-size:13px;border-top:1px solid #e2e8f0;padding-top:12px}
table{border-collapse:collapse;width:100%;font-size:13px}th,td{text-align:left;padding:8px;border-bottom:1px solid #e2e8f0}
a{color:#245ad1}.tablewrap{overflow-x:auto}h2{font-size:19px}
@media(max-width:700px){main{margin:0;padding:12px}h1{font-size:21px}}</style></head><body><main>
<h1>AGN + GW candidates</h1><p class="muted">Equatorial Mollweide projection · RA increases to the left · Center: RA 180°, Dec 0°</p>
<p>Hover for object details. Scroll to zoom, drag to pan, and double-click the map to reset.
Click a legend entry to toggle its layer; double-click it to isolate the layer.</p>
'''
    html += f'<p class="status">{status}</p>{plot}'
    if published:
        html += ('<h2>Candidate events and sources</h2><p>Curated literature sample: the nine O3 events '
                 'in <a href="https://arxiv.org/html/2209.13004">Graham et al. (2023), Table 3</a>, '
                 'plus S230922g and its favored candidate AT2023aagj from '
                 '<a href="https://arxiv.org/html/2407.10698">Cabrera et al. (2024)</a>. '
                 'This is not an exhaustive GW catalogue or a new association search.</p>'
                 '<p><a href="gw_events.csv" download>GW events CSV</a> · '
                 '<a href="gw_events.txt" download>Event names TXT</a> · '
                 '<a href="gw_agn_candidate_pairs.csv" download>Candidate pairings CSV</a></p>'
                 '<div class="tablewrap"><table><thead><tr><th>GW event</th><th>Proposed counterparts</th>'
                 '<th>90% area (deg²)</th><th>Map / waveform</th><th>GW status</th></tr></thead><tbody>')
        for row in gw.itertuples(index=False):
            html += (f'<tr><td><a href="{escape(row.source_url)}">{escape(row.event)}</a></td>'
                     f'<td>{escape(row.counterparts)}</td><td>{row.sky_area_90_deg2:.1f}</td>'
                     f'<td><a href="{escape(row.skymap_url)}">{escape(row.map_release)}</a> / '
                     f'{escape(row.waveform)}</td><td>{escape(row.gw_status)}</td></tr>')
        html += ('</tbody></table></div><p class="muted">Contours are display approximations on a 0.25° grid; '
                 'areas are computed from the native HEALPix probability maps. These selected waveform maps '
                 'can differ from the maps used in the counterpart papers. A published counterpart need not '
                 'lie within the selected map’s 90% region. The overlay tests neither distance nor timing. '
                 'See <a href="https://gwosc.org/GWTC-2/">GWOSC</a> for the downgraded event status and '
                 '<a href="https://doi.org/10.1093/mnras/stae2787">Veronesi et al.</a> for a reassessment '
                 'of the O3 associations.</p>')
    html += (f'<footer>{coarse_count:,} coarse AGN catalogue rows, including {refined_count:,} refined rows. '
             'The blue layer excludes the refined objects to avoid plotting them twice. '
             'AGN sources: AGNFCC.csv and AGNFRC.csv. Coordinates are treated as ICRS degrees. '
             'Sky overlap alone does not establish a physical association. '
             'This file includes its plotting library and works offline.</footer></main></body></html>')
    return html


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gw-csv', type=Path, default=FOLDER/'gw_events.csv',
                        help='GW positions: event, RA, DEC in degrees; default: gw_events.csv')
    parser.add_argument('--pairs-csv', type=Path, default=FOLDER/'gw_agn_candidate_pairs.csv')
    parser.add_argument('--contours-json', type=Path, default=FOLDER/'gw_data/gw_contours.json')
    parser.add_argument('--gwtc-positions', type=Path, default=FOLDER/'gw_data/gwtc_skymap_peaks.csv',
                        help='Map-derived GWTC sky peaks; omitted if file is absent')
    parser.add_argument('--output', type=Path, default=FOLDER/'AGN_GW_mollweide_interactive.html')
    args = parser.parse_args()
    coarse = read_positions(FOLDER/'AGNFCC.csv', 'AGN_name')
    refined = read_positions(FOLDER/'AGNFRC.csv', 'AGN_name')
    gw = read_positions(args.gw_csv, 'event') if args.gw_csv.exists() else None
    published = (gw is not None and 'skymap_file' in gw and
                 args.pairs_csv.exists() and args.contours_json.exists())
    fig = build_figure(coarse, refined, None if published else gw)
    gwtc = read_positions(args.gwtc_positions, 'event') if args.gwtc_positions.exists() else None
    if gwtc is not None:
        if gw is not None:
            gwtc = gwtc.loc[~gwtc.event.isin(gw.event)]
        add_catalogue(fig, gwtc, 'GWTC sky-map peaks', '#a52a5b', 9, 'event', gw=True)
    if published:
        pairs = pd.read_csv(args.pairs_csv)
        add_published_candidates(fig, gw, pairs, json.loads(args.contours_json.read_text()))
        status = (f'{len(gw)} curated GW events · {len(gwtc) if gwtc is not None else 0} additional GWTC sky-map peaks · '
                  f'{pairs.counterpart.nunique()} proposed optical counterparts · '
                  f'{len(pairs)} published pairings. Select an event with the menu above the map.<br>'
                  '<b>Diamonds:</b> GW probability-density peaks. <b>Stars:</b> proposed optical counterparts. '
                  '<b>Curated-event contours:</b> solid 50% and dotted 90% sky probability; '
                  'no contours are shown for the additional GWTC peaks. '
                  'Gray events were downgraded in GWTC-2.1. Associations are unconfirmed.')
    elif gw is not None:
        status = (f'{len(gw):,} GW candidates · Source: {escape(args.gw_csv.name)} · '
                  'Diamonds show supplied point estimates; localization uncertainty is not shown.')
    else:
        status = 'GW layer pending: no GW sky coordinates or localization maps were supplied.'
    plot = pio.to_html(fig, full_html=False, include_plotlyjs=True, config={
        'scrollZoom': True, 'displaylogo': False, 'responsive': True,
        'toImageButtonOptions': {'format': 'png', 'filename': 'AGN_GW_mollweide', 'scale': 3}})
    html = render_page(plot, status, gw, published, len(coarse), len(refined))
    args.output.write_text(html, encoding='utf-8')
    print(f'Saved {args.output} ({args.output.stat().st_size/1e6:.1f} MB)')
    print(status)


if __name__ == '__main__':
    main()
