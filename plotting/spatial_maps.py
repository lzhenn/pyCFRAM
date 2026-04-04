"""Spatial map plots for CFRAM decomposition results.

Reproduces Fig.3 of Wu et al. (2025): multi-panel maps of partial
temperature changes from each physical process.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False

from .style import (LEVELS_DT, LEVELS_AER, CMAP_DIVERGE,
                    FONTSIZE_PANEL, setup_style, make_norm)


def plot_decomposition_maps(data_dict, lats, lons, title='',
                            levels=None, figsize=(12, 16),
                            key_areas=None, save_path=None):
    """Plot multi-panel spatial maps of CFRAM decomposition.

    Args:
        data_dict: OrderedDict or dict of {label: (nlat, nlon) array}
            Each entry is one panel. Typical order:
            WV, Cloud, Aerosol, Surface, Atmos.dyn, Sum
        lats: 1D latitude array
        lons: 1D longitude array
        title: figure super title
        levels: contour levels (default: LEVELS_DT)
        figsize: figure size
        key_areas: list of dicts with 'lon' and 'lat' tuples for box overlay
        save_path: if provided, save figure to this path
    """
    setup_style()
    if levels is None:
        levels = LEVELS_DT
    norm = make_norm(levels)

    n_panels = len(data_dict)
    ncols = 2
    nrows = (n_panels + 1) // 2

    projection = ccrs.PlateCarree() if HAS_CARTOPY else None

    fig, axes = plt.subplots(
        nrows, ncols, figsize=figsize,
        subplot_kw={'projection': projection} if HAS_CARTOPY else {},
    )
    axes = axes.flatten()

    lon2d, lat2d = np.meshgrid(lons, lats)

    for idx, (label, field) in enumerate(data_dict.items()):
        ax = axes[idx]

        if HAS_CARTOPY:
            ax.set_extent([lons.min(), lons.max(), lats.min(), lats.max()],
                          crs=ccrs.PlateCarree())
            ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
            ax.add_feature(cfeature.BORDERS, linewidth=0.3)
            gl = ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)
            gl.top_labels = False
            gl.right_labels = False

            cf = ax.contourf(lon2d, lat2d, field, levels=levels,
                             cmap=CMAP_DIVERGE, norm=norm, extend='both',
                             transform=ccrs.PlateCarree())
        else:
            cf = ax.contourf(lon2d, lat2d, field, levels=levels,
                             cmap=CMAP_DIVERGE, norm=norm, extend='both')
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')

        ax.set_title(f'({chr(97+idx)}) {label}', fontsize=FONTSIZE_PANEL,
                     loc='left', fontweight='bold')

        # Overlay key area box
        if key_areas:
            for ka in key_areas:
                lon0, lon1 = ka['lon']
                lat0, lat1 = ka['lat']
                rect_lons = [lon0, lon1, lon1, lon0, lon0]
                rect_lats = [lat0, lat0, lat1, lat1, lat0]
                if HAS_CARTOPY:
                    ax.plot(rect_lons, rect_lats, 'b-', linewidth=1.5,
                            transform=ccrs.PlateCarree())
                else:
                    ax.plot(rect_lons, rect_lats, 'b-', linewidth=1.5)

    # Hide unused axes
    for idx in range(n_panels, len(axes)):
        axes[idx].set_visible(False)

    # Colorbar
    fig.subplots_adjust(bottom=0.08, hspace=0.25, wspace=0.15)
    cbar_ax = fig.add_axes([0.15, 0.03, 0.7, 0.015])
    cbar = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal',
                        ticks=levels)
    cbar.set_label('K', fontsize=FONTSIZE_PANEL)

    if title:
        fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig


def plot_eh13_eh22_comparison(results_eh13, results_eh22, lats, lons,
                               key_area_13=None, key_area_22=None,
                               save_path=None):
    """Plot side-by-side EH13 vs EH22 decomposition (Fig.3 style).

    Args:
        results_eh13: dict with 'dT_q', 'dT_cloud', etc. as (nlat, nlon) arrays
        results_eh22: same for EH22
        lats, lons: coordinate arrays
        save_path: output file path
    """
    terms = [
        ('Water vapor', 'dT_q'),
        ('Cloud', 'dT_cloud'),
        ('Aerosol', 'dT_aerosol'),
        ('Surface', 'dT_ts'),
        ('Atmos.dyn', 'dT_atd'),
        ('Sum', 'dT_sum'),
    ]

    setup_style()
    norm = make_norm(LEVELS_DT)

    projection = ccrs.PlateCarree() if HAS_CARTOPY else None
    fig, axes = plt.subplots(
        6, 2, figsize=(10, 18),
        subplot_kw={'projection': projection} if HAS_CARTOPY else {},
    )

    lon2d, lat2d = np.meshgrid(lons, lats)

    for row, (label, key) in enumerate(terms):
        for col, (results, yr, ka) in enumerate([
            (results_eh13, '2013', key_area_13),
            (results_eh22, '2022', key_area_22),
        ]):
            ax = axes[row, col]
            field = results.get(key, np.zeros((len(lats), len(lons))))

            if HAS_CARTOPY:
                ax.set_extent([lons.min(), lons.max(), lats.min(), lats.max()],
                              crs=ccrs.PlateCarree())
                ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
                cf = ax.contourf(lon2d, lat2d, field, levels=LEVELS_DT,
                                 cmap=CMAP_DIVERGE, norm=norm, extend='both',
                                 transform=ccrs.PlateCarree())
            else:
                cf = ax.contourf(lon2d, lat2d, field, levels=LEVELS_DT,
                                 cmap=CMAP_DIVERGE, norm=norm, extend='both')

            panel_idx = row * 2 + col
            ax.set_title(f'({chr(97+panel_idx)}) {label}',
                         fontsize=9, loc='left')

            if row == 0:
                ax.set_title(f'{yr}\n({chr(97+panel_idx)}) {label}',
                             fontsize=9, loc='left')

            if ka:
                rect_lons = [ka['lon'][0], ka['lon'][1], ka['lon'][1],
                             ka['lon'][0], ka['lon'][0]]
                rect_lats = [ka['lat'][0], ka['lat'][0], ka['lat'][1],
                             ka['lat'][1], ka['lat'][0]]
                if HAS_CARTOPY:
                    ax.plot(rect_lons, rect_lats, 'b-', lw=1.5,
                            transform=ccrs.PlateCarree())
                else:
                    ax.plot(rect_lons, rect_lats, 'b-', lw=1.5)

    fig.subplots_adjust(bottom=0.05, hspace=0.2, wspace=0.1)
    cbar_ax = fig.add_axes([0.15, 0.02, 0.7, 0.012])
    fig.colorbar(cf, cax=cbar_ax, orientation='horizontal',
                 ticks=LEVELS_DT, label='K')

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig
