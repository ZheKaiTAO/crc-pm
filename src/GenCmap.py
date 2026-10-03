import numpy as np
import matplotlib.pyplot as plt
from numbers import Integral
from matplotlib.colors import LinearSegmentedColormap, rgb_to_hsv, hsv_to_rgb
from matplotlib.patches import Circle


# Hwang et al., Nature 2025 (user-supplied HEX transcription).
nature_hwang = [
    "#dd6572", "#e896a0", "#d96573", "#b9d696", "#5154a3", "#afc7e8",
    "#a3ca69", "#b35154", "#7ea254", "#7f7f7f", "#ffed6f", "#bc80bb",
]

# Butler et al., Nature 2025 (user-supplied HEX transcription).
nature_butler = [
    "#0f82bf", "#6ac6e9", "#3d4092", "#e92633", "#e4852b", "#fae41e",
    "#0a8648", "#83bd55", "#b96497",
]

# Jeffries et al., Nature 2025 (user-supplied HEX transcription).
nature_jeffries = [
    # red / pink
    "#cb2426", "#ea3b2c", "#aa1a7d", "#d83890", "#ed6b9f",
    # orange / yellow
    "#d64b23", "#f08c44", "#ca9424", "#f6be2a", "#ffd4af",
    # purple
    "#4c54a0", "#68559d", "#7f7cb6", "#9f9ac4", "#bebed8",
    # blue
    "#2573b4", "#02779b", "#4392c4", "#6baed5", "#9bc9dd",
    # green
    "#035830", "#148843", "#3bab5a", "#76c277", "#a2d59b",
]


_PALETTE_NAMES = (
    "nature_hwang",
    "nature_butler",
    "nature_jeffries",
    "saturated",
)


def lst(n=8):
    """Show every registered palette, wrapping at ``n`` colors per row.

    Return the Matplotlib ``(fig, ax)`` after displaying the preview.
    ``n`` must be a positive integer (booleans are not accepted).
    The continuous ``saturated()`` colormap is previewed with eight evenly
    spaced samples; changing ``n`` only changes the number of columns.
    """
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, Integral) or n <= 0:
        raise ValueError("n must be a positive integer")
    n = int(n)

    palettes = []
    for name in _PALETTE_NAMES:
        palette = globals()[name]
        if name == "saturated":
            palettes.append(("saturated()", palette()(np.linspace(0, 1, 8))))
        else:
            palettes.append((name, palette))
    columns = min(n, max(len(colors) for _, colors in palettes))
    row_counts = [(len(colors) + n - 1) // n for _, colors in palettes]
    group_gap = 0.6
    height = sum(row_counts) + group_gap * (len(palettes) - 1)
    label_width = 4.2

    fig, ax = plt.subplots(figsize=((label_width + columns) * 0.5, height * 0.5))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_xlim(-label_width, columns)
    ax.set_ylim(height - 0.5, -0.5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()

    row_start = 0
    for (name, colors), rows in zip(palettes, row_counts):
        ax.text(-label_width + 0.3, row_start, name, va="center", fontsize=10)
        for index, color in enumerate(colors):
            row, column = divmod(index, n)
            ax.add_patch(Circle((column + 0.5, row_start + row), radius=0.32,
                                facecolor=color, edgecolor="black", linewidth=0.8))
        row_start += rows + group_gap

    plt.show()
    return fig, ax

def saturated(cmap="Spectral_r", saturation=2, minval=0.08, maxval=0.92, n=256):
    cmap = plt.get_cmap(cmap)

    colors = cmap(np.linspace(minval, maxval, n))
    rgb = colors[:, :3]

    hsv = rgb_to_hsv(rgb)
    hsv[:, 1] = np.clip(hsv[:, 1] * saturation, 0, 1)

    rgb_sat = hsv_to_rgb(hsv)

    return LinearSegmentedColormap.from_list(
        f"{cmap.name}_sat_{saturation}",
        rgb_sat,
        N=n,
    )
