import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, rgb_to_hsv, hsv_to_rgb

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