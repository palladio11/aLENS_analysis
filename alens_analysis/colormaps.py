#!/usr/bin/env python

"""@package docstring
File: colomaps.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

# Setup new colormap


def create_cmaps(cmap_data):
    """
    Generate a dict mapping standard colormap names to standard colormaps, as
    well as the reversed colormaps.
    """
    cmap_d = dict()
    lutsize = mpl.rcParams['image.lut']
    for name, spec in cmap_data.items():
        cmap_d[name] = (
            mpl.colors.LinearSegmentedColormap(name, spec, lutsize)
            if 'red' in spec else
            mpl.colors.ListedColormap(spec['listed'], name)
            if 'listed' in spec else
            mpl.colors.LinearSegmentedColormap.from_list(name, spec, lutsize))

    # Generate reversed cmaps.
    for cmap in list(cmap_d.values()):
        rmap = cmap.reversed()
        cmap._global = True
        rmap._global = True
        cmap_d[rmap.name] = rmap
    for k, v in cmap_d.items():
        plt.register_cmap(name=k, cmap=v)


# Based on https://doi.org/10.1016/j.molcel.2018.05.032
emct8_data = {'red': ((0.00, 1.00, 1.00),
                      (0.33, 0.26953125, 0.26953125),
                      (0.66, 0.40625, 0.40625),
                      (0.95, 0.91015625, 0.91015625),
                      (1.00, 0.9453125, 0.9453125)),

              'green': ((0.00, 1.00, 1.00),
                        (0.33, 0.33203125, 0.33203125),
                        (0.66, 0.078125, 0.078125),
                        (0.95, 0.83203125, 0.83203125),
                        (1.00, 0.8828125, 0.8828125)),

              'blue': ((0.00, 1.00, 1.00),
                       (0.33, 0.6328125, 0.6328125),
                       (0.66, 0.0703125, 0.0703125),
                       (0.95, 0.3359375, 0.3359375),
                       (1.00, 0.6484375, 0.6484375))}


#plt.rcParams['image.cmap'] = 'emct8'
#plt.rcParams['image.cmap'] = 'warm'


def register_cmaps():
    """Register useful colormaps
    @return: None

    """
    create_cmaps({"emct8": emct8_data})
    cm = mpl.cm.coolwarm
    cdict = {"red": [], "green": [], "blue": [], "alpha": []}

    # Get top range of coolwarm and interpolate between.
    for i in range(127, 256):
        r, g, b, a = cm(i)
        si = ((i - 127) * 2) / 256
        cdict["red"].append((si, r, r))
        cdict["green"].append((si, g, g))
        cdict["blue"].append((si, b, b))
        cdict["alpha"].append((si, a, a))
    warm = mpl.colors.LinearSegmentedColormap("warm", cdict)
    plt.register_cmap(cmap=warm)


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
