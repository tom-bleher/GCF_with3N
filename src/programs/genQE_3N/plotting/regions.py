"""Line-bounded angular regions in (theta12, theta23) [deg].

Vectorized port of the ``is_in_region2`` lambda that
``analysis/SRC_analysis_3N.cpp`` contained before commit fe37fce (the
version at ``fe37fce^``); the current analysis writes a TTree and no longer
has it. Region division visualized at
https://www.desmos.com/calculator/x8p0w0wbjt.

Constants: a = 54 deg, b = a/2 = 27 deg, c = 160 deg, d = 7c/8 = 140 deg.

Main regions (per the C++ comment "0,2,3,6 are main regions"):
    0  head rocket          (lead rocket near (180, 0))
    2  tail rocket          (recoil rocket near (180, 180))
    3  tail rocket          (recoil rocket near (0, 180))
    6  star                 (symmetric, around (120, 120))
Transitions 1, 4, 5 lie between the main regions.
"""
import numpy as np

A = 54.0
B = A / 2.0     # 27
C = 160.0
D = 7.0 / 8.0 * C  # 140


def in_region(t12, t23, iregion):
    """Return a boolean ndarray selecting events in region ``iregion``.

    ``t12`` and ``t23`` are arrays of angles in degrees.
    """
    t12 = np.asarray(t12)
    t23 = np.asarray(t23)
    if iregion == 0:                                  # head rocket
        return t23 < 0.5 * (t12 - B)
    if iregion == 1:                                  # head <-> star
        return (t23 > 0.5 * (t12 - B)) & (t23 < 0.5 * (t12 + A))
    if iregion == 2:                                  # tail rocket (180, 180)
        return t23 > -(t12 - C) + C
    if iregion == 3:                                  # tail rocket (0, 180)
        return t23 > 2.0 * t12 + B
    if iregion == 4:                                  # (180,180) <-> star
        return (t23 > -(t12 - D) + D) & (t23 < -(t12 - C) + C)
    if iregion == 5:                                  # (0,180) <-> star
        return (t23 < 2.0 * t12 + B) & (t23 > 2.0 * t12 - A)
    if iregion == 6:                                  # star
        return ((t23 < 2.0 * t12 - A)
                & (t23 < -(t12 - D) + D)
                & (t23 > 0.5 * (t12 + A)))
    raise ValueError(f'Unknown region index: {iregion}')


# Convenience names for the main regions
REGION_NAMES = {
    0: 'Lead Rocket',     # head rocket (near (180, 0))
    2: 'Recoil Rocket',   # tail rocket (near (180, 180))
    3: 'Recoil Rocket',   # tail rocket (near (0, 180))   (alternate)
    6: 'Star',
}
