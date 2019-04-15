import numpy as np
import matplotlib.pyplot as plt

def redundantly_populate(angles, values, refine_angle):
    """
    Parameters
    ----------
    angles : np.ndarray
    values : np.ndarray
    refine_angle : int
        Maximum angle width in radians

    Assumes evenly spaced azimuthal bins
    """
    n_angles = len(angles)
    bin_width = angles[1] - angles[0]
    refine_level = int(np.ceil(bin_width/refine_angle))
    # Check if refinement is required
    if refine_level < 2:
        return angles, values

    angle_increments = np.linspace(0, bin_width, refine_level)
    # new_angles = np.zeros(refine_level*(n_angles-1) + 1)
    new_angles = np.zeros(refine_level*(n_angles))
    new_values = np.zeros_like(new_angles)

    for i in range(n_angles):
        new_angles[i*refine_level:i*refine_level+refine_level] = angles[i] + angle_increments
        new_values[i*refine_level:i*refine_level+refine_level] = [values[i]]*refine_level
    # Might need to set first or last value of bins
    # new_angles[1] = angles[1]
    # new_values[1] = values[1]

    return new_angles, new_values

