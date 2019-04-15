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
    n_repeat = int(np.ceil(bin_width/refine_angle))
    # Check if refinement is required, return original data if not
    if n_repeat < 2:
        return angles, values

    angle_increments = np.linspace(0, bin_width, n_repeat)
    new_angles = np.zeros(n_repeat*(n_angles))
    new_values = np.zeros_like(new_angles)

    for i in range(n_angles):
        new_angles[i*n_repeat:i*n_repeat+n_repeat] = angles[i] + angle_increments
        new_values[i*n_repeat:i*n_repeat+n_repeat] = [values[i]]*n_repeat

    return new_angles, new_values

