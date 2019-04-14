import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle

from radial_plot import redundantly_populate


class Viewer(object):
    def __init__(self, df):
        # Initialize data attribues.
        self.df = df.copy()
        self.energy_edges = np.array([])
        self.azim_edges = np.array([])
        self.coords = [0, 0, 0.0, 0.025]

        # Initialize the figure and the axes.
        self.fig = plt.figure()
        self.pin_ax = self.fig.add_axes((0.05, 0.05, 0.40, 0.90),
             aspect='equal')
        self.spec_ax = self.fig.add_axes((0.55, 0.55, 0.40, 0.40))
        self.spec_ax.set_xscale('log')
        self.azim_ax = self.fig.add_axes((0.55, 0.05, 0.40, 0.40), 
                                         projection='polar')

        # Format the data.
        self.format_dataframe()

        # Draw the pincell and register event handlers.
        self.draw_pincell()
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)

        # Initialize plot actors.
        self.cur_square = None
        self.energy_line = None
        self.azim_line = None
        self.energy_marker = None
        self.azim_marker = None
        self.azim_ray = None

        # Draw the plots with the starting data.
        self.update_plot()

    def format_dataframe(self):
        # Get rid of useless columns.
        del self.df['nuclide']
        del self.df['score']

        # Replace the mesh index labels.
        self.df['i_x'] = self.df[('mesh 1', 'x')] - 1
        self.df['i_y'] = self.df[('mesh 1', 'y')] - 1
        del self.df[('mesh 1', 'x')]
        del self.df[('mesh 1', 'y')]
        del self.df[('mesh 1', 'z')]

        # Remove rows for the unused mesh bins above the diagonal.
        self.df = self.df[self.df['i_x'] >= self.df['i_y']]

        # Extract the energy and azimuthal bin edges.
        self.energy_edges = np.unique(self.df['energy low [MeV]'])
        self.energy_edges = np.concatenate((self.energy_edges,
                                          [self.df['energy high [MeV]'].max()]))
        self.azim_edges = np.unique(self.df['azimuthal low'])
        self.azim_edges = np.concatenate((self.azim_edges,
                                          [self.df['azimuthal high'].max()]))

        # Replace energy bounds with bin numbers.
        bins = self.df['energy low [MeV]']
        bins = [np.argmin(np.abs(self.energy_edges - a)) for a in bins]
        self.df['energy'] = bins
        del self.df['energy low [MeV]']
        del self.df['energy high [MeV]']

        # Replace azimuthal bounds with bin numbers.
        bins = self.df['azimuthal low']
        bins = [np.argmin(np.abs(self.azim_edges - a)) for a in bins]
        self.df['azim'] = bins
        del self.df['azimuthal low']
        del self.df['azimuthal high']

        # Convert MeV to eV.
        self.energy_edges *= 1e6

        # Double fluxes on boundary elements to account for area difference.
        diagonal = self.df['i_x'] == self.df['i_y']
        bottom = self.df['i_y'] == 0
        right = self.df['i_x'] == self.df['i_x'].max()
        corners = diagonal & (bottom | right)
        for rows in (diagonal, bottom, right, corners):
            for col in ('mean', 'std. dev.'):
                self.df.loc[rows, col] *= 2

    def draw_pincell(self):
        # Draw patches for the materials.
        fuel_color = np.array([255, 50, 50]) / 255.0
        gap_color = np.array([220, 220, 220]) / 255.0
        clad_color = np.array([150, 150, 150]) / 255.0
        mod_color = np.array([100, 200, 200]) / 255.0
        pitch = 0.62992 * 2.0
        mod_patch = matplotlib.patches.Rectangle((-pitch/2.0, -pitch/2.0),
             pitch, pitch, facecolor=mod_color)
        clad_patch = matplotlib.patches.Circle((0.0, 0.0), radius=0.45720,
             facecolor=clad_color)
        gap_patch = matplotlib.patches.Circle((0.0, 0.0), radius=0.40005,
             facecolor=gap_color)
        fuel_patch = matplotlib.patches.Circle((0.0, 0.0), radius=0.39218,
             facecolor=fuel_color)
        self.pin_ax.add_patch(mod_patch)
        self.pin_ax.add_patch(clad_patch)
        self.pin_ax.add_patch(gap_patch)
        self.pin_ax.add_patch(fuel_patch)

        # Draw centerlines.
        self.pin_ax.plot((-pitch/2.0, pitch/2.0), (0.0, 0.0), c='black',
                         linestyle='-.')
        self.pin_ax.plot((0.0, 0.0), (-pitch/2.0, pitch/2.0), c='black',
                         linestyle='-.')

        # Draw mesh gridlines.
        delta = pitch / 2.0 / (10 - 1)
        d0 = -delta / 2.0
        self.pin_ax.plot((d0, pitch/2.0), (d0, d0), c='black')
        self.pin_ax.plot((d0, d0), (d0, d0+delta), c='black')
        for i in range(10):
            d1 = d0 + i*delta
            d2 = d1 + delta
            self.pin_ax.plot((d1, pitch/2.0), (d2, d2), c='black')
            self.pin_ax.plot((d1, d1), (d0, d2), c='black')
        self.pin_ax.set_xlim((-pitch/2.0, pitch/2.0))
        self.pin_ax.set_ylim((-pitch/2.0, pitch/2.0))

    def onclick(self, event):
        if event.inaxes == self.pin_ax:
            # Index the selected bin.
            pitch = 0.62992 * 2.0
            delta = pitch / 2.0 / 9
            x = event.xdata + delta/2.0
            y = event.ydata + delta/2.0
            i_x = int(x / delta)
            i_y = int(y / delta)

            # Ignore the click if it was outside of the mesh.
            if i_x < 0 or i_y < 0 or i_y > i_x: return

            # Update the coordinates.
            self.coords[0] = i_x
            self.coords[1] = i_y

            # Update the plot.
            self.update_plot()

        elif event.inaxes == self.spec_ax:
            self.coords[3] = event.xdata
            self.update_plot()

        elif event.inaxes == self.azim_ax:
            self.coords[2] = event.xdata
            self.update_plot()

    def update_plot(self):
        # Unpack the current phase-space coordinates.
        i_x, i_y, azim, energy = self.coords

        # Find the azimuthal and energy bin indices.
        a_bin = np.searchsorted(self.azim_edges, azim) - 1
        e_bin = np.searchsorted(self.energy_edges, energy) - 1

        # Highlight the selected spatial bin.
        pitch = 0.62992 * 2.0
        delta = pitch / 2.0 / 9
        if self.cur_square is not None: self.cur_square.remove()
        self.cur_square = matplotlib.patches.Rectangle(
             ((i_x-0.5)*delta, (i_y-0.5)*delta), delta, delta,
             facecolor='green')
        self.pin_ax.add_patch(self.cur_square)

        # Plot the energy spectrum.
        df = self.df[(self.df['i_x'] == i_x) & (self.df['i_y'] == i_y)
                     & (self.df['azim'] == a_bin)]
        flux = np.concatenate(([df['mean'].values[0]], df['mean'].values))
        if self.energy_line is not None: self.energy_line.remove()
        self.energy_line, = self.spec_ax.step(self.energy_edges, flux,
                                              c='green')
        if self.energy_marker is not None: self.energy_marker.remove()
        self.energy_marker = self.spec_ax.axvline(energy)

        # Plot the azimuthal spectrum.
        df = self.df[(self.df['i_x'] == i_x) & (self.df['i_y'] == i_y)
                     & (self.df['energy'] == e_bin)]
        flux = np.concatenate(([df['mean'].values[0]], df['mean'].values))
        if self.azim_line is not None: self.azim_line.remove()
        # self.azim_line, = self.azim_ax.step(self.azim_edges, flux, c='green')
        angle_refine = np.pi/360
        angles, flux = redundantly_populate(self.azim_edges[:-1], flux[1:], angle_refine)
        self.azim_line, = self.azim_ax.plot(angles, flux, c='green')
        ### TODO fix marker
        if self.azim_marker is not None: self.azim_marker.remove()
        self.azim_marker = self.azim_ax.axvline(azim)

        # Adjust the azimuthal spectrum yscale.
        pad = 0.05 * (flux.max() - flux.min())
        pad = max(pad, 0.1*flux.max())
        self.azim_ax.set_rmax(flux.max() + pad)
        self.azim_ax.set_rmin(0)
        self.azim_ax.set_theta_zero_location("W")

        # Draw an azimuthal ray on the spatial plot.
        x0 = i_x*delta
        y0 = i_y*delta
        x_target = pitch/2.0 if abs(azim) > np.pi / 2.0 else -pitch/2.0
        y_target = pitch/2.0 if azim < 0.0 else -pitch/2.0
        d0 = (x_target - x0) / np.cos(azim)
        d1 = (y_target - y0) / np.sin(azim)
        if d0 < d1:
            x1 = x_target
            y1 = y0 + d0*np.sin(azim)
        else:
            x1 = x0 + d1*np.cos(azim)
            y1 = y_target
        if self.azim_ray is not None: self.azim_ray.remove()
        self.azim_ray, = self.pin_ax.plot((x0, x1), (y0, y1), c='black',
                                          linewidth=2)

        # Redraw the image.
        self.fig.canvas.draw()


if __name__ == '__main__':
    # Get the path to the data directory.
    my_path = os.path.abspath(__file__)
    my_dir = os.path.dirname(my_path)
    root_path = os.path.join(my_dir, os.pardir)
    data_path = os.path.join(root_path, 'data')

    # Pick an appropriate data file.
    data_file = 'fine.p'
    data_path = os.path.join(data_path, data_file)

    # Read the tally data and start the gui.
    with open(data_path, 'rb') as fh:
        df = pickle.load(fh)
    viewer = Viewer(df)
    plt.show()
