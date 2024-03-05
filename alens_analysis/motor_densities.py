#!/usr/bin/env python

import numpy as np
import scipy
import matplotlib.pyplot as plt
import h5py
import json
import yaml
import pathlib
from collections import defaultdict


class Empirical_Motor_Density_Constructor:

    def __init__(self, simulation_data_path, simulation_config_file, retrieve_num_T_steps="all", autosave=False, **kwargs):

        self.simulation_data_path = simulation_data_path
        self.simulation_config_file = simulation_config_file

        # Declare self.P_data, self.S_data, self.num_T_steps
        self.retrieve_simulation_data(retrieve_num_T_steps)

        # Declare self.box_side_length, self.syl_L (rod/sylinder length)
        self.retrieve_simulation_config_data()

        # Construct self.empirical_motor_density
        self.construct_empirical_motor_density()

        if autosave:
            self.save_empirical_motor_density(kwargs["save_file_name"])

    def save_empirical_motor_density(self, save_file_name):

        with save_file_name.open("w") as f:
            json.dump(self.empirical_motor_density, f)

    def retrieve_simulation_data(self, retrieve_num_T_steps="all"):
        """
        Removes first two timesteps where there are no attachments.
        Removes P and S data which isn't coordinates or int indices.
        """
        with h5py.File(self.simulation_data_path, 'r+') as h5_data:

            if retrieve_num_T_steps == "all": retrieve_num_T_steps = None
            else: retrieve_num_T_steps += 2

            self.P_data = h5_data['raw_data']['proteins'][:, 2:, 2:retrieve_num_T_steps]
            self.S_data = h5_data['raw_data']['sylinders'][:, 2:-1, 2:retrieve_num_T_steps]
            self.num_T_steps = len(h5_data['time'][2:retrieve_num_T_steps])     

    def retrieve_simulation_config_data(self):

        with open(self.simulation_config_file, 'r') as file:
            RunConfig_data = yaml.safe_load(file)
        self.box_side_length = np.float64(RunConfig_data['simBoxHigh'][0] * 2)
        self.syl_L = np.float64(RunConfig_data['sylinderLength'])

    def calculate_attachment_pos(self, p_attachment_coords, s_data, tol=1e-7, periodic_box=True):
        """
        Calculates the attachment position s \in [-syl_L/2, syl_L/2] for protein p and sylinder s.
        To deal with periodicity, we test the triangle inequality. 
        If it fails, then we simply run through each coordinate of p_attachment_coords to see which one(s)
        we need to add/subtract box_side_length to. If processing needs to be sped up, it's worth thinking
        about the more efficient way of doing this.
        """
        s_end_0_coords = s_data[:3]
        s_end_1_coords = s_data[3:]

        l0 = np.linalg.norm(s_end_0_coords - p_attachment_coords)
        l1 = np.linalg.norm(s_end_1_coords - p_attachment_coords)

        if periodic_box:
            if np.abs(l0 + l1 - self.syl_L) > tol:
                A = np.abs(s_end_0_coords - p_attachment_coords) > self.syl_L
                B = np.abs(s_end_1_coords - p_attachment_coords) > self.syl_L
                C = np.argwhere(A | B)
                for bad_arg in C: 
                    if p_attachment_coords[bad_arg] < 0:
                        p_attachment_coords[bad_arg] += self.box_side_length
                    else:
                        p_attachment_coords[bad_arg] -= self.box_side_length
                l0, l1 = np.linalg.norm(p_attachment_coords - s_end_0_coords), np.linalg.norm(p_attachment_coords - s_end_1_coords)
                if np.abs(l0 + l1 - self.syl_L) > tol:
                    raise ValueError(f"Periodicity calculation broke l0 + l1 - syl_L = {l0 + l1 - self.syl_L}: Syl numbers = {s_end_0}, {s_end_1}. P number = {p_num}.")
        
        s = l0 - self.syl_L / 2
        if np.abs(np.abs(s) - self.syl_L / 2 > tol):
            raise ValueError(f"Periodicity calculation broke. Maybe strange corner case. Syl numbers = {s_end_0}, {s_end_1}. P number = {p_num}.")

        return s

    def differentiate_motors_by_bound_state(self):
        """
        Returns dictionary of proportions of unbound, singly bound, doubly bound motors
        """
        motors_by_bound_state = {}
        
        num_Ps = self.P_data.shape[0]
        pdata = self.P_data[:, -2, :] * self.P_data[:, -1, :]
        motors_by_bound_state["unbound"] = sum(pdata==1) / num_Ps
        motors_by_bound_state["singly bound"] = sum((pdata < 0)) / num_Ps
        motors_by_bound_state["doubly bound"] = sum((pdata >= 0) & (pdata!=1)) / num_Ps

        return motors_by_bound_state
    
    def construct_empirical_motor_density(self):

        empirical_motor_density = []

        for t_step in range(self.num_T_steps):
            if (t_step + 1) % 10 == 0:
                print(f"{t_step + 1} / {self.num_T_steps}")

            empirical_motor_density.append(defaultdict(list))

            for p_data in self.P_data[:, :, t_step]:
                s_end_0, s_end_1 = p_data[-2:].astype(int)

                if s_end_0 > -1 and s_end_1 > -1:
                    si = self.calculate_attachment_pos(p_data[:3], self.S_data[s_end_0, :, t_step], tol=1e-5)
                    sj = self.calculate_attachment_pos(p_data[3:6], self.S_data[s_end_1, :, t_step], tol=1e-5)
                    
                    if s_end_0 < s_end_1:
                        empirical_motor_density[-1][f"{s_end_0},{s_end_1}"].append([si, sj])
                    else:
                        empirical_motor_density[-1][f"{s_end_1},{s_end_0}"].append([sj, si])

            empirical_motor_density[-1] = dict(empirical_motor_density[-1])
            
        self.empirical_motor_density = empirical_motor_density


class Empirical_Motor_Density_Smoother:

    def __init__(self, empirical_motor_density, simulation_config_file, std_param=0.04, smoothing_type="full", num_bins=16):

        self.load_empirical_motor_density(empirical_motor_density)
        self.simulation_config_file = simulation_config_file
        self.num_bins = num_bins
        self.num_T_steps = len(self.empirical_motor_density)

        # Declare self.box_side_length, self.syl_L (rod/sylinder length)
        self.retrieve_simulation_config_data()

        # Construct self.discrete_motor_density and corresponding self.x, self.y
        self.construct_binned_empirical_motor_density()

        # Smooth out self.discrete_motor_density
        self.smooth_discrete_motor_density(std_param=std_param, smoothing_type=smoothing_type)

    def load_empirical_motor_density(self, empirical_motor_density):

        if isinstance(empirical_motor_density, pathlib.PosixPath):
            with empirical_motor_density.open("r") as f:
                self.empirical_motor_density = json.load(f)
        else:
            self.empirical_motor_density = empirical_motor_density

    def construct_bin_edges(self):

        half_syl_L = self.syl_L / 2
        dummy_locs = np.zeros((1, 2))
        dbin = self.syl_L / (2 * self.num_bins)

        # Bin empirical data
        _, x_edges, y_edges = np.histogram2d(dummy_locs[:, 0], dummy_locs[:, 1], bins=self.num_bins, range=[[-half_syl_L, half_syl_L], [-half_syl_L, half_syl_L]])

        # Construct centre locations of bins with edges: xedges, yedges
        x, y = np.meshgrid(x_edges[:-1] + dbin, y_edges[:-1] + dbin, indexing="ij")

        return x_edges, y_edges, x, y

    def bin_empirical_motor_density_at_t_step(self, locs, x_edges, y_edges):
        """ 
        This bins the empirical locations si, sj.
        """

        z, x_edges, y_edges = np.histogram2d(locs[:, 0], locs[:, 1], bins=(x_edges, y_edges))
        z /= np.max(z)

        return z

    def retrieve_simulation_config_data(self):

        with open(self.simulation_config_file, 'r') as file:
            RunConfig_data = yaml.safe_load(file)
        self.box_side_length = np.float64(RunConfig_data['simBoxHigh'][0] * 2)
        self.syl_L = np.float64(RunConfig_data['sylinderLength'])

    def construct_binned_empirical_motor_density(self):

        discrete_motor_density = []
        x_edges, y_edges, x, y = self.construct_bin_edges()

        for t_step in range(self.num_T_steps):
            locs = np.vstack((list(self.empirical_motor_density[t_step].values())))

            z = self.bin_empirical_motor_density_at_t_step(locs, x_edges, y_edges)
            discrete_motor_density.append(z)

        self.discrete_motor_density, self.x, self.y = np.array(discrete_motor_density), x, y

    def smooth_discrete_motor_density(self, std_param=0.02, smoothing_type=None):
        """
        std_param is a % applied to length of data. See scipy.ndimage.gaussian_filter docs
        """ 

        if not smoothing_type:
            smoothed_motor_density = self.discrete_motor_density.copy()

        elif smoothing_type == "full":
            sigmas = np.array(self.discrete_motor_density.shape) * std_param #std is 2%
            smoothed_motor_density = scipy.ndimage.gaussian_filter(self.discrete_motor_density, sigma=sigmas, mode='reflect') #mode='reflect' is scipy default

        elif smoothing_type == "in time":
            smoothed_motor_density = np.zeros_like(self.discrete_motor_density)
            sigma = self.discrete_motor_density.shape[0] * std_param

            for i in range(self.num_bins):
                for j in range(self.num_bins):
                    smoothed_motor_density[:, i, j] = scipy.ndimage.gaussian_filter1d(self.discrete_motor_density[:, i, j], sigma=sigma, mode='reflect')

        smoothed_motor_density /= np.max(smoothed_motor_density, axis=(1, 2), keepdims=True)
        self.smoothed_motor_density = smoothed_motor_density
    
    def plot_smoothed_motor_density(self, plot_smoothed_motor_density_dir):

        for t_step, z in enumerate(self.smoothed_motor_density):

            if (t_step + 1) % 10 == 0:
                print(f"{t_step + 1} / {self.num_T_steps}")

            fig = plt.figure(figsize=(15, 15))
            ax = fig.add_subplot(projection='3d')
            ax.plot_surface(self.x, self.y, z, rstride=1, cstride=1, cmap="plasma", antialiased=True)
            ax.set_zlim(0, 1)
            plt.savefig(f"{plot_smoothed_motor_density_dir}/T{t_step}.png", dpi=50)
            plt.close()