# -*- coding: utf-8 -*-

"""Top-level package for AMSOS analysis."""

__author__ = """Adam Reay Lamson"""
__email__ = 'alamson@flatironinstitute.org'
__version__ = '0.1.0'

from .chrom_analysis import (gauss_weighted_contact,
                             log_gauss_weighted_contact,
                             get_link_energy_arrays,
                             get_sep_hist,
                             get_sep_dist_mat,
                             get_overlap_arrs,
                             )
