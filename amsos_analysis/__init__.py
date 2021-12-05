# -*- coding: utf-8 -*-

"""Top-level package for AMSOS analysis."""

__author__ = """Adam Reay Lamson"""
__email__ = 'alamson@flatironinstitute.org'
__version__ = '0.1.0'

from .chrom_analysis import (gauss_weighted_contact,
                             log_gauss_weighted_contact,
                             get_link_energy_arrays,
                             get_link_tension,
                             get_sep_hist,
                             get_sep_dist_mat,
                             get_overlap_arrs,
                             autocorr_bead_pos,
                             distr_hists,
                             total_distr_hists,
                             cart_distr_hists,
                             cylin_distr_hists,
                             rad_distr_hists,
                             get_all_rog_stats,
                             get_time_avg_contact_mat,
                             get_end_end_distance,
                             calc_rad_of_gyration,
                             find_neighbors
                             )

from .rouse_mode_analysis import (get_rouse_modes_at_t,
                                  get_rouse_modes,
                                  get_rouse_mode_corr,
                                  get_rouse_mode_corr_fast,)

from .chrom_graph_funcs import (graph_link_energy_vs_time,
                                make_total_distr_plots,
                                make_min_distr_plots,
                                make_segment_distr_graphs,
                                make_rog_vs_time_graph,
                                plot_rog_vs_time_graph,
                                make_hic_plot,
                                make_summed_contact_kymo_graph
                                )

from .physical_scales import (get_drag_coeff,
                              get_char_time,
                              get_poly_diffuse_dist,
                              get_rouse_time,
                              get_link_relax_time,
                              get_pfract
                              )

from .colormaps import register_cmaps
