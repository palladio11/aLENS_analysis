# -*- coding: utf-8 -*-

"""Top-level package for aLENS analysis."""

__author__ = """Adam Reay Lamson"""
__email__ = 'alamson@flatironinstitute.org'
__version__ = '0.1.0'


from .physical_scales import (get_drag_coeff,
                              get_char_time,
                              get_poly_diffuse_dist,
                              get_rouse_time,
                              get_link_relax_time,
                              get_pfract
                              )

from .colormaps import register_cmaps

from .controller_funcs import (TYPE_FUNC_DICT,
                               seed_analysis)

from .rouse_mode_analysis import (get_rouse_modes_at_t,
                                  get_rouse_modes,
                                  get_rouse_mode_corr,
                                  get_rouse_mode_corr_fast,)

from . import chromatin
from .scripts import Util
from . import hydro
from . import helpers
from . import nematic_order
#######################################################################
#                        Soon to be depricated                        #
#######################################################################


from .chromatin.chrom_condensate_analysis import (Condensate,
                                                  gen_condensate_track_info,
                                                  get_max_and_total_cond_size,
                                                  extract_condensates,
                                                  )

from .chromatin.chrom_analysis import (gauss_weighted_contact,
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
                                       get_contact_mat_analysis,
                                       get_end_end_distance,
                                       calc_rad_of_gyration,
                                       find_neighbors,
                                       get_pos_kymo_data,
                                       get_pos_cond_data,
                                       smooth_kymo_mat,
                                       get_contact_cond_data,
                                       )

from .chromatin.chrom_seed_scan_analysis import(get_scan_cond_data,
                                                get_scan_avg_contact_mat,
                                                get_scan_avg_kymo)


from .chromatin.chrom_graph_funcs import (make_total_distr_plots,
                                          make_min_distr_plots,
                                          make_segment_distr_graphs,
                                          make_rog_vs_time_graph,
                                          make_hic_plot,
                                          make_summed_contact_kymo_graph,
                                          plot_link_energy_vs_time,
                                          plot_rog_vs_time_graph,
                                          plot_condensate_kymo,
                                          plot_contact_kymo,
                                          plot_pos_kymo,
                                          plot_condensate_characterize,
                                          plot_condensate_tracks,
                                          plot_condensate_size_vs_time,
                                          plot_condensate_avg_contact_vs_time,
                                          )

from .chromatin.chrom_seed_scan_graph_funcs import (plot_condensate_num_sd_scan,
                                                    plot_condensate_size_sd_scan,
                                                    sd_num)
