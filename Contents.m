% Stochastic Integration Modeling Toolbox.
% Version 2.9 (29-Mar-2020)
% J. Ditterich, Center for Neuroscience, Univ. of California, Davis
%
% First passage time problems.
%   msprt_fb_3d_3b_sim            - Feedback implementation of MSPRT, 3 integrators, 3 boundaries, simulation.
%   ou_1d_1b_mar                  - 1D Ornstein-Uhlenbeck process, range limit, 1 boundary, Markov chain approximation.
%   ou_1d_2b_mar                  - 1D Ornstein-Uhlenbeck process, 2 boundaries, Markov chain approximation.
%   ou_1d_2b_num                  - 1D Ornstein-Uhlenbeck process, 2 boundaries, numerical solution.
%   ou_2d_2b_mar                  - 2D Ornstein-Uhlenbeck process, 2 boundaries, Markov chain approximation.
%   ou_2d_2b_sim                  - 2D Ornstein-Uhlenbeck process, 2 boundaries, simulation.
%   ou_2d_3b_mar                  - 2D Ornstein-Uhlenbeck process, 3 boundaries, Markov chain approximation.
%   ou_2d_3b_sim                  - 2D Ornstein-Uhlenbeck process, 3 boundaries, simulation.
%   ou_2d_3b_sim_sc               - 2D Ornstein-Uhlenbeck process, 3 boundaries, simulation,
%                                   also reports second choices (second guesses) based on the states of the integrators
%                                   at threshold crossing.
%   ou_2d_3b_sim_sc_add_time      - 2D Ornstein-Uhlenbeck process, 3 boundaries, simulation,
%                                   also reports second choices (second guesses) based on the states of the integrators
%                                   a fixed amount of time after the threshold crossing.
%   ou_2d_3b_two_cross_sim        - 2D Ornstein-Uhlenbeck prosess, 3 boundaries, simulation,
%                                   the process waits for a second threshold crossing, which is then
%                                   reported as the decision time.
%   ou_2d_3b_1d_2b_sim_sc         - 2D Ornstein-Uhlenbeck process, 3 boundaries, simulation;
%                                   once the first threshold crossing occurs, a new 1D OU process with two boundaries is
%                                   started to decide between the two remaining options as a second choice; the
%                                   decision time is given by the second threshold crossing.
%   ou_2d_3b_1d_fixed_time_sim_sc - 2D Ornstein-Uhlenbeck process, 3 boundaries, simulation;
%                                   once the first threshold crossing occurs, a new 1D OU process is started
%                                   to decide between the two remaining options as a second choice; the
%                                   decision is made after a fixed amount of time based on the sign of the
%                                   current state of the process.
%   ou_3d_3b_sim                  - 3D Ornstein-Uhlenbeck process, 3 boundaries, simulation.
%   ou_3d_inh_fb_3b_sim           - 3D Ornstein-Uhlenbeck process with inhibitory feedback, 3 boundaries, simulation.
%   ou_dropout_1d_2b_mar          - 1D Ornstein-Uhlenbeck process, dropout rate, 2 boundaries, Markov chain approximation.
%   ou_dropout_vd_1d_2b_mar       - 1D Ornstein-Uhlenbeck process, variable drift, dropout rate, 2 boundaries,
%                                   Markov chain approximation.
%   ou_timelim_1d_2b_mar          - 1D Ornstein-Uhlenbeck process, time limit, 2 boundaries, Markov chain approximation.
%   ou_vd_1d_2b_num               - 1D Ornstein-Uhlenbeck process, variable drift, 2 boundaries, numerical solution.
%   probsum_1d_2b_dis             - probability summation model, 1D normal process, 2 thresholds, discrete solution.
%   probsum_uncor_2d_2b_dis       - probability summation model, 2 normal processes, uncorrelated noise, discrete solution.
%   probsum_vs_1d_2b_dis          - probability summation model, 1D normal process, variable mean of the signal,
%                                   2 thresholds, discrete solution.
%   wiener_1d_1b_lin_ana          - 1D Wiener process, 1 linear boundary, analytical solution.
%   wiener_1d_2b_num              - 1D Wiener process, 2 boundaries, numerical solution.
%   wiener_vb_1d_2b_num           - 1D Wiener process, 2 variable boundaries, numerical solution.
%   wiener_vd_1d_2b_num           - 1D Wiener process, variable drift, 2 boundaries, numerical solution.
%   wiener_vi_1d_2b_num           - 1D Wiener process, variable initial value, 2 boundaries, numerical solution.
%
% Probability density functions (pdf).
%   normpdf2                      - Bivariate normal (Gaussian) density.
%
% Trajectories.
%   traj_aet_ou_1d_2b_del_sim     - 1D Ornstein-Uhlenbeck process, 2 boundaries, variable delays between
%                                   start of trial and integration onset and between first boundary crossing
%                                   and end of trial, simulation, aligned with end of trial.
%   traj_aet_ou_1d_2b_lim_del_sim - 1D Ornstein-Uhlenbeck process, range limit, 2 boundaries, variable delays
%                                   between start of trial and integration onset and between first boundary
%                                   crossing and end of trial, simulation, aligned with end of trial.
%   traj_aet_ou_2d_2b_del_sim     - 2D Ornstein-Uhlenbeck process, 2 boundaries, variable delays between
%                                   start of trial and integration onset and between first boundary crossing
%                                   and end of trial, simulation, aligned with end of trial.
%   traj_aet_ou_vd_1d_2b_del_sim  - 1D Ornstein-Uhlenbeck process, variable drift, 2 boundaries,
%                                   variable delays between start of trial and integration onset and between
%                                   first boundary crossing and end of trial, simulation,
%                                   aligned with end of trial.
%   traj_aet_ou_vd_2d_2b_del_sim  - 2D Ornstein-Uhlenbeck process, variable drift, 2 boundaries,
%                                   variable delays between start of trial and integration onset and between
%                                   first boundary crossing and end of trial, simulation,
%                                   aligned with end of trial.
%   traj_afp_msprt_fb_3d_3b_sim   - Feedback implementation of MSPRT, 3 integrators, 3 boundaries, simulation,
%                                   aligned with respect to first passage.
%   traj_afp_ou_1d_2b_sim         - 1D Ornstein-Uhlenbeck process, 2 boundaries, simulation,
%                                   aligned with respect to first passage.
%   traj_afp_ou_2d_2b_sim         - 2D Ornstein-Uhlenbeck process, 2 boundaries, simulation,
%                                   aligned with respect to first passage.
%   traj_afp_ou_2d_3b_sim         - 2D Ornstein-Uhlenbeck process, 3 boundaries, simulation,
%                                   aligned with respect to first passage.
%   traj_afp_ou_3d_inh_fb_3b_sim  - 3D Ornstein-Uhlenbeck process with (or also without) inhibitory feedback,
%                                   3 boundaries, simulation, aligned with respect to first passage.
%   traj_afp_ou_vd_1d_2b_sim      - 1D Ornstein-Uhlenbeck process, variable drift, 2 boundaries,
%                                   simulation, aligned with respect to first passage.
%   traj_msprt_fb_3d_3b_sim       - Feedback implementation of MSPRT, 3 integrators, 3 boundaries, simulation,
%                                   aligned with respect to t=0.
%   traj_ou_1d_2b_del_sim         - 1D Ornstein-Uhlenbeck process, 2 boundaries, variable delays between
%                                   start of trial and integration onset and between first boundary crossing
%                                   and end of trial, simulation, aligned with start of trial.
%   traj_ou_1d_2b_lim_del_sim     - 1D Ornstein-Uhlenbeck process, range limit, 2 boundaries, variable delays
%                                   between start of trial and integration onset and between first boundary
%                                   crossing and end of trial, simulation, aligned with start of trial.
%   traj_ou_2d_2b_del_sim         - 2D Ornstein-Uhlenbeck process, 2 boundaries, variable delays between
%                                   start of trial and integration onset and between first boundary crossing
%                                   and end of trial, simulation, aligned with start of trial.
%   traj_ou_1d_2b_sim             - 1D Ornstein-Uhlenbeck process, 2 boundaries, simulation,
%                                   aligned with respect to t=0.
%   traj_ou_2d_2b_sim             - 2D Ornstein-Uhlenbeck process, 2 boundaries, simulation,
%                                   aligned with respect to t=0.
%   traj_ou_2d_3b_sim             - 2D Ornstein-Uhlenbeck process, 3 boundaries, simulation,
%                                   aligned with respect to t=0.
%   traj_ou_3d_inh_fb_3b_sim      - 3D Ornstein-Uhlenbeck process with (or also without) inhibitory feedback,
%                                   3 boundaries, simulation, aligned with respect to t=0.
%   traj_ou_vd_1d_2b_del_sim      - 1D Ornstein-Uhlenbeck process, variable drift, 2 boundaries,
%                                   variable delays between start of trial and integration onset and between
%                                   first boundary crossing and end of trial, simulation,
%                                   aligned with start of trial.
%   traj_ou_vd_2d_2b_del_sim      - 2D Ornstein-Uhlenbeck process, variable drift, 2 boundaries,
%                                   variable delays between start of trial and integration onset and between
%                                   first boundary crossing and end of trial, simulation,
%                                   aligned with start of trial.
%   traj_ou_vd_1d_2b_sim          - 1D Ornstein-Uhlenbeck process, variable drift, 2 boundaries,
%                                   simulation, aligned with respect to t=0.
%
% Comments:
% The analytical solution is the fastest and most accurate solution.
% The numerical solution of the OU process doesn't accept 0 for the leakage of the integrator.
% This is why there is a separate function for the Wiener process. The temporal resolution
% doesn't seem to be very critical. I haven't seen any big systematic errors.
% The Markov chain approximation should return results, which are identical to the numerical
% solution, if the chosen temporal and spatial resolutions are high enough. Apparently,
% they also have to match. The distributions are too flat when a too high spatial resolution is chosen.
% This is why I have built a heuristic method for automatically choosing an optimal number of states
% into the functions.
