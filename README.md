# StochInt
Stochastic Integration Modeling Toolbox

A MATLAB toolbox for modeling decision processes based on stochastic integration. In addition to standard drift-diffusion models it also supports leaky integration, 2- and 3-dimensional stochastic processes, second guesses, etc.

Some additional information (e.g., publications in which the toolbox has been used) can be found here: https://www.peractionlab.org/software?id=5

The individual functions provide more information about what they do and their input and output arguments.

Content of Contents.m:

`% Stochastic Integration Modeling Toolbox.`<br/>
`% Version 2.9 (29-Mar-2020)`<br/>
`% J. Ditterich, Center for Neuroscience, Univ. of California, Davis`<br/>
`%`<br/>
`% First passage time problems.`<br/>
`%   msprt_fb_3d_3b_sim            - Feedback implementation of MSPRT, 3 integrators, 3 boundaries, simulation.`<br/>
`%   ou_1d_1b_mar                  - 1D Ornstein-Uhlenbeck process, range limit, 1 boundary, Markov chain approximation.`<br/>
`%   ou_1d_2b_mar                  - 1D Ornstein-Uhlenbeck process, 2 boundaries, Markov chain approximation.`<br/>
`%   ou_1d_2b_num                  - 1D Ornstein-Uhlenbeck process, 2 boundaries, numerical solution.`<br/>
`%   ou_2d_2b_mar                  - 2D Ornstein-Uhlenbeck process, 2 boundaries, Markov chain approximation.`<br/>
`%   ou_2d_2b_sim                  - 2D Ornstein-Uhlenbeck process, 2 boundaries, simulation.`<br/>
`%   ou_2d_3b_mar                  - 2D Ornstein-Uhlenbeck process, 3 boundaries, Markov chain approximation.`<br/>
`%   ou_2d_3b_sim                  - 2D Ornstein-Uhlenbeck process, 3 boundaries, simulation.`<br/>
`%   ou_2d_3b_sim_sc               - 2D Ornstein-Uhlenbeck process, 3 boundaries, simulation, also reports second choices (second guesses).`<br/>
`%   ou_2d_3b_sim_sc_add_time      - 2D Ornstein-Uhlenbeck process, 3 boundaries, simulation, also reports second choices (second guesses) based on the states of the integrators a fixed amount of time after the threshold crossing.`<br/>
`%   ou_2d_3b_two_cross_sim        - 2D Ornstein-Uhlenbeck prosess, 3 boundaries, simulation, the process waits for a second threshold crossing, which is then reported as the decision time.`<br/>
`%   ou_2d_3b_1d_2b_sim_sc         - 2D Ornstein-Uhlenbeck process, 3 boundaries, simulation; once the first threshold crossing occurs, a new 1D OU process with two boundaries is started to decide between the two remaining options as a second choice; the decision time is given by the second threshold crossing.`<br/>
`%   ou_2d_3b_1d_fixed_time_sim_sc - 2D Ornstein-Uhlenbeck process, 3 boundaries, simulation; once the first threshold crossing occurs, a new 1D OU process is started to decide between the two remaining options as a second choice; the decision is made after a fixed amount of time based on the sign of the current state of the process.`<br/>
`%   ou_3d_3b_sim                  - 3D Ornstein-Uhlenbeck process, 3 boundaries, simulation.`<br/>
`%   ou_3d_inh_fb_3b_sim           - 3D Ornstein-Uhlenbeck process with inhibitory feedback, 3 boundaries, simulation.`<br/>
`%   ou_dropout_1d_2b_mar          - 1D Ornstein-Uhlenbeck process, dropout rate, 2 boundaries, Markov chain approximation.`<br/>
`%   ou_dropout_vd_1d_2b_mar       - 1D Ornstein-Uhlenbeck process, variable drift, dropout rate, 2 boundaries, Markov chain approximation.`<br/>
`%   ou_timelim_1d_2b_mar          - 1D Ornstein-Uhlenbeck process, time limit, 2 boundaries, Markov chain approximation.`<br/>
`%   ou_vd_1d_2b_num               - 1D Ornstein-Uhlenbeck process, variable drift, 2 boundaries, numerical solution.`<br/>
`%   probsum_1d_2b_dis             - probability summation model, 1D normal process, 2 thresholds, discrete solution.`<br/>
`%   probsum_uncor_2d_2b_dis       - probability summation model, 2 normal processes, uncorrelated noise, discrete solution.`<br/>
`%   probsum_vs_1d_2b_dis          - probability summation model, 1D normal process, variable mean of the signal, 2 thresholds, discrete solution.`<br/>
`%   wiener_1d_1b_lin_ana          - 1D Wiener process, 1 linear boundary, analytical solution.`<br/>
`%   wiener_1d_2b_num              - 1D Wiener process, 2 boundaries, numerical solution.`<br/>
`%   wiener_vb_1d_2b_num           - 1D Wiener process, 2 variable boundaries, numerical solution.`<br/>
`%   wiener_vd_1d_2b_num           - 1D Wiener process, variable drift, 2 boundaries, numerical solution.`<br/>
`%   wiener_vi_1d_2b_num           - 1D Wiener process, variable initial value, 2 boundaries, numerical solution.`<br/>
`%`<br/>
`% Probability density functions (pdf).`<br/>
`%   normpdf2                      - Bivariate normal (Gaussian) density.`<br/>
`%`<br/>
`% Trajectories.`<br/>
`%   traj_aet_ou_1d_2b_del_sim     - 1D Ornstein-Uhlenbeck process, 2 boundaries, variable delays between start of trial and integration onset and between first boundary crossing and end of trial, simulation, aligned with end of trial.`<br/>
`%   traj_aet_ou_1d_2b_lim_del_sim - 1D Ornstein-Uhlenbeck process, range limit, 2 boundaries, variable delays between start of trial and integration onset and between first boundary crossing and end of trial, simulation, aligned with end of trial.`<br/>
`%   traj_aet_ou_2d_2b_del_sim     - 2D Ornstein-Uhlenbeck process, 2 boundaries, variable delays between start of trial and integration onset and between first boundary crossing and end of trial, simulation, aligned with end of trial.`<br/>
`%   traj_aet_ou_vd_1d_2b_del_sim  - 1D Ornstein-Uhlenbeck process, variable drift, 2 boundaries, variable delays between start of trial and integration onset and between first boundary crossing and end of trial, simulation, aligned with end of trial.`<br/>
`%   traj_aet_ou_vd_2d_2b_del_sim  - 2D Ornstein-Uhlenbeck process, variable drift, 2 boundaries, variable delays between start of trial and integration onset and between first boundary crossing and end of trial, simulation, aligned with end of trial.`<br/>
`%   traj_afp_msprt_fb_3d_3b_sim   - Feedback implementation of MSPRT, 3 integrators, 3 boundaries, simulation, aligned with respect to first passage.`<br/>
`%   traj_afp_ou_1d_2b_sim         - 1D Ornstein-Uhlenbeck process, 2 boundaries, simulation, aligned with respect to first passage.`<br/>
`%   traj_afp_ou_2d_2b_sim         - 2D Ornstein-Uhlenbeck process, 2 boundaries, simulation, aligned with respect to first passage.`<br/>
`%   traj_afp_ou_2d_3b_sim         - 2D Ornstein-Uhlenbeck process, 3 boundaries, simulation, aligned with respect to first passage.`<br/>
`%   traj_afp_ou_3d_inh_fb_3b_sim  - 3D Ornstein-Uhlenbeck process with (or also without) inhibitory feedback, 3 boundaries, simulation, aligned with respect to first passage.`<br/>
`%   traj_afp_ou_vd_1d_2b_sim      - 1D Ornstein-Uhlenbeck process, variable drift, 2 boundaries, simulation, aligned with respect to first passage.`<br/>
`%   traj_msprt_fb_3d_3b_sim       - Feedback implementation of MSPRT, 3 integrators, 3 boundaries, simulation, aligned with respect to t=0.`<br/>
`%   traj_ou_1d_2b_del_sim         - 1D Ornstein-Uhlenbeck process, 2 boundaries, variable delays between start of trial and integration onset and between first boundary crossing and end of trial, simulation, aligned with start of trial.`<br/>
`%   traj_ou_1d_2b_lim_del_sim     - 1D Ornstein-Uhlenbeck process, range limit, 2 boundaries, variable delays between start of trial and integration onset and between first boundary crossing and end of trial, simulation, aligned with start of trial.`<br/>
`%   traj_ou_2d_2b_del_sim         - 2D Ornstein-Uhlenbeck process, 2 boundaries, variable delays between start of trial and integration onset and between first boundary crossing and end of trial, simulation, aligned with start of trial.`<br/>
`%   traj_ou_1d_2b_sim             - 1D Ornstein-Uhlenbeck process, 2 boundaries, simulation, aligned with respect to t=0.`<br/>
`%   traj_ou_2d_2b_sim             - 2D Ornstein-Uhlenbeck process, 2 boundaries, simulation, aligned with respect to t=0.`<br/>
`%   traj_ou_2d_3b_sim             - 2D Ornstein-Uhlenbeck process, 3 boundaries, simulation, aligned with respect to t=0.`<br/>
`%   traj_ou_3d_inh_fb_3b_sim      - 3D Ornstein-Uhlenbeck process with (or also without) inhibitory feedback, 3 boundaries, simulation, aligned with respect to t=0.`<br/>
`%   traj_ou_vd_1d_2b_del_sim      - 1D Ornstein-Uhlenbeck process, variable drift, 2 boundaries, variable delays between start of trial and integration onset and between first boundary crossing and end of trial, simulation, aligned with start of trial.`<br/>
`%   traj_ou_vd_2d_2b_del_sim      - 2D Ornstein-Uhlenbeck process, variable drift, 2 boundaries, variable delays between start of trial and integration onset and between first boundary crossing and end of trial, simulation, aligned with start of trial.`<br/>
`%   traj_ou_vd_1d_2b_sim          - 1D Ornstein-Uhlenbeck process, variable drift, 2 boundaries, simulation, aligned with respect to t=0.`<br/>
`%`<br/>
`% Comments:`<br/>
`% The analytical solution is the fastest and most accurate solution.`<br/>
`% The numerical solution of the OU process doesn't accept 0 for the leakage of the integrator.`<br/>
`% This is why there is a separate function for the Wiener process. The temporal resolution`<br/>
`% doesn't seem to be very critical. I haven't seen any big systematic errors.`<br/>
`% The Markov chain approximation should return results, which are identical to the numerical`<br/>
`% solution, if the chosen temporal and spatial resolutions are high enough. Apparently,`<br/>
`% they also have to match. The distributions are too flat when a too high spatial resolution is chosen.`<br/>
`% This is why I have built a heuristic method for automatically choosing an optimal number of states`<br/>
`% into the functions.`<br/>