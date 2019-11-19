function [g_upper,g_lower,t_vec]=wiener_vd_1d_2b_num(wie_drift_mean,wie_drift_sd,wie_var,wie_init,a_upper,a_lower,delta_t,stop_time,num_calls,wie_drift_range)
% Numerical solution of the first passage time problem
% for a 1D Wiener process with 2 (time-variant) boundaries.
% The drift varies from trial to trial. This function calls
% WIENER_1D_2B_NUM several times in order to get an approximate
% solution.
%
% J. Ditterich, 10/01
%
% [g_upper,g_lower,t_vec] = wiener_vd_1d_2b_num (wie_drift_mean,wie_drift_sd,wie_var,wie_init,a_upper,a_lower,
%                                                delta_t,stop_time[,num_calls[,wie_drift_range]])
%
% g_upper is the first passage time density for the upper boundary multiplied
%         by the probability of hitting the upper boundary first, evaluated
%         at the times given in t_vec.
% g_lower is the first passage time density for the lower boundary multiplied
%         by the probability of hitting the lower boundary first, evaluated
%         at the times given in t_vec.
%
% wie_drift_mean is the mean drift of the Wiener process.
% wie_drift_sd is the (between trial) standard deviation of the drift. The drift is constant
%              during a trial.
% wie_var is the variance of the Wiener process.
% wie_init is the initial value of the Wiener process.
% a_upper defines the upper absorbing boundary. a_upper can either be a constant or
%         the name of a function, which must return the location of the boundary
%         and its temporal derivative when called with the time as the argument.
%         In case an additional argument needs to be passed to the function,
%         a_upper can also be a cell array with two elements. In this case the first cell has
%         to be the name of the function, the second cell is passed as
%         the second function argument (after the time).
% a_lower defines the lower absorbing boundary. See a_upper for the format.
% delta_t is the temporal resolution used by the algorithm.
% stop_time: The densities are calculated for times up to stop_time.
% num_calls is an optional argument, which defines, how often WIENER_1D_2B_NUM
%           is called. The default value is 9.
% wie_drift_range is an optional argument, which defines the range used for the
%                 calculation. The range is WIE_DRIFT_MEAN +/- WIE_DRIFT_RANGE *
%                 WIE_DRIFT_SD. The default value is 2.5.

% History:
% released on 8/23/01 as part of toolbox V 0.5 Beta
% bug regarding input arguments fixed on 10/12/01
% new argument checks added on 10/22/01
% a_upper and a_lower adjusted for Parallel Computing Toolbox on 12/12/08

if nargin<10 % wie_drift_range not given?
    wie_drift_range=2.5; % default value
end;

if nargin<9 % num_calls not given?
    num_calls=9; % default value
end;

num_calls=round(num_calls);

% Some checks
if ~isnumeric(wie_drift_mean)
    error('WIENER_VD_1D_2B_NUM: WIE_DRIFT_MEAN must be a number!');
end;

if wie_drift_sd<=0
    error('WIENER_VD_1D_2B_NUM: WIE_DRIFT_SD must be a positive number!');
end;

if wie_var<=0
    error('WIENER_VD_1D_2B_NUM: The variance must be a positive number!');
end;

if delta_t<=0
    error('WIENER_VD_1D_2B_NUM: The time step must be a positive number!');
end;

if stop_time<=0
    error('WIENER_VD_1D_2B_NUM: STOP_TIME must be a positive number!');
end;

if num_calls<2
    error('WIENER_VD_1D_2B_NUM: NUM_CALLS must be at least 2!');
end;

if wie_drift_range<=0
    error('WIENER_VD_1D_2B_NUM: WIE_DRIFT_RANGE must be a positive number!');
end;

% Initialization
temp=[wie_drift_mean-wie_drift_range*wie_drift_sd+wie_drift_range*wie_drift_sd/(num_calls-1): ...
        2*wie_drift_range*wie_drift_sd/(num_calls-1): ...
        wie_drift_mean+wie_drift_range*wie_drift_sd-wie_drift_range*wie_drift_sd/(num_calls-1)];
temp2=[0 normcdf(temp,wie_drift_mean,wie_drift_sd) 1];
weights=diff(temp2); % These are the weights, which have to be applied to the results of the single calls.
drifts=[wie_drift_mean-wie_drift_range*wie_drift_sd:2*wie_drift_range*wie_drift_sd/(num_calls-1): ...
        wie_drift_mean+wie_drift_range*wie_drift_sd]; % These are the drifts, which have to be used in the single calls.

% Loop
for i=1:num_calls
    [g_upper_temp,g_lower_temp,t_temp]=wiener_1d_2b_num(drifts(i),wie_var,wie_init,a_upper,a_lower,delta_t,1,stop_time);
    
    if i==1 % first call?
        g_upper=g_upper_temp*weights(i);
        g_lower=g_lower_temp*weights(i);
        t_vec=t_temp;
    else
        g_upper=g_upper+g_upper_temp*weights(i);
        g_lower=g_lower+g_lower_temp*weights(i);
    end;
end;
