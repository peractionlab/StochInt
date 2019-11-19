function [g_upper,g_lower,t_vec]=wiener_vi_1d_2b_num(wie_drift,wie_var,wie_init_mean,wie_init_sd,a_upper,a_lower,delta_t,stop_time,num_calls,wie_init_range)
% Numerical solution of the first passage time problem
% for a 1D Wiener process with 2 (time-variant) boundaries.
% The initial value varies from trial to trial. This function calls
% WIENER_1D_2B_NUM several times in order to get an approximate
% solution.
%
% J. Ditterich, 3/02
%
% [g_upper,g_lower,t_vec] = wiener_vi_1d_2b_num (wie_drift,wie_var,wie_init_mean,wie_init_sd,a_upper,a_lower,
%                                                delta_t,stop_time[,num_calls[,wie_init_range]])
%
% g_upper is the first passage time density for the upper boundary multiplied
%         by the probability of hitting the upper boundary first, evaluated
%         at the times given in t_vec.
% g_lower is the first passage time density for the lower boundary multiplied
%         by the probability of hitting the lower boundary first, evaluated
%         at the times given in t_vec.
%
% wie_drift is the drift of the Wiener process. It can either be a constant or,
%           for the time-variant case, the name of a function, which must return
%           the drift when called with the time as the argument.
% wie_var is the variance of the Wiener process.
% wie_init_mean is the mean initial value of the Wiener process.
% wie_init_sd is the standard deviation of the initial value.
% a_upper defines the upper absorbing boundary. a_upper can either be a constant or
%         the name of a function, which must return the location of the boundary
%         and its temporal derivative when called with the time as the argument.
% a_lower defines the lower absorbing boundary. See a_upper for the format.
% delta_t is the temporal resolution used by the algorithm.
% stop_time: The densities are calculated for times up to stop_time.
% num_calls is an optional argument, which defines, how often WIENER_1D_2B_NUM
%           is called. The default value is 9.
% wie_init_range is an optional argument, which defines the range used for the
%                calculation. The range is WIE_INIT_MEAN +/- WIE_INIT_RANGE *
%                WIE_INIT_SD. The default value is 2.5.

% History:
% released on 6/6/02 as part of toolbox V 2.0

if nargin<10 % wie_init_range not given?
    wie_init_range=2.5; % default value
end;

if nargin<9 % num_calls not given?
    num_calls=9; % default value
end;

num_calls=round(num_calls);

% Some checks
if wie_var<=0
    error('WIENER_VI_1D_2B_NUM: The variance must be a positive number!');
end;

if wie_init_sd<=0
    error('WIENER_VI_1D_2B_NUM: WIE_INIT_SD must be a positive number!');
end;

if delta_t<=0
    error('WIENER_VI_1D_2B_NUM: The time step must be a positive number!');
end;

if stop_time<=0
    error('WIENER_VI_1D_2B_NUM: STOP_TIME must be a positive number!');
end;

if num_calls<2
    error('WIENER_VI_1D_2B_NUM: NUM_CALLS must be at least 2!');
end;

if wie_init_range<=0
    error('WIENER_VI_1D_2B_NUM: WIE_INIT_RANGE must be a positive number!');
end;

% Initialization
temp=[wie_init_mean-wie_init_range*wie_init_sd+wie_init_range*wie_init_sd/(num_calls-1): ...
        2*wie_init_range*wie_init_sd/(num_calls-1): ...
        wie_init_mean+wie_init_range*wie_init_sd-wie_init_range*wie_init_sd/(num_calls-1)];
temp2=[0 normcdf(temp,wie_init_mean,wie_init_sd) 1];
weights=diff(temp2); % These are the weights, which have to be applied to the results of the single calls.
inits=[wie_init_mean-wie_init_range*wie_init_sd:2*wie_init_range*wie_init_sd/(num_calls-1): ...
        wie_init_mean+wie_init_range*wie_init_sd]; % These are the initial values, which have to be used in the single calls.

% Loop
for i=1:num_calls
    [g_upper_temp,g_lower_temp,t_temp]=wiener_1d_2b_num(wie_drift,wie_var,inits(i),a_upper,a_lower,delta_t,1,stop_time);
    
    if i==1 % first call?
        g_upper=g_upper_temp*weights(i);
        g_lower=g_lower_temp*weights(i);
        t_vec=t_temp;
    else
        g_upper=g_upper+g_upper_temp*weights(i);
        g_lower=g_lower+g_lower_temp*weights(i);
    end;
end;
