function [g_upper,g_lower,t_vec]=wiener_vb_1d_2b_num(wie_drift,wie_var,wie_init,a_upper_mean,a_lower_mean,bound_sd,delta_t,stop_time,num_calls,wie_bound_range)
% Numerical solution of the first passage time problem
% for a 1D Wiener process with 2 (fixed) boundaries.
% The boundaries vary from trial to trial, but this variability
% is fully anti-correlated which means that both boundaries move
% either closer together or away from each other. This function calls
% WIENER_1D_2B_NUM several times in order to get an approximate
% solution.
%
% J. Ditterich, 3/02
%
% [g_upper,g_lower,t_vec] = wiener_vb_1d_2b_num (wie_drift,wie_var,wie_init,a_upper_mean,a_lower_mean,
%                                                bound_sd,delta_t,stop_time[,num_calls[,bound_range]])
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
% wie_init is the initial value of the Wiener process.
% a_upper_mean defines the mean upper absorbing boundary. a_upper_mean must be
%              a constant.
% a_lower_mean defines the mean lower absorbing boundary. a_lower_mean must be
%              a constant.
% bound_sd is the (between trial) standard deviation of the boundaries. The
%          boundaries are constant during a trial. The variations in the boundaries
%          are fully anti-correlated which means that both boundaries move either
%          closer together or away from each other.
% delta_t is the temporal resolution used by the algorithm.
% stop_time: The densities are calculated for times up to stop_time.
% num_calls is an optional argument, which defines, how often WIENER_1D_2B_NUM
%           is called. The default value is 9.
% wie_bound_range is an optional argument, which defines the ranges used for the
%                 calculation. The ranges are A_UPPER_MEAN +/- BOUND_RANGE *
%                 BOUND_SD and A_LOWER_MEAN +/- BOUND_RANGE * BOUND_SD. The default
%                 value is 2.5.

% History:
% released on 6/6/02 as part of toolbox V 2.0

if nargin<10 % wie_init_range not given?
    bound_range=2.5; % default value
end;

if nargin<9 % num_calls not given?
    num_calls=9; % default value
end;

num_calls=round(num_calls);

% Some checks
if wie_var<=0
    error('WIENER_VB_1D_2B_NUM: The variance must be a positive number!');
end;

if ~isnumeric(a_upper_mean)
    error('WIENER_VD_1D_2B_NUM: A_UPPER_MEAN must be a number!');
end;

if ~isnumeric(a_lower_mean)
    error('WIENER_VD_1D_2B_NUM: A_LOWER_MEAN must be a number!');
end;

if bound_sd<=0
    error('WIENER_VB_1D_2B_NUM: BOUND_SD must be a positive number!');
end;

if delta_t<=0
    error('WIENER_VB_1D_2B_NUM: The time step must be a positive number!');
end;

if stop_time<=0
    error('WIENER_VB_1D_2B_NUM: STOP_TIME must be a positive number!');
end;

if num_calls<2
    error('WIENER_VB_1D_2B_NUM: NUM_CALLS must be at least 2!');
end;

if bound_range<=0
    error('WIENER_VB_1D_2B_NUM: BOUND_RANGE must be a positive number!');
end;

% Initialization
temp=[-bound_range*bound_sd+bound_range*bound_sd/(num_calls-1): ...
        2*bound_range*bound_sd/(num_calls-1): ...
        bound_range*bound_sd-bound_range*bound_sd/(num_calls-1)];
temp2=[0 normcdf(temp,0,bound_sd) 1];
weights=diff(temp2); % These are the weights, which have to be applied to the results of the single calls.
deltas=[-bound_range*bound_sd:2*bound_range*bound_sd/(num_calls-1): ...
        bound_range*bound_sd]; % These are the deviations, which have to be used in the single calls.
uppers=a_upper_mean+deltas; % These are the upper bounds, which have to be used in the single calls.
lowers=a_lower_mean-deltas; % These are the lower bounds, which have to be used in the single calls.
                            % The negative sign determines the anti-correlation!

% Loop
for i=1:num_calls
    [g_upper_temp,g_lower_temp,t_temp]=wiener_1d_2b_num(wie_drift,wie_var,wie_init,uppers(i),lowers(i),delta_t,1,stop_time);
    
    if i==1 % first call?
        g_upper=g_upper_temp*weights(i);
        g_lower=g_lower_temp*weights(i);
        t_vec=t_temp;
    else
        g_upper=g_upper+g_upper_temp*weights(i);
        g_lower=g_lower+g_lower_temp*weights(i);
    end;
end;
