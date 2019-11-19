function [g_upper,g_lower,t_vec]=ou_vd_1d_2b_num(ou_drift_mean,ou_drift_sd,ou_leak,ou_var,ou_init,a_upper,a_lower,delta_t,stop_time,num_calls,ou_drift_range)
% Numerical solution of the first passage time problem
% for a 1D Ornstein-Uhlenbeck process with 2 (time-variant) boundaries.
% The drift varies from trial to trial. This function calls
% OU_1D_2B_NUM several times in order to get an approximate
% solution.
%
% J. Ditterich, 8/01
%
% [g_upper,g_lower,t_vec] = ou_vd_1d_2b_num (ou_drift_mean,ou_drift_sd,ou_leak,ou_var,ou_init,a_upper,a_lower,
%                                            delta_t,stop_time[,num_calls[,ou_drift_range]])
%
% g_upper is the first passage time density for the upper boundary multiplied
%         by the probability of hitting the upper boundary first, evaluated
%         at the times given in t_vec.
% g_lower is the first passage time density for the lower boundary multiplied
%         by the probability of hitting the lower boundary first, evaluated
%         at the times given in t_vec.
%
% ou_drift_mean is the mean drift of the OU process.
% ou_drift_sd is the (between trial) standard deviation of the drift. The drift is constant
%              during a trial.
% ou_leak defines the "leakiness" of the integrator. The deterministic part
%         of the stochastic differential equation is given by
%         ou_drift - ou_leak * current_value. ou_leak must not be 0! Please use
%         WIENER_VD_1D_2B_NUM in this case.
% ou_var is the variance of the OU process.
% ou_init is the initial value of the OU process.
% a_upper defines the upper absorbing boundary. a_upper can either be a constant or
%         the name of a function, which must return the location of the boundary
%         and its temporal derivative when called with the time as the argument.
% a_lower defines the lower absorbing boundary. See a_upper for the format.
% delta_t is the temporal resolution used by the algorithm.
% stop_time: The densities are calculated for times up to stop_time.
% num_calls is an optional argument, which defines, how often OU_1D_2B_NUM
%           is called. The default value is 9.
% ou_drift_range is an optional argument, which defines the range used for the
%                 calculation. The range is OU_DRIFT_MEAN +/- OU_DRIFT_RANGE *
%                 OU_DRIFT_SD. The default value is 2.5.

% History:
% released on 8/23/01 as part of toolbox V 0.5 Beta
% bug in argument checks fixed on 10/16/01
% new argument checks added on 10/22/01
% documentation completed on 04/08/02

if nargin<11 % ou_drift_range not given?
    ou_drift_range=2.5; % default value
end;

if nargin<10 % num_calls not given?
    num_calls=9; % default value
end;

num_calls=round(num_calls);

% Some checks
if ~isnumeric(ou_drift_mean)
    error('OU_VD_1D_2B_NUM: OU_DRIFT_MEAN must be a number!');
end;

if ou_drift_sd<=0
    error('OU_VD_1D_2B_NUM: OU_DRIFT_SD must be a positive number!');
end;

if ou_leak<0
    error('OU_VD_1D_2B_NUM: OU_LEAK must not be negative!');
end;

if ou_leak==0
    error('OU_VD_1B_2D_NUM: OU_LEAK must not be 0! Please use WIENER_VD_1B_2D_NUM instead.');
end;

if ou_var<=0
    error('OU_VD_1D_2B_NUM: The variance must be a positive number!');
end;

if delta_t<=0
    error('OU_VD_1D_2B_NUM: The time step must be a positive number!');
end;

if stop_time<=0
    error('OU_VD_1D_2B_NUM: STOP_TIME must be a positive number!');
end;

if num_calls<2
    error('OU_VD_1D_2B_NUM: NUM_CALLS must be at least 2!');
end;

if ou_drift_range<=0
    error('OU_VD_1D_2B_NUM: OU_DRIFT_RANGE must be a positive number!');
end;

% Initialization
temp=[ou_drift_mean-ou_drift_range*ou_drift_sd+ou_drift_range*ou_drift_sd/(num_calls-1): ...
        2*ou_drift_range*ou_drift_sd/(num_calls-1): ...
        ou_drift_mean+ou_drift_range*ou_drift_sd-ou_drift_range*ou_drift_sd/(num_calls-1)];
temp2=[0 normcdf(temp,ou_drift_mean,ou_drift_sd) 1];
weights=diff(temp2); % These are the weights, which have to be applied to the results of the single calls.
drifts=[ou_drift_mean-ou_drift_range*ou_drift_sd:2*ou_drift_range*ou_drift_sd/(num_calls-1): ...
        ou_drift_mean+ou_drift_range*ou_drift_sd]; % These are the drifts, which have to be used in the single calls.

% Loop
for i=1:num_calls
    [g_upper_temp,g_lower_temp,t_temp]=ou_1d_2b_num(drifts(i),ou_leak,ou_var,ou_init,a_upper,a_lower,delta_t,1,stop_time);
    
    if i==1 % first call?
        g_upper=g_upper_temp*weights(i);
        g_lower=g_lower_temp*weights(i);
        t_vec=t_temp;
    else
        g_upper=g_upper+g_upper_temp*weights(i);
        g_lower=g_lower+g_lower_temp*weights(i);
    end;
end;
