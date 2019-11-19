function [g_upper,g_lower,t_vec,g_dropout]=ou_dropout_vd_1d_2b_mar(ou_drift_mean,ou_drift_sd,ou_leak,ou_var,ou_init,a_upper,a_lower,delta_t,covered_range,num_tr_states,stop_time,dropout_rate,num_calls,ou_drift_range)
% Markov chain approximation of the first passage time problem
% for a 1D Ornstein-Uhlenbeck process with 2 (time-variant) boundaries.
% The drift varies from trial to trial. There is a risk of aborting
% the process before hitting a boundary given by a dropout rate.
% This function calls OU_DROPOUT_1D_2B_MAR several times in order to
% get an approximate solution.
%
% J. Ditterich, 1/03
%
% [g_upper,g_lower,t_vec,g_dropout] = ou_dropout_vd_1d_2b_mar (ou_drift_mean,ou_drift_sd,ou_leak,ou_var,ou_init,a_upper,a_lower,
%                                                              delta_t,covered_range,num_tr_states,stop_time,
%                                                              dropout_rate[,num_calls[,ou_drift_range]])
%
% g_upper is the first passage time density for the upper boundary multiplied
%         by the probability of hitting the upper boundary first, evaluated
%         at the times given in t_vec.
% g_lower is the first passage time density for the lower boundary multiplied
%         by the probability of hitting the lower boundary first, evaluated
%         at the times given in t_vec.
% g_dropout is the PDF for dropout times multiplied by the probability of dropping
%           out before hitting a boundary, evaluated at the times given in t_vec.
%
% ou_drift_mean is the mean drift of the OU process.
% ou_drift_sd is the (between trial) standard deviation of the drift. The drift is constant
%              during a trial.
% ou_leak defines the "leakiness" of the integrator. The deterministic part
%         of the stochastic differential equation is given by
%         ou_drift - ou_leak * current_value. A Wiener process can be studied by setting
%         ou_leak to 0.
% ou_var is the variance of the OU process.
% ou_init is the initial value of the OU process.
% a_upper defines the upper absorbing boundary. a_upper can either be a constant or
%         the name of a function, which must return the location of the boundary
%         when called with the time as the argument.
% a_lower defines the lower absorbing boundary. See a_upper for the format.
% delta_t is the temporal resolution used by the algorithm.
% covered_range defines the range, which should be covered by the transient states.
%               It has to be a 2 element vector. Make sure that you never define boundaries
%               outside the specified range! (In this case the algorithm doesn't solve the
%               problem you are interested in.) In the case of constant boundaries the best
%               results are obtained when the boundaries are identical to
%               covered_range(1) and covered_range(2).
% num_tr_states defines the number of transient states used for the discrete approximation
%               of the problem. When using a symmetrical range you should specify
%               an odd number of states for an unbiased representation of 0.
%               When passing 0 a built-in mechanism is used for choosing an optimal
%               number of states (based on a heuristic method).
% stop_time: The densities are calculated for times up to stop_time.
% dropout_rate is the rate defining the risk of dropping out of the process
%              at any time. It can be either a constant or, for the time-variant
%              case, the name of a function which must return the rate when
%              called with the time as the argument.
% num_calls is an optional argument, which defines, how often OU_1D_2B_NUM
%           is called. The default value is 9.
% ou_drift_range is an optional argument, which defines the range used for the
%                 calculation. The range is OU_DRIFT_MEAN +/- OU_DRIFT_RANGE *
%                 OU_DRIFT_SD. The default value is 2.5.

% History:
% released on 3/17/03 as part of toolbox V 2.5

if nargin<14 % ou_drift_range not given?
    ou_drift_range=2.5; % default value
end;

if nargin<13 % num_calls not given?
    num_calls=9; % default value
end;

num_calls=round(num_calls);

% Some checks
if ~isnumeric(ou_drift_mean)
    error('OU_DROPOUT_VD_1D_2B_MAR: OU_DRIFT_MEAN must be a number!');
end;

if ou_drift_sd<=0
    error('OU_DROPOUT_VD_1D_2B_MAR: OU_DRIFT_SD must be a positive number!');
end;

if ou_leak<0
    error('OU_DROPOUT_VD_1D_2B_MAR: OU_LEAK must not be negative!');
end;

if ou_var<=0
    error('OU_DROPOUT_VD_1D_2B_MAR: The variance must be a positive number!');
end;

if delta_t<=0
    error('OU_DROPOUT_VD_1D_2B_MAR: The time step must be a positive number!');
end;

if diff(covered_range)<0 % screwed up range?
    error('OU_DROPOUT_VD_1D_2B_MAR: Invalid range!');
end;

if stop_time<=0
    error('OU_DROPOUT_VD_1D_2B_MAR: STOP_TIME must be a positive number!');
end;

if num_calls<2
    error('OU_DROPOUT_VD_1D_2B_MAR: NUM_CALLS must be at least 2!');
end;

if ou_drift_range<=0
    error('OU_DROPOUT_VD_1D_2B_MAR: OU_DRIFT_RANGE must be a positive number!');
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
    [g_upper_temp,g_lower_temp,t_temp,dummy,g_dropout_temp]=ou_dropout_1d_2b_mar(drifts(i),ou_leak,ou_var,ou_init,a_upper,a_lower,delta_t,covered_range,num_tr_states,1,stop_time,dropout_rate);
    
    if i==1 % first call?
        g_upper=g_upper_temp*weights(i);
        g_lower=g_lower_temp*weights(i);
        g_dropout=g_dropout_temp*weights(i);
        t_vec=t_temp;
    else
        g_upper=g_upper+g_upper_temp*weights(i);
        g_lower=g_lower+g_lower_temp*weights(i);
        g_dropout=g_dropout+g_dropout_temp*weights(i);
    end;
end;
