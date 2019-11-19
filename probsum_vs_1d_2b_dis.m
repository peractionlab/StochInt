function [g_upper,g_lower,t_vec]=probsum_vs_1d_2b_dis(sig_mean_mean,sig_mean_sd,sig_var,thr_upper,thr_lower,stop_time,num_calls,sig_mean_range)
% Discrete solution of the first passage time problem for
% a probability summation model assuming a
% normal process and two thresholds. The mean of the signal
% varies from trial to trial. Delta_t is set to 1
% to make sure that the transition probabilities are
% identical to the transition rates.
% This function calls PROBSUM_1D_2B_DIS several times in order
% to get an approximate solution.
%
% J. Ditterich, 5/02
%
% [g_upper,g_lower,t_vec] = probsum_vs_1d_2b_dis (sig_mean_mean,sig_mean_sd,sig_var,thr_upper,thr_lower,
%                                                 stop_time[,num_calls[,sig_mean_range]])
%
% g_upper is the first passage time density for the upper boundary multiplied
%         by the probability of hitting the upper boundary first, evaluated
%         at the times given in t_vec.
% g_lower is the first passage time density for the lower boundary multiplied
%         by the probability of hitting the lower boundary first, evaluated
%         at the times given in t_vec.
%
% sig_mean_mean is the expected value of the mean of the normal process.
% sig_mean_sd is the (between trial) standard deviation of the mean of the process.
% sig_var is the variance of the normal process. It can either be a constant or,
%        for the time-variant case, the name of a function, which must return
%        the variance when called with the time as the argument.
% thr_upper defines the upper threshold. thr_upper can either be a constant
%         or the name of a function, which must return the location of the threshold
%         when called with the time as the argument.
% thr_lower defines the lower threshold. See thr_upper for the format.
% stop_time: The densities are calculated for times up to stop_time.
% num_calls is an optional argument, which defines, how often PROBSUM_1D_2B_DIS
%           is called. The default value is 9.
% sig_mean_range is an optional argument, which defines the range used for the
%                calculation. The range is SIG_MEAN_MEAN +/- SIG_MEAN_RANGE *
%                 SIG_MEAN_SD. The default value is 2.5.

% History:
% released on 6/6/02 as part of toolbox V 2.0
% bug fixed on 7/25/02; erroneous argument check removed

if nargin<8 % sig_mean_range not given?
    sig_mean_range=2.5; % default value
end;

if nargin<7 % num_calls not given?
    num_calls=9; % default value
end;

num_calls=round(num_calls);

% Some checks
if ~isnumeric(sig_mean_mean)
    error('PROBSUM_VS_1D_2B_DIS: SIG_MEAN_MEAN must be a number!');
end;

if sig_mean_sd<=0
    error('PROBSUM_VS_1D_2B_DIS: SIG_MEAN_SD must be a positive number!');
end;

if stop_time<=0
    error('PROBSUM_VS_1D_2B_DIS: STOP_TIME must be a positive number!');
end;

if num_calls<2
    error('PROBSUM_VS_1D_2B_DIS: NUM_CALLS must be at least 2!');
end;

if sig_mean_range<=0
    error('PROBSUM_VS_1D_2B_DIS: SIG_MEAN_RANGE must be a positive number!');
end;

% Initialization
temp=[sig_mean_mean-sig_mean_range*sig_mean_sd+sig_mean_range*sig_mean_sd/(num_calls-1): ...
        2*sig_mean_range*sig_mean_sd/(num_calls-1): ...
        sig_mean_mean+sig_mean_range*sig_mean_sd-sig_mean_range*sig_mean_sd/(num_calls-1)];
temp2=[0 normcdf(temp,sig_mean_mean,sig_mean_sd) 1];
weights=diff(temp2); % These are the weights, which have to be applied to the results of the single calls.
sig_means=[sig_mean_mean-sig_mean_range*sig_mean_sd:2*sig_mean_range*sig_mean_sd/(num_calls-1): ...
        sig_mean_mean+sig_mean_range*sig_mean_sd]; % These are the means of the signal, which have to be used in the single calls.

% Loop
for i=1:num_calls
    [g_upper_temp,g_lower_temp,t_temp]=probsum_1d_2b_dis(sig_means(i),sig_var,thr_upper,thr_lower,1,stop_time);
    
    if i==1 % first call?
        g_upper=g_upper_temp*weights(i);
        g_lower=g_lower_temp*weights(i);
        t_vec=t_temp;
    else
        g_upper=g_upper+g_upper_temp*weights(i);
        g_lower=g_lower+g_lower_temp*weights(i);
    end;
end;
