function [g_upper,g_lower,t_vec]=probsum_1d_2b_dis(sig_mean,sig_var,thr_upper,thr_lower,stop_type,stop_val)
% Discrete solution of the first passage time problem for
% a probability summation model assuming a
% normal process and two thresholds. Delta_t is set to 1
% to make sure that the transition probabilities are
% identical to the transition rates.
%
% J. Ditterich, 4/02
%
% [g_upper,g_lower,t_vec] = probsum_1d_2b_dis (sig_mean,sig_var,thr_upper,thr_lower,
%                                              stop_type,stop_val)
%
% g_upper is the first passage time density for the upper boundary multiplied
%         by the probability of hitting the upper boundary first, evaluated
%         at the times given in t_vec.
% g_lower is the first passage time density for the lower boundary multiplied
%         by the probability of hitting the lower boundary first, evaluated
%         at the times given in t_vec.
%
% sig_mean is the mean of the normal process. It can either be a constant or,
%          for the time-variant case, the name of a function, which must return
%          the mean when called with the time as the argument.
% sig_var is the variance of the normal process. It can either be a constant or,
%        for the time-variant case, the name of a function, which must return
%        the variance when called with the time as the argument.
% thr_upper defines the upper threshold. thr_upper can either be a constant
%         or the name of a function, which must return the location of the threshold
%         when called with the time as the argument.
% thr_lower defines the lower threshold. See thr_upper for the format.
% stop_type defines, what determines when the algorithm stops.
%           1 = time: The densities are calculated for times up to stop_val.
%           2 = deviation of the distribution from 1:
%               The densities are calculated until the deviation of the sum
%               of the integrals over the densities from 1 is smaller than
%               stop_val.
% stop_val defines the time or the error, which determines when the algorithm
%          will stop.

% History:
% released on 6/6/02 as part of toolbox V 2.0

% Compiler flag:
%#realonly

% Some checks
if (stop_type~=1)&(stop_type~=2)
    error('PROBSUM_1D_2B_DIS: STOP_TYPE must be either 1 or 2!');
end;

if stop_val<=0
    error('PROBSUM_1D_2B_DIS: STOP_VAL must be a positive number!');
end;

% Initialization
g_lower=[];
g_upper=[];
t_vec=[];
t=1;
last_thr_lower=nan;
last_thr_upper=nan;
last_mean=nan;
last_sd=nan;

if isnumeric(thr_lower) % Is the lower threshold time-invariant?
    thr_lower_const=1;
    thr_lower_cur=thr_lower;
else
    thr_lower_const=0;
end;

if isnumeric(thr_upper) % Is the upper threshold time-invariant?
    thr_upper_const=1;
    thr_upper_cur=thr_upper;
else
    thr_upper_const=0;
end;

if isnumeric(sig_mean) % Is the mean time-invariant?
    mean_const=1;
    mean_cur=sig_mean;
else
    mean_const=0;
end;

if isnumeric(sig_var) % Is the variance time-invariant?
    var_const=1;
    
    if sig_var<=0
        error('PROBSUM_1D_2B_DIS: The variance must be a positive number!');
    end;   

    sd_cur=sqrt(sig_var);
else
    var_const=0;
end;

p_not_crossed=1; % initial value for the probability of not having crossed the thresholds so far

% Loop
while (1)
    if (((stop_type==1)&(t>stop_val))|((stop_type==2)&(1-(sum(g_lower')+sum(g_upper'))<stop_val))) % enough?
        break;
    end;
    
    if ~thr_lower_const % Do we have to update the lower threshold?
        thr_lower_cur=feval(thr_lower,t);
    end;
    
    if ~thr_upper_const % Do we have to update the upper threshold?
        thr_upper_cur=feval(thr_upper,t);
    end;
    
    if thr_upper_cur<=thr_lower_cur % crossing thresholds?
        warning('PROBSUM_1D_2B_DIS: Algorithm stopped due to crossing thresholds!');
        break;
    end;
    
    if ~mean_const % Do we have to update the current mean?
        mean_cur=feval(sig_mean,t);
    end;
    
    if ~var_const % Do we have to update the current variance?
        temp=feval(sig_var,t);
        
        if temp<=0
            warning('PROBSUM_1D_2B_DIS: Algorithm stopped due to a non-positive variance!');
            break;
        end;      

        sd_cur=sqrt(temp);        
    end;
    
    % Do we have to recalculate the transition probabilities?
    if (thr_lower_cur~=last_thr_lower)|(thr_upper_cur~=last_thr_upper)|(mean_cur~=last_mean)|(sd_cur~=last_sd)
        last_thr_lower=thr_lower_cur;
        last_thr_upper=thr_upper_cur;
        last_mean=mean_cur;
        last_sd=sd_cur;
        p_lower=normcdf(thr_lower_cur,mean_cur,sd_cur); % transition probability for the lower threshold
        p_upper=1-normcdf(thr_upper_cur,mean_cur,sd_cur); % transition probability for the upper threshold
    end;
    
    g_lower=[g_lower p_not_crossed*p_lower];
    g_upper=[g_upper p_not_crossed*p_upper];
    t_vec=[t_vec t];
    
    p_not_crossed=p_not_crossed*(1-p_lower-p_upper);
    t=t+1;
end; % loop
