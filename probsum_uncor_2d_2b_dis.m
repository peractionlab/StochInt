function [g_1,g_2,t_vec]=probsum_uncor_2d_2b_dis(sig_mean,sig_var,thr_1,thr_2,stop_type,stop_val)
% Discrete solution of the first passage time problem for
% a probability summation model assuming a race between 2 independent
% normal processes with separate thresholds. Delta_t is set to 1
% to make sure that the transition probabilities are
% identical to the transition rates.
%
% J. Ditterich, 10/02
%
% [g_1,g_2,t_vec] = probsum_uncor_2d_2b_dis (sig_mean,sig_var,thr_1,thr_2,
%                                            stop_type,stop_val)
%
% g_1 is the first passage time density for the first boundary multiplied
%     by the probability of hitting the first boundary first, evaluated
%     at the times given in t_vec.
% g_2 is the first passage time density for the second boundary multiplied
%     by the probability of hitting the second boundary first, evaluated
%     at the times given in t_vec.
%
% sig_mean defines the means of the normal processes. It can either be a vector of
%          length 2 or, for the time-variant case, the name of a function, which
%          must return a vector of length 2 when called with the time as the argument.
% sig_var defines the variances of the normal processes. It can either be a vector of
%         length 2 or, for the time-variant case, the name of a function, which must
%         return a vector of length 2 when called with the time as the argument.
% thr_1 defines the threshold for the first process. thr_1 can either be a constant
%       or the name of a function, which must return the location of the threshold
%       when called with the time as the argument.
% thr_2 defines the threshold for the second process. See thr_1 for the format.
% stop_type defines, what determines when the algorithm stops.
%           1 = time: The densities are calculated for times up to stop_val.
%           2 = deviation of the distribution from 1:
%               The densities are calculated until the deviation of the sum
%               of the integrals over the densities from 1 is smaller than
%               stop_val.
% stop_val defines the time or the error, which determines when the algorithm
%          will stop.

% History:
% released on 10/3/02 as part of toolbox V 2.2

% Compiler flag:
%#realonly

% Some checks
if (stop_type~=1)&(stop_type~=2)
    error('PROBSUM_UNCOR_2D_2B_DIS: STOP_TYPE must be either 1 or 2!');
end;

if stop_val<=0
    error('PROBSUM_UNCOR_2D_2B_DIS: STOP_VAL must be a positive number!');
end;

% Initialization
g_1=[];
g_2=[];
t_vec=[];
t=1;
last_thr_1=nan;
last_thr_2=nan;
last_mean=nan;
last_sd=nan;

if isnumeric(thr_1) % Is the first threshold time-invariant?
    thr_1_const=1;
    thr_1_cur=thr_1;
else
    thr_1_const=0;
end;

if isnumeric(thr_2) % Is the second threshold time-invariant?
    thr_2_const=1;
    thr_2_cur=thr_2;
else
    thr_2_const=0;
end;

if isnumeric(sig_mean) % Is the mean time-invariant?
    if length(sig_mean)~=2
        error('PROBSUM_UNCOR_2D_2B_DIS: Two means have to be given!');
    end;
    
    mean_const=1;
    mean_cur=sig_mean;
else
    mean_const=0;
end;

if isnumeric(sig_var) % Is the variance time-invariant?
    if length(sig_var)~=2
        error('PROBSUM_UNCOR_2D_2B_DIS: Two variances have to be given!');
    end;
    
    var_const=1;
    
    if ~isempty(find(sig_var<=0))
        error('PROBSUM_UNCOR_2D_2B_DIS: The variances must be positive numbers!');
    end;   

    sd_cur=sqrt(sig_var);
else
    var_const=0;
end;

p_not_crossed=1; % initial value for the probability of not having crossed the thresholds so far

% Loop
while (1)
    if (((stop_type==1)&(t>stop_val))|((stop_type==2)&(1-(sum(g_1')+sum(g_2'))<stop_val))) % enough?
        break;
    end;
    
    if ~thr_1_const % Do we have to update the first threshold?
        thr_1_cur=feval(thr_1,t);
    end;
    
    if ~thr_2_const % Do we have to update the second threshold?
        thr_2_cur=feval(thr_2,t);
    end;
    
    if ~mean_const % Do we have to update the current mean?
        mean_cur=feval(sig_mean,t);
        
        if length(mean_cur)~=2
            error('PROBSUM_UNCOR_2D_2B_DIS: The function calculating time-variant means has to return a vector of length 2!');
        end;
    end;
    
    if ~var_const % Do we have to update the current variance?
        temp=feval(sig_var,t);
        
        if length(temp)~=2
            error('PROBSUM_UNCOR_2D_2B_DIS: The function calculating time-variant variances has to return a vector of length 2!');
        end;
        
        if ~isempty(find(temp<=0))
            warning('PROBSUM_UNCOR_2D_2B_DIS: Algorithm stopped due to a non-positive variance!');
            break;
        end;      

        sd_cur=sqrt(temp);        
    end;
    
    % Do we have to recalculate the transition probabilities?
    if (thr_1_cur~=last_thr_1)|(thr_2_cur~=last_thr_2)|(mean_cur~=last_mean)|(sd_cur~=last_sd)
        last_thr_1=thr_1_cur;
        last_thr_2=thr_2_cur;
        last_mean=mean_cur;
        last_sd=sd_cur;
        p_1=1-normcdf(thr_1_cur,mean_cur(1),sd_cur(1)); % transition probability for the first threshold
        p_2=1-normcdf(thr_2_cur,mean_cur(2),sd_cur(2)); % transition probability for the second threshold
    end;
    
    g_1=[g_1 p_not_crossed*p_1];
    g_2=[g_2 p_not_crossed*p_2];
    t_vec=[t_vec t];
    
    p_not_crossed=p_not_crossed*(1-p_1)*(1-p_2);
    t=t+1;
end; % loop
