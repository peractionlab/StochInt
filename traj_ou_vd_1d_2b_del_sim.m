function [exp_traj,exp_var,t_vec,num_sim]=traj_ou_vd_1d_2b_del_sim(ou_drift_mean,ou_drift_sd,ou_leak,ou_var,ou_init,a_upper,a_lower,del_1_mean,del_1_sd,del_2_mean,del_2_sd,delta_t,t_min,sel,stop_type,stop_val,runmed_width,align_timevar)
% Expected trajectory and expected variance of a bounded
% (2 (time-variant) boundaries), (time-variant) 1D Ornstein-
% Uhlenbeck process. The drift varies from trial to trial.
% This version allows variable delays between start of a trial
% and integration onset and between the first passage and the
% end of the trial. If a negative random delay should be picked
% a new value will be drawn. It is assumed that the integration
% process continues undisturbed between the first passage and
% the end of the trial. The calculation is based on a simulation.
%
% J. Ditterich, 7/02
%
% [exp_traj,exp_var,t_vec,num_sim] = traj_ou_vd_1d_2b_del_sim (ou_drift_mean,ou_drift_sd,ou_leak,ou_var,ou_init,a_upper,a_lower,
%                                                              del_1_mean,del_1_sd,del_2_mean,del_2_sd,delta_t,
%                                                              t_min,sel,stop_type,stop_val[,runmed_width[,align_timevar]])
%
% exp_traj is the expected trajectory as a function of time, evaluated at
%          the times given in t_vec. exp_traj is only valid if num_sim is
%          at least 1.
% exp_var is the expected (trial-to-trial) variance as a function of time,
%         evaluated at the times given in t_vec. exp_var is only valid if
%         num_sim is at least 2.
% num_sim is the number of simulations, which have effectively
%         contributed to the result.
%
% ou_drift_mean is the mean drift of the OU process.
% ou_drift_sd is the (between trial) standard deviation of the drift.
%             The drift is constant during a trial.
% ou_leak defines the "leakiness" of the integrator. The deterministic part
%         of the stochastic differential equation is given by
%         ou_drift - ou_leak * current_value. A Wiener process can be studied
%         by setting ou_leak to 0.
% ou_var is the variance of the OU process. It can either be a constant or,
%        for the time-variant case, the name of a function, which must return
%        the variance when called with the time as the argument.
% ou_init is the initial value of the OU process.
% a_upper defines the upper absorbing boundary. a_upper can either be a constant
%         or the name of a function, which must return the location of the boundary
%         when called with the time as the argument.
% a_lower defines the lower absorbing boundary. See a_upper for the format.
% del_1_mean defines the mean of a random delay between start of a trial and
%            integration onset.
% del_1_sd defines the standard deviation of the random delay between start of
%          a trial and integration onset.
% del_2_mean defines the mean of a random delay between the first boundary crossing
%            and the end of a trial.
% del_2_sd defines the standard deviation of the random delay between the first
%          boundary crossing and the end of the trial.
% delta_t is the temporal step size.
% t_min defines the temporal interval for studying the trajectory. Only trials
%       with a minimum RT of t_min contribute to the result!
% sel defines the selection criterion.
%     0 = All trials with a minimum RT of t_min contribute to the result.
%     1 = Only trials, which will eventually cross the upper boundary first,
%         contribute to the result.
%     2 = Only trials, which will eventually cross the lower boundary first,
%         contribute to the result.
% stop_type defines, what determines when the algorithm stops.
%           1 = total number of simulated trials;
%           2 = number of simulated trials contributing to the result
% stop_val defines the number of trials, which determines when the algorithm
%          will stop.
% runmed_width is an optional parameter, which defines the width of a running median filter
%              applied to the output. It has to be an odd number. 0 deactivates the filter.
%              The default value is 0.
% align_timevar is an optional parameter, which defines whether for time-variant
%               parameters t=0 should be aligned with
%               1 = the start of the trial or
%               2 = the integration onset.
%               By default t=0 is aligned with integration onset (2).

% History:
% released on 7/17/02 as part of toolbox V 2.1

% Compiler flag:
%#realonly

if nargin<18 % align_timevar not given?
    align_timevar=2;
end;

if nargin<17 % runmed_width not given?
    runmed_width=0;
end;

runmed_width=round(runmed_width);
stop_val=round(stop_val);

% Some checks
if ~isnumeric(ou_drift_mean)
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: OU_DRIFT_MEAN must be a number!');
end;

if ou_drift_sd<=0
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: OU_DRIFT_SD must be a positive number!');
end;

if ou_leak<0
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: OU_LEAK must be a non-negative number!');
end;

if (size(ou_init,1)~=1)|(size(ou_init,2)~=1)
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: OU_INIT must be a scalar!');
end;

if del_1_mean<0
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: DEL_1_MEAN must not be negative!');
end;

if del_1_sd<0
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: DEL_1_SD must not be negative!');
end;

if del_2_mean<0
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: DEL_2_MEAN must not be negative!');
end;

if del_2_sd<0
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: DEL_2_SD must not be negative!');
end;

if delta_t<=0
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: The time step must be a positive number!');
end;

if t_min<=0
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: T_MIN must be a positive number!');
end;

if (sel~=0)&(sel~=1)&(sel~=2)
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: SEL must be either 0, 1, or 2!');
end;

if (stop_type~=1)&(stop_type~=2)
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: STOP_TYPE must be either 1 or 2!');
end;

if stop_val<=0
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: STOP_VAL must be a positive number!');
end;

if runmed_width<0
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: RUNMED_WIDTH must not be negative!');
end;

if (runmed_width>0)&(~mod(runmed_width,2))
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: RUNMED_WIDTH must be an odd number!');
end;

if (align_timevar~=1)&(align_timevar~=2)
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: ALIGN_TIMEVAR must be either 1 or 2!');
end;

% Initialization
vec_length=floor(t_min/delta_t);
traj_mat=[];
exp_traj=[];
exp_var=[];
t_vec=[delta_t:delta_t:vec_length*delta_t];
num_sim=0;
total_sim=0;

if isnumeric(a_lower) % Is the first boundary time-invariant?
    a_lower_const=1;
    a_lower_cur=a_lower;
    
    if ou_init<=a_lower
        error('TRAJ_OU_VD_1D_2B_DEL_SIM: The initial value is out of range!');
    end;
else
    a_lower_const=0;
end;

if isnumeric(a_upper) % Is the second boundary time-invariant?
    a_upper_const=1;
    a_upper_cur=a_upper;
    
    if ou_init>=a_upper
        error('TRAJ_OU_VD_1D_2B_DEL_SIM: The initial value is out of range!');
    end;
else
    a_upper_const=0;
end;

if isnumeric(ou_var) % Is the variance time-invariant?
    var_const=1;
    
    if ou_var<=0
        error('TRAJ_OU_VD_1D_2B_DEL_SIM: The variance must be a positive number!');
    end;
    
    sd_cur=sqrt(ou_var*delta_t);    
else
    var_const=0;
end;

% Further checks
if runmed_width>vec_length
    error('TRAJ_OU_VD_1D_2B_DEL_SIM: Inappropriate filter width!');
end;

% Loop
while (1)
    % pick random delays
    del_1=random('norm',del_1_mean,del_1_sd);
    
    while del_1<0 % no negative delays
        del_1=random('norm',del_1_mean,del_1_sd);
    end;
    
    del_1_dis=round(del_1/delta_t);
    
    if align_timevar==1 % align time-variant parameters with start of trial
        time_corr=0;
    else % align time-variant parameters with integration onset
        time_corr=del_1_dis;
    end;
    
    del_2=random('norm',del_2_mean,del_2_sd);
    
    while del_2<0 % no negative delays
        del_2=random('norm',del_2_mean,del_2_sd);
    end;
    
    del_2_dis=round(del_2/delta_t);
    
    % draw a new drift rate
    drift_cur=random('norm',ou_drift_mean,ou_drift_sd)*delta_t;
    
    % create a trajectory with a length of vec_length
    if var_const&(ou_leak==0) % In this case we can do it in a single step.
        rand_vec=random('norm',drift_cur,sd_cur,vec_length,1); % drift & noise
        
        if del_1_dis % initial delay?
            rand_vec(1:min(del_1_dis,vec_length))=0;
        end;
        
        cur_traj=(ou_init+tril(ones(vec_length,vec_length))*rand_vec)'; % integration
    else % separate random calls necessary
        cur_val=ou_init; % start with the initial value
        cur_traj=[];
        
        for i=1:vec_length
            if i>del_1_dis % initial delay?
                if ou_leak>0 % OU process?
                    cur_val=cur_val*(1-ou_leak*delta_t); % leaky integrator part
                end;
                
                cur_val=cur_val+drift_cur; % drift part
                
                if var_const % time-invariant variance?
                    cur_val=cur_val+random('norm',0,sd_cur); % noise part
                else
                    temp=feval(ou_var,(i-time_corr)*delta_t);
                    
                    if temp<=0
                        error('TRAJ_OU_VD_1D_2B_DEL_SIM: The variance must be positive!');
                    end;
                    
                    sd_cur=sqrt(temp*delta_t);
                    cur_val=cur_val+random('norm',0,sd_cur);
                end;
            end;
            
            cur_traj(i)=cur_val;                
        end;
    end;
    
    % check for boundary crossings
    if a_upper_const % time-invariant upper boundary?
        upper_ind=find(cur_traj>=a_upper_cur); % find upper boundary crossings
    else
        upper_ind=[];
        
        for i=1:vec_length
            if i>del_1_dis
                a_upper_cur=feval(a_upper,(i-time_corr)*delta_t);
            else
                a_upper_cur=feval(a_upper,0);
            end;
            
            if cur_traj(i)>=a_upper_cur
                upper_ind=i;
                break; % We only need the first crossing.
            end;
        end;
    end;
    
    if (~isempty(upper_ind))&(upper_ind(1)+del_2_dis<=vec_length) % Upper boundary has been crossed. => forget trial
        total_sim=total_sim+1;
        
        if (stop_type==1)&(total_sim==stop_val) % We are done ...
            break;
        end;
        
        continue; % next trial
    end;
    
    if a_lower_const % time-invariant lower boundary?
        lower_ind=find(cur_traj<=a_lower_cur); % find lower boundary crossings
    else
        lower_ind=[];
        
        for i=1:vec_length
            if i>del_1_dis
                a_lower_cur=feval(a_lower,(i-time_corr)*delta_t);
            else
                a_lower_cur=feval(a_lower,0);
            end;
            
            if cur_traj(i)<=a_lower_cur
                lower_ind=i;
                break; % We only need the first crossing.
            end;
        end;
    end;
    
    if (~isempty(lower_ind))&(lower_ind(1)+del_2_dis<=vec_length) % Lower boundary has been crossed. => forget trial
        total_sim=total_sim+1;
        
        if (stop_type==1)&(total_sim==stop_val) % We are done ...
            break;
        end;
        
        continue; % next trial
    end;
    
    if sel==0 % valid trial?
        total_sim=total_sim+1;
        num_sim=num_sim+1;
        traj_mat(num_sim,:)=cur_traj; % store the trial
    else % selection requested
        if (~isempty(upper_ind))|(~isempty(lower_ind)) % Has the first passage already happened?
            if isempty(upper_ind)
                upper_ind=vec_length+1;
            end;
            
            if isempty(lower_ind)
                lower_ind=vec_length+1;
            end;
            
            if upper_ind(1)<=lower_ind(1) % upper boundary crossing
                if sel==1 % good trial?
                    total_sim=total_sim+1;
                    num_sim=num_sim+1;
                    traj_mat(num_sim,:)=cur_traj; % store the trial
                else % forget trial
                    total_sim=total_sim+1;
                end;
            else % lower boundary crossing
                if sel==2 % good trial?
                    total_sim=total_sim+1;
                    num_sim=num_sim+1;
                    traj_mat(num_sim,:)=cur_traj; % store the trial
                else % forget trial
                    total_sim=total_sim+1;
                end;
            end;
        else % We have to continue the simulation until the first boundary crossing.
            i=vec_length;
            cur_val=cur_traj(vec_length);
            
            while 1 % inner loop
                i=i+1; % update time
                
                if i>del_1_dis % initial delay over?
                    if ou_leak>0 % OU process?
                        cur_val=cur_val*(1-ou_leak*delta_t); % leaky integrator part
                    end;
                    
                    cur_val=cur_val+drift_cur; % drift part
                    
                    if var_const % time-invariant variance?
                        cur_val=cur_val+random('norm',0,sd_cur); % noise part
                    else
                        temp=feval(ou_var,(i-time_corr)*delta_t);
                        
                        if temp<=0
                            error('TRAJ_OU_VD_1D_2B_DEL_SIM: The variance must be positive!');
                        end;
                        
                        sd_cur=sqrt(temp*delta_t);
                        cur_val=cur_val+random('norm',0,sd_cur);
                    end;
                end;
                
                % check for boundary crossings
                if ~a_upper_const % time-variant upper boundary?
                    if i>del_1_dis
                        a_upper_cur=feval(a_upper,(i-time_corr)*delta_t);
                    else
                        a_upper_cur=feval(a_upper,0);
                    end;
                end;
                
                if cur_val>=a_upper_cur % upper boundary crossed?
                    if sel==1 % good trial?
                        total_sim=total_sim+1;
                        num_sim=num_sim+1;
                        traj_mat(num_sim,:)=cur_traj; % store the trial
                    else % forget trial
                        total_sim=total_sim+1;
                    end;
                    
                    break; % We are done with this trial.
                end;
                
                if ~a_lower_const % time-variant lower boundary?
                    if i>del_1_dis
                        a_lower_cur=feval(a_lower,(i-time_corr)*delta_t);
                    else
                        a_lower_cur=feval(a_lower,0);
                    end;
                end;
                
                if cur_val<=a_lower_cur % lower boundary crossed?
                    if sel==2 % good trial?
                        total_sim=total_sim+1;
                        num_sim=num_sim+1;
                        traj_mat(num_sim,:)=cur_traj; % store the trial
                    else % forget trial
                        total_sim=total_sim+1;
                    end;
                    
                    break; % We are done with this trial.
                end;            
            end; % inner loop
        end; % Has the first passage already happened?
    end; % if sel==0
    
    if ((stop_type==1)&(total_sim==stop_val))|((stop_type==2)&(num_sim==stop_val)) % Are we done?
        break; % leave the loop
    end;
end;

if num_sim % at least 1 effective trial?
    exp_traj=mean(traj_mat);
end;

if num_sim>=2 % at least 2 effective trials?
    exp_var=var(traj_mat);
end;

if runmed_width % filtering?
    exp_traj=runmed(exp_traj,runmed_width,1,0);
    exp_var=runmed(exp_var,runmed_width,1,0);
end;
