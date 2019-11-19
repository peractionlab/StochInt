function [exp_traj,exp_var,t_vec,num_sim]=traj_afp_ou_vd_1d_2b_sim(ou_drift_mean,ou_drift_sd,ou_leak,ou_var,ou_init,a_upper,a_lower,delta_t,t_min,sel,stop_type,stop_val,runmed_width)
% Expected trajectory and expected variance (aligned
% with respect to the first passage) of a bounded
% (2 (time-variant) boundaries), (time-variant) 1D Ornstein-
% Uhlenbeck process. The drift varies from trial to trial.
% The calculation is based on a simulation.
%
% J. Ditterich, 10/01
%
% [exp_traj,exp_var,t_vec,num_sim] = traj_afp_ou_vd_1d_2b_sim (ou_drift_mean,ou_drift_sd,ou_leak,ou_var,ou_init,a_upper,a_lower,
%                                                              delta_t,t_min,sel,stop_type,stop_val[,runmed_width])
%
% exp_traj is the expected trajectory as a function of time, evaluated at
%          the times given in t_vec (with respect to the first passage).
%          exp_traj is only valid if num_sim is at least 1.
% exp_var is the expected (trial-to-trial) variance as a function of time,
%         evaluated at the times given in t_vec (with respect to the first
%         passage). exp_var is only valid if num_sim is at least 2.
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

% History:
% released on 10/22/01 as part of toolbox V 1.2
% new implementation on 10/31/01

% Compiler flag:
%#realonly

if nargin<13 % runmed_width not given?
    runmed_width=0;
end;

runmed_width=round(runmed_width);
stop_val=round(stop_val);

% Some checks
if ~isnumeric(ou_drift_mean)
    error('TRAJ_AFP_OU_VD_1D_2B_SIM: OU_DRIFT_MEAN must be a number!');
end;

if ou_drift_sd<=0
    error('TRAJ_AFP_OU_VD_1D_2B_SIM: OU_DRIFT_SD must be a positive number!');
end;

if ou_leak<0
    error('TRAJ_AFP_OU_VD_1D_2B_SIM: OU_LEAK must be a non-negative number!');
end;

if (size(ou_init,1)~=1)|(size(ou_init,2)~=1)
    error('TRAJ_AFP_OU_VD_1D_2B_SIM: OU_INIT must be a scalar!');
end;

if delta_t<=0
    error('TRAJ_AFP_OU_VD_1D_2B_SIM: The time step must be a positive number!');
end;

if t_min<=0
    error('TRAJ_AFP_OU_VD_1D_2B_SIM: T_MIN must be a positive number!');
end;

if (sel~=0)&(sel~=1)&(sel~=2)
    error('TRAJ_AFP_OU_VD_1D_2B_SIM: SEL must be either 0, 1, or 2!');
end;

if (stop_type~=1)&(stop_type~=2)
    error('TRAJ_AFP_OU_VD_1D_2B_SIM: STOP_TYPE must be either 1 or 2!');
end;

if stop_val<=0
    error('TRAJ_AFP_OU_VD_1D_2B_SIM: STOP_VAL must be a positive number!');
end;

if runmed_width<0
    error('TRAJ_AFP_OU_VD_1D_2B_SIM: RUNMED_WIDTH must not be negative!');
end;

if (runmed_width>0)&(~mod(runmed_width,2))
    error('TRAJ_AFP_OU_VD_1D_2B_SIM: RUNMED_WIDTH must be an odd number!');
end;

% Initialization
vec_length=floor(t_min/delta_t);
traj_mat=[];
exp_traj=[];
exp_var=[];
t_vec=[-(vec_length-1)*delta_t:delta_t:0];
num_sim=0;
total_sim=0;

if isnumeric(a_lower) % Is the first boundary time-invariant?
    a_lower_const=1;
    a_lower_cur=a_lower;
    
    if ou_init<=a_lower
        error('TRAJ_AFP_OU_VD_1D_2B_SIM: The initial value is out of range!');
    end;
else
    a_lower_const=0;
end;

if isnumeric(a_upper) % Is the second boundary time-invariant?
    a_upper_const=1;
    a_upper_cur=a_upper;
    
    if ou_init>=a_upper
        error('TRAJ_AFP_OU_VD_1D_2B_SIM: The initial value is out of range!');
    end;
else
    a_upper_const=0;
end;

if isnumeric(ou_var) % Is the variance time-invariant?
    var_const=1;
    
    if ou_var<=0
        error('TRAJ_AFP_OU_VD_1D_2B_SIM: The variance must be a positive number!');
    end;
    
    sd_cur=sqrt(ou_var*delta_t);    
else
    var_const=0;
end;

% Further checks
if runmed_width>vec_length
    error('TRAJ_AFP_OU_VD_1D_2B_SIM: Inappropriate filter width!');
end;

% Loop
while (1)
    % draw a new drift rate
    drift_cur=random('norm',ou_drift_mean,ou_drift_sd)*delta_t;
    
    % create a trajectory with a length of vec_length
    if var_const&(ou_leak==0) % In this case we can do it in a single step.
        rand_vec=random('norm',drift_cur,sd_cur,vec_length,1); % drift & noise
        cur_traj=(ou_init+tril(ones(vec_length,vec_length))*rand_vec)'; % integration
    else % separate random calls necessary
        cur_val=ou_init; % start with the initial value
        cur_traj=[];
        
        for i=1:vec_length
            if ou_leak>0 % OU process?
                cur_val=cur_val*(1-ou_leak*delta_t); % leaky integrator part
            end;
            
            cur_val=cur_val+drift_cur; % drift part
            
            if var_const % time-invariant variance?
                cur_val=cur_val+random('norm',0,sd_cur); % noise part
            else
                temp=feval(ou_var,i*delta_t);
                
                if temp<=0
                    error('TRAJ_AFP_OU_VD_1D_2B_SIM: The variance must be positive!');
                end;
                
                sd_cur=sqrt(temp*delta_t);
                cur_val=cur_val+random('norm',0,sd_cur);
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
            a_upper_cur=feval(a_upper,i*delta_t);
            
            if cur_traj(i)>=a_upper_cur
                upper_ind=i;
                break; % We only need the first crossing.
            end;
        end;
    end;
    
    if ~isempty(upper_ind) % Upper boundary has been crossed. => forget trial
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
            a_lower_cur=feval(a_lower,i*delta_t);
            
            if cur_traj(i)<=a_lower_cur
                lower_ind=i;
                break; % We only need the first crossing.
            end;
        end;
    end;
    
    if ~isempty(lower_ind) % Lower boundary has been crossed. => forget trial
        total_sim=total_sim+1;
        
        if (stop_type==1)&(total_sim==stop_val) % We are done ...
            break;
        end;
        
        continue; % next trial
    end;
    
    % no boundaries crossed
    % We have to continue the simulation until the first boundary crossing.
    i=vec_length;
    cur_val=cur_traj(vec_length);
    
    while 1 % inner loop
        i=i+1; % update time
        
        if ou_leak>0 % OU process?
            cur_val=cur_val*(1-ou_leak*delta_t); % leaky integrator part
        end;
        
        cur_val=cur_val+drift_cur; % drift part
        
        if var_const % time-invariant variance?
            cur_val=cur_val+random('norm',0,sd_cur); % noise part
        else
            temp=feval(ou_var,i*delta_t);
            
            if temp<=0
                error('TRAJ_AFP_OU_VD_1D_2B_SIM: The variance must be positive!');
            end;
            
            sd_cur=sqrt(temp*delta_t);
            cur_val=cur_val+random('norm',0,sd_cur);
        end;
        
        cur_traj(i)=cur_val;
        
        % check for boundary crossings
        if ~a_upper_const % time-variant upper boundary?
            a_upper_cur=feval(a_upper,i*delta_t);
        end;
        
        if cur_val>=a_upper_cur % upper boundary crossed?
            if (sel==0)|(sel==1) % good trial?
                total_sim=total_sim+1;
                num_sim=num_sim+1;
                traj_mat(num_sim,:)=cur_traj(i-vec_length+1:i); % store the trial
            else % forget trial
                total_sim=total_sim+1;
            end;
            
            break; % We are done with this trial.
        end;
        
        if ~a_lower_const % time-variant lower boundary?
            a_lower_cur=feval(a_lower,i*delta_t);
        end;
        
        if cur_val<=a_lower_cur % lower boundary crossed?
            if (sel==0)|(sel==2) % good trial?
                total_sim=total_sim+1;
                num_sim=num_sim+1;
                traj_mat(num_sim,:)=cur_traj(i-vec_length+1:i); % store the trial
            else % forget trial
                total_sim=total_sim+1;
            end;
            
            break; % We are done with this trial.
        end;            
    end; % inner loop
    
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
