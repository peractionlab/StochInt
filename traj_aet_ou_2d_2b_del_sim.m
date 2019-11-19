function [exp_traj,exp_var,t_vec,num_sim]=traj_aet_ou_2d_2b_del_sim(ou_drift,ou_leak,ou_cov,ou_init,a_1,a_2,del_1_mean,del_1_sd,del_2_mean,del_2_sd,delta_t,process_range,t_min,sel,stop_type,stop_val,runmed_width,align_timevar)
% Expected trajectory and expected variance (aligned
% with respect to the end of the trial) of a bounded
% (2 (time-variant) boundaries), (time-variant) 2D Ornstein-
% Uhlenbeck process. This version allows variable delays
% between start of a trial and integration onset and between
% the first passage and the end of the trial. If a negative
% random delay should be picked a new value will be drawn.
% It is assumed that the integration process continues
% undisturbed between the first passage and the end of the trial.
% The calculation is based on a simulation.
%
% J. Ditterich, 10/02
%
% [exp_traj,exp_var,t_vec,num_sim] = traj_aet_ou_2d_2b_del_sim (ou_drift,ou_leak,ou_cov,ou_init,a_1,a_2,
%                                                               del_1_mean,del_1_sd,del_2_mean,del_2_sd,
%                                                               delta_t,process_range,t_min,sel,stop_type,
%                                                               stop_val[,runmed_width[,align_timevar]])
%
% exp_traj is the expected trajectory as a function of time, evaluated at
%          the times given in t_vec (with respect to the end of the trial).
%          The first row contains the coordinates of the first dimension,
%          the second row the coordinates of the second dimension.
%          exp_traj is only valid if num_sim is at least 1.
% exp_var is the expected (trial-to-trial) variance as a function of time,
%         evaluated at the times given in t_vec (with respect to the end of
%         the trial). The first row contains the variance in the first dimension,
%         the second row the variance in the second dimension.
%         exp_var is only valid if num_sim is at least 2.
% num_sim is the number of simulations, which have effectively contributed
%         to the result.
%
% ou_drift is the drift vector of the OU process. It can either be a vector of length 2 or,
%          for the time-variant case, the name of a function, which must return
%          the drift vector when called with the time as the argument.
% ou_leak defines the "leakiness" of the integrator(s) and has to be a scalar.
%         The deterministic part of the stochastic differential equation is given by
%         ou_drift - ou_leak * current_value. A Wiener process can be studied
%         by setting ou_leak to 0.
% ou_cov is the covariance matrix of the OU process. It can either be a 2-by-2 matrix or,
%        for the time-variant case, the name of a function, which must return
%        the covariance matrix when called with the time as the argument.
%        The absolute value of the correlation coefficient must be smaller than 1.
%        Please use the 1D function for fully correlated processes.
% ou_init is the initial vector of the OU process. It must be a vector of length 2.
% a_1 defines the first absorbing boundary. a_1 is the name of a function,
%     which must return 1, if a certain location is located on or outside the boundary,
%     and 0, if a certain location is located inside the boundary, when called
%     with a 1-by-2 vector defining the location as the first and time as the
%     second argument.
% a_2 defines the second absorbing boundary. See a_1 for the format. The boundaries
%     should be defined in such a way that both "boundary crossed" regions do not
%     overlap. Since the algorithm checks the first boundary first, a crossing of
%     both boundaries in the same time step will be registered as a crossing of
%     the first boundary.
% del_1_mean defines the mean of a random delay between start of a trial and
%            integration onset.
% del_1_sd defines the standard deviation of the random delay between start of
%          a trial and integration onset.
% del_2_mean defines the mean of a random delay between the first boundary crossing
%            and the end of a trial.
% del_2_sd defines the standard deviation of the random delay between the first
%          boundary crossing and the end of the trial.
% delta_t is the temporal step size.
% process_range defines the valid process range. It normally has to be a 2-by-2 matrix.
%               The first row defines the lower and the upper limit of the first dimension,
%               the second row the lower and the upper limit of the second dimension.
%               Make sure that you define the boundaries in such a way that they are
%               located within this range. Otherwise the algorithm will block.
%               Limiting the process range allows to study the development of the variance
%               of processes with natural limits. When passing 0 the process range is
%               unlimited.
% t_min defines the temporal interval for studying the trajectory. Only trials
%       with a minimum RT of t_min contribute to the result!
% sel defines the selection criterion.
%     0 = All trials with a minimum RT of t_min contribute to the result.
%     1 = Only trials, which will eventually cross the first boundary first,
%         contribute to the result.
%     2 = Only trials, which will eventually cross the second boundary first,
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
% released on 10/30/02 as part of toolbox V 2.3

% Compiler flag:
%#realonly

if nargin<18 % align_timevar not given?
    align_timevar=2; % default value
end;

if nargin<17 % runmed_width not given?
    runmed_width=0; % default value
end;

stop_val=round(stop_val);
runmed_width=round(runmed_width);

% Some checks
if isnumeric(ou_drift)&(length(ou_drift)~=2)
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: OU_DRIFT must be either a vector of length 2 or the name of a function!');
end;

if isnumeric(ou_drift)&(size(ou_drift,1)==2) % wrong orientation?
    ou_drift=ou_drift'; % transpose it
end;

if ou_leak<0
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: OU_LEAK must be a non-negative number!');
end;

if isnumeric(ou_cov)&((size(ou_cov,1)~=2)|(size(ou_cov,2)~=2))
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: OU_COV must either be a 2-by-2 matrix or the name of a function!');
end;

if isnumeric(ou_cov)&((ou_cov(1,1)<=0)|(ou_cov(2,2)<=0))
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: The main diagonal elements of OU_COV must be positive!');
end;

if isnumeric(ou_cov)&(det(ou_cov)==0)
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: The covariance matrix must not be singular!');
end;

if isnumeric(ou_cov)&((det(ou_cov)<0)|(ou_cov(1,2)~=ou_cov(2,1)))
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: Invalid covariance matrix!');
end;

if ~((size(ou_init,1)==1)&(size(ou_init,2)==2))&~((size(ou_init,1)==2)&(size(ou_init,2)==1))
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: OU_INIT must be a vector of length 2!');
end;

if size(ou_init,1)==2 % wrong orientation?
    ou_init=ou_init'; % transpose it
end;

if del_1_mean<0
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: DEL_1_MEAN must not be negative!');
end;

if del_1_sd<0
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: DEL_1_SD must not be negative!');
end;

if del_2_mean<0
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: DEL_2_MEAN must not be negative!');
end;

if del_2_sd<0
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: DEL_2_SD must not be negative!');
end;

if delta_t<=0
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: The time step must be a positive number!');
end;

if (size(process_range,1)~=2)|(size(process_range,2)~=2)
    if (size(process_range,1)~=1)|(size(process_range,2)~=1)|(process_range~=0)
        error('TRAJ_AET_OU_2D_2B_DEL_SIM: PROCESS_RANGE must either be a 2-by-2 matrix or 0!');
    end;
end;

limited_range=(size(process_range,1)==2); % limited range?

if limited_range
    if (diff(process_range(1,:))<0)|(diff(process_range(2,:))<0) % screwed up range?
        error('TRAJ_AET_OU_2D_2B_DEL_SIM: Invalid range!');
    end;
end;

if limited_range
    if (ou_init(1)<process_range(1,1))|(ou_init(1)>process_range(1,2))|(ou_init(2)<process_range(2,1))|(ou_init(2)>process_range(2,2))
        error('TRAJ_AET_OU_2D_2B_DEL_SIM: Initial value out of range!');
    end;
end;

if t_min<=0
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: T_MIN must be a positive number!');
end;

if (sel~=0)&(sel~=1)&(sel~=2)
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: SEL must be either 0, 1, or 2!');
end;

if (stop_type~=1)&(stop_type~=2)
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: STOP_TYPE must be either 1 or 2!');
end;

if stop_val<=0
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: STOP_VAL must be a positive number!');
end;

if runmed_width<0
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: RUNMED_WIDTH must not be negative!');
end;

if runmed_width&(~mod(runmed_width,2))
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: RUNMED_WIDTH must be an odd number!');
end;

if (align_timevar~=1)&(align_timevar~=2)
    error('TRAJ_AET_OU_2D_2B_DEL_SIM: ALIGN_TIMEVAR must be either 1 or 2!');
end;

% Initialization
vec_length=floor(t_min/delta_t);
traj_mat_1=[];
traj_mat_2=[];
exp_traj=[];
exp_var=[];
t_vec=[-(vec_length-1)*delta_t:delta_t:0];
num_sim=0;
total_sim=0;

if isnumeric(ou_drift) % Is the drift time-invariant?
    drift_const=1;
    drift_cur=ou_drift*delta_t;
else
    drift_const=0;
end;

if isnumeric(ou_cov) % Is the covariance matrix time-invariant?
    cov_const=1;
    sqrtm_cov_cur=sqrtm(ou_cov*delta_t);
else
    cov_const=0;
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
    
    boundary_1_crossed=0;
    boundary_2_crossed=0;
    boundary_tested=0;
    
    % create a trajectory with a length of vec_length
    if cov_const&drift_const&(ou_leak==0)&(limited_range==0) % In this case we can do it in a single step.
        rand_vec=random('norm',0,1,2,vec_length); % independent noise
        rand_vec=repmat(drift_cur',1,vec_length)+sqrtm_cov_cur*rand_vec; % drift & correlated noise

        if del_1_dis % initial delay?
            temp=min(del_1_dis,vec_length);
            rand_vec(:,1:temp)=zeros(2,temp);
        end;
        
        cur_traj=(repmat(ou_init,vec_length,1)+tril(ones(vec_length,vec_length))*rand_vec')'; % integration
    elseif cov_const&(ou_leak==0)&(limited_range==0) % This case is not much more difficult ...
        rand_vec=random('norm',0,1,2,vec_length); % independent noise
        rand_vec=sqrtm_cov_cur*rand_vec; % correlated noise
        
        for i=1:vec_length
            if i>del_1_dis % initial delay over?
                temp=feval(ou_drift,(i-time_corr)*delta_t); % get current drift
                
                if ~((size(temp,1)==1)&(size(temp,2)==2))&~((size(temp,1)==2)&(size(temp,2)==1))
                    error('TRAJ_AET_OU_2D_2B_DEL_SIM: The drift returned by a function must be a vector of length 2!');
                end;
                
                if size(temp,2)==2 % wrong orientation?
                    temp=temp'; % transpose it
                end;
                
                rand_vec(:,i)=rand_vec(:,i)+temp*delta_t; % drift part
            else
                rand_vec(:,i)=[0;0];
            end;
        end;
        
        cur_traj=(repmat(ou_init,vec_length,1)+tril(ones(vec_length,vec_length))*rand_vec')'; % integration
    else % separate random calls necessary
        boundary_tested=1;
        cur_val=ou_init; % start with the initial value
        cur_traj=[];
        
        for i=1:vec_length
            if i>del_1_dis % initial delay over?
                if ou_leak>0 % OU process?
                    cur_val=cur_val*(1-ou_leak*delta_t); % leaky integrator part
                end;
                
                if drift_const % time-invariant drift?
                    cur_val=cur_val+drift_cur; % drift part
                else
                    temp=feval(ou_drift,(i-time_corr)*delta_t); % get current drift
                    
                    if ~((size(temp,1)==1)&(size(temp,2)==2))&~((size(temp,1)==2)&(size(temp,2)==1))
                        error('TRAJ_AET_OU_2D_2B_DEL_SIM: The drift returned by a function must be a vector of length 2!');
                    end;
                    
                    if size(temp,1)==2 % wrong orientation?
                        temp=temp'; % transpose it
                    end;
                    
                    cur_val=cur_val+temp*delta_t;
                end;
                
                if cov_const % time-invariant covariance matrix?
                    rand_vec=random('norm',0,1,2,1); % independent noise
                    rand_vec=sqrtm_cov_cur*rand_vec; % correlated noise
                    cur_val=cur_val+rand_vec'; % noise part
                else
                    temp=feval(ou_cov,(i-time_corr)*delta_t); % get current covariance matrix
                    
                    if (size(temp,1)~=2)|(size(temp,2)~=2)
                        error('TRAJ_AET_OU_2D_2B_DEL_SIM: The covariance matrix returned by a function must be a 2-by-2 matrix!');
                    end;
                    
                    if (temp(1,1)<=0)|(temp(2,2)<=0)
                        error('TRAJ_AET_OU_2D_2B_DEL_SIM: Algorithm stopped due to a non-positive variance!');
                    end;
                    
                    if det(temp)==0
                        error('TRAJ_AET_OU_2D_2B_DEL_SIM: Algorithm stopped due to a singular covariance matrix!');
                    end;
                    
                    if (det(temp)<0)|(temp(1,2)~=temp(2,1))
                        error('TRAJ_AET_OU_2D_2B_DEL_SIM: Algorithm stopped due to an invalid covariance matrix!');
                    end;
                    
                    sqrtm_cov_cur=sqrtm(temp*delta_t);
                    rand_vec=random('norm',0,1,2,1); % independent noise
                    rand_vec=sqrtm_cov_cur*rand_vec; % correlated noise
                    cur_val=cur_val+rand_vec'; % noise part
                end;
                
                if limited_range % Do we have to test the range?
                    if cur_val(1)<process_range(1,1)
                        cur_val(1)=process_range(1,1);
                    end;
                    
                    if cur_val(1)>process_range(1,2)
                        cur_val(1)=process_range(1,2);
                    end;
                    
                    if cur_val(2)<process_range(2,1)
                        cur_val(2)=process_range(2,1);
                    end;
                    
                    if cur_val(2)>process_range(2,2)
                        cur_val(2)=process_range(2,2);
                    end;
                end;
                
                if feval(a_1,cur_val,(i-time_corr)*delta_t) % first boundary crossed?
                    if ~boundary_1_crossed&~boundary_2_crossed
                        boundary_1_crossed=1;
                        boundary_crossed=i;
                    end;
                end;

                if feval(a_2,cur_val,(i-time_corr)*delta_t) % second boundary crossed?
                    if ~boundary_1_crossed&~boundary_2_crossed
                        boundary_2_crossed=1;
                        boundary_crossed=i;
                    end;
                end;
            end; % if i>del_1_dis
            
            cur_traj(:,i)=cur_val';                
        end; % for i
    end;
    
    % check for boundary crossing
    if ~boundary_tested    
        for i=1:vec_length
            if i>del_1_dis
                if feval(a_1,cur_traj(:,i)',(i-time_corr)*delta_t) % first boundary crossed?
                    boundary_1_crossed=1;
                    boundary_crossed=i;
                    break;
                end;
                
                if feval(a_2,cur_traj(:,i)',(i-time_corr)*delta_t) % second boundary crossed?
                    boundary_2_crossed=1;
                    boundary_crossed=i;
                    break;
                end;
            end;
        end;
    end;
            
    if (boundary_1_crossed|boundary_2_crossed)&(boundary_crossed+del_2_dis<=vec_length) % boundary crossed?
        total_sim=total_sim+1;
        
        if (stop_type==1)&(total_sim==stop_val) % We are done ...
            break;
        end;
        
        continue; % next trial
    end;
    
    stop_ind=nan; % We don't know when to stop the simulation yet ...
    
    if boundary_1_crossed % Has the first boundary already been crossed?
        if (sel==0)|(sel==1) % good trial?
            stop_ind=boundary_crossed+del_2_dis; % This is when we have to stop the simulation.
        else % forget trial
            total_sim=total_sim+1;
            
            if (stop_type==1)&(total_sim==stop_val) % We are done ...
                break;
            end;
            
            continue; % next trial
        end;
    elseif boundary_2_crossed % Has the second boundary already been crossed?
        if (sel==0)|(sel==2) % good trial?
            stop_ind=boundary_crossed+del_2_dis; % This is when we have to stop the simulation.
        else % forget trial
            total_sim=total_sim+1;
            
            if (stop_type==1)&(total_sim==stop_val) % We are done ...
                break;
            end;
            
            continue; % next trial
        end;
    end;
    
    % We have to continue the simulation until the end of the trial.
    i=vec_length;
    cur_val=cur_traj(:,vec_length)';
    
    while 1 % inner loop
        i=i+1; % update time
        
        if i>del_1_dis % initial delay over?
            if ou_leak>0 % OU process?
                cur_val=cur_val*(1-ou_leak*delta_t); % leaky integrator part
            end;
            
            if drift_const % time-invariant drift?
                cur_val=cur_val+drift_cur; % drift part
            else
                temp=feval(ou_drift,(i-time_corr)*delta_t); % get current drift
                
                if ~((size(temp,1)==1)&(size(temp,2)==2))&~((size(temp,1)==2)&(size(temp,2)==1))
                    error('TRAJ_AET_OU_2D_2B_DEL_SIM: The drift returned by a function must be a vector of length 2!');
                end;
                
                if size(temp,1)==2 % wrong orientation?
                    temp=temp'; % transpose it
                end;
                
                cur_val=cur_val+temp*delta_t;
            end;
            
            if cov_const % time-invariant covariance matrix?
                rand_vec=random('norm',0,1,2,1); % independent noise
                rand_vec=sqrtm_cov_cur*rand_vec; % correlated noise
                cur_val=cur_val+rand_vec'; % noise part
            else
                temp=feval(ou_cov,(i-time_corr)*delta_t); % get current covariance matrix
                
                if (size(temp,1)~=2)|(size(temp,2)~=2)
                    error('TRAJ_AET_OU_2D_2B_DEL_SIM: The covariance matrix returned by a function must be a 2-by-2 matrix!');
                end;
                
                if (temp(1,1)<=0)|(temp(2,2)<=0)
                    error('TRAJ_AET_OU_2D_2B_DEL_SIM: Algorithm stopped due to a non-positive variance!');
                end;
                
                if det(temp)==0
                    error('TRAJ_AET_OU_2D_2B_DEL_SIM: Algorithm stopped due to a singular covariance matrix!');
                end;
                
                if (det(temp)<0)|(temp(1,2)~=temp(2,1))
                    error('TRAJ_AET_OU_2D_2B_DEL_SIM: Algorithm stopped due to an invalid covariance matrix!');
                end;
                
                sqrtm_cov_cur=sqrtm(temp*delta_t);
                rand_vec=random('norm',0,1,2,1); % independent noise
                rand_vec=sqrtm_cov_cur*rand_vec; % correlated noise
                cur_val=cur_val+rand_vec'; % noise part
            end;
            
            if limited_range % Do we have to test the range?
                if cur_val(1)<process_range(1,1)
                    cur_val(1)=process_range(1,1);
                end;
                
                if cur_val(1)>process_range(1,2)
                    cur_val(1)=process_range(1,2);
                end;
                
                if cur_val(2)<process_range(2,1)
                    cur_val(2)=process_range(2,1);
                end;
                
                if cur_val(2)>process_range(2,2)
                    cur_val(2)=process_range(2,2);
                end;
            end;
        end; % if i>del_1_dis
        
        cur_traj(:,i)=cur_val';                
        
        if isnan(stop_ind) % no boundary crossing detected yet?
            % check for boundary crossings
            if feval(a_1,cur_val,(i-time_corr)*delta_t) % first boundary crossed?
                if (sel==0)|(sel==1) % good trial?
                    stop_ind=i+del_2_dis; % This is when we have to stop the simulation.
                else % forget trial
                    total_sim=total_sim+1;
                    break; % We are done with this trial.
                end;
            end;
            
            if feval(a_2,cur_val,(i-time_corr)*delta_t) % second boundary crossed?                
                if (sel==0)|(sel==2) % good trial?
                    stop_ind=i+del_2_dis; % This is when we have to stop the simulation.
                else % forget trial
                    total_sim=total_sim+1;
                    break; % We are done with this trial.
                end;
            end;
        end; % if isnan(stop_ind)

        if ~isnan(stop_ind) % Do we already know the stop time?
            if i==stop_ind % Are we done with this trial?
                total_sim=total_sim+1;
                num_sim=num_sim+1;
                traj_mat_1(num_sim,:)=cur_traj(1,i-vec_length+1:i); % store the trial
                traj_mat_2(num_sim,:)=cur_traj(2,i-vec_length+1:i);
                break; % We are done with this trial.
            end;
        end;
    end; % inner loop
    
    if ((stop_type==1)&(total_sim==stop_val))|((stop_type==2)&(num_sim==stop_val)) % Are we done?
        break; % leave the loop
    end;
end;
    
if num_sim % at least 1 effective trial?
    exp_traj(1,:)=mean(traj_mat_1);
    exp_traj(2,:)=mean(traj_mat_2);
end;

if num_sim>=2 % at least 2 effective trials?
    exp_var(1,:)=var(traj_mat_1);
    exp_var(2,:)=var(traj_mat_2);
end;

if runmed_width % filtering?
    exp_traj(1,:)=runmed(exp_traj(1,:),runmed_width,1,0);
    exp_traj(2,:)=runmed(exp_traj(2,:),runmed_width,1,0);
    exp_var(1,:)=runmed(exp_var(1,:),runmed_width,1,0);
    exp_var(2,:)=runmed(exp_var(2,:),runmed_width,1,0);
end;
