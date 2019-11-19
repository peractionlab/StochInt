function [exp_traj,exp_var,t_vec,num_sim]=traj_ou_2d_2b_sim(ou_drift,ou_leak,ou_cov,ou_init,a_1,a_2,delta_t,process_range,t_min,sel,stop_type,stop_val,runmed_width)
% Expected trajectory and expected variance of a bounded
% (2 (time-variant) boundaries), (time-variant) 2D Ornstein-
% Uhlenbeck process. The calculation is based on a simulation.
%
% J. Ditterich, 10/01
%
% [exp_traj,exp_var,t_vec,num_sim] = traj_ou_2d_2b_sim (ou_drift,ou_leak,ou_cov,ou_init,a_1,a_2,delta_t,
%                                                       process_range,t_min,sel,stop_type,stop_val
%                                                       [,runmed_width])
%
% exp_traj is the expected trajectory as a function of time, evaluated at
%          the times given in t_vec. The first row contains the coordinates
%          of the first dimension, the second row the coordinates of the
%          second dimension. exp_traj is only valid if num_sim is at least 1.
% exp_var is the expected (trial-to-trial) variance as a function of time,
%         evaluated at the times given in t_vec. The first row contains the
%         variance in the first dimension, the second row the variance in
%         the second dimension. exp_var is only valid if num_sim is at least 2.
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

% History:
% released on 10/11/01 as part of toolbox V 1.0
% initialization bug fixed on 10/29/01

% Compiler flag:
%#realonly

if nargin<13 % runmed_width not given?
    runmed_width=0; % default value
end;

stop_val=round(stop_val);
runmed_width=round(runmed_width);

% Some checks
if isnumeric(ou_drift)&(length(ou_drift)~=2)
    error('TRAJ_OU_2D_2B_SIM: OU_DRIFT must be either a vector of length 2 or the name of a function!');
end;

if isnumeric(ou_drift)&(size(ou_drift,1)==2) % wrong orientation?
    ou_drift=ou_drift'; % transpose it
end;

if ou_leak<0
    error('TRAJ_OU_2D_2B_SIM: OU_LEAK must be a non-negative number!');
end;

if isnumeric(ou_cov)&((size(ou_cov,1)~=2)|(size(ou_cov,2)~=2))
    error('TRAJ_OU_2D_2B_SIM: OU_COV must either be a 2-by-2 matrix or the name of a function!');
end;

if isnumeric(ou_cov)&((ou_cov(1,1)<=0)|(ou_cov(2,2)<=0))
    error('TRAJ_OU_2D_2B_SIM: The main diagonal elements of OU_COV must be positive!');
end;

if isnumeric(ou_cov)&(det(ou_cov)==0)
    error('TRAJ_OU_2D_2B_SIM: The covariance matrix must not be singular!');
end;

if isnumeric(ou_cov)&((det(ou_cov)<0)|(ou_cov(1,2)~=ou_cov(2,1)))
    error('TRAJ_OU_2D_2B_SIM: Invalid covariance matrix!');
end;

if ~((size(ou_init,1)==1)&(size(ou_init,2)==2))&~((size(ou_init,1)==2)&(size(ou_init,2)==1))
    error('TRAJ_OU_2D_2B_SIM: OU_INIT must be a vector of length 2!');
end;

if size(ou_init,1)==2 % wrong orientation?
    ou_init=ou_init'; % transpose it
end;

if delta_t<=0
    error('TRAJ_OU_2D_2B_SIM: The time step must be a positive number!');
end;

if (size(process_range,1)~=2)|(size(process_range,2)~=2)
    if (size(process_range,1)~=1)|(size(process_range,2)~=1)|(process_range~=0)
        error('TRAJ_OU_2D_2B_SIM: PROCESS_RANGE must either be a 2-by-2 matrix or 0!');
    end;
end;

limited_range=(size(process_range,1)==2); % limited range?

if limited_range
    if (diff(process_range(1,:))<0)|(diff(process_range(2,:))<0) % screwed up range?
        error('TRAJ_OU_2D_2B_SIM: Invalid range!');
    end;
end;

if limited_range
    if (ou_init(1)<process_range(1,1))|(ou_init(1)>process_range(1,2))|(ou_init(2)<process_range(2,1))|(ou_init(2)>process_range(2,2))
        error('TRAJ_OU_2D_2B_SIM: Initial value out of range!');
    end;
end;

if t_min<=0
    error('TRAJ_OU_2D_2B_SIM: T_MIN must be a positive number!');
end;

if (sel~=0)&(sel~=1)&(sel~=2)
    error('TRAJ_OU_2D_2B_SIM: SEL must be either 0, 1, or 2!');
end;

if (stop_type~=1)&(stop_type~=2)
    error('TRAJ_OU_2D_2B_SIM: STOP_TYPE must be either 1 or 2!');
end;

if stop_val<=0
    error('TRAJ_OU_2D_2B_SIM: STOP_VAL must be a positive number!');
end;

if runmed_width<0
    error('TRAJ_OU_2D_2B_SIM: RUNMED_WIDTH must not be negative!');
end;

if runmed_width&(~mod(runmed_width,2))
    error('TRAJ_OU_2D_2B_SIM: RUNMED_WIDTH must be an odd number!');
end;

% Initialization
vec_length=floor(t_min/delta_t);
traj_mat_1=[];
traj_mat_2=[];
exp_traj=[];
exp_var=[];
t_vec=[delta_t:delta_t:vec_length*delta_t];
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
    boundary_crossed=0;
    boundary_tested=0;
    
    % create a trajectory with a length of vec_length
    if cov_const&drift_const&(ou_leak==0)&(limited_range==0) % In this case we can do it in a single step.
        rand_vec=random('norm',0,1,2,vec_length); % independent noise
        rand_vec=repmat(drift_cur',1,vec_length)+sqrtm_cov_cur*rand_vec; % drift & correlated noise
        cur_traj=(repmat(ou_init,vec_length,1)+tril(ones(vec_length,vec_length))*rand_vec')'; % integration
    elseif cov_const&(ou_leak==0)&(limited_range==0) % This case is not much more difficult ...
        rand_vec=random('norm',0,1,2,vec_length); % independent noise
        rand_vec=sqrtm_cov_cur*rand_vec; % correlated noise
        
        for i=1:vec_length
            temp=feval(ou_drift,i*delta_t); % get current drift
            
            if ~((size(temp,1)==1)&(size(temp,2)==2))&~((size(temp,1)==2)&(size(temp,2)==1))
                error('TRAJ_OU_2D_2B_SIM: The drift returned by a function must be a vector of length 2!');
            end;
            
            if size(temp,2)==2 % wrong orientation?
                temp=temp'; % transpose it
            end;
            
            rand_vec(:,i)=rand_vec(:,i)+temp*delta_t; % drift part
        end;
        
        cur_traj=(repmat(ou_init,vec_length,1)+tril(ones(vec_length,vec_length))*rand_vec')'; % integration
    else % separate random calls necessary
        boundary_tested=1;
        cur_val=ou_init; % start with the initial value
        cur_traj=[];
        
        for i=1:vec_length
            if ou_leak>0 % OU process?
                cur_val=cur_val*(1-ou_leak*delta_t); % leaky integrator part
            end;
            
            if drift_const % time-invariant drift?
                cur_val=cur_val+drift_cur; % drift part
            else
                temp=feval(ou_drift,i*delta_t); % get current drift
                
                if ~((size(temp,1)==1)&(size(temp,2)==2))&~((size(temp,1)==2)&(size(temp,2)==1))
                    error('TRAJ_OU_2D_2B_SIM: The drift returned by a function must be a vector of length 2!');
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
                temp=feval(ou_cov,i*delta_t); % get current covariance matrix
                
                if (size(temp,1)~=2)|(size(temp,2)~=2)
                    error('TRAJ_OU_2D_2B_SIM: The covariance matrix returned by a function must be a 2-by-2 matrix!');
                end;
                
                if (temp(1,1)<=0)|(temp(2,2)<=0)
                    error('TRAJ_OU_2D_2B_SIM: Algorithm stopped due to a non-positive variance!');
                end;
                
                if det(temp)==0
                    error('TRAJ_OU_2D_2B_SIM: Algorithm stopped due to a singular covariance matrix!');
                end;
                
                if (det(temp)<0)|(temp(1,2)~=temp(2,1))
                    error('TRAJ_OU_2D_2B_SIM: Algorithm stopped due to an invalid covariance matrix!');
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
            
            if feval(a_1,cur_val,i*delta_t)|feval(a_2,cur_val,i*delta_t) % boundary crossed?
                boundary_crossed=1;
                break; % We no longer have to calculate the rest of the trajectory.
            end;
            
            cur_traj(:,i)=cur_val';                
        end; % for i
    end;
    
    % check for boundary crossing
    if ~boundary_tested    
        for i=1:vec_length
            if feval(a_1,cur_traj(:,i)',i*delta_t)|feval(a_2,cur_traj(:,i)',i*delta_t) % boundary crossed?
                boundary_crossed=1;
                break;
            end;
        end;
    end;
            
    if boundary_crossed % boundary crossed?
        total_sim=total_sim+1;
        
        if (stop_type==1)&(total_sim==stop_val) % We are done ...
            break;
        end;
        
        continue; % next trial
    end;
    
    % no boundaries crossed
    if sel==0 % valid trial?
        total_sim=total_sim+1;
        num_sim=num_sim+1;
        traj_mat_1(num_sim,:)=cur_traj(1,:); % store the trial
        traj_mat_2(num_sim,:)=cur_traj(2,:);
    else % We have to continue the simulation until the first boundary crossing.
        i=vec_length;
        cur_val=cur_traj(:,vec_length)';
        
        while 1 % inner loop
            i=i+1; % update time
            
            if ou_leak>0 % OU process?
                cur_val=cur_val*(1-ou_leak*delta_t); % leaky integrator part
            end;
            
            if drift_const % time-invariant drift?
                cur_val=cur_val+drift_cur; % drift part
            else
                temp=feval(ou_drift,i*delta_t); % get current drift
                
                if ~((size(temp,1)==1)&(size(temp,2)==2))&~((size(temp,1)==2)&(size(temp,2)==1))
                    error('TRAJ_OU_2D_2B_SIM: The drift returned by a function must be a vector of length 2!');
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
                temp=feval(ou_cov,i*delta_t); % get current covariance matrix
                
                if (size(temp,1)~=2)|(size(temp,2)~=2)
                    error('TRAJ_OU_2D_2B_SIM: The covariance matrix returned by a function must be a 2-by-2 matrix!');
                end;
                
                if (temp(1,1)<=0)|(temp(2,2)<=0)
                    error('TRAJ_OU_2D_2B_SIM: Algorithm stopped due to a non-positive variance!');
                end;
                
                if det(temp)==0
                    error('TRAJ_OU_2D_2B_SIM: Algorithm stopped due to a singular covariance matrix!');
                end;
                
                if (det(temp)<0)|(temp(1,2)~=temp(2,1))
                    error('TRAJ_OU_2D_2B_SIM: Algorithm stopped due to an invalid covariance matrix!');
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
            
            % check for boundary crossings
            if feval(a_1,cur_val,i*delta_t) % first boundary crossed?
                if sel==1 % good trial?
                    total_sim=total_sim+1;
                    num_sim=num_sim+1;
                    traj_mat_1(num_sim,:)=cur_traj(1,:); % store the trial
                    traj_mat_2(num_sim,:)=cur_traj(2,:);
                else % forget trial
                    total_sim=total_sim+1;
                end;
                
                break; % We are done with this trial.
            end;
                
            if feval(a_2,cur_val,i*delta_t) % second boundary crossed?                
                if sel==2 % good trial?
                    total_sim=total_sim+1;
                    num_sim=num_sim+1;
                    traj_mat_1(num_sim,:)=cur_traj(1,:); % store the trial
                    traj_mat_2(num_sim,:)=cur_traj(2,:);
                else % forget trial
                    total_sim=total_sim+1;
                end;
                
                break; % We are done with this trial.
            end;
        end; % inner loop
    end; % if sel==0
    
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
