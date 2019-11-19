function [exp_traj,exp_var,t_vec,num_sim]=traj_msprt_fb_3d_3b_sim(inp_mean,inp_cov,msprt_gain,a_1,a_2,a_3,delta_t,process_range,t_min,sel,stop_type,stop_val,runmed_width)
% Expected trajectory and expected variance of a bounded
% (3 (time-variant) boundaries), (time-variant) 3D MSPRT
% process with the MSPRT normalization being provided by inhibitory
% feedback. The calculation is based on a simulation.
%
% J. Ditterich, 6/10
%
% [exp_traj,exp_var,t_vec,num_sim] = traj_msprt_fb_3d_3b_sim (inp_mean,inp_cov,msprt_gain,a_1,a_2,a_3,
%                                                              delta_t,process_range,t_min,sel,
%                                                              stop_type,stop_val[,runmed_width])
%
% exp_traj is the expected trajectory as a function of time, evaluated at
%          the times given in t_vec. The first row contains the coordinates
%          of the first dimension, the second row the coordinates of the
%          second dimension, and the third row the coordinates of the third
%          dimension. exp_traj is only valid if num_sim is at least 1.
% exp_var is the expected (trial-to-trial) variance as a function of time,
%         evaluated at the times given in t_vec. The first row contains the
%         variance in the first dimension, the second row the variance in
%         the second dimension, and the third row the variance in the third
%         dimension. exp_var is only valid if num_sim is at least 2.
% num_sim is the number of simulations, which have effectively contributed
%         to the result.
%
% inp_mean is the mean of the input vector to the MSPRT process (per unit time).
%          It can either be a vector of length 3 or, for the time-variant case,
%          the name of a function, which must return the vector of means when
%          called with the time as the argument. Add a fixed positive offset
%          to each of the means for moving the threshold crossing into the
%          positive range. E.g., if you want the process to stop when the
%          log posterior probability exceeds -0.4, but you want this crossing
%          to happen when an integrator value exceeds 1.0, you have to apply
%          an offset of 1.4 (the difference between the two) divided by delta_t.
%          (The log of the posterior probability does not change when the same
%          offset is applied to all integrators.)
% inp_cov is the covariance matrix of the input to the MSPRT process (per unit time).
%         It can either be a 3-by-3 matrix or, for the time-variant case,
%         the name of a function, which must return the covariance matrix
%         when called with the time as the argument. The absolute value of
%         the correlation coefficients must be smaller than 1.
% msprt_gain defines the gain of the MSPRT normalization. Each integrator
%            gets normalized to msprt_gain times its value minus the log of
%            the sum of the exponentials of each integrator value
%            multiplied by msprt_gain. Note that a gain different from 1
%            leads to either leaky integration or an unstable process.
% a_1 defines the first absorbing boundary. a_1 is the name of a function,
%     which must return 1, if a certain location is located on or outside the boundary,
%     and 0, if a certain location is located inside the boundary, when called
%     with a 1-by-3 vector defining the location as the first and time as the
%     second argument.
% a_2 defines the second absorbing boundary. See a_1 for the format. The boundaries
%     should be defined in such a way that "boundary crossed" regions do not
%     overlap. Since the algorithm checks the first boundary first, a crossing of
%     multiple boundaries in the same time step will be registered as a crossing of
%     the first boundary.
% a_3 defines the third absorbing boundary. See a_1 for the format. The boundaries
%     should be defined in such a way that "boundary crossed" regions do not
%     overlap. Since the algorithm checks the first boundary first, a crossing of
%     multiple boundaries in the same time step will be registered as a crossing of
%     the first boundary.
% delta_t is the temporal step size.
% process_range defines the valid process range. It normally has to be a 3-by-2 matrix.
%               The first row defines the lower and the upper limit of the first dimension,
%               the second row the lower and the upper limit of the second dimension,
%               and the third row the lower and upper limit of the third dimension.
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
%     3 = Only trials, which will eventually cross the third boundary
%         first, contribute to the result.
% stop_type defines, what determines when the algorithm stops.
%           1 = total number of simulated trials;
%           2 = number of simulated trials contributing to the result
% stop_val defines the number of trials, which determines when the algorithm
%          will stop.
% runmed_width is an optional parameter, which defines the width of a running median filter
%              applied to the output. It has to be an odd number. 0 deactivates the filter.
%              The default value is 0.

% History:
% released on 8/13/10 as part of toolbox V 2.7

if nargin<13 % runmed_width not given?
    runmed_width=0; % default value
end;

stop_val=round(stop_val);
runmed_width=round(runmed_width);

% Some checks
if isnumeric(inp_mean)&&(length(inp_mean)~=3)
    error('TRAJ_MSPRT_FB_3D_3B_SIM: inp_mean must be either a vector of length 3 or the name of a function!');
end;

if isnumeric(inp_mean)&&(size(inp_mean,1)==3) % wrong orientation?
    inp_mean=inp_mean'; % transpose it
end;

if isnumeric(inp_cov)&&((size(inp_cov,1)~=3)||(size(inp_cov,2)~=3))
    error('TRAJ_MSPRT_FB_3D_3B_SIM: inp_cov must either be a 3-by-3 matrix or the name of a function!');
end;

if isnumeric(inp_cov)&&((inp_cov(1,1)<=0)||(inp_cov(2,2)<=0)||(inp_cov(3,3)<=0))
    error('TRAJ_MSPRT_FB_3D_3B_SIM: The main diagonal elements of inp_cov must be positive!');
end;

if isnumeric(inp_cov)&&(det(inp_cov)==0)
    error('TRAJ_MSPRT_FB_3D_3B_SIM: The covariance matrix must not be singular!');
end;

if isnumeric(inp_cov)&&((det(inp_cov)<0)||(inp_cov(1,2)~=inp_cov(2,1))||(inp_cov(1,3)~=inp_cov(3,1))||(inp_cov(2,3)~=inp_cov(3,2)))
    error('TRAJ_MSPRT_FB_3D_3B_SIM: Invalid covariance matrix!');
end;

if msprt_gain<=0
    error('TRAJ_MSPRT_FB_3D_3B_SIM: The MSPRT gain must be a positive number!');
end;

if delta_t<=0
    error('TRAJ_MSPRT_FB_3D_3B_SIM: The time step must be a positive number!');
end;

if (size(process_range,1)~=3)||(size(process_range,2)~=2)
    if (size(process_range,1)~=1)||(size(process_range,2)~=1)||(process_range~=0)
        error('TRAJ_MSPRT_FB_3D_3B_SIM: PROCESS_RANGE must either be a 3-by-2 matrix or 0!');
    end;
end;

limited_range=(size(process_range,1)==3); % limited range?

if limited_range
    if (diff(process_range(1,:))<0)||(diff(process_range(2,:))<0||(diff(process_range(3,:))<0)) % screwed up range?
        error('TRAJ_MSPRT_FB_3D_3B_SIM: Invalid range!');
    end;
end;

if t_min<=0
    error('TRAJ_MSPRT_FB_3D_3B_SIM: T_MIN must be a positive number!');
end;

if (sel~=0)&&(sel~=1)&&(sel~=2)&&(sel~=3)
    error('TRAJ_MSPRT_FB_3D_3B_SIM: SEL must be either 0, 1, 2, or 3!');
end;

if (stop_type~=1)&&(stop_type~=2)
    error('TRAJ_MSPRT_FB_3D_3B_SIM: STOP_TYPE must be either 1 or 2!');
end;

if stop_val<=0
    error('TRAJ_MSPRT_FB_3D_3B_SIM: STOP_VAL must be a positive number!');
end;

if runmed_width<0
    error('TRAJ_MSPRT_FB_3D_3B_SIM: RUNMED_WIDTH must not be negative!');
end;

if runmed_width&&(~mod(runmed_width,2))
    error('TRAJ_MSPRT_FB_3D_3B_SIM: RUNMED_WIDTH must be an odd number!');
end;

% Initialization
vec_length=floor(t_min/delta_t);
traj_mat_1=[];
traj_mat_2=[];
traj_mat_3=[];
exp_traj=[];
exp_var=[];
t_vec=delta_t:delta_t:vec_length*delta_t;
num_sim=0;
total_sim=0;

if isnumeric(inp_mean) % Is the drift time-invariant?
    drift_const=1;
    drift_cur=inp_mean*delta_t;
else
    drift_const=0;
end;

if isnumeric(inp_cov) % Is the covariance matrix time-invariant?
    cov_const=1;
    sqrtm_cov_cur=sqrtm(inp_cov*delta_t);
else
    cov_const=0;
end;

% Loop
while (1)
    boundary_crossed=0;
    
    % create a trajectory with a length of vec_length
    cur_val=[0 0 0];
    cur_traj=[];
    
    for i=1:vec_length
        cur_val=msprt_gain*cur_val-log(sum(exp(msprt_gain*cur_val))); % MSPRT normalization
        
        if drift_const % time-invariant drift?
            cur_val=cur_val+drift_cur; % drift part
        else
            temp=feval(inp_mean,i*delta_t); % get current drift
            
            if ~((size(temp,1)==1)&&(size(temp,2)==3))&&~((size(temp,1)==3)&&(size(temp,2)==1))
                error('TRAJ_MSPRT_FB_3D_3B_SIM: The drift returned by a function must be a vector of length 3!');
            end;
            
            if size(temp,1)==3 % wrong orientation?
                temp=temp'; % transpose it
            end;
            
            cur_val=cur_val+temp*delta_t;
        end;
        
        if cov_const % time-invariant covariance matrix?
            rand_vec=random('norm',0,1,3,1); % independent noise
            rand_vec=sqrtm_cov_cur*rand_vec; % correlated noise
            cur_val=cur_val+rand_vec'; % noise part
        else
            temp=feval(inp_cov,i*delta_t); % get current covariance matrix
            
            if (size(temp,1)~=3)||(size(temp,2)~=3)
                error('TRAJ_MSPRT_FB_3D_3B_SIM: The covariance matrix returned by a function must be a 3-by-3 matrix!');
            end;
            
            if (temp(1,1)<=0)||(temp(2,2)<=0)||(temp(3,3)<=0)
                error('TRAJ_MSPRT_FB_3D_3B_SIM: Algorithm stopped due to a non-positive variance!');
            end;
            
            if det(temp)==0
                error('TRAJ_MSPRT_FB_3D_3B_SIM: Algorithm stopped due to a singular covariance matrix!');
            end;
            
            if (det(temp)<0)||(temp(1,2)~=temp(2,1))||(temp(1,3)~=temp(3,1))||(temp(2,3)~=temp(3,2))
                error('TRAJ_MSPRT_FB_3D_3B_SIM: Algorithm stopped due to an invalid covariance matrix!');
            end;
            
            sqrtm_cov_cur=sqrtm(temp*delta_t);
            rand_vec=random('norm',0,1,3,1); % independent noise
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
            
            if cur_val(3)<process_range(3,1)
                cur_val(3)=process_range(3,1);
            end;
            
            if cur_val(3)>process_range(3,2)
                cur_val(3)=process_range(3,2);
            end;
        end;
        
        if feval(a_1,cur_val,i*delta_t)||feval(a_2,cur_val,i*delta_t)||feval(a_3,cur_val,i*delta_t) % boundary crossed?
            boundary_crossed=1;
            break; % We no longer have to calculate the rest of the trajectory.
        end;
        
        cur_traj(:,i)=cur_val';
    end; % for i
    
    if boundary_crossed % boundary crossed?
        total_sim=total_sim+1;
        
        if (stop_type==1)&&(total_sim==stop_val) % We are done ...
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
        traj_mat_3(num_sim,:)=cur_traj(3,:);
    else % We have to continue the simulation until the first boundary crossing.
        i=vec_length;
        cur_val=cur_traj(:,vec_length)';
        
        while 1 % inner loop
            i=i+1; % update time
            cur_val=msprt_gain*cur_val-log(sum(exp(msprt_gain*cur_val))); % MSPRT normalization
            
            if drift_const % time-invariant drift?
                cur_val=cur_val+drift_cur; % drift part
            else
                temp=feval(inp_mean,i*delta_t); % get current drift
                
                if ~((size(temp,1)==1)&&(size(temp,2)==3))&&~((size(temp,1)==3)&&(size(temp,2)==1))
                    error('TRAJ_MSPRT_FB_3D_3B_SIM: The drift returned by a function must be a vector of length 3!');
                end;
                
                if size(temp,1)==3 % wrong orientation?
                    temp=temp'; % transpose it
                end;
                
                cur_val=cur_val+temp*delta_t;
            end;
            
            if cov_const % time-invariant covariance matrix?
                rand_vec=random('norm',0,1,3,1); % independent noise
                rand_vec=sqrtm_cov_cur*rand_vec; % correlated noise
                cur_val=cur_val+rand_vec'; % noise part
            else
                temp=feval(inp_cov,i*delta_t); % get current covariance matrix
                
                if (size(temp,1)~=3)||(size(temp,2)~=3)
                    error('TRAJ_MSPRT_FB_3D_3B_SIM: The covariance matrix returned by a function must be a 3-by-3 matrix!');
                end;
                
                if (temp(1,1)<=0)||(temp(2,2)<=0)||(temp(3,3)<=0)
                    error('TRAJ_MSPRT_FB_3D_3B_SIM: Algorithm stopped due to a non-positive variance!');
                end;
                
                if det(temp)==0
                    error('TRAJ_MSPRT_FB_3D_3B_SIM: Algorithm stopped due to a singular covariance matrix!');
                end;
                
                if (det(temp)<0)||(temp(1,2)~=temp(2,1))||(temp(1,3)~=temp(3,1))||(temp(2,3)~=temp(3,2))
                    error('TRAJ_MSPRT_FB_3D_3B_SIM: Algorithm stopped due to an invalid covariance matrix!');
                end;
                
                sqrtm_cov_cur=sqrtm(temp*delta_t);
                rand_vec=random('norm',0,1,3,1); % independent noise
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
                
                if cur_val(3)<process_range(3,1)
                    cur_val(3)=process_range(3,1);
                end;
                
                if cur_val(3)>process_range(3,2)
                    cur_val(3)=process_range(3,2);
                end;
            end;
            
            % check for boundary crossings
            if feval(a_1,cur_val,i*delta_t) % first boundary crossed?
                if sel==1 % good trial?
                    total_sim=total_sim+1;
                    num_sim=num_sim+1;
                    traj_mat_1(num_sim,:)=cur_traj(1,:); % store the trial
                    traj_mat_2(num_sim,:)=cur_traj(2,:);
                    traj_mat_3(num_sim,:)=cur_traj(3,:);
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
                    traj_mat_3(num_sim,:)=cur_traj(3,:);
                else % forget trial
                    total_sim=total_sim+1;
                end;
                
                break; % We are done with this trial.
            end;
            
            if feval(a_3,cur_val,i*delta_t) % third boundary crossed?
                if sel==3 % good trial?
                    total_sim=total_sim+1;
                    num_sim=num_sim+1;
                    traj_mat_1(num_sim,:)=cur_traj(1,:); % store the trial
                    traj_mat_2(num_sim,:)=cur_traj(2,:);
                    traj_mat_3(num_sim,:)=cur_traj(3,:);
                else % forget trial
                    total_sim=total_sim+1;
                end;
                
                break; % We are done with this trial.
            end;
        end; % inner loop
    end; % if sel==0
    
    if ((stop_type==1)&&(total_sim==stop_val))||((stop_type==2)&&(num_sim==stop_val)) % Are we done?
        break; % leave the loop
    end;
end;

if num_sim % at least 1 effective trial?
    exp_traj(1,:)=mean(traj_mat_1);
    exp_traj(2,:)=mean(traj_mat_2);
    exp_traj(3,:)=mean(traj_mat_3);
end;

if num_sim>=2 % at least 2 effective trials?
    exp_var(1,:)=var(traj_mat_1);
    exp_var(2,:)=var(traj_mat_2);
    exp_var(3,:)=var(traj_mat_3);
end;

if runmed_width % filtering?
    exp_traj(1,:)=runmed(exp_traj(1,:),runmed_width,1,0);
    exp_traj(2,:)=runmed(exp_traj(2,:),runmed_width,1,0);
    exp_traj(3,:)=runmed(exp_traj(3,:),runmed_width,1,0);
    exp_var(1,:)=runmed(exp_var(1,:),runmed_width,1,0);
    exp_var(2,:)=runmed(exp_var(2,:),runmed_width,1,0);
    exp_var(3,:)=runmed(exp_var(3,:),runmed_width,1,0);
end;
