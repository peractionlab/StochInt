function [g_1,g_2,g_3,t_vec,sc]=ou_2d_3b_1d_2b_sim_sc(ou_means,ou_leak,ou_variances,a_1,a_2,a_3,b_1,b_2,delta_t,max_time,num_sim,runmed_width,report_incomplete)
% Estimate of the first passage time densities for a (time-variant)
% 2D Ornstein-Uhlenbeck process with 3 (time-variant) boundaries
% based on a simulation. Second choices are also reported based on starting
% a new 1D Ornstein-Uhlenbeck process with 2 (time-variant) boundaries to
% decide between the two remaining options. The OU processes are
% constructed from 3 independent Gaussian random processes x_1, x_2, and x_3:
% 
% The first dimension of the 2D OU process is given by integrating, potentially with leak,
% x_1-0.5*(x_2+x_3).
% The scond dimension is given by integrating x_2-0.5*(x_1+x_3).
%
% Once the first threshold crossing has occurred, corresponding to identifying
% either x_1, x_2, or x_3 as the process with the largest mean, a new 1D OU process is
% initiated as the integral, again potentially with leak, of the two
% remaining Gaussian random processes, i.e., either x_2-x_3, x_1-x_3, or
% x_1-x_2.
%
% J. Ditterich, 3/20
%
% [g_1,g_2,g_3,t_vec,sc] = ou_2d_3b_1d_2b_sim_sc (ou_means,ou_leak,ou_variances,a_1,a_2,a_3,b_1,b_2,delta_t,
%                                                 max_time,num_sim[,runmed_width[,report_incomplete]])
%
% g_1 is the first passage time density for the first boundary multiplied
%     by the probability of hitting the first boundary first, evaluated
%     at the times given in t_vec.
% g_2 is the first passage time density for the second boundary multiplied
%     by the probability of hitting the second boundary first, evaluated
%     at the times given in t_vec.
% g_3 is the first passage time density for the third boundary multiplied
%     by the probability of hitting the third boundary first, evaluated
%     at the times given in t_vec.
% sc is a 3x3 matrix of the probabilities of all possible combinations of
%    first and second choices. The row determines the first choice, the column
%    the second choice.
%
% ou_means is the vector of the means of 3 random processes feeding into the OU processes.
% ou_leak defines the "leakiness" of the integrator(s) and has to be a scalar.
%         The deterministic part of the stochastic differential equation is given by
%         ou_drift - ou_leak * current_value. A Wiener process can be studied
%         by setting ou_leak to 0.
% ou_variances is the vector of the variances of 3 random processes feeding into the OU processes.
% a_1 defines the first absorbing boundary for the 2D OU process. a_1 is the name of a function,
%     which must return 1, if a certain location is located on or outside the boundary,
%     and 0, if a certain location is located inside the boundary, when called
%     with a 1-by-2 vector defining the location as the first and time as the
%     second argument.
% a_2 defines the second absorbing boundary for the 2D OU process. See a_1 for the format. The boundaries
%     should be defined in such a way that "boundary crossed" regions do not
%     overlap. Since the algorithm checks the first boundary first, a crossing of
%     multiple boundaries in the same time step will be registered as a crossing of
%     the first boundary.
% a_3 defines the third absorbing boundary for the 2D OU process. See a_1 for the format. The boundaries
%     should be defined in such a way that "boundary crossed" regions do not
%     overlap. Since the algorithm checks the first boundary first, a crossing of
%     multiple boundaries in the same time step will be registered as a crossing of
%     the first boundary.
% b_1 defines the first absorbing boundary for the 1D OU process. b_1 is
%     the name of a function, which must return 1, if a certain location is
%     located on or outside the boundary, and 0, if a certain location is
%     located inside the boundary, when called with the location as the
%     first and time as the second argument. Alternatively, if additional
%     parameters need to be passed into the function, b_1 can be a cell
%     array with the first element being the name of the function and the
%     second element being the additional parameters, which are passed into
%     the function as a third argument.
% b_2 defines the second absorbing boundary for the 1D OU process. See b_1
%     for the format.
% delta_t is the temporal step size.
% max_time defines the maximum first passage time taken into account by the algorithm.
% num_sim is the number of simulations used for calculating the result.
% runmed_width is an optional parameter, which defines the width of a running median filter
%              applied to the output. It has to be an odd number. 0 deactivates the filter.
%              The default value is 0.
% report_incomplete is an optional parameter, which defines what happens
%                   when max_time is reached and no second threshold crossing has been
%                   registered yet. 0 means that no first threshold crossing is reported. 1 means that
%                   the first threshold crossing is reported with a decision time of
%                   max_time. The default value is 1. A second choice,
%                   obviously, cannot be reported in either case.

% History:
% derived on 3/3/20 from OU_2D_3B_SIM
% released on 3/29/20 as part of toolbox V 2.9

if nargin<13 % report_incomplete not given?
    report_incomplete=1; % default value
end;

if nargin<12 % runmed_width not given?
    runmed_width=0; % default value
end;

num_sim=round(num_sim);
runmed_width=round(runmed_width);

% Some checks
if length(ou_means)~=3
    error('OU_2D_3B_1D_2B_SIM_SC: OU_MEANS must be a vector of length 3!');
end;

if ou_leak<0
    error('OU_2D_3B_1D_2B_SIM_SC: OU_LEAK must be a non-negative number!');
end;

if length(ou_variances)~=3
    error('OU_2D_3B_1D_2B_SIM_SC: OU_VARIANCES must be a vector of length 3!');
end;

if delta_t<=0
    error('OU_2D_3B_1D_2B_SIM_SC: The time step must be a positive number!');
end;

if max_time<=0
    error('OU_2D_3B_1D_2B_SIM_SC: MAX_TIME must be a positive number!');
end;

if num_sim<1
    error('OU_2D_3B_1D_2B_SIM_SC: A minimum of 1 simulation is required!');
end;

if runmed_width<0
    error('OU_2D_3B_1D_2B_SIM_SC: RUNMED_WIDTH must not be negative!');
end;

if runmed_width&(~mod(runmed_width,2))
    error('OU_2D_3B_1D_2B_SIM_SC: RUNMED_WIDTH must be an odd number!');
end;

if (report_incomplete~=0)&(report_incomplete~=1)
    error('OU_2D_3B_1D_2B_SIM_SC: REPORT_INCOMPLETE has to be either 0 or 1!');
end;

% Initialization
vec_length=floor(max_time/delta_t);
t_vec=[delta_t:delta_t:vec_length*delta_t];
g_1=zeros(1,vec_length);
g_2=zeros(1,vec_length);
g_3=zeros(1,vec_length);
sc=zeros(3,3);

ou_drift=[ou_means(1)-.5*ou_means(2)-.5*ou_means(3) ou_means(2)-.5*ou_means(1)-.5*ou_means(3)];
my_cov=-.5*ou_variances(1)-.5*ou_variances(2)+.25*ou_variances(3);
ou_cov=[ou_variances(1)+.25*ou_variances(2)+.25*ou_variances(3) my_cov;my_cov .25*ou_variances(1)+ou_variances(2)+.25*ou_variances(3)];
drift_cur=ou_drift*delta_t;
sqrtm_cov_cur=sqrtm(ou_cov*delta_t);

% Loop
for k=1:num_sim
    boundary_tested=0;
    first_crossing=0;
    index_first_crossing=0;
    
    % create a trajectory with a length of vec_length
    if ou_leak==0 % In this case we can do it in a single step.
        rand_vec=random('norm',0,1,2,vec_length); % independent noise
        rand_vec=repmat(drift_cur',1,vec_length)+sqrtm_cov_cur*rand_vec; % drift & correlated noise
        cur_traj=(tril(ones(vec_length,vec_length))*rand_vec')'; % integration
    else % separate random calls necessary
        boundary_tested=1;
        cur_val=[0 0]; % start with the initial value
        cur_traj=[];
        
        for i=1:vec_length
            cur_val=cur_val*(1-ou_leak*delta_t); % leaky integrator part
            cur_val=cur_val+drift_cur; % drift part
            rand_vec=random('norm',0,1,2,1); % independent noise
            rand_vec=sqrtm_cov_cur*rand_vec; % correlated noise
            cur_val=cur_val+rand_vec'; % noise part
                        
            if feval(a_1,cur_val,i*delta_t) % first boundary crossed?
                first_crossing=1; % register boundary crossing
                index_first_crossing=i;
                break; % We no longer have to calculate the rest of the trajectory.
            end;
            
            if feval(a_2,cur_val,i*delta_t) % second boundary crossed?
                first_crossing=2; % register boundary crossing
                index_first_crossing=i;
                break; % We no longer have to calculate the rest of the trajectory.
            end;
            
            if feval(a_3,cur_val,i*delta_t) % third boundary crossed?
                first_crossing=3; % register boundary crossing\
                index_first_crossing=i;
                break; % We no longer have to calculate the rest of the trajectory.
            end;
        end; % for i
    end;
    
    % check for boundary crossing
    if ~boundary_tested    
        for i=1:vec_length
            if feval(a_1,cur_traj(:,i)',i*delta_t) % first boundary crossed?
                first_crossing=1; % register boundary crossing
                index_first_crossing=i;
                break;
            end;

            if feval(a_2,cur_traj(:,i)',i*delta_t) % second boundary crossed?
                first_crossing=2; % register boundary crossing
                index_first_crossing=i;
                break;
            end;

            if feval(a_3,cur_traj(:,i)',i*delta_t) % third boundary crossed?
                first_crossing=3; % register boundary crossing
                index_first_crossing=i;
                break;
            end;
        end;
    end;

    % 1D OU process for second crossing
    if first_crossing % make sure that we had a first crossing
        remaining_vec_length=vec_length-index_first_crossing;
        
        switch first_crossing
            case 1
                ou_drift_2=ou_means(2)-ou_means(3);
                ou_var_2=ou_variances(2)+ou_variances(3);
            case 2
                ou_drift_2=ou_means(1)-ou_means(3);
                ou_var_2=ou_variances(1)+ou_variances(3);
            case 3
                ou_drift_2=ou_means(1)-ou_means(2);
                ou_var_2=ou_variances(1)+ou_variances(2);
        end;
        
        drift_cur_2=ou_drift_2*delta_t;
        std_cur_2=sqrt(ou_var_2*delta_t);
        
        boundary_tested=0;
        second_crossing=0;
        
        % create a trajectory with a length of remaining_vec_length
        if ou_leak==0 % In this case we can do it in a single step.
            rand_vec=random('norm',0,1,1,remaining_vec_length); % independent noise
            rand_vec=repmat(drift_cur_2,1,remaining_vec_length)+std_cur_2*rand_vec; % drift & correlated noise
            cur_traj=(tril(ones(remaining_vec_length,remaining_vec_length))*rand_vec')'; % integration
        else % separate random calls necessary
            boundary_tested=1;
            cur_val=0; % start with the initial value
            cur_traj=[];
            
            for i=1:remaining_vec_length
                cur_val=cur_val*(1-ou_leak*delta_t); % leaky integrator part
                cur_val=cur_val+drift_cur_2; % drift part
                cur_val=cur_val+std_cur_2*random('norm',0,1,1,1); % noise part
                
                if iscell(b_1) % Is b_1 a cell array?
                    if feval(b_1{1},cur_val,i*delta_t,b_1{2}) % first boundary crossed?
                        total_index=index_first_crossing+i;
                        
                        switch first_crossing
                            case 1
                                second_crossing=2;
                                g_1(total_index)=g_1(total_index)+1;
                                sc(1,2)=sc(1,2)+1;
                            case 2
                                second_crossing=1;
                                g_2(total_index)=g_2(total_index)+1;
                                sc(2,1)=sc(2,1)+1;
                            case 3
                                second_crossing=1;
                                g_3(total_index)=g_3(total_index)+1;
                                sc(3,1)=sc(3,1)+1;
                        end;
                        
                        break; % We no longer have to calculate the rest of the trajectory.
                    end;
                else
                    if feval(b_1,cur_val,i*delta_t) % first boundary crossed?
                        total_index=index_first_crossing+i;
                        
                        switch first_crossing
                            case 1
                                second_crossing=2;
                                g_1(total_index)=g_1(total_index)+1;
                                sc(1,2)=sc(1,2)+1;
                            case 2
                                second_crossing=1;
                                g_2(total_index)=g_2(total_index)+1;
                                sc(2,1)=sc(2,1)+1;
                            case 3
                                second_crossing=1;
                                g_3(total_index)=g_3(total_index)+1;
                                sc(3,1)=sc(3,1)+1;
                        end;
                        
                        break; % We no longer have to calculate the rest of the trajectory.
                    end;
                end;
                
                if iscell(b_2) % Is b_2 a cell array?
                    if feval(b_2{1},cur_val,i*delta_t,b_2{2}) % second boundary crossed?
                        total_index=index_first_crossing+i;
                        
                        switch first_crossing
                            case 1
                                second_crossing=3;
                                g_1(total_index)=g_1(total_index)+1;
                                sc(1,3)=sc(1,3)+1;
                            case 2
                                second_crossing=3;
                                g_2(total_index)=g_2(total_index)+1;
                                sc(2,3)=sc(2,3)+1;
                            case 3
                                second_crossing=2;
                                g_3(total_index)=g_3(total_index)+1;
                                sc(3,2)=sc(3,2)+1;
                        end;
                        
                        break; % We no longer have to calculate the rest of the trajectory.
                    end;
                else
                    if feval(b_2,cur_val,i*delta_t) % second boundary crossed?
                        total_index=index_first_crossing+i;
                        
                        switch first_crossing
                            case 1
                                second_crossing=3;
                                g_1(total_index)=g_1(total_index)+1;
                                sc(1,3)=sc(1,3)+1;
                            case 2
                                second_crossing=3;
                                g_2(total_index)=g_2(total_index)+1;
                                sc(2,3)=sc(2,3)+1;
                            case 3
                                second_crossing=2;
                                g_3(total_index)=g_3(total_index)+1;
                                sc(3,2)=sc(3,2)+1;
                        end;
                        
                        break; % We no longer have to calculate the rest of the trajectory.
                    end;
                end;
            end; % for i
        end;
        
        % check for boundary crossing
        if ~boundary_tested
            for i=1:remaining_vec_length
                if iscell(b_1) % Is b_1 a cell array?
                    if feval(b_1{1},cur_traj(:,i)',i*delta_t,b_1{2}) % first boundary crossed?
                        total_index=index_first_crossing+i;
                        
                        switch first_crossing
                            case 1
                                second_crossing=2;
                                g_1(total_index)=g_1(total_index)+1;
                                sc(1,2)=sc(1,2)+1;
                            case 2
                                second_crossing=1;
                                g_2(total_index)=g_2(total_index)+1;
                                sc(2,1)=sc(2,1)+1;
                            case 3
                                second_crossing=1;
                                g_3(total_index)=g_3(total_index)+1;
                                sc(3,1)=sc(3,1)+1;
                        end;
                        
                        break;
                    end;
                else
                    if feval(b_1,cur_traj(:,i)',i*delta_t) % first boundary crossed?
                        total_index=index_first_crossing+i;
                        
                        switch first_crossing
                            case 1
                                second_crossing=2;
                                g_1(total_index)=g_1(total_index)+1;
                                sc(1,2)=sc(1,2)+1;
                            case 2
                                second_crossing=1;
                                g_2(total_index)=g_2(total_index)+1;
                                sc(2,1)=sc(2,1)+1;
                            case 3
                                second_crossing=1;
                                g_3(total_index)=g_3(total_index)+1;
                                sc(3,1)=sc(3,1)+1;
                        end;
                        
                        break;
                    end;
                end;
                
                if iscell(b_2) % Is b_2 a cell array?
                    if feval(b_2{1},cur_traj(:,i)',i*delta_t,b_2{2}) % second boundary crossed?
                        total_index=index_first_crossing+i;
                        
                        switch first_crossing
                            case 1
                                second_crossing=3;
                                g_1(total_index)=g_1(total_index)+1;
                                sc(1,3)=sc(1,3)+1;
                            case 2
                                second_crossing=3;
                                g_2(total_index)=g_2(total_index)+1;
                                sc(2,3)=sc(2,3)+1;
                            case 3
                                second_crossing=2;
                                g_3(total_index)=g_3(total_index)+1;
                                sc(3,2)=sc(3,2)+1;
                        end;
                        
                        break;
                    end;
                else
                    if feval(b_2,cur_traj(:,i)',i*delta_t) % second boundary crossed?
                        total_index=index_first_crossing+i;
                        
                        switch first_crossing
                            case 1
                                second_crossing=3;
                                g_1(total_index)=g_1(total_index)+1;
                                sc(1,3)=sc(1,3)+1;
                            case 2
                                second_crossing=3;
                                g_2(total_index)=g_2(total_index)+1;
                                sc(2,3)=sc(2,3)+1;
                            case 3
                                second_crossing=2;
                                g_3(total_index)=g_3(total_index)+1;
                                sc(3,2)=sc(3,2)+1;
                        end;
                        
                        break;
                    end;
                end;
            end;
        end;
        
        if (~second_crossing)&&report_incomplete % no second crossing? report it?
            switch first_crossing
                case 1
                    g_1(vec_length)=g_1(vec_length)+1;
                case 2
                    g_2(vec_length)=g_2(vec_length)+1;
                case 3
                    g_3(vec_length)=g_3(vec_length)+1;
            end;
        end;
    end;
end;

g_1=g_1/num_sim/delta_t;
g_2=g_2/num_sim/delta_t;
g_3=g_3/num_sim/delta_t;
    
if runmed_width % filtering?
    g_1=runmed(g_1,runmed_width,1,0);
    g_2=runmed(g_2,runmed_width,1,0);
    g_3=runmed(g_3,runmed_width,1,0);
end;

sc=sc/sum(sum(sc));
