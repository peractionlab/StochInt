function [g_1,g_2,g_3,t_vec,ret_num_tr_states,g_oor]=ou_2d_3b_mar(ou_drift,ou_leak,ou_cov,ou_init,a_1,a_2,a_3,delta_t,covered_range,num_tr_states,stop_type,stop_val,oor_type,mapping_range,runmed_width)
% Markov chain approximation of the first passage time problem
% for a (time-variant) 2D Ornstein-Uhlenbeck process with 3
% (time-variant) boundaries. The calculation is based on a
% Random Walk with Gaussian Increments.
%
% J. Ditterich, 8/05
%
% [g_1,g_2,g_3,t_vec,ret_num_tr_states,g_oor] = ou_2d_3b_mar (ou_drift,ou_leak,ou_cov,ou_init,a_1,a_2,a_3,delta_t,
%                                                             covered_range,num_tr_states,stop_type,stop_val,
%                                                             oor_type[,mapping_range[,runmed_width]])
%
% g_1 is the first passage time density for the first boundary multiplied
%     by the probability of hitting the first boundary first, evaluated
%     at the times given in t_vec. If ou_init is a matrix, g_1 is a matrix.
%     Each row represents the FPT density for one of the initial values.
% g_2 is the first passage time density for the second boundary multiplied
%     by the probability of hitting the second boundary first, evaluated
%     at the times given in t_vec. If ou_init is a matrix, g_2 is a matrix.
%     Each row represents the FPT density for one of the initial values.
% g_3 is the first passage time density for the third boundary multiplied
%     by the probability of hitting the third boundary first, evaluated
%     at the times given in t_vec. If ou_init is a matrix, g_3 is a matrix.
%     Each row represents the FPT density for one of the initial values.
% ret_num_tr_states is the number of transient states. This is for the case
%                   when you use the built-in mechanism for choosing an optimal
%                   number of states and you would like to know what spatial
%                   resolution was chosen. ret_num_tr_states is a 1-by-2 vector.
% g_oor is the first passage time density for leaving the range defined by
%       covered_range without hitting a boundary. It has the same format as
%       g_1 and is only calculated if oor_type is set to 2. NaN is
%       returned for oor_type 1.
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
%        Please use one of the 1D functions for fully correlated processes.
% ou_init is the initial vector of the OU process. ou_init must have 2 columns.
%         A matrix can be passed, if you are interested in studying several inital vectors
%         at once. This doesn't cause additional computational overhead.
% a_1 defines the first absorbing boundary. a_1 is the name of a function,
%     which must return 1, if a certain location is located on or outside the boundary,
%     and 0, if a certain location is located inside the boundary, when called
%     with a 1-by-2 vector defining the location as the first and time as the
%     second argument.
% a_2 defines the second absorbing boundary. See a_1 for the format. The boundaries
%     should be defined in such a way that "boundary crossed" regions do not
%     overlap. If the algorithm detects that multiple boundaries have been crossed
%     in the same time step it assumes that each of the boundaries has been crossed
%     first in an equal number of cases.
% a_3 defines the third absorbing boundary. See a_1 for the format. The boundaries
%     should be defined in such a way that "boundary crossed" regions do not
%     overlap. If the algorithm detects that multiple boundaries have been crossed
%     in the same time step it assumes that each of the boundaries has been crossed
%     first in an equal number of cases.
% delta_t is the temporal step size.
% covered_range defines the range, which should be covered by the transient states.
%               It has to be a 2-by-2 matrix. The first row defines the lower and
%               the upper limit of the first dimension, the second row the lower and
%               the upper limit of the second dimension. Make sure that you define the
%               boundaries in such a way that there are always transient states outside
%               the boundaries. Otherwise the algorithm cannot determine the location
%               of the boundaries. On the other hand, the results are more accurate
%               if the boundaries are not too far from the edges of the range covered
%               by the transient states. The optimum would be a thin chain of transient
%               states in between the boundaries and the edges of the covered range.
% num_tr_states defines the number of transient states for each dimension used for the
%               discrete approximation of the problem. num_tr_states normally has to be
%               a vector of length 2. When passing 0 a built-in mechanism is used for
%               choosing an optimal number of states (based on a heuristic method).
% stop_type defines, what determines when the algorithm stops.
%           1 = time: The densities are calculated for times up to stop_val.
%           2 = deviation of the distribution from 1:
%               The densities are calculated until the deviation of the sum
%               of the integrals over the densities from 1 is smaller than
%               stop_val.
% stop_val defines the time or the error, which determines when the algorithm
%          will stop.
% oor_type defines what happens when the process leaves the range of values
%          covered by the transient states without having hit a boundary.
%          If oor_type is set to 1, the process cannot leave the area. Any
%          change component, which would lead to this situation, is reduced
%          accordingly.
%          If oor_type is set to 2, any attempt to leave the area results in
%          a transition to a fourth absorbing state, which means that you have
%          control over how often this happens.
% mapping_range is an optional parameter, which defines how big the matrices
%               used for constructing the transition matrix are and therefore
%               how sparse the transition matrix will be. The mapped range for
%               each dimension is MEAN +/- mapping_range * SD. The default value is 3.
% runmed_width is an optional parameter, which defines the width of a running median filter
%              applied to the output. It has to be an odd number. 0 deactivates the filter.
%              The default value is 0. Such a filter may be applied in the case of
%              time-variant boundaries to remove spikes. 15 would be reasonable value
%              to start with in this case.

% History:
% released on 7/24/07 as part of toolbox V 2.6

% State layout:
%
%      row #:
%                                                                                    -- covered_range(2,2)
% num_tr_states(2)                                               prod(num_tr_states)
%
%        .                                              .                                ^
%        .                                              .                                |
%        .                                              .                                | second dimension
%                                                                                        |                                                                             
%        2         num_tr_states(1)+1                  ...        2*num_tr_states(1)
%
%        1                 1                 2         ...         num_tr_states(1)
%                                                                                    -- covered_range(2,1)
%                 |                                                                 |
%          covered_range(1,1)                                                covered_range(1,2)
%
% column #:                1                 2         ...         num_tr_states(1)
%
%                                 --------------->
%                                 first dimension
%
% State # i is located in row # r=floor((i-1)/num_tr_states(1))+1 and column # c=i-(r-1)*num_tr_states(1).
% It is centered around the coordinates covered_range(1,1)+(c-.5)*diff(covered_range(1,:))/num_tr_states(1)
% and covered_range(2,1)+(r-.5)*diff(covered_range(2,:))/num_tr_states(2).

if nargin<15 % runmed_width not given?
    runmed_width=0; % default value
end;

if nargin<14 % mapping_range not given?
    mapping_range=3; % default value
end;

num_tr_states=round(num_tr_states);
runmed_width=round(runmed_width);

% Heuristic method for choosing an optimal number of states
if length(num_tr_states)==1 % scalar passed?
    if num_tr_states==0 % Should the number of states be determined automatically?
        if isnumeric(ou_cov) % constant covariance matrix?
            temp_var=diag(ou_cov);
        else % time-variant covariance matrix
            temp_var=diag(feval(ou_cov,0));
        end;
        
        temp1=6+diff(covered_range(1,:))^2/(45*temp_var(1)*delta_t); % optimal number of states for the first dimension
        temp2=6+diff(covered_range(2,:))^2/(45*temp_var(2)*delta_t); % optimal number of states for the second dimension
        temp=[temp1 temp2];
        % There is no theoretical foundation for this formula.
        % It was determined by playing around with the parameters a lot and comparing the result
        % to the numerical solution.
        
        % choose the next odd numbers
        temp2=round(temp);
        
        if round(temp2(1)/2)*2==temp2(1) % even?
            poss1=temp2(1)-1; % first possibility
            poss2=temp2(1)+1; % second possibility
            
            if abs(poss2-temp(1))<abs(poss1-temp(1)) % Which one is closer to the original guess?
                num_tr_states(1)=poss2; % choose the larger number
            else
                num_tr_states(1)=poss1; % choose the smaller number
            end;
        else % odd
            num_tr_states(1)=temp2(1);
        end;
        
        if round(temp2(2)/2)*2==temp2(2) % even?
            poss1=temp2(2)-1; % first possibility
            poss2=temp2(2)+1; % second possibility
            
            if abs(poss2-temp(2))<abs(poss1-temp(2)) % Which one is closer to the original guess?
                num_tr_states(2)=poss2; % choose the larger number
            else
                num_tr_states(2)=poss1; % choose the smaller number
            end;
        else % odd
            num_tr_states(2)=temp2(2);
        end;
        
        disp(['OU_2D_3B_MAR: The automatically chosen number of transient states for the first dimension is ' num2str(num_tr_states(1)) '.']);
        disp(['OU_2D_3B_MAR: The automatically chosen number of transient states for the second dimension is ' num2str(num_tr_states(2)) '.']);
        
        if (num_tr_states(1)<9)|(num_tr_states(2)<9) % number of states too small?
            error('OU_2D_3B_MAR: The results would be likely to be inaccurate. Please choose a smaller DELTA_T!');
        end;   
    else % number of states was given by the user
        if length(num_tr_states)~=2 % wrong vector length?
            error('OU_2D_3B_MAR: The number of transient states must be either a vector of length 2 or the scalar 0!');
        end;
        
        if size(num_tr_states,1)==2 % wrong orientation?
            num_tr_states=num_tr_states'; % transpose it
        end;
        
        if (num_tr_states(1)<0)|(num_tr_states(2)<0) % negative number of states?
            error('OU_2D_3B_MAR: There must be at least one transient state for both dimensions!');
        end;
        
        if (num_tr_states(1)<9)|(num_tr_states(2)<9)
            warning('OU_2D_3B_MAR: You should choose at least 9 transient states per dimension. Otherwise the results are very likely to be inaccurate!');
        end;   
    end;
end;

ret_num_tr_states=num_tr_states;

% Some checks
if isnumeric(ou_drift)&~((size(ou_drift,1)==1)&(size(ou_drift,2)==2))&~((size(ou_drift,1)==2)&(size(ou_drift,2)==1))
    error('OU_2D_3B_MAR: OU_DRIFT must be either a vector of length 2 or the name of a function!');
end;

if isnumeric(ou_drift)&(size(ou_drift,1)==2) % wrong orientation?
    ou_drift=ou_drift'; % transpose it
end;

if ou_leak<0
    error('OU_2D_3B_MAR: OU_LEAK must be a non-negative number!');
end;

if isnumeric(ou_cov)&((size(ou_cov,1)~=2)|(size(ou_cov,2)~=2))
    error('OU_2D_3B_MAR: OU_COV must either be a 2-by-2 matrix or the name of a function!');
end;

if isnumeric(ou_cov)&((ou_cov(1,1)<=0)|(ou_cov(2,2)<=0))
    error('OU_2D_3B_MAR: The main diagonal elements of OU_COV must be positive!');
end;

if isnumeric(ou_cov)&(det(ou_cov)==0)
    error('OU_2D_3B_MAR: The covariance matrix must not be singular!');
end;

if isnumeric(ou_cov)&((det(ou_cov)<0)|(ou_cov(1,2)~=ou_cov(2,1)))
    error('OU_2D_3B_MAR: Invalid covariance matrix!');
end;

if size(ou_init,2)~=2
    error('OU_2D_3B_MAR: OU_INIT must be a n-by-2 matrix!');
end;

if delta_t<=0
    error('OU_2D_3B_MAR: The time step must be a positive number!');
end;

if (size(covered_range,1)~=2)|(size(covered_range,2)~=2)
    error('OU_2D_3B_MAR: COVERED_RANGE must be a 2-by-2 matrix!');
end;

if (diff(covered_range(1,:))<0)|(diff(covered_range(2,:))<0) % screwed up range?
    error('OU_2D_3B_MAR: Invalid range!');
end;

if (stop_type~=1)&(stop_type~=2)
    error('OU_2D_3B_MAR: STOP_TYPE must be either 1 or 2!');
end;

if stop_val<=0
    error('OU_2D_3B_MAR: STOP_VAL must be a positive number!');
end;

if (oor_type~=1)&(oor_type~=2)
    error('OU_2D_3B_MAR: OOR_TYPE must be either 1 or 2!');
end;

if mapping_range<=0
    error('OU_2D_3B_MAR: MAPPING_RANGE must be a positive number!');
end;

if runmed_width<0
    error('OU_2D_3B_MAR: RUNMED_WIDTH must not be negative!');
end;

if runmed_width&(~mod(runmed_width,2))
    error('OU_2D_3B_MAR: RUNMED_WIDTH must be an odd number!');
end;

% Initialization
g_1=[];
g_2=[];
g_3=[];

if oor_type==1
    g_oor=nan;
else
    g_oor=[];
end;

t_vec=[];
k=1;
last_a_1=nan;
last_a_2=nan;
last_a_3=nan;
last_drift=nan;
last_cov=nan;

if isnumeric(ou_drift) % Is the drift time-invariant?
    drift_const=1;
    drift_cur=incval2stateinc(ou_drift*delta_t,covered_range,num_tr_states); % We need the increment in the state space per time step.
else
    drift_const=0;
end;

if isnumeric(ou_cov) % Is the covariance matrix time-invariant?
    cov_const=1;
    cov_cur=cov_scale(ou_cov*delta_t,covered_range,num_tr_states); % scale covariance matrix for Markov coordinate system
else
    cov_const=0;
end;

ou_init=val2state(ou_init,covered_range,num_tr_states); % We need the initial values in the state space.

current_mat=speye(prod(num_tr_states)); % start with a sparse identity matrix

% Further checks
ind=find(ou_init==0);

if ~isempty(ind) % Are there initial values, which do not map to transient states?
    error('OU_2D_3B_MAR: All initial values must map to transient states!');
end;

% Loop
while (1)
    t=k*delta_t; % current time
    
    if oor_type==1 % no fourth absorbing state
        if (((stop_type==1)&(t>stop_val))|((stop_type==2)&((1-(sum(g_1')+sum(g_2')+sum(g_3'))*delta_t)<stop_val))) % enough?
            break;
        end;
    else % fourth absorbing state
        if (((stop_type==1)&(t>stop_val))|((stop_type==2)&((1-(sum(g_1')+sum(g_2')+sum(g_3')+sum(g_oor'))*delta_t)<stop_val))) % enough?
            break;
        end;
    end;
    
    % check first boundary
    a_1_cur=sparse(1,prod(num_tr_states));
    
    for i=1:prod(num_tr_states)
        a_1_cur(i)=feval(a_1,state2val(i,covered_range,num_tr_states),t); % check if location of the transient state is inside or outside the boundary
    end;
    
    if nnz(a_1_cur)==0 % no transient states outside or on boundary?
        warning('OU_2D_3B_MAR: First boundary out of range!');
    end;
    
    % check second boundary
    a_2_cur=sparse(1,prod(num_tr_states));
    
    for i=1:prod(num_tr_states)
        a_2_cur(i)=feval(a_2,state2val(i,covered_range,num_tr_states),t); % check if location of the transient state is inside or outside the boundary
    end;
    
    if nnz(a_2_cur)==0 % no transient states outside or on boundary?
        warning('OU_2D_3B_MAR: Second boundary out of range!');
    end;
    
    % check third boundary
    a_3_cur=sparse(1,prod(num_tr_states));
    
    for i=1:prod(num_tr_states)
        a_3_cur(i)=feval(a_3,state2val(i,covered_range,num_tr_states),t); % check if location of the transient state is inside or outside the boundary
    end;
    
    if nnz(a_3_cur)==0 % no transient states outside or on boundary?
        warning('OU_2D_3B_MAR: Third boundary out of range!');
    end;
    
    if ~drift_const % Do we have to update the current drift?
        temp=feval(ou_drift,t); % get current drift
        
        if ~((size(temp,1)==1)&(size(temp,2)==2))&~((size(temp,1)==2)&(size(temp,2)==1))
            error('OU_2D_3B_MAR: The drift returned by a function must be a vector of length 2!');
        end;
        
        if size(temp,1)==2 % wrong orientation?
            temp=temp'; % transpose it
        end;
        
        drift_cur=incval2stateinc(temp*delta_t,covered_range,num_tr_states);
    end;
    
    if ~cov_const % Do we have to update the current variance?
        temp=feval(ou_cov,t);
        
        if (size(temp,1)~=2)|(size(temp,2)~=2)
            error('OU_2D_3B_MAR: The covariance matrix returned by a function must be a 2-by-2 matrix!');
        end;
        
        if (temp(1,1)<=0)|(temp(2,2)<=0)
            warning('OU_2D_3B_MAR: Algorithm stopped due to a non-positive variance!');
            break;
        end;
        
        if det(temp)==0
            warning('OU_2D_3B_MAR: Algorithm stopped due to a singular covariance matrix!');
        end;
        
        if (det(temp)<0)|(temp(1,2)~=temp(2,1))
            warning('OU_2D_3B_MAR: Algorithm stopped due to an invalid covariance matrix!');
        end;
        
        cov_cur=cov_scale(temp*delta_t,covered_range,num_tr_states);
    end;
    
    % Do we have to construct a new transition matrix?
    if nnz(a_1_cur~=last_a_1)|nnz(a_2_cur~=last_a_2)|nnz(a_3_cur~=last_a_3)|sum(drift_cur~=last_drift)|sum(sum(cov_cur~=last_cov))
        last_a_1=a_1_cur;
        last_a_2=a_2_cur;
        last_a_3=a_3_cur;
        last_drift=drift_cur;
        last_cov=cov_cur;

        if ou_leak==0 % Wiener process? -> The transition "block" has to be calculated only once.
            m_cur=drift_cur; % The drift is the only deterministic component.
            lower_bound_1=floor(m_cur(1)-mapping_range*sqrt(cov_cur(1,1)));
            upper_bound_1=ceil(m_cur(1)+mapping_range*sqrt(cov_cur(1,1)));
            lower_bound_2=floor(m_cur(2)-mapping_range*sqrt(cov_cur(2,2)));
            upper_bound_2=ceil(m_cur(2)+mapping_range*sqrt(cov_cur(2,2)));
            x1=[lower_bound_1:upper_bound_1];
            x2=[lower_bound_2:upper_bound_2];
            x1_=[lower_bound_1-.25:.5:upper_bound_1+.25]; % used for getting a better approximation
            x2_=[lower_bound_2-.25:.5:upper_bound_2+.25]; % of the CDF
            [x1_,x2_]=meshgrid(x1_,x2_); % construct the grid for evaluating the PDF
            y_=normpdf2(x1_,x2_,m_cur,cov_cur); % evaluate the PDF
            y_=y_/sum(sum(y_)); % normalize the transition probabilities
            
            for i=1:length(x2)
                for j=1:length(x1)
                    y(i,j)=sum(sum(y_([(i-1)*2+1:i*2],[(j-1)*2+1:j*2]))); % 4 points in the calculated PDF contribute to each value in the transition "block".
                end;
            end;
        end;      
        
        Q=sparse(prod(num_tr_states),prod(num_tr_states)); % This will be the transition matrix.
        
        if oor_type==1 % 3 absorbing states
            R=sparse(prod(num_tr_states),3); % These will be the transition vectors to the absorbing states.
        else % 4 absorbing states
            R=sparse(prod(num_tr_states),4);
        end;
        
        for from_state=1:prod(num_tr_states) % This is the "from" state.
            if ou_leak~=0 % OU process? -> The transition "block" has to be recalculated for each row.
                m_cur=drift_cur+leak2stateinc(from_state,ou_leak,delta_t,covered_range,num_tr_states); % Take the "leakage" of the integrator into account!
                lower_bound_1=floor(m_cur(1)-mapping_range*sqrt(cov_cur(1,1)));
                upper_bound_1=ceil(m_cur(1)+mapping_range*sqrt(cov_cur(1,1)));
                lower_bound_2=floor(m_cur(2)-mapping_range*sqrt(cov_cur(2,2)));
                upper_bound_2=ceil(m_cur(2)+mapping_range*sqrt(cov_cur(2,2)));
                x1=[lower_bound_1:upper_bound_1];
                x2=[lower_bound_2:upper_bound_2];
                x1_=[lower_bound_1-.25:.5:upper_bound_1+.25]; % used for getting a better approximation
                x2_=[lower_bound_2-.25:.5:upper_bound_2+.25]; % of the CDF
                [x1_,x2_]=meshgrid(x1_,x2_); % construct the grid for evaluating the PDF
                y_=normpdf2(x1_,x2_,m_cur,cov_cur); % evaluate the PDF
                y_=y_/sum(sum(y_)); % normalize the transition probabilities
                
                for i=1:length(x2)
                    for j=1:length(x1)
                        y(i,j)=sum(sum(y_([(i-1)*2+1:i*2],[(j-1)*2+1:j*2]))); % 4 points in the calculated PDF contribute to each value in the transition "block".
                    end;
                end;
            end;
            
            [from_row,from_col]=state2rowcol(from_state,num_tr_states); % get the location of the "from" state on the grid of transient states         
            
            for row_ind=1:length(x2) % walk through the transition "block"
                for col_ind=1:length(x1)
                    to_row=from_row+x2(row_ind); % This is the row of the "to" state.
                    to_col=from_col+x1(col_ind); % This is the column of the "to" state.
                    to_state=rowcol2state(to_row,to_col,num_tr_states); % This is the "to" state.
                    
                    if ~to_state % out of range?
                        b1_crossed=feval(a_1,rowcol2val(to_row,to_col,covered_range,num_tr_states),t); % check if the destination location is outside the first boundary
                        b2_crossed=feval(a_2,rowcol2val(to_row,to_col,covered_range,num_tr_states),t); % check if the destination location is outside the second boundary
                        b3_crossed=feval(a_3,rowcol2val(to_row,to_col,covered_range,num_tr_states),t); % check if the destination location is outside the third boundary
                        if b1_crossed&~b2_crossed&~b3_crossed % first boundary will be crossed, but not second or third boundary
                            R(from_state,1)=R(from_state,1)+y(row_ind,col_ind); % update the probability of crossing the first boundary
                        elseif ~b1_crossed&b2_crossed&~b3_crossed % second boundary will be crossed, but no first or third boundary
                            R(from_state,2)=R(from_state,2)+y(row_ind,col_ind); % update the probability of crossing the second boundary
                        elseif ~b1_crossed&~b2_crossed&b3_crossed % third boundary will be crossed, but no first or second boundary
                            R(from_state,3)=R(from_state,3)+y(row_ind,col_ind); % update the probability of crossing the third boundary
                        elseif b1_crossed&b2_crossed&~b3_crossed % first and second boundary would be crossed, but not the third boundary
                            R(from_state,1)=R(from_state,1)+.5*y(row_ind,col_ind); % This is the critical case.
                            % I have decided to add half of the probability to the probability of crossing the first boundary
                            % and the other half to the probability of crossing the second boundary.
                            R(from_state,2)=R(from_state,2)+.5*y(row_ind,col_ind);
                        elseif b1_crossed&~b2_crossed&b3_crossed % first and third boundary would be crossed, but not the second boundary
                            R(from_state,1)=R(from_state,1)+.5*y(row_ind,col_ind);
                            R(from_state,3)=R(from_state,3)+.5*y(row_ind,col_ind);
                        elseif ~b1_crossed&b2_crossed&b3_crossed % second and third boundary would be crossed, but not the first boundary
                            R(from_state,2)=R(from_state,2)+.5*y(row_ind,col_ind);
                            R(from_state,3)=R(from_state,3)+.5*y(row_ind,col_ind);
                        elseif b1_crossed&b2_crossed&b3_crossed % all boundaries would be crossed
                            R(from_state,1)=R(from_state,1)+y(row_ind,col_ind)/3;
                            R(from_state,2)=R(from_state,2)+y(row_ind,col_ind)/3;
                            R(from_state,3)=R(from_state,3)+y(row_ind,col_ind)/3;
                        else % We are leaving the range without crossing any boundary.                     
                            if oor_type==2 % fourth absorbing state?
                                R(from_state,4)=R(from_state,4)+y(row_ind,col_ind); % update the probability of leaving the range without crossing a boundary
                            else % only 3 absorbing states
                                temp_state=nearest_state(to_row,to_col,covered_range,num_tr_states); % stay inside the range
                                Q(from_state,temp_state)=Q(from_state,temp_state)+y(row_ind,col_ind); % update the probability of jumping to the nearest state inside the range
                            end;
                        end;
                    else % within range
                        if a_1_cur(to_state)&~a_2_cur(to_state)&~a_3_cur(to_state) % first boundary will be crossed, but not second or third boundary
                            R(from_state,1)=R(from_state,1)+y(row_ind,col_ind); % update the probability of crossing the first boundary
                        elseif ~a_1_cur(to_state)&a_2_cur(to_state)&~a_3_cur(to_state) % second boundary will be crossed, but not first or third boundary
                            R(from_state,2)=R(from_state,2)+y(row_ind,col_ind); % update the probability of crossing the second boundary
                        elseif ~a_1_cur(to_state)&~a_2_cur(to_state)&a_3_cur(to_state) % third boundary will be crossed, but not first or second boundary
                            R(from_state,3)=R(from_state,3)+y(row_ind,col_ind); % update the probability of crossing the third boundary
                        elseif a_1_cur(to_state)&a_2_cur(to_state)&~a_3_cur(to_state) % first and second boundary would be crossed, but not the third boundary
                            R(from_state,1)=R(from_state,1)+.5*y(row_ind,col_ind); % This is the critical case.
                            % I have decided to add half of the probability to the probability of crossing the first boundary
                            % and the other half to the probability of crossing the second boundary.
                            R(from_state,2)=R(from_state,2)+.5*y(row_ind,col_ind);
                        elseif a_1_cur(to_state)&~a_2_cur(to_state)&a_3_cur(to_state) % first and third boundary would be crossed, but not the second boundary
                            R(from_state,1)=R(from_state,1)+.5*y(row_ind,col_ind);
                            R(from_state,3)=R(from_state,3)+.5*y(row_ind,col_ind);
                        elseif ~a_1_cur(to_state)&a_2_cur(to_state)&a_3_cur(to_state) % second and third boundary would be crossed, but not the first boundary
                            R(from_state,2)=R(from_state,2)+.5*y(row_ind,col_ind);
                            R(from_state,3)=R(from_state,3)+.5*y(row_ind,col_ind);
                        elseif a_1_cur(to_state)&a_2_cur(to_state)&a_3_cur(to_state) % all boundaries would be crossed
                            R(from_state,1)=R(from_state,1)+y(row_ind,col_ind)/3;
                            R(from_state,2)=R(from_state,2)+y(row_ind,col_ind)/3;
                            R(from_state,3)=R(from_state,3)+y(row_ind,col_ind)/3;
                        else % no boundary will be crossed
                            Q(from_state,to_state)=Q(from_state,to_state)+y(row_ind,col_ind); % update transition matrix
                        end;                  
                    end;               
                end; % for col_ind
            end; % for row_ind
        end; % for all "from" states             
    end; % new transition matrix?
    
    current_probs=current_mat*R;
    
    g_1=[g_1 current_probs(ou_init,1)/delta_t];
    g_2=[g_2 current_probs(ou_init,2)/delta_t];
    g_3=[g_3 current_probs(ou_init,3)/delta_t];
    
    if oor_type==2 % 4 absorbing states?
        g_oor=[g_oor current_probs(ou_init,4)/delta_t];
    end;
    
    t_vec=[t_vec t];
    current_mat=current_mat*Q;
    k=k+1;
end; % loop

if runmed_width % filtering requested?
    % Run the output through a running median filter to remove spikes.
    if length(g_1(1,:))<runmed_width % inappropriate filter width?
        error('OU_2D_3B_MAR: Inappropriate filter width!');
    end;
    
    for i=1:length(ou_init)
        g_1(i,:)=runmed(g_1(i,:),runmed_width,1,0);
        g_2(i,:)=runmed(g_2(i,:),runmed_width,1,0);
        g_3(i,:)=runmed(g_3(i,:),runmed_width,1,0);
    end;
end;

% Local functions
function state=val2state(val,covered_range,num_tr_states)
% val must be a n-by-2 matrix. state is a column vector.
% The transient states are numbered from 1 to prod(num_tr_states).
% 0 is returned if the coordinates are outside the range of the transient states.
col=(val(:,1)-covered_range(1,1)+diff(covered_range(1,:))/(2*num_tr_states(1)))/diff(covered_range(1,:))*num_tr_states(1);
row=(val(:,2)-covered_range(2,1)+diff(covered_range(2,:))/(2*num_tr_states(2)))/diff(covered_range(2,:))*num_tr_states(2);
ind=find(col<.5);
col(ind)=0;
ind=find(col>num_tr_states(1)+.5);
col(ind)=0;
ind=find(row<.5);
row(ind)=0;
ind=find(row>num_tr_states(2)+.5);
row(ind)=0;
col=round(col);
row=round(row);

if (col==0)|(row==0) % out of range?
    state=0;
else
    state=(row-1)*num_tr_states(1)+col;
end;

function val=state2val(state,covered_range,num_tr_states)
% state must be a column vector. val is a n-by-2 matrix.
row=floor((state-1)/num_tr_states(1))+1;
col=state-(row-1)*num_tr_states(1);
val1=covered_range(1,1)+(col-.5)*diff(covered_range(1,:))/num_tr_states(1);
val2=covered_range(2,1)+(row-.5)*diff(covered_range(2,:))/num_tr_states(2);
val=[val1 val2];

function state_inc=incval2stateinc(inc_val,covered_range,num_tr_states)
% inc_val must be a n-by-2 matrix. state_inc is a n-by-2 matrix.
delta_col=inc_val(:,1)/diff(covered_range(1,:))*num_tr_states(1);
delta_row=inc_val(:,2)/diff(covered_range(2,:))*num_tr_states(2);
state_inc=[delta_col delta_row];

function state_inc=leak2stateinc(state,ou_leak,delta_t,covered_range,num_tr_states)
% state must be a column vector. state_inc is a n-by-2 matrix.
state_inc=incval2stateinc(-state2val(state,covered_range,num_tr_states)*ou_leak*delta_t,covered_range,num_tr_states);

function new_cov=cov_scale(covm,covered_range,num_tr_states)
% scale a covariance matrix for the Markov coordinate system
new_cov=[];
div1=diff(covered_range(1,:))/num_tr_states(1);
div2=diff(covered_range(2,:))/num_tr_states(2);
new_cov(1,1)=covm(1,1)/div1^2;
new_cov(1,2)=covm(1,2)/(div1*div2);
new_cov(2,1)=covm(2,1)/(div1*div2);
new_cov(2,2)=covm(2,2)/div2^2;

function [row,col]=state2rowcol(state,num_tr_states)
% state, row, and col are column vectors.
row=floor((state-1)/num_tr_states(1))+1;
col=state-(row-1)*num_tr_states(1);

function state=rowcol2state(row,col,num_tr_states)
% row, col, and state are column vectors.
% 0 is returned if the location is out of range.
state=(row-1)*num_tr_states(1)+col;

if (col<1)|(col>num_tr_states(1))|(row<1)|(row>num_tr_states(2)) % out of range?
    state=0;
end;

function val=rowcol2val(row,col,covered_range,num_tr_states)
% row and col must be column vectors. val is a n-by-2 matrix.
val1=covered_range(1,1)+(col-.5)*diff(covered_range(1,:))/num_tr_states(1);
val2=covered_range(2,1)+(row-.5)*diff(covered_range(2,:))/num_tr_states(2);
val=[val1 val2];

function state=nearest_state(row,col,covered_range,num_tr_states)
% get the nearest transient state for an out-of-range location
if (row>num_tr_states(2))&(col<1) % upper left quadrant
    state=rowcol2state(num_tr_states(2),1,num_tr_states);
elseif (row>num_tr_states(2))&(col>=1)&(col<=num_tr_states(1)) % upper part
    state=rowcol2state(num_tr_states(2),col,num_tr_states);
elseif (row>num_tr_states(2))&(col>num_tr_states(1)) % upper right quadrant
    state=rowcol2state(num_tr_states(2),num_tr_states(1),num_tr_states);
elseif (row>=1)&(row<=num_tr_states(2))&(col>num_tr_states(1)) % right part
    state=rowcol2state(row,num_tr_states(1),num_tr_states);
elseif (row<1)&(col>num_tr_states(1)) % lower right quadrant
    state=num_tr_states(1);
elseif (row<1)&(col>=1)&(col<=num_tr_states(1)) % lower part
    state=col;
elseif (row<1)&(col<1) % lower left quadrant
    state=1;
else % left part
    state=rowcol2state(row,1,num_tr_states);
end;
