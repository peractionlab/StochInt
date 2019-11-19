function [g_upper,g_lower,t_vec,ret_num_tr_states,g_dropout]=ou_dropout_1d_2b_mar(ou_drift,ou_leak,ou_var,ou_init,a_upper,a_lower,delta_t,covered_range,num_tr_states,stop_type,stop_val,dropout_rate,mapping_range,runmed_width)
% Markov chain approximation of the first passage time problem
% for a (time-variant) 1D Ornstein-Uhlenbeck process with 2
% (time-variant) boundaries. The risk of aborting the process
% without hitting a boundary is given by a dropout rate.
% The calculation is based on a Random Walk with Gaussian Increments.
%
% J. Ditterich, 1/03
%
% [g_upper,g_lower,t_vec,ret_num_tr_states,g_dropout] = ou_dropout_1d_2b_mar (ou_drift,ou_leak,ou_var,ou_init,a_upper,a_lower,delta_t,
%                                                                             covered_range,num_tr_states,stop_type,stop_val,
%                                                                             dropout_rate[,mapping_range[,runmed_width]])
%
% g_upper is the first passage time density for the upper boundary multiplied
%         by the probability of hitting the upper boundary first, evaluated
%         at the times given in t_vec. If ou_init is a vector, g_upper is a matrix.
%         Each row represents the FPT density for one of the initial values.
% g_lower is the first passage time density for the lower boundary multiplied
%         by the probability of hitting the lower boundary first, evaluated
%         at the times given in t_vec. If ou_init is a vector, g_lower is a matrix.
%         Each row represents the FPT density for one of the initial values.
% ret_num_tr_states is the number of transient states. This is for the case
%                   when you use the built-in mechanism for choosing an optimal
%                   number of states and you would like to know what spatial
%                   resolution was chosen.
% g_dropout is the PDF for dropout times multiplied by the probability of dropping
%           out before hitting a boundary, evaluated at the times given in
%           t_vec. If ou_init is a vector, g_dropout is a matrix. Each row
%           represents the distribution for one of the initial values.
%
% ou_drift is the drift of the OU process. It can either be a constant or,
%          for the time-variant case, the name of a function, which must return
%          the drift when called with the time as the argument.
% ou_leak defines the "leakiness" of the integrator. The deterministic part
%         of the stochastic differential equation is given by
%         ou_drift - ou_leak * current_value. A Wiener process can be studied
%         by setting ou_leak to 0.
% ou_var is the variance of the OU process. It can either be a constant or,
%        for the time-variant case, the name of a function, which must return
%        the variance when called with the time as the argument.
% ou_init is the initial value of the OU process. A vector can be passed,
%         if you are interested in studying several inital values at once.
%         This doesn't cause additional computational overhead.
% a_upper defines the upper absorbing boundary. a_upper can either be a constant
%         or the name of a function, which must return the location of the boundary
%         when called with the time as the argument.
% a_lower defines the lower absorbing boundary. See a_upper for the format.
% delta_t is the temporal step size.
% covered_range defines the range, which should be covered by the transient states.
%               It has to be a 2 element vector. Make sure that you never define boundaries
%               outside the specified range! (In this case the algorithm doesn't solve the
%               problem you are interested in.) In the case of constant boundaries the best
%               results are obtained when the boundaries are identical to
%               covered_range(1) and covered_range(2).
% num_tr_states defines the number of transient states used for the discrete approximation
%               of the problem. When using a symmetrical range you should specify
%               an odd number of states for an unbiased representation of 0.
%               When passing 0 a built-in mechanism is used for choosing an optimal
%               number of states (based on a heuristic method).
% stop_type defines, what determines when the algorithm stops.
%           1 = time: The densities are calculated for times up to stop_val.
%           2 = deviation of the distribution from 1:
%               The densities are calculated until the deviation of the sum
%               of the integrals over the densities from 1 is smaller than
%               stop_val.
% stop_val defines the time or the error, which determines when the algorithm
%          will stop.
% dropout_rate is the rate defining the risk of dropping out of the process
%              at any time. It can be either a constant or, for the time-variant
%              case, the name of a function which must return the rate when
%              called with the time as the argument.
% mapping_range is an optional parameter, which defines how long the vectors
%               used for constructing the transition matrix are and therefore
%               how sparse the transition matrix will be. The mapped range is
%               MEAN +/- mapping_range * SD. The default value is 3.
% runmed_width is an optional parameter, which defines the width of a running median filter
%              applied to the output in case of time-variant boundaries to remove spikes.
%              It has to be an odd number. The default value is 15. 0 deactivates the filter.

% History:
% released on 3/17/03 as part of toolbox V 2.5

% State layout:
%
% absorbing                        transient                                        absorbing
%     0          1          2         ...         num_tr_states         num_tr_states+1   num_tr_states+2
%           |                                                       |                     (dropout state)
%    covered_range(1)                                        covered_range(2)
%
% State # i is centered around covered_range(1)+(i-.5)*diff(covered_range)/num_tr_states. 

% Compiler flag:
%#realonly

if nargin<14 % runmed_width not given?
    runmed_width=15; % default value
end;

if nargin<13 % mapping_range not given?
    mapping_range=3; % default value
end;

num_tr_states=round(num_tr_states);
num_states=num_tr_states+3; % The total number of states is the number of transient states plus the three absorbing states.
runmed_width=round(runmed_width);

% Heuristic method for choosing an optimal number of states
if num_tr_states==0 % Should the number of states be determined automatically?
    if isnumeric(ou_var) % constant variance?
        temp_var=ou_var;
    else % time-variant variance
        temp_var=feval(ou_var,0);
    end;
    
    temp=8+diff(covered_range)^2/(45*temp_var*delta_t); % optimal total number of states (excluding dropout state)
    % There is no theoretical foundation for this formula.
    % It was determined by playing around with the parameters a lot and comparing the result
    % to the numerical solution.
    
    % choose the next odd number
    temp2=round(temp);
    
    if round(temp2/2)*2==temp2 % even?
        poss1=temp2-1; % first possibility
        poss2=temp2+1; % second possibility
        
        if abs(poss2-temp)<abs(poss1-temp) % Which one is closer to the original guess?
            num_states=poss2; % choose the larger number
        else
            num_states=poss1; % choose the smaller number
        end;
    else % odd
        num_states=temp2;
    end;
    
    num_states=num_states+1; % add the dropout state
    num_tr_states=num_states-3;
    disp(['OU_DROPOUT_1D_2B_MAR: The automatically chosen number of transient states is ' num2str(num_tr_states) '.']);
    
    if num_tr_states<9 % number of states too small?
        error('OU_DROPOUT_1D_2B_MAR: The results would be likely to be inaccurate. Please choose a smaller DELTA_T!');
    end;   
else % number of states was given by the user
    if num_tr_states<0 % negative number of states?
        error('OU_DROPOUT_1D_2B_MAR: There must be at least one transient state!');
    end;
    
    if num_tr_states<9
        warning('OU_DROPOUT_1D_2B_MAR: You should choose at least 9 transient states. Otherwise the results are very likely to be inaccurate!');
    end;   
end;

ret_num_tr_states=num_tr_states;

% Some checks
if ou_leak<0
    error('OU_DROPOUT_1D_2B_MAR: OU_LEAK must be a non-negative number!');
end;

if delta_t<=0
    error('OU_DROPOUT_1D_2B_MAR: The time step must be a positive number!');
end;

if diff(covered_range)<0 % screwed up range?
    error('OU_DROPOUT_1D_2B_MAR: Invalid range!');
end;

if (stop_type~=1)&(stop_type~=2)
    error('OU_DROPOUT_1D_2B_MAR: STOP_TYPE must be either 1 or 2!');
end;

if stop_val<=0
    error('OU_DROPOUT_1D_2B_MAR: STOP_VAL must be a positive number!');
end;

if mapping_range<=0
    error('OU_DROPOUT_1D_2B_MAR: MAPPING_RANGE must be a positive number!');
end;

if runmed_width<0
    error('OU_DROPOUT_1D_2B_MAR: RUNMED_WIDTH must not be negative!');
end;

if runmed_width&(~mod(runmed_width,2))
    error('OU_DROPOUT_1D_2B_MAR: RUNMED_WIDTH must be an odd number!');
end;

% Initialization
g_lower=[];
g_upper=[];
g_dropout=[];
t_vec=[];
k=1;
last_a_lower=nan;
last_a_upper=nan;
last_drift=nan;
last_sd=nan;
last_p_dropout=nan;

if isnumeric(a_lower) % Is the first boundary time-invariant?
    a_lower_const=1;
    
    if (a_lower<covered_range(1))|(a_lower>covered_range(2)) % boundary out of range?
        warning('OU_DROPOUT_1D_2B_MAR: Boundary out of range!');
    end;
    
    a_lower_cur=lb2state(a_lower,covered_range,num_states-1);
else
    a_lower_const=0;
end;

if isnumeric(a_upper) % Is the second boundary time-invariant?
    a_upper_const=1;
    
    if (a_upper<covered_range(1))|(a_upper>covered_range(2)) % boundary out of range?
        warning('OU_DROPOUT_1D_2B_MAR: Boundary out of range!');
    end;
    
    a_upper_cur=ub2state(a_upper,covered_range,num_states-1);
else
    a_upper_const=0;
end;

if isnumeric(ou_drift) % Is the drift time-invariant?
    drift_const=1;
    drift_cur=incval2stateinc(ou_drift*delta_t,covered_range,num_states-1); % We need the increment in the state space per time step.
else
    drift_const=0;
end;

if isnumeric(ou_var) % Is the variance time-invariant?
    var_const=1;
    
    if ou_var<=0
        error('OU_DROPOUT_1D_2B_MAR: The variance must be a positive number!');
    end;   

    sd_cur=incval2stateinc(sqrt(ou_var*delta_t),covered_range,num_states-1); % We need the standard deviation in the state space per time step.
else
    var_const=0;
end;

ou_init=val2state(ou_init,covered_range,num_states-1); % We need the initial values in the state space.

if isnumeric(dropout_rate) % Is the dropout rate time-invariant?
    dropout_const=1;
    
    if ((dropout_rate*delta_t)<0)|((dropout_rate*delta_t)>1)
        error('OU_DROPOUT_1D_2B_MAR: The dropout probability must be between 0 and 1!');
    end;
    
    p_dropout_cur=dropout_rate*delta_t;
else
    dropout_const=0;
end;

current_mat=diag(ones(1,num_tr_states)); % start with a unity matrix

% Further checks
ind1=find(ou_init<1);
ind2=find(ou_init>num_tr_states);

if ~isempty(ind1)|~isempty(ind2) % Are there initial values, which map to absorbing states?
    error('OU_DROPOUT_1D_2B_MAR: All initial values must map to transient states!');
end;

% Loop
while (1)
    t=k*delta_t; % current time
    
    if (((stop_type==1)&(t>stop_val))|((stop_type==2)&((1-(sum(g_lower')+sum(g_upper')+sum(g_dropout'))*delta_t)<stop_val))) % enough?
        break;
    end;
    
    if ~a_lower_const % Do we have to update the first boundary?
        temp=feval(a_lower,t);
        
        if (temp<covered_range(1))|(temp>covered_range(2)) % boundary out of range?
            warning('OU_DROPOUT_1D_2B_MAR: Boundary out of range!');
        end;
        
        a_lower_cur=lb2state(temp,covered_range,num_states-1);
    end;
    
    if ~a_upper_const % Do we have to update the second boundary?
        temp=feval(a_upper,t);
        
        if (temp<covered_range(1))|(temp>covered_range(2)) % boundary out of range?
            warning('OU_DROPOUT_1D_2B_MAR: Boundary out of range!');
        end;
        
        a_upper_cur=ub2state(temp,covered_range,num_states-1);
    end;
    
    if a_upper_cur<=a_lower_cur % crossing boundaries?
        warning('OU_DROPOUT_1D_2B_MAR: Algorithm stopped due to crossing boundaries!');
        break;
    end;
    
    if ~drift_const % Do we have to update the current drift?
        drift_cur=incval2stateinc(feval(ou_drift,t)*delta_t,covered_range,num_states-1);
    end;
    
    if ~var_const % Do we have to update the current variance?
        temp=feval(ou_var,t);
        
        if temp<=0
            warning('OU_DROPOUT_1D_2B_MAR: Algorithm stopped due to a non-positive variance!');
            break;
        end;      

        sd_cur=incval2stateinc(sqrt(temp*delta_t),covered_range,num_states-1);        
    end;
    
    if ~dropout_const % Do we have to update the current dropout probability?
        temp=feval(dropout_rate,t);
        
        if ((temp*delta_t)<0)|((temp*delta_t)>1)
            warning('OU_DROPOUT_1D_2B_MAR: Algorithm stopped because the dropout probability was out of range!');
            break;
        end;
        
        p_dropout_cur=temp*delta_t;
    end;
    
    % Do we have to construct a new transition matrix?
    if (a_lower_cur~=last_a_lower)|(a_upper_cur~=last_a_upper)|(drift_cur~=last_drift)|(sd_cur~=last_sd)|(p_dropout_cur~=last_p_dropout)
        last_a_lower=a_lower_cur;
        last_a_upper=a_upper_cur;
        last_drift=drift_cur;
        last_sd=sd_cur;
        last_p_dropout=p_dropout_cur;

        if ou_leak==0 % Wiener process? -> The transition vector has to be calculated only once.
            m_cur=drift_cur; % The drift is the only deterministic component.
            lower_bound=floor(m_cur-mapping_range*sd_cur);
            upper_bound=ceil(m_cur+mapping_range*sd_cur);
            x=[lower_bound+.5:upper_bound-.5];
            y=diff([0 normcdf(x,m_cur,sd_cur) 1]); % get the transition probabilities
            y=y*(1-p_dropout_cur); % take dropout probability into account
        end;      
        
        helpm=zeros(num_tr_states,num_states); % This is a temporary matrix for constructing the transition matrix.
        
        for row=1:num_tr_states % This is the "from" state.
            if ou_leak~=0 % OU process? -> The transition vector has to be recalculated for each row.
                m_cur=drift_cur+leak2stateinc(row,ou_leak,delta_t,covered_range,num_states-1); % Take the "leakage" of the integrator into account!
                lower_bound=floor(m_cur-mapping_range*sd_cur);
                upper_bound=ceil(m_cur+mapping_range*sd_cur);
                x=[lower_bound+.5:upper_bound-.5];
                y=diff([0 normcdf(x,m_cur,sd_cur) 1]); % get the transition probabilities
                y=y*(1-p_dropout_cur); % take dropout probability into account
            end;  
            
            index=-row;
            
            for col=1:num_states % This is the "to" state (+1).
                if (col==1) % first col is special (lower absorbing boundary)
                    if (row<=a_lower_cur) % source state below (or part of) lower boundary?
                        helpm(row,col)=1; % is absorbed by lower boundary
                    elseif (row>=a_upper_cur) % source state above (or part of) upper boundary?
                        %helpm(row,col)=0; % is absorbed by UPPER boundary
                    elseif (-row+a_lower_cur>=lower_bound) % entry necessary ?
                        if (-row+a_lower_cur>upper_bound)
                            helpm(row,col)=1;
                        else
                            helpm(row,col)=sum(y(1:-row+a_lower_cur-lower_bound+1));
                        end;               
                    end;               
                    
                elseif (col==(num_states-1)) % second last col is special (upper absorbing boundary)
                    if (row<=a_lower_cur) % source state below (or part of) lower boundary?
                        %helpm(row,col)=0; % is absorbed by LOWER boundary
                    elseif (row>=a_upper_cur) % source state above (or part of) upper boundary?
                        helpm(row,col)=1; % is absorbed by upper boundary
                    elseif (-row+a_upper_cur<=upper_bound) % entry necessary ?
                        if (-row+a_upper_cur<lower_bound)
                            helpm(row,col)=1;
                        else
                            helpm(row,col)=sum(y(-row+a_upper_cur-lower_bound+1:size(y,2)));
                        end;               
                    end;
                    
                elseif (col==num_states) % last col is special (dropout state)
                    helpm(row,col)=p_dropout_cur;
                    
                else % standard col
                    if (row<=a_lower_cur) % source state below (or part of) lower boundary?
                        %helpm(row,col)=0; % is absorbed by LOWER boundary
                    elseif (row>=a_upper_cur) % source state above (or part of) upper boundary?
                        %helpm(row,col)=0; % is absorbed by UPPER boundary
                    elseif (col-1<=a_lower_cur | col-1>=a_upper_cur) % destination state outside (or part of) the boundary range?
                        %helpm(row,col)=0; % is absorbed by the respective boundary
                    elseif ((index>=lower_bound)&(index<=upper_bound)) % entry necessary ?
                        helpm(row,col)=y(index-lower_bound+1);
                    end;            
                end;
                
                index=index+1;   
            end; % for col         
        end; % for row
        
        % construct R
        R=[helpm(:,1) helpm(:,num_states-1) helpm(:,num_states)];
        
        % construct Q
        Q=helpm(:,2:num_states-2);
        
    end; % new transition matrix?
    
    current_probs=current_mat*R;
    
    g_lower=[g_lower current_probs(ou_init,1)/delta_t];
    g_upper=[g_upper current_probs(ou_init,2)/delta_t];
    g_dropout=[g_dropout current_probs(ou_init,3)/delta_t];
    t_vec=[t_vec t];
    
    current_mat=current_mat*Q;
    k=k+1;
end; % loop

if (~a_lower_const|~a_upper_const)&(runmed_width>0) % Are there time-variant boundaries?
    % Run the output through a running median filter to remove spikes.
    if length(g_lower(1,:))<runmed_width % inappropriate filter width?
        error('OU_DROPOUT_1D_2B_MAR: Inappropriate filter width!');
    end;
    
    for i=1:length(ou_init)
        g_lower(i,:)=runmed(g_lower(i,:),runmed_width,1,0);
        g_upper(i,:)=runmed(g_upper(i,:),runmed_width,1,0);
    end;
end;

% Local functions
function state=val2state(val,covered_range,num_states)
% The states are numbered from 0 to num_states-1. These are the 2 absorbing states.
% The transient states are numbered from 1 to num_states-2.
temp=(val-covered_range(1)+diff(covered_range)/(2*(num_states-2)))/diff(covered_range)*(num_states-2);
ind=find(temp<.5);
temp(ind)=0;
ind=find(temp>num_states-1.5);
temp(ind)=num_states-1;
state=round(temp);

function val=state2val(state,covered_range,num_states)
val=covered_range(1)+(state-.5)*diff(covered_range)/(num_states-2);

function state_inc=incval2stateinc(inc_val,covered_range,num_states)
state_inc=inc_val/diff(covered_range)*(num_states-2);

function state_inc=leak2stateinc(state,ou_leak,delta_t,covered_range,num_states)
state_inc=incval2stateinc(-state2val(state,covered_range,num_states)*ou_leak*delta_t,covered_range,num_states);

function state=lb2state(val,covered_range,num_states)
% Transform lower boundary into the number of the virtual absorbing state.
% The boundary should be located at the transition from one state to another,
% not centered on a state. This is why val2state cannot be used.
temp=(val-covered_range(1))/diff(covered_range)*(num_states-2);
ind=find(temp<.5);
temp(ind)=0;
state=round(temp);

function state=ub2state(val,covered_range,num_states)
% Transform upper boundary into the number of the virtual absorbing state.
% The boundary should be located at the transition from one state to another,
% not centered on a state. This is why val2state cannot be used.
temp=(val-covered_range(1))/diff(covered_range)*(num_states-2)+1;
ind=find(temp>num_states-1.5);
temp(ind)=num_states-1;
state=round(temp);
