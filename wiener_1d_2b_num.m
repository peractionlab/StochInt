function [g_upper,g_lower,t_vec]=wiener_1d_2b_num(wie_drift,wie_var,wie_init,a_upper,a_lower,delta_t,stop_type,stop_val)
% Numerical solution of the first passage time problem
% for a 1D Wiener process with 2 (time-variant) boundaries.
% The algorithm can be found in Smith, P. L.: Stochastic Dynamic
% Models of Response Time and Accuracy: A Foundational Primer.
% Journal of Mathematical Psychology, 44 (2000) 408-463.
%
% J. Ditterich, 8/01
%
% [g_upper,g_lower,t_vec] = wiener_1d_2b_num (wie_drift,wie_var,wie_init,a_upper,a_lower,delta_t,stop_type,stop_val)
%
% g_upper is the first passage time density for the upper boundary multiplied
%         by the probability of hitting the upper boundary first, evaluated
%         at the times given in t_vec.
% g_lower is the first passage time density for the lower boundary multiplied
%         by the probability of hitting the lower boundary first, evaluated
%         at the times given in t_vec.
%
% wie_drift is the drift of the Wiener process. It can either be a constant or,
%           for the time-variant case, the name of a function, which must return
%           the drift when called with the time as the argument.
% wie_var is the variance of the Wiener process.
% wie_init is the initial value of the Wiener process.
% a_upper defines the upper absorbing boundary. a_upper can either be a constant or
%         the name of a function, which must return the location of the boundary
%         and its temporal derivative when called with the time as the argument.
%         In case an additional argument needs to be passed to the function,
%         a_upper can also be a cell array with two elements. In this case the first cell has
%         to be the name of the function, the second cell is passed as
%         the second function argument (after the time).
% a_lower defines the lower absorbing boundary. See a_upper for the format.
% delta_t is the temporal resolution used by the algorithm.
% stop_type defines, what determines when the algorithm stops.
%           1 = time: The densities are calculated for times up to stop_val.
%           2 = deviation of the distribution from 1:
%               The densities are calculated until the deviation of the sum
%               of the integrals over the densities from 1 is smaller than
%               stop_val.
% stop_val defines the time or the error, which determines when the algorithm
%          will stop.

% History:
% released on 8/23/01 as part of toolbox V 0.5 Beta
% a_upper and a_lower adjusted for Parallel Computing Toolbox on 12/12/08

% Initialization
g_upper=[];
g_lower=[];
t_vec=[];
a_upper_vec=[];
a_lower_vec=[];
drift_vec=[];
k=1;

if isnumeric(a_upper) % Is the first boundary time-invariant?
    a_upper_const=1;
    a_upper_cur=a_upper;
    a_upper_deriv_cur=0;
else
    a_upper_const=0;
end;

if isnumeric(a_lower) % Is the second boundary time-invariant?
   a_lower_const=1;
   a_lower_cur=a_lower;
   a_lower_deriv_cur=0;
else
   a_lower_const=0;
end;

if isnumeric(wie_drift) % Is the drift time-invariant?
   drift_const=1;
   drift_cur=wie_drift;
else
   drift_const=0;
end;

% Some checks
if wie_var<=0
   error('WIENER_1D_2B_NUM: The variance must be a positive number!');
end;

if (iscell(a_upper))&&(length(a_upper)~=2)
    error('WIENER_1D_2B_NUM: A_UPPER must have 2 elements when given as a cell array!');
end;

if (iscell(a_lower))&&(length(a_lower)~=2)
    error('WIENER_1D_2B_NUM: A_LOWER must have 2 elements when given as a cell array!');
end;

if delta_t<=0
   error('WIENER_1D_2B_NUM: The time step must be a positive number!');
end;

if (stop_type~=1)&(stop_type~=2)
   error('WIENER_1D_2B_NUM: STOP_TYPE must be either 1 or 2!');
end;

if stop_val<=0
   error('WIENER_1D_2B_NUM: STOP_VAL must be a positive number!');
end;

% Loop
while (1)
   t=k*delta_t; % current time
   
   if (((stop_type==1)&(t>stop_val))|((stop_type==2)&((1-(sum(g_upper)+sum(g_lower))*delta_t)<stop_val))) % enough?
      break;
   end;
      
   if ~a_upper_const % Do we have to update the first boundary?
       if iscell(a_upper) % Do we need to pass a second function argument?
           [a_upper_cur,a_upper_deriv_cur]=feval(a_upper{1},t,a_upper{2});
       else % no second function argument
           [a_upper_cur,a_upper_deriv_cur]=feval(a_upper,t);
       end;
   end;
   
   if ~a_lower_const % Do we have to update the second boundary?
       if iscell(a_lower) % Do we need to pass a second function argument?
           [a_lower_cur,a_lower_deriv_cur]=feval(a_lower{1},t,a_lower{2});
       else % no second function argument
           [a_lower_cur,a_lower_deriv_cur]=feval(a_lower,t);
       end;
   end;
   
   if a_upper_cur<=a_lower_cur % crossing boundaries?
      warning('WIENER_1D_2B_NUM: Algorithm stopped due to crossing boundaries!');
      break;
   end;
   
   if (k==1)&((wie_init>=a_upper_cur)|(wie_init<=a_lower_cur)) % initial value out of range?
      error('WIENER_1D_2B_NUM: The initial value must be in between the two boundaries!');
   end;
      
   if ~drift_const % Do we have to update the current drift?
      drift_cur=feval(wie_drift,t);
      drift_vec=[drift_vec feval(wie_drift,(k-.5)*delta_t)]; % new value for numerical integration
      drift_int=sum(drift_vec)*delta_t; % approx. of the integral over the time-variant drift
   else % constant drift
      drift_int=wie_drift*t;
   end;
         
   g_upper_cur=-2*psi(a_upper_cur,a_upper_deriv_cur,t,wie_init,0,drift_cur,wie_var,drift_int);
   g_lower_cur=2*psi(a_lower_cur,a_lower_deriv_cur,t,wie_init,0,drift_cur,wie_var,drift_int);
      
   if (k>1) % not first step?
      for j=1:k-1
         
         if ~drift_const % time-variant drift
            drift_int=sum(drift_vec(j+1:k))*delta_t; % approx. of the integral over the time-variant drift
         else % constant_drift
            drift_int=wie_drift*(k-j)*delta_t;
         end;
                     
         g_upper_cur=g_upper_cur+2*delta_t*(g_upper(j)*psi(a_upper_cur,a_upper_deriv_cur,t,a_upper_vec(j),j*delta_t,drift_cur,wie_var,drift_int) ...
            +g_lower(j)*psi(a_upper_cur,a_upper_deriv_cur,t,a_lower_vec(j),j*delta_t,drift_cur,wie_var,drift_int));
         g_lower_cur=g_lower_cur-2*delta_t*(g_upper(j)*psi(a_lower_cur,a_lower_deriv_cur,t,a_upper_vec(j),j*delta_t,drift_cur,wie_var,drift_int) ...
            +g_lower(j)*psi(a_lower_cur,a_lower_deriv_cur,t,a_lower_vec(j),j*delta_t,drift_cur,wie_var,drift_int));
      end;      
   end;
   
   g_upper=[g_upper g_upper_cur];
   g_lower=[g_lower g_lower_cur];
   t_vec=[t_vec t];
   a_upper_vec=[a_upper_vec a_upper_cur];
   a_lower_vec=[a_lower_vec a_lower_cur];
   k=k+1;
end;

% Local functions
function val=psi(a,a_deriv,t,y,tau,wie_drift,wie_var,drift_int)
val=f(a,t,y,tau,wie_var,drift_int)/2*(a_deriv-wie_drift-(a-y-drift_int)/(t-tau));

function val=f(x,t,y,tau,wie_var,drift_int)
val=1/sqrt(2*pi*wie_var*(t-tau))*exp(-(x-y-drift_int)^2/(2*wie_var*(t-tau)));
