function [g_upper,g_lower,t_vec]=ou_1d_2b_num(ou_drift,ou_leak,ou_var,ou_init,a_upper,a_lower,delta_t,stop_type,stop_val)
% Numerical solution of the first passage time problem
% for a 1D Ornstein-Uhlenbeck process with 2 (time-variant) boundaries.
% The algorithm can be found in Smith, P. L.: Stochastic Dynamic
% Models of Response Time and Accuracy: A Foundational Primer.
% Journal of Mathematical Psychology, 44 (2000) 408-463.
%
% J. Ditterich, 8/01
%
% [g_upper,g_lower,t_vec] = ou_1d_2b_num (ou_drift,ou_leak,ou_var,ou_init,a_upper,a_lower,delta_t,stop_type,stop_val)
%
% g_upper is the first passage time density for the upper boundary multiplied
%         by the probability of hitting the upper boundary first, evaluated
%         at the times given in t_vec.
% g_lower is the first passage time density for the lower boundary multiplied
%         by the probability of hitting the lower boundary first, evaluated
%         at the times given in t_vec.
%
% ou_drift is the drift of the OU process. It can either be a constant or,
%          for the time-variant case, the name of a function, which must return
%          the drift when called with the time as the argument.
% ou_leak defines the "leakiness" of the integrator. The deterministic part
%         of the stochastic differential equation is given by
%         ou_drift - ou_leak * current_value. ou_leak must not be 0! Please use
%         WIENER_1D_2B_NUM in this case.
% ou_var is the variance of the OU process.
% ou_init is the initial value of the OU process.
% a_upper defines the upper absorbing boundary. a_upper can either be a constant or
%         the name of a function, which must return the location of the boundary
%         and its temporal derivative when called with the time as the argument.
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

% Compiler flag:
%#realonly

% Some checks
if (ou_leak==0) % Wiener Process?
   error('OU_1D_2B_NUM: OU_LEAK must not be 0! Please use WIENER_1D_2B_NUM.');
end;

if ou_leak<0
   error('OU_1D_2B_NUM: OU_LEAK must be a positive number!');
end;

if ou_var<=0
   error('OU_1D_2B_NUM: The variance must be a positive number!');
end;

if delta_t<=0
   error('OU_1D_2B_NUM: The time step must be a positive number!');
end;

if (stop_type~=1)&(stop_type~=2)
   error('OU_1D_2B_NUM: STOP_TYPE must be either 1 or 2!');
end;

if stop_val<=0
   error('OU_1D_2B_NUM: STOP_VAL must be a positive number!');
end;

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

if isnumeric(ou_drift) % Is the drift time-invariant?
   drift_const=1;
   drift_cur=ou_drift;
else
   drift_const=0;
end;

% Loop
while (1)
   t=k*delta_t; % current time
   
   if (((stop_type==1)&(t>stop_val))|((stop_type==2)&((1-(sum(g_upper)+sum(g_lower))*delta_t)<stop_val))) % enough?
      break;
   end;
      
   if ~a_upper_const % Do we have to update the first boundary?
      [a_upper_cur,a_upper_deriv_cur]=feval(a_upper,t);
   end;
   
   if ~a_lower_const % Do we have to update the second boundary?
      [a_lower_cur,a_lower_deriv_cur]=feval(a_lower,t);
   end;
   
   if a_upper_cur<=a_lower_cur % crossing boundaries?
      warning('OU_1D_2B_NUM: Algorithm stopped due to crossing boundaries!');
      break;
   end;
   
   if (k==1)&((ou_init>=a_upper_cur)|(ou_init<=a_lower_cur)) % initial value out of range?
      error('OU_1D_2B_NUM: The initial value must be in between the two boundaries!');
   end;
   
   if ~drift_const % Do we have to update the current drift?
      drift_cur=feval(ou_drift,t);
      drift_vec=[drift_vec feval(ou_drift,(k-.5)*delta_t)]; % new value for numerical integration
      drift_weights=exp(([-length(drift_vec):-1]+.5)*ou_leak*delta_t);
      drift_int=(drift_vec*drift_weights')*delta_t; % approx. of the integral over the time-variant drift
   else % constant drift
      drift_int=(1-exp(-ou_leak*t))/ou_leak*ou_drift;
   end;
         
   g_upper_cur=-2*psi(a_upper_cur,a_upper_deriv_cur,t,ou_init,0,drift_cur,ou_leak,ou_var,drift_int);
   g_lower_cur=2*psi(a_lower_cur,a_lower_deriv_cur,t,ou_init,0,drift_cur,ou_leak,ou_var,drift_int);
      
   if (k>1) % not first step?
      for j=1:k-1
         
         if ~drift_const % time-variant drift
            drift_weights=exp(([-(k-j):-1]+.5)*ou_leak*delta_t);
            drift_int=(drift_vec(j+1:k)*drift_weights')*delta_t; % approx. of the integral over the time-variant drift
         else % constant_drift
            drift_int=(1-exp(-ou_leak*(t-j*delta_t)))/ou_leak*ou_drift;
         end;
                     
         g_upper_cur=g_upper_cur+2*delta_t*(g_upper(j)*psi(a_upper_cur,a_upper_deriv_cur,t,a_upper_vec(j),j*delta_t,drift_cur,ou_leak,ou_var,drift_int) ...
            +g_lower(j)*psi(a_upper_cur,a_upper_deriv_cur,t,a_lower_vec(j),j*delta_t,drift_cur,ou_leak,ou_var,drift_int));
         g_lower_cur=g_lower_cur-2*delta_t*(g_upper(j)*psi(a_lower_cur,a_lower_deriv_cur,t,a_upper_vec(j),j*delta_t,drift_cur,ou_leak,ou_var,drift_int) ...
            +g_lower(j)*psi(a_lower_cur,a_lower_deriv_cur,t,a_lower_vec(j),j*delta_t,drift_cur,ou_leak,ou_var,drift_int));
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
function val=psi(a,a_deriv,t,y,tau,ou_drift,ou_leak,ou_var,drift_int)
val=f(a,t,y,tau,ou_leak,ou_var,drift_int)/2*(a_deriv+ou_leak*a-ou_drift-(2*ou_leak)/(1-exp(-2*ou_leak*(t-tau))) ...
   *(a-exp(-ou_leak*(t-tau))*y-drift_int));

function val=f(x,t,y,tau,ou_leak,ou_var,drift_int)
val=sqrt(ou_leak/(pi*ou_var*(1-exp(-2*ou_leak*(t-tau)))))*exp((-ou_leak*(x-y*exp(-ou_leak*(t-tau))-drift_int)^2) ...
   /(ou_var*(1-exp(-2*ou_leak*(t-tau)))));
