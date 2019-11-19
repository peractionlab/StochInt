function [g,t_vec]=wiener_1d_1b_lin_ana(wie_drift,wie_var,wie_init,b_slope,b0,delta_t,stop_type,stop_val)
% Analytical solution of the first passage time problem
% for a 1D Wiener process with a linear boundary.
% The formula can be found in Buonocore, A., Nobile, A. G.,
% Ricciardi, L. M.: A new integral equation for the
% evaluation of first-passage-time probability
% densities. Adv. Appl. Prob., 19 (1987) 784-800.
%
% J. Ditterich, 8/01
%
% [g,t_vec] = wiener_1d_1b_lin_ana (wie_drift,wie_var,wie_init,b_slope,b0,delta_t,stop_type,stop_val)
%
% g is the first passage time density, evaluated at the times
%   given in t_vec.
%
% wie_drift is the drift of the Wiener process.
% wie_var is the variance of the Wiener process.
% wie_init is the initial value of the Wiener process.
% b_slope is the slope of the boundary.
% b0 is the inital value of the boundary.
% delta_t is the temporal resolution.
% stop_type defines, what determines when the algorithm stops.
%           1 = time: The density is calculated for times up to stop_val.
%           2 = deviation of the distribution from 1:
%               The density is calculated until the deviation of the
%               integral over the density from 1 is smaller than
%               stop_val.
% stop_val defines the time or the error, which determines when the algorithm
%          will stop.

% History:
% released on 8/23/01 as part of toolbox V 0.5 Beta

% Compiler flag:
%#realonly

% Initialization
g=[];
t_vec=[];
k=1;

% Some checks
if wie_var<=0
   error('WIENER_1D_1B_LIN_ANA: The variance must be a positive number!');
end;

if delta_t<=0
   error('WIENER_1D_1B_LIN_ANA: The time step must be a positive number!');
end;

if (stop_type~=1)&(stop_type~=2)
   error('WIENER_1D_1B_LIN_ANA: STOP_TYPE must be either 1 or 2!');
end;

if stop_val<=0
   error('WIENER_1D_1B_LIN_ANA: STOP_VAL must be a positive number!');
end;

% Calculation
if (stop_type==1) % time criterion?
   t_max=floor(stop_val/delta_t)*delta_t;
   t_vec=[delta_t:delta_t:t_max];
   g=g_local(wie_drift,wie_var,wie_init,b_slope,b0,t_vec);
else % error criterion
   while (1)
      t=k*delta_t;
      
      if ((1-sum(g)*delta_t)<stop_val) % enough?
         break;
      end;
      
      g=[g g_local(wie_drift,wie_var,wie_init,b_slope,b0,t)];
      t_vec=[t_vec t];
      k=k+1;
   end;   
end;

% Local function
function val=g_local(wie_drift,wie_var,wie_init,b_slope,b0,t)
val=abs(b0-wie_init)./sqrt(2*pi*wie_var*t.^3).*exp(-(b0-wie_init+(b_slope-wie_drift)*t).^2./(2*wie_var*t));
