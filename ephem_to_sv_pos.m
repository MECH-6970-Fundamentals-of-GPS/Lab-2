function [svPos,svVel,svB, svD] = ephem_to_sv_pos(ephem, t_trans, t_trav)
%sv_positions.m Calculates SV state from ephemeris data
%   Processes ephemeris data from a single SV to give position, velocity,
%   bias, and drift.
% Inputs:
%    ephem  : vector of ephemeris data for a single SV
%    t_trans: time of transmission [s]
%    t_trav : length of time to travel to receiver [s]
%
% Outputs:
%    svPos  : a matrix of all visible SV positions [x,y,z]
%    svVel  : a matrix of all visible SV velocities [vx,vy,vz]
%    svB    : a vector of all visible SV clock biases
%    svD    : a vector of all visible SV clock drifts
%
% Author: Michaela Barksdale, Walter Livingston, Seth Porter
% Modified from work done by David Bevly

T_GD=ephem.T_GD;       %ephem_data(1)
t_oc=ephem.t_oc;       %ephem_data(2)
a_f2=ephem.a_f2;       %ephem_data(3)
a_f1=ephem.a_f1;
a_f0=ephem.a_f0;
C_rc=ephem.C_rc;
C_rs=ephem.C_rs;
C_uc=ephem.C_uc;
C_us=ephem.C_us;
C_ic=ephem.C_ic;
C_is=ephem.C_is;
Delta_n=ephem.deltan;
M_0=ephem.M_0;
e=ephem.e;
sqrt_A=sqrt(ephem.A); % Note this code assumes the data is A!!
t_oe=ephem.t_oe;
Omega_0=ephem.omega_0;
i_0=ephem.i_0;
omega=ephem.omega;
dot_Omega=ephem.omegaDot;
Idot=ephem.iDot;
% GPS constants

gpsPi = 3.1415926535898;  % Pi used in the GPS coordinate system

%--- Constants for satellite position calculation -------------------------
Omegae_dot = 7.2921151467e-5;  % Earth rotation rate, [rad/s]
GM = 3.986005e14;      % Earth's universal [m^3/s^2]
F = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]

%--- Find time difference ---------------------------------------------
dt = check_t(t_trans - t_oc);

%--- Calculate clock correction ---------------------------------------
svB = (a_f2 * dt + a_f1) * dt + a_f0 - T_GD;

time = t_trans - svB;

% Find satellite's position ----------------------------------------------

%Restore semi-major axis (this is not needed if ephemeris is already given
% as "a" - In that case simply use:  a=semiMajAx);
a   = sqrt_A*sqrt_A;

%Time correction
tk  = check_t(time - t_oe);

%Initial mean motion
n0  = sqrt(GM / a^3);
%Mean motion
n   = n0 + Delta_n;

%Mean anomaly
M   = M_0 + n * tk;

%Reduce mean anomaly to between 0 and 360 deg
M   = rem(M + 2*gpsPi, 2*gpsPi);

%Initial guess of eccentric anomaly
E   = M;

%--- Iteratively compute eccentric anomaly ----------------------------
for ii = 1:10
	E_old = E;
	E = M + e * sin(E);
	dE = rem(E - E_old, 2*gpsPi);

	if abs(dE) < 1.e-12
		% Necessary precision is reached, exit from the loop
		break;
	end
end

%Reduce eccentric anomaly to between 0 and 360 deg
E = rem(E + 2*gpsPi, 2*gpsPi);

%Compute relativistic correction term
dtr = F * e * sqrt_A * sin(E);

%Calculate the true anomaly
nu = atan2(sqrt(1 - e^2) * sin(E), cos(E)-e);

%Compute angle phi
phi = nu + omega;

%Reduce phi to between 0 and 360 deg
phi = rem(phi, 2*gpsPi);

%Correct argument of latitude
u = phi + C_uc * cos(2*phi) + C_us * sin(2*phi);
%Correct radius
r = a * (1 - e*cos(E)) + C_rc * cos(2*phi) + C_rs * sin(2*phi);
%Correct inclination
i = i_0 + Idot * tk + C_ic * cos(2*phi) + C_is * sin(2*phi);

%Compute the angle between the ascending node and the Greenwich meridian
Omega = Omega_0 + (dot_Omega - Omegae_dot)*tk - Omegae_dot * t_oe-Omegae_dot *t_trav;
%Reduce to between 0 and 360 deg
Omega = rem(Omega + 2*gpsPi, 2*gpsPi);

X = r*cos(u);
Y = r*sin(u);

%--- Compute satellite coordinates ------------------------------------
x = X*cos(Omega) - Y*cos(i)*sin(Omega);
y = X*sin(Omega) + Y*cos(i)*cos(Omega);
z = Y*sin(i);

%% SV velocity calculations
Edot = (n0 + Delta_n)/(1-e*cos(E));

phidot = (sqrt(1-e^2)/(1-e*cos(E)))*Edot;

udot = (1+2*C_us*cos(2*phi)-2*C_uc*sin(2*phi))*phidot;

rdot = 2*(C_rs*cos(2*phi)-C_rc*sin(2*phi))*phidot + a*e*sin(E)*Edot;

idot = 2*(C_is*cos(2*phi)-C_ic*sin(2*phi))*phidot + Idot;

Xdot = rdot*cos(u) - r*sin(u)*udot;

Ydot = rdot*sin(u) + r*cos(u)*udot;

Omegadot = dot_Omega - Omegae_dot;

xdot = Xdot*cos(Omega) - Ydot*cos(i)*sin(Omega) + Y*sin(i)*sin(Omega)*idot - y*Omegadot;

ydot = Xdot*sin(Omega) + Ydot*cos(i)*cos(Omega) - Y*sin(i)*cos(Omega)*idot + x*Omegadot;

zdot = Ydot*sin(i) + Y*cos(i)*idot;

%%


%--- Compute satellite coordinates ------------------------------------
svPos(1) = x;
svPos(2) = y;
svPos(3) = z;

svVel(1) = xdot;
svVel(2) = ydot;
svVel(3) = zdot;


% Include relativistic correction in clock correction -----------------
svB = (a_f2 * dt + a_f1) * dt + a_f0 - T_GD + dtr;

%% Drift Calculations
svD = a_f1 + 2*a_f2*dt + (n*F*e*sqrt(a)*cos(E))/(1-e*cos(E));

end

function corrTime = check_t(time)
%CHECK_T accounting for beginning or end of week crossover.
%
%corrTime = check_t(time);
%
%   Inputs:
%       time        - time in seconds
%
%   Outputs:
%       corrTime    - corrected time (seconds)

%Kai Borre 04-01-96
%Copyright (c) by Kai Borre
%
% CVS record:
% $Id: check_t.m,v 1.1.1.1.2.4 2006/08/22 13:45:59 dpl Exp $
%==================================================================
    half_week = 302400;     % seconds
    
    corrTime = time;
    
    if time > half_week
	    corrTime = time - 2*half_week;
    elseif time < -half_week
	    corrTime = time + 2*half_week;
    end
end
%%%%%%% end check_t.m  %%%%%%%%%%%%%%%%%


