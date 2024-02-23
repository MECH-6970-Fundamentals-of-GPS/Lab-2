function [svPos, svVel, svB, svD, svPrns] = sv_positions(ephem, psr, t)
%sv_positions.m Provides SV states given ephemeris data and pseudoranges
%   Given ephemeris data, pseudoranges, and GPS time this function returns
%   SV position, velocity, clock bias, drift, and SV prns
% Inputs:
%    ephem  : a matrix of ephemeris data from all SVs
%    psr    : a vector of pseudoranges from all SVs
%    t      : the current time
%
% Outputs:
%    svPos  : a matrix of all visible SV positions [x,y,z]
%    svVel  : a matrix of all visible SV velocities [vx,vy,vz]
%    svB    : a vector of all visible SV clock biases
%    svD    : a vector of all visible SV clock drifts
%    svPrns : a vector of all visible SV PRNs

    % CONSTANTS
    c = 299792458;                              % Speed of Light [m/s]
    
    % FUNCTION VARIABLES
    prnList = 1:length(psr);                    % List of all SV PRNs
    svPrns = prnList(~isnan(psr));              % List of currently visible SV PRNs
    N = length(svPrns);                         % Number of currently visible SVs
    svPos = zeros(N,3);                         % SV Positions [m]
    svVel = zeros(N,3);                         % SV Velocities [m/s]
    svB = zeros(N,1);                           % SV Clock Bias [s]
    svD = zeros(N,1);                           % SV Clock Drift [s/s]

    % CALCULATIONS
    % Iterating over all visible SVs
    for i = 1:N
        tau = psr(svPrns(i)) / c;               % Message travel time [s]
        % Calculate SV state from ephemeris data
        [svPos(i, :), svVel(i, :), svB(i), svD(i)] = ephem_to_sv_pos(ephem(svPrns(i)), t - tau, tau);
    end
end