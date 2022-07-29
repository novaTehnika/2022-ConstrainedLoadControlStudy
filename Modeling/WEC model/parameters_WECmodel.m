function par = parameters_WECmodel(par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters_WECmodel.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/29/2021
%
% PURPOSE/DESCRIPTION:
% This function loads prepares parameters used in the hydrodynamcic WEC
% model. Parameters are loaded from two files. One contains hydrodynamic
% data that are outputs of boundary element method code and the other
% contains general parameters and coefficients of a radiation damping model
% calculated seperately from the BEM code results. This funciton calculates
% the wave specrum based on the significant wave height and peak period,
% and prepares the matrices for the radiation damping, state-space model.
%
% FILE DEPENDENCY: 
% nemohResults_vantHoff2009_20180802.mat
% vantHoffTFCoeff.mat
%
% UPDATES:
% 6/29/2021 - created from buoyTorque.m developed in 2018 
% by the author. The file is rewritten to take a 
% structure of parameters as an input instead of using hard coded 
% variables.
% 7/9/2021 - implemented the equal-area (equal-energy) method of 
% generating wave elevation timeseries and excitation force. This method 
% involves interpolation of the excitation force. Additions include the use
% of nested functions for the wave spectrum and obj. func.
%
% Copyright (C) 2022  Jeremy W. Simmons II
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program. If not, see <https://www.gnu.org/licenses/>.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load saved WEC parameters
load('nemohResults_vantHoff2009_20180802.mat','OSWC_vantHoff2009_simmons')
par.imprt.WEC.w = OSWC_vantHoff2009_simmons.w;          % [rad/s] frequency
par.imprt.WEC.F_amp = OSWC_vantHoff2009_simmons.F_amp;  % [Nm/m] frequency dependent excitation torque amplitude coefficient
par.imprt.WEC.phi_e = OSWC_vantHoff2009_simmons.phi_e;  % [rad] frequency dependent phase of excitation torque
par.imprt.WEC.I_inf = OSWC_vantHoff2009_simmons.A_inf(1);  % [kgm^2] added mass at infinite frequency
par.WEC.I_inf = par.imprt.WEC.I_inf;
load('vantHoffTFCoeff_v2.mat','KradDen','KradNum')
par.WEC.a = KradDen';     % transfer funciton denominator coefficients for raditation damping derived to approximate the hydrodynamic response
par.WEC.b = KradNum';     % transfer funciton numerator coefficients for raditation damping derived to approximate the hydrodynamic response
par.WEC.L_flap = 11;      % [m] length of flap (hinge to top edge)
par.WEC.T = 2;            % [m] flap thickness
par.WEC.h = 8.9;          % [m] mean water depth
par.WEC.W = 18;           % [m] flap width
par.WEC.m = 127000;       % [kg] mass of flap
par.WEC.I = 1.85e6+0*par.WEC.m*3.9^2;     % [kgm^2] mom. of inertia of flap
par.WEC.rho = 1025;       % [kg/m^3] density of sea water
par.WEC.g = 9.81;         % [m/s^2] gravitational acceleration at sea level

% construct matrices for radiation damping state-space model
 % number of states describing radiation damping memory
par.WEC.ny_rad = length(par.WEC.a)-1; 

% Cast linear radiation model to Observable Cononical Form
A1 = [ eye(par.WEC.ny_rad-1); ...
       zeros(1,par.WEC.ny_rad-1)] ; 
par.WEC.A_rad  = [-par.WEC.a(2:end) A1];
par.WEC.B_rad = [zeros(par.WEC.ny_rad-length(par.WEC.b),1); ...
                 par.WEC.b];
             
             
% determine new freqencies based on equal energy method
 % determine the area under the curve of the wave spectrum for
 % discretization
  % high density grid of frequencies evenly spaced
w = linspace(par.imprt.WEC.w(1),par.imprt.WEC.w(end),1000); 

  % wave spectrum calculated for high-density even grid
S_w = PiersonSpec(w,par);

  % area for each bin
a = trapz(w,S_w)/(par.WEC.nw+1);

  % calculate the cumulative area under the curve
A = cumtrapz(w,S_w);
        
 % determine the frequencies of the equal area grid using optimization of 
 % an obj. function involving the error between the target area under the
 % curve and the area calculated from the tested frequency
par.WEC.w = zeros(par.WEC.nw,1);
par.WEC.dw = zeros(par.WEC.nw,1);
wmax = par.imprt.WEC.w(end);
for iw = 1:par.WEC.nw
    if iw == 1 % set wmin to min frequency in range considered
        wmin = par.imprt.WEC.w(1);
    else % set wmin to right bound of previous bin
        wmin = par.WEC.w(iw-1) + par.WEC.dw(iw-1)/2;
    end
    Atarget = a*iw; % target cumlative area under the curve (integral)
    
    % solve for the root for the difference between Atarget and A
    par.WEC.w(iw) = ...
        fzero(@(w_ea) Atarget - interp1(w,A,w_ea,'spline'),...
        (wmax - wmin)/10);
    
    % calculate the width of the bin based on the target area and the value
    % for the spectrum density function at this frequency 
	par.WEC.dw(iw) = a/PiersonSpec(par.WEC.w(iw),par);

end

% calculate the wave spectrum for the equal area frequencies
par.wave.S_w = PiersonSpec(par.WEC.w,par);

% interpolate frequency dependent coefficients
 % frequency dependent excitation torque amplitude coefficient
par.WEC.F_amp = interp1(par.imprt.WEC.w,par.imprt.WEC.F_amp,par.WEC.w,'spline');

 % frequency dependent phase of excitation torque
par.WEC.phi_e = interp1(par.imprt.WEC.w,par.imprt.WEC.phi_e,par.WEC.w,'spline');

% generate random phases for each frequency component for the wave elevation
 % seed the random number generator
rng(par.wave.rngSeedPhase); 

 % generate uniformly distributed random numbers [-pi,pi]
par.wave.phi = 2*pi*(rand(par.WEC.nw,1)-0.5);


% regular waves (supersedes previous calculations)
if par.WEC.nw == 1
    par.wave.S_w = 0.5*par.wave.Hs^2;
    par.WEC.dw = 1;
    par.WEC.w = 2*pi/par.wave.Tp;
    par.WEC.F_amp = interp1(par.imprt.WEC.w,par.imprt.WEC.F_amp,par.WEC.w,'spline');
    par.WEC.phi_e = interp1(par.imprt.WEC.w,par.imprt.WEC.phi_e,par.WEC.w,'spline');
    par.WEC.phi = 0;
end

%% %%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%
%     function Aerr = areaErrFun(w,A,Atarget,w_ea,wmin,wmax)
%         if w_ea < wmin 
%             Aerr = abs(w_ea-wmin)*1e1 + ((Atarget-A(1))/Atarget)^2;
%         elseif w_ea > wmax
%             Aerr = abs(w_ea-wmax)*1e1 + ((Atarget-A(end))/Atarget)^2;
%         else
%             Aerr = ((Atarget-interp1(w,A,w_ea,'spline'))/Atarget)^2;
%         end
%     end

    function S_w = PiersonSpec(w,par)
        % Based on Falnes (2002) "Ocean Waves and Oscillating Systems:..."
        S_w = 10*pi^5*par.wave.Hs^2/par.wave.Tp^4./w.^5 ... 
            .*exp(-20*pi^4/par.wave.Tp^4./w.^4);
    end

end