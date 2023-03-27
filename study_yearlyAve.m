% study_yearlyAve.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 07/7/2022
%
% PURPOSE/DESCRIPTION:
% This function evaluate WEC load scheduling performance as a function of 
% sea state, upper load limit and lower fractional load limit using the 
% function m-file modelPredictiveLoadScheduling.m .
%
% This script is called as a function with the arguements
% - variable-set index (linear index within an n-dimentional grid)
%
% FILE DEPENDENCY:
% modelPredictiveLoadScheduling.m
% model_OLloadSchedule.m
% parameters_OLloadSchedule.m
% SSdata_HumboltBay_1D.mat
% data_coulombPTO_dampingStudy_20220927_slim.mat
%
% UPDATES:
% 01/10/2023 - Adapted from study_loadScheduleConstraints.m.
%
% Copyright (C) 2023  Jeremy W. Simmons II
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
function [] = study_yearlyAve(iVar)
path_models = ['.' filesep 'Modeling'];
addpath([path_models filesep 'Solvers'])
addpath([path_models filesep 'Sea States'])
addpath([path_models filesep 'WEC model']) 
addpath([path_models filesep 'WEC model' filesep 'WECdata']) 
addpath([path_models filesep 'Open-loop load schedule PTO']) 

%% %%%%%%%%%%%%   Study Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_max = logspace(log10(1e6),log10(8e6),10);
lbFrac = [0 0.25 0.5 0.75];
[NDmeshVar.SS,NDmeshVar.T_max,NDmeshVar.lbFrac] = ...
                                                ndgrid(1:114,T_max,lbFrac);

%% %%%%%%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation length
par.Tramp = 250; % [s] excitation force ramp period
par.tstart = 0; %[s] start time of simulation
par.tend = 2000; %[s] end time of simulation

% Solver parameters
par.odeSolverRelTol = 1e-9; % Rel. error tolerance parameter for ODE solver
par.odeSolverAbsTol = 1e-9; % Abs. error tolerance parameter for ODE solver
par.MaxStep = 1e-2; % Step size in the case of fixed step size solver

% Wave construction parameters
par.WEC.nw = 1000; % num. of frequency components for harmonic superposition 
par.wave.rngSeedPhase = 3; % seed for the random number generator

% Sea State specification
load('SSdata_HumboltBay_1D.mat','Hs','Tp');
iSS = NDmeshVar.SS(iVar);
par.wave.Hs = Hs(iSS);
par.wave.Tp = Tp(iSS);

% load remaining parameters
par = parameters_OLloadSchedule(par,...
    'nemohResults_vantHoff2009_20180802.mat','vantHoffTFCoeff.mat');

% Load Coulomb damping study results for limits on torque in each Sea State
load("data_coulombPTO_dampingStudy_20220927_slim.mat","T_c_data");

% Define intial conditions
y0 = [  0, ...
        0, ...
        zeros(1,par.WEC.ny_rad)];

%% %%%% Model Predictive Load Scheduling (MPLS) Parameters %%%%%%%%%%%%%%%%
% Scheduling parameters
par.dt_ctrl = 2;            % interval between control updates
par.MPLS.tc = 10;                  % control horizon
par.MPLS.tp = par.MPLS.tc + 1.5*par.dt_ctrl;      % prediction horizon

% Load control constraints
par.T_max = NDmeshVar.T_max(iVar);        % [Nm] max PTO reaction torque
deltat_Tmax = 3;
par.dTdt_max = par.T_max/deltat_Tmax;  % [Nm/s] max rate of change in WEC load
par.T_min = par.T_max*NDmeshVar.lbFrac(iVar);           % [Nm] min PTO reaction torque

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic_1 = tic; % record start time of optimization

if par.T_min > max(T_c_data(iSS,:))
    tMPLS = NaN;
    Tpto = NaN;
    PP = 0;
    dur = toc(tic_1);
    save(['data_study_yearlyAve_',num2str(iVar),'.mat'])
    return
elseif par.T_max > max(T_c_data(iSS,:))
    par.T_max = max(T_c_data(iSS,:));
end

% Update initial conditions post ramp
Tpto = 0.5*(par.T_max+par.T_min);
temp = model_OLloadSchedule([],y0,Tpto,[],par,4);
y0 = temp(end,:); clearvars temp

% Perform optimization
[tMPLS,Tpto] = modelPredictiveLoadScheduling(y0,par);

% record the average power absorption
PP = model_OLloadSchedule(tMPLS(1),y0,Tpto,tMPLS(end),par,1);

% record the time the optimization took
dur = toc(tic_1);

% Save results
save(['data_study_yearlyAve_',num2str(iVar),'.mat'])

end