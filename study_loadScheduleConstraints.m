% study_loadScheduleConstraints.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 07/7/2022
%
% PURPOSE/DESCRIPTION:
% This script serves as a shell for performing parameter variation studies
% for WEC load scheduling performance using the function m-file 
% modelPredictiveLoadScheduling.m. The parameters to be varied are the 
% upper limit and the limit to the rate of change. The lower limit will be 
% fixed as a fraction of the upper limit for each run of the study.
%
% This script is called as a function with the arguements
% - variable-set index
% - sea-state
% - the lower load limit as a fraction of the upper limit as 
%
% FILE DEPENDENCY:
% modelPredictiveLoadScheduling.m
% model_OLloadSchedule.m
% parameters_OLloadSchedule.m
% SSdata_HumboltBay_1D.mat
%
% UPDATES:
% 07/07/2022 - Created.
% 07/12/2022 - Changed study variable from maximum rate of change in load 
% to the minimum time taken to go from zero to max load (dTdt_max to 
% deltat_Tmax) and updated the calculation of dTdt_max accordingly.
% 05/08/2023 - added specification for a ramp period for the MPLS algorithm
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
function [] = study_loadScheduleConstraints(iVar,iSS)
path_models = ['.' filesep 'Modeling'];
addpath([path_models filesep 'Solvers'])
addpath([path_models filesep 'Sea States'])
addpath([path_models filesep 'WEC model']) 
addpath([path_models filesep 'WEC model' filesep 'WECdata']) 
addpath([path_models filesep 'Open-loop load schedule PTO']) 

%% %%%%%%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation length
par.Tramp = 150; % [s] excitation force ramp period
par.TrampMPLS = 200; % [s] MPLS ramp period
par.tstart = 0; %[s] start time of simulation for results
par.tend = 2000; %[s] end time of simulation for results

% Solver parameters
par.odeSolverRelTol = 1e-9; % Rel. error tolerance parameter for ODE solver
par.odeSolverAbsTol = 1e-9; % Abs. error tolerance parameter for ODE solver
par.MaxStep = 1e-2; % Step size in the case of fixed step size solver

% Wave construction parameters
par.WEC.nw = 1000; % num. of frequency components for harmonic superposition 
par.wave.rngSeedPhase = 3; % seed for the random number generator

% Sea State specification
% Tp = [7.31 9.86 11.52 12.71 15.23 16.50];
% Hs = [2.34 2.64 5.36 2.05 5.84 3.25];
load('SSdata_HumboltBay_1D.mat','Hs','Tp');
% nSS = length(Tp);
par.wave.Hs = Hs(iSS);
par.wave.Tp = Tp(iSS);

% load remaining parameters
par = parameters_OLloadSchedule(par,...
    'nemohResults_vantHoff2009_20180802.mat','vantHoffTFCoeff.mat');

% Define initial conditions
y0 = [  0, ...
        0, ...
        zeros(1,par.WEC.ny_rad)];

%% %%%% Model Predictive Load Scheduling (MPLS) Parameters %%%%%%%%%%%%%%%%
par.dt_ctrl = 2;            % interval between control updates
par.MPLS.tc = 10;                  % control horizon
par.MPLS.tp = par.MPLS.tc + 1.5*par.dt_ctrl;      % prediction horizon

%% %%%%%%%%%%%%   Study Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_max = logspace(log10(1e5),log10(7e6),25);
deltat_Tmax = logspace(log10(0.3),log10(300),4);
lbFrac = [0 0.25 0.5 0.75];
[meshVar.T_max,meshVar.deltat_Tmax,meshVar.lbFrac] = ...
                                        meshgrid(T_max,deltat_Tmax,lbFrac);

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic_1 = tic; % record start time of optimization

% Load control parameters
par.T_max = meshVar.T_max(iVar);        % [Nm] max PTO reaction torque
par.dTdt_max = par.T_max/meshVar.deltat_Tmax(iVar);  % [Nm/s] max rate of change in WEC load
par.T_min = par.T_max*meshVar.lbFrac(iVar);           % [Nm] min PTO reaction torque

% Update initial conditions post ramp
Tpto = 0.5*(par.T_max+par.T_min);
temp = model_OLloadSchedule([],y0,Tpto,[],par,4);
y0 = temp(end,:); clearvars temp

% Perform optimization
[tMPLS,Tpto] = modelPredictiveLoadScheduling(y0,par);

% record the average power absorption post MPC ramp
it_start = find(tMPLS >= par.tstart,1,'first');
PP = model_OLloadSchedule(tMPLS(it_start),y0,Tpto(it_start:end),tMPLS(end),par,1);

% record the time the optimization took
dur = toc(tic_1);

% Save results
save(['data_study_loadScheduleConstraints_SS',num2str(iSS),'_',num2str(iVar),'.mat'])

end