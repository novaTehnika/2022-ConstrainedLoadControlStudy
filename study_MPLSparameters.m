% study_MPLSparameters.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 06/14/2022
%
% PURPOSE/DESCRIPTION:
% This script serves as a shell for performing parameter variation studies
% using simulations executed by the simPTO_VXXxXX.m function m-file.
% Parameter and function initiallization fuctions are called within this
% script before the sim_coulombPTO_VXXxXX script is called. 
%
% FILE DEPENDENCY:
% modelPredictiveLoadScheduling.m
% model_OLloadSchedule.m
% parameters_OLloadSchedule.m
% SSdata_HumboltBay_1D.mat
%
% UPDATES:
% 06/14/2022 - Created.
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
function [] = study_MPLSparameters(iVar,iSS)
path_models = ['.' filesep 'Modeling'];
addpath([path_models filesep 'Solvers'])
addpath([path_models filesep 'Sea States'])
addpath([path_models filesep 'WEC model'])
addpath([path_models filesep 'WEC model' filesep 'WECdata'])
addpath([path_models filesep 'Open-loop load schedule PTO'])

%% %%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation length
par.Tramp = 250; % [s] excitation force ramp period
par.tstart = 0; %[s] start time of simulation
par.tend = 500; %[s] end time of simulation

% Solver parameters
par.odeSolverRelTol = 1e-9; % Rel. error tolerance parameter for ODE solver
par.odeSolverAbsTol = 1e-9; % Abs. error tolerance parameter for ODE solver
par.MaxStep = 1e-2;

% Wave construction parameters
par.WEC.nw = 1000;
par.wave.rngSeedPhase = 3; % set the seed for the random number generator

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

% Define intial conditions
y0 = [  0, ...
        0, ...
        zeros(1,par.WEC.ny_rad)];

%% %%%%%%%%%%%%   Study Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt_ctrl = logspace(log10(0.5),log10(3),10);
tc = logspace(log10(3),log10(30),10);
[meshVar.dt_ctrl,meshVar.tc] = meshgrid(dt_ctrl,tc);

% Load control parameters
par.T_max = 5e6;    % [Nm] PTO reaction torque
par.dTdt_max = par.T_max/5;
par.T_min = par.T_max*0;           % [Nm] min PTO reaction torque

% Model predictive control parameters
par.dt_ctrl = meshVar.dt_ctrl(iVar);            % interval between control updates
par.MPLS.tc = meshVar.tc(iVar);                  % control horizon
par.MPLS.tp = par.MPLS.tc + 1.5*par.dt_ctrl;      % prediction horizon

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic_1 = tic; % record start time of optimization

% Update initial conditions post ramp
Tpto = 0.5*(par.T_max+par.T_min);
temp = model_OLloadSchedule([],y0,Tpto,[],par,4);
y0 = temp(end,:); clearvars temp

% Perform optimization
[tMPLS,Tpto] = modelPredictiveLoadScheduling(y0,par);

% record the average power absorption
PP = model_OLloadcontrol(tMPLS(1),y0,Tpto,tMPLS(end),par,1);

% record the time the optimization took
dur = toc(tic_1);

% Save results
save(['data_test_MPLSparameters',num2str(iVar),'.mat'])

end