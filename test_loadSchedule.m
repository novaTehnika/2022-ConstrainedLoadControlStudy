% test_loadSchedule.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 06/14/2022
%
% PURPOSE/DESCRIPTION:
% This script was used in developing and testing the a model predicitve 
% load scheduling algorithm.
%
% FILE DEPENDENCY:
% sim_coulombPTO_V02x00.m
% parameters_coulombPTO_V01x00.m
%
% UPDATES:
% 06/14/2022 - (V01x00) - Created from studyPTO_V02x03 found in the 2021Q2 
% pipeline modeling project folder.
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
clear
clc
path_models = ['.' filesep 'Modeling'];
addpath([path_models filesep 'Sea States']) 
addpath([path_models filesep 'Solvers']) 
addpath([path_models filesep 'WEC model']) 
addpath([path_models filesep 'WEC model' filesep 'WECdata']) 
addpath([path_models filesep 'Open-loop load schedule PTO']) 
%% %%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Simulation Parameters
par.tramp = 100; % [s] excitation force ramp period
par.tstart = 0; %[s] start time of simulation
par.tend = 100; %[s] end time of simulation

par.odeSolverRelTol = 1e-9; % Rel. error tolerance parameter for ODE solver
par.odeSolverAbsTol = 1e-9; % Abs. error tolerance parameter for ODE solver
par.MaxStep = 1e-2;

par.wave.rngSeedPhase = 3; % set the seed for the random number generator

 % set the number of frequency components used to calculate the wave 
 % elevation and excitation force using harmonic superposition 
par.WEC.nw = 100;  

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Tp = [7.31 9.86 11.52 12.71 15.23 16.50];
% Hs = [2.34 2.64 5.36 2.05 5.84 3.25];
load('SSdata_HumboltBay_1D.mat')
nSS = length(Tp);

nVar1 = 1;
Tcoulomb = 1e6;%1e6*logspace(log10(0.1),log10(20),nVar1);% [Nm] PTO reaction torque

saveSimData = 1;
for iSS = 7%:nSS
    par.wave.Hs = Hs(iSS);
    par.wave.Tp = Tp(iSS);
    
    % load parameters
    par = parameters_OLloadSchedule(par);
    
    for iVar1 = 1:nVar1
        param = par;
        
        %% test model_OLloadSchedule.m
        if 0
        % define start and end times
        tstart = param.tstart;
        tend = param.tend;
        tspan = [tstart; tend];

        % Define intial conditions
        y0 = [  0, ...
                0, ...
                zeros(1,param.WEC.ny_rad)];

        
        param.dt_ctrl = 2;
        param.T_max = Tcoulomb(iVar1);
        param.dTdt_max = param.T_max/5;
        Tbar = Tcoulomb(iVar1)*ones(10,1);
        tc = param.dt_ctrl*length(Tbar); % control horizon
        tp = tend - tstart; % prediction horizon
        tic
        PPmean = model_OLloadSchedule(tstart,y0,Tbar,[tc tp],param,1);
        toc

        tic
        [c, ceq] = model_OLloadSchedule(tstart,y0,Tbar,[tc tp],param,2);
        toc

        tic
        yend = model_OLloadSchedule(tstart,y0,Tbar,[tc tp],param,3);
        toc
        end
        %% test sim_OLloadSchedule with constant 
        if 0
        % define anonymous function for open-loop load control
        Tpto = @(t) Tcoulomb(iVar1);

        % run simulation
        tic
        out = sim_OLloadSchedule(tspan,y0,Tpto,param);
        toc
        % Post-process
        PPmean(iSS,iVar1) = mean(out.power.P_WEC);%mean(out.theta_dot.*out.T_pto);

        if saveSimData
            simOut(iSS,iVar1) = out;
        end

        end
        %% test modelPredictiveLoadScheduling.m
        if 1
        % Load control parameters
        lbFrac = 0.25;
        param.T_max = Tcoulomb(iVar1);
        param.dTdt_max = 1e6;
        param.T_min = param.T_max*lbFrac;           % [Nm] min PTO reaction torque
        % Define intial conditions
        y0 = [  0, ...
                0, ...
                zeros(1,param.WEC.ny_rad)];

        % Model predictive control parameters
        param.dt_ctrl = 2;          % interval between control updates
        param.MPLS.tc = 10;               % control horizon
        param.MPLS.tp = param.MPLS.tc + 1.5*param.dt_ctrl;  % prediction horizon
        
        ticMPLS = tic;
        [tMPLS,Tpto] = modelPredictiveLoadScheduling(y0,param);
        toc(ticMPLS)
        
        figure; plot(tMPLS,Tpto)
        end
    end
end

% save('data_coulombPTO_dampingStudy_20220601.mat')
return
