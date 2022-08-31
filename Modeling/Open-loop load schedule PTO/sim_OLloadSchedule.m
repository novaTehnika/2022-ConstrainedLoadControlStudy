% sim_OLloadcontrol.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 06/14/2022
%
% PURPOSE/DESCRIPTION:
% This script simulates a simple wave energy PTO with open-loop control of
% the load on the WEC. Tpto is an anonymous function giving the reaction
% torque from the PTO as a function of time.
%
% FILE DEPENDENCY:
% sys_OLloadSchedule.m
%
% UPDATES:
% 06/14/2022 - First version creation. Adapted from sim_coulombPTO.m
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
function out = sim_OLloadSchedule(tspan,y0,Tpto,par)
%% %%%%%%%%%%%%   SOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve for states of the dynamic system
 % Set-up solver

 % Solver options
    options = odeset('RelTol',par.odeSolverRelTol,...
                     'AbsTol',par.odeSolverAbsTol,...
                     'MaxOrder',3);
    options.MaxStep = par.MaxStep; % min(0.1*par.WEC.w(:));
    
    
 % Run solver    
    dt = par.MaxStep;
    y = ode1(@(t,y) sys(t,y,Tpto,par)',tspan(1),dt,tspan(2),y0(:));
    t = tspan(1) : dt : tspan(2);
%% %%%%%%%%%%%%   POST-PROCESS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters
    out.par = par;

    % determine number of time steps
    nt = length(t); % number of time steps

    % Extract system states from simulation results
%     out.t = t;
    out.y = y;
    
%     out.theta = y(:,1); % [rad] position of the WEC
    out.theta_dot = y(:,2); % [rad/s] angular velocity of the WEC
    
    % Post-process non-state variables and state derivatives
    syspost = @(t,y,par) sysPost(t,y,Tpto,par);
    
    parfor it = 1:nt 
        [~, nonState(it)] = syspost(t(it),y(it,:),par); 
        
        % Move WEC torque results up a level in nonState stucture so that
        % they can be used like arrays in assiging to the output structure
%         temp(it).T_hydroStatic = nonState(it).torqueWEC.hydroStatic;
%         temp(it).T_wave = nonState(it).torqueWEC.wave;
%         temp(it).T_rad = nonState(it).torqueWEC.radiation;  
                    
    end

    % State derivatives
%     out.dydt = dydt;
    
     % System input: wave elevation
%     out.waveElev = [nonState(:).waveElev]';
    
     % forces on WEC
    out.T_pto = [nonState(:).T_pto]';
%     out.T_hydroStatic = [temp(:).T_hydroStatic]';
%     out.T_wave = [temp(:).T_wave]';
%     out.T_rad = [temp(:).T_rad]';
        
%% Post-process analysis
    %% Energy analysis
    out.power.P_WEC = -out.T_pto.*out.theta_dot;

    %% Energy Balance
    
    %% Mass Balance


%% %%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dydt = sys(t,y,Tpto,par)
    	[dydt, ~ ] = sys_OLloadSchedule(t,y,Tpto,par);
    end

    function [dydt, nonState] = sysPost(t,y,Tpto,par)
    	[dydt, nonState] = sys_OLloadSchedule(t,y,Tpto,par);
    end

end
