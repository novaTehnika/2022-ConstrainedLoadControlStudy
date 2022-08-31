function [dydt, nonState] = sys_OLloadSchedule(t,y,Tpto,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys_OLloadSchedule.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/14/2022
%
% PURPOSE/DESCRIPTION:
% Calculate the state derivatives for a simple wave energy PTO with a 
% specified load on the from the PTO.
%
% FILE DEPENDENCY: 
%
% UPDATES:
% 06/14/2022 - First version creation. Adapted from sys_coulombPTO.m
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
iytheta = 1;
iytheta_dot = 2;
iyrad = (1:par.WEC.ny_rad) + iytheta_dot; % state vector indices for radiation damping states for WEC model
iWEC = [iytheta iytheta_dot iyrad];

par.Tcoulomb = Tpto(t);
nonState = nonStateVars(t,y,par); % recover nonstate variables

[dydt_WEC, nonState.torqueWEC, nonState.waveElev] = ...
                    flapModel(t,y(iWEC),nonState.T_pto,par);

dydt = zeros(2+par.WEC.ny_rad,1);
dydt(iytheta) = dydt_WEC(1); % angular velocity
dydt(iytheta_dot) = dydt_WEC(2); % angular acceleration
dydt(iyrad) = dydt_WEC(3:end); % radiation damping states for WEC model

%% %%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function nonState = nonStateVars(t,y,par)
        nonState.T_pto = (abs(y(iytheta_dot)) > 5e-3)* ...
                        -sign(y(iytheta_dot))*par.Tcoulomb;
    end
end
