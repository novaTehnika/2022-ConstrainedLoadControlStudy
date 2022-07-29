% model_OLloadSchedule.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/15/2022
%
% PURPOSE:
% This function implements three functions for the model predictive load
% scheduling algorithm based on the input "outputConfig".
% outputConfig: 
%   1 - perform a simulation over the prediction horizon for the model 
%   predictive loadsensing algorithm and output the average power captured
%   2 - perform a simulation for a state update through the control update
%   step
%   3 - evaluate constraints on the load schedule
% The inputs include all variables to define a simulation and include a
% specified load schedule.
%
% FUNCTION M-FILES
%
% UPDATES
% 6/15/2022 - created.
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

function varargout = model_OLloadcontrol(t,y,Tbar,tp,par,outputConfig)
    switch outputConfig
        case {1 3}
            % run simulation
            switch outputConfig
                case 1; tspan = [t,t+tp];
                case 3; tspan = [t,t+par.dt_ctrl];
            end
            
            Tpto = @(tprime) Tpto_ramp(tprime-t,Tbar,par.dt_ctrl);
            out = sim_OLloadSchedule(tspan,y,Tpto,par);

            switch outputConfig
                case 1; varargout = {mean(out.power.P_WEC)};
                case 3; varargout = {out.y};
            end

        case 2
            % constraints c <= 0, ceq = 0
            N = length(Tbar) - 1;
            c = zeros(N,1);
            if N > 1
                for n = 2:N
                    c(n) = abs(Tbar(n+1)-Tbar(n)) ...
                           - par.dTdt_max*par.dt_ctrl;
                end
            else
                c(1) = 0;
            end
            ceq = [];
            varargout = {c,ceq};

    end

end