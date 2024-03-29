% modelPredictiveLoadScheduling.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/15/2022
%
% PURPOSE:
% This function implements a model predictive algorithm for scheduling the
% load on a wave energy converter. It performs an optimization of the load
% schedule over a control and prediction horizon at each control update 
% step.
% The optimization starts at a specified time ("[par.tstart par.tend]") and
% with specified initial conditions ("y0").
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
function [t, Tpto] = modelPredictiveLoadScheduling(y0,par)
% Options for optimization routine fmincon()
% options = optimset('Display','iter','PlotFcns',@optimplotfval);
options = optimoptions('fmincon');
options.Display = 'iter';
options.UseParallel = true;
% options.Algorithm = 'sqp';
% options.ScaleProblem = true;
options.MaxFunEvals = 1e3;
options.MaxIter = 1e2;
%     options.TolFun = 1e-17;
%     options.TolCon = 1e-17;
options.OptimalityTolerance = 1e-20;
options.StepTolerance = 1e-10;

% Model predictive control parameters
dt_ctrl = par.dt_ctrl;  % Update period/step size
tc = par.MPLS.tc;       % Control horizon
tp = par.MPLS.tp;       % Prediction horizon
N = ceil(tc/dt_ctrl);   % number of control updates in control horizon

% Linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Upper and lower bounds of WEC load
lb = par.T_min*ones(1,N);
ub = par.T_max*ones(1,N);

t = (par.tstart-par.TrampMPLS):dt_ctrl:par.tend;
nt = length(t);
Tpto = zeros(nt,1);
Tpto(1) = 0.5*(par.T_max+par.T_min);
x0 = Tpto(1)*ones(size(lb));
y = y0;
for it = 1:nt-1
    
    % Build the objective function for the current case
    obj = @(Tbar) -model_OLloadSchedule(t(it),y,[Tpto(it) Tbar],tp,par,1);
    
    % Build the contraint function for the current case
    nonlcon = @(Tbar) model_OLloadSchedule(t(it),y,[Tpto(it) Tbar],tp,par,2);
    
    % Execute optimization
    tic
    x = fmincon(obj,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    toc

    Tpto(it+1) = x(1);      % selected control update
    x0 = [x(2:end) x(end)]; % next intial guess based on previous results

    % Update state
    temp = model_OLloadSchedule(t(it),y,[Tpto(it) x(1)],[],par,3);
    y = temp(end,:);
    
    % display progress of the optimization
    display([num2str(t(it)),'seconds of ',num2str(t(end))])
    display([num2str(x)])

end
