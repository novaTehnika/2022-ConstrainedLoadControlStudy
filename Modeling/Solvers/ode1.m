function yout = ode1(F,t0,dt,tfinal,y0)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ode1.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 
%
% PURPOSE/DESCRIPTION:
% Created from ode1() written by MathWorks
%
% FILE DEPENDENCY: 
% parameters_WECmodel_V02x00.m
%
% UPDATES:
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
   
    t = t0:dt:tfinal;
    nt = length(t);
    yout = zeros(nt,length(y0));
    yout(1,:) = y0;
    
    for it = 2:nt
        dy = F(t(it),yout(it-1,:));
        yout(it,:) = yout(it-1,:) + dy.*dt;
    end
      
end    