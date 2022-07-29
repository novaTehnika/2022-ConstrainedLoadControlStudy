% Tpto_ramp.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/15/2022
%
% PURPOSE:
% This function provides a function of time defined by values specified at   
% a fixed interval. The value of the function is linearly interpolated 
% between the fixed intervals and, beyond the range of specified values, 
% takes the value of the last value specified.
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

function Tpto = Tpto_ramp(t,Tbar,dt)
    N = length(Tbar);
    n = floor(t/dt)+1;
    if n < N
        d = mod(t,dt)/dt;
        Tpto = d*(Tbar(n+1) - Tbar(n)) + Tbar(n);
    else
        Tpto = Tbar(N);
    end
end