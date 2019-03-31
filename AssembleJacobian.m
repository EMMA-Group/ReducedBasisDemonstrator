%% Assembles the Jacobi matrix of the RB method
% cf. equation (29) in the paper

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% AssembleJacobian.m
% Copyright (C) 2019 by Oliver Kunc
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% (the full license is distributed together with the software
% in a file name LICENSE)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This software package is related to the research article
%
% Authors: Oliver Kunc and Felix Fritzen
% Title:   Finite strain homogenization using a reduced basis and efficient sampling
% Journal: Mathematical and Computational Applications
%          Special Issue "Machine Learning, Low-Rank Approximations and Reduced Order Modeling in Computational Mechanics"
% Year:    2019
% Volume:  24
% Number:  1
% DOI   ...
% URL   dx.doi.org/...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Jacobian] = AssembleJacobian()                                   % cf. equation (29)
    global N N_qp B C w
    Jacobian = zeros(N,N);
    for i_qp=1:N_qp
        Jacobian = Jacobian + ...
            B((i_qp-1)*9+1:i_qp*9,:)' * ...
            reshape(C((i_qp-1)*81+1:i_qp*81),9,9)' * ...
            B((i_qp-1)*9+1:i_qp*9,:) * w(i_qp);
    end
end