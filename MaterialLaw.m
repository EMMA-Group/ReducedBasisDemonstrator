%% Saint Venant-Kirchhoff model w.r.t. deformation gradient F
% This is modular and not actually considered to be part of the code
% (i.e. lines are not counted)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% MaterialLaw.m
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
% Number:  2
% DOI   :  10.3390/mca24020056
% URL   :  https://dx.doi.org/10.3390/mca24020056
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = MaterialLaw(N_qp_phase1)
%% initialize variables
    global F W P C N_qp
    persistent C_local S_local
    if isempty(C_local)
        C_local = zeros(81,1);
    end
    if isempty(S_local)
        S_local = zeros(3,3);
    end
    CreateIndexmatrix;                                                     % create helpful index matrix
%% material parameters and stiffness tensor w.r.t. Green's strain in vector notation, for both materials
    E=[10,70000]; nu=[0.49,0.3];
    K = E./(3.*(1.-2.*nu));
    mu = E./(2.*(1.+nu));
    lambda = K - 2.*mu/3.;
    C_EGreen_matrix = zeros(9,9,2);
    C_EGreen_matrix(:,:,1) = ...
        [ lambda(1) + 2*mu(1),     0,     0,     0,           lambda(1),     0,     0,     0,           lambda(1)
                            0, mu(1),     0, mu(1),                   0,     0,     0,     0,                   0
                            0,     0, mu(1),     0,                   0,     0, mu(1),     0,                   0
                            0, mu(1),     0, mu(1),                   0,     0,     0,     0,                   0
                    lambda(1),     0,     0,     0, lambda(1) + 2*mu(1),     0,     0,     0,           lambda(1)
                            0,     0,     0,     0,                   0, mu(1),     0, mu(1),                   0
                            0,     0, mu(1),     0,                   0,     0, mu(1),     0,                   0
                            0,     0,     0,     0,                   0, mu(1),     0, mu(1),                   0
                    lambda(1),     0,     0,     0,           lambda(1),     0,     0,     0, lambda(1) + 2*mu(1)];
    C_EGreen_matrix(:,:,2) = ...
        [ lambda(2) + 2*mu(2),     0,     0,     0,           lambda(2),     0,     0,     0,           lambda(2)
                            0, mu(2),     0, mu(2),                   0,     0,     0,     0,                   0
                            0,     0, mu(2),     0,                   0,     0, mu(2),     0,                   0
                            0, mu(2),     0, mu(2),                   0,     0,     0,     0,                   0
                    lambda(2),     0,     0,     0, lambda(2) + 2*mu(2),     0,     0,     0,           lambda(2)
                            0,     0,     0,     0,                   0, mu(2),     0, mu(2),                   0
                            0,     0, mu(2),     0,                   0,     0, mu(2),     0,                   0
                            0,     0,     0,     0,                   0, mu(2),     0, mu(2),                   0
                    lambda(2),     0,     0,     0,           lambda(2),     0,     0,     0, lambda(2) + 2*mu(2)];
%% loop over quadrature points
    for i_qp=1:N_qp
        %% set up auxiliary index sets and kinematic quantities
        idx_9 = (i_qp-1)*9+1:i_qp*9;
        idx_81= (i_qp-1)*81+1:i_qp*81;
        F_qp = reshape(F(idx_9),3,3)'; % deformation gradient at this quadrature point
        E_Green = .5*(F_qp'*F_qp-eye(3)); % Green's strain tensor
        %% determine which material we are in -> determines material constants
        if i_qp<=N_qp_phase1
            i_phase = 1;
        else
            i_phase = 2;
        end
        %% energy
        W(i_qp) = lambda(i_phase)/2*trace(E_Green)^2 + mu(i_phase)*norm(E_Green,'fro')^2;
        %% stress
        S_local  = lambda(i_phase)*trace(E_Green)*eye(3)+2*mu(i_phase)*E_Green; % second Piola-Kirchhoff stress
        P(idx_9) = reshape( (F_qp*S_local)',9,1 );                              % first Piola-Kirchhoff stress
        %% stiffness
        for i_comp=1:81 % loop components of 9x9-matrix in row-major format
            % get indices [i,j,k,l] of a 9x9 matrix that's actually a 3x3x3x3 tensor
            i = indexmatrix(i_comp,1);
            j = indexmatrix(i_comp,2);
            k = indexmatrix(i_comp,3);
            l = indexmatrix(i_comp,4);
            
            % "geometric stiffness"
            if i==k
                C_local(i_comp) = S_local(l,j);
            else
                C_local(i_comp) = 0;
            end
            % Rayleigh-like transformation of the material stiffness
            for m=1:3
                for p=1:3
                    C_local(i_comp) = C_local(i_comp) + ...
                        F_qp((i-1)*3+m) * C_EGreen_matrix((m-1)*3+j, (p-1)*3+l, i_phase) * F_qp((k-1)*3+p);
                end
            end
        end
        C(idx_81)= C_local;
    end
end