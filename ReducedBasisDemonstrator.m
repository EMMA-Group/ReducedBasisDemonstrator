%% Main program of the Reduced Basis Demonstrator

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% ReducedBasisDemonstrator.m
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
%% erase workspace and set global variables
clear all
global N B w N_qp xi W P C F                                               % for better readability, most quantities are global
%% read Reduced Basis matrix and the weight vector
N = 10;                                                                    % number of basis elements
B = [dlmread('data/example_Bmatrix_material1.txt');dlmread('data/example_Bmatrix_material2.txt')]; % read the RB from files, one for each material
B = B(:,1:N);                                                              % adjust the size of the RB. now size(B) = [9*N_qp, N]
w =  dlmread('data/example_weights_material1.txt');                        % read the quadrature weights for the first material
N_qp_phase1 = length(w);                                                   % store the number of quadrature points of this material
w = [w; dlmread('data/example_weights_material2.txt')];                    % read the remaining quadrature weights. now sum(w)=1.
N_qp        = length(w);                                                   % total number of quadrature points
w9 = reshape(repmat(w,1,9)',9*N_qp,1);                                     % repeat weigts 9 times at every quadrature point, store as single column of length 9*N_qp
%% read deviatoric directions and choose one (row)
directions = dlmread('data/directions_16.txt');                            % reach row is a vector in R^5 with norm 1
direction = directions(1,:);                                               % choose any of these directions, e.g. the first row
%% set deviatoric basis (cf. equation (A10))
X = { sqrt(1/6)*[2,0,0;0,-1,0;0,0,-1], sqrt(0.5)*[0,0,0;0,1,0;0,0,-1], sqrt(0.5)*[0,1,0;1,0,0;0,0,0], sqrt(0.5)*[0,0,1;0,0,0;1,0,0], sqrt(0.5)*[0,0,0;0,0,1;0,1,0] };
%% set boundary condition Ubar (cf. equation (45))
Jbar = 1.001; t = 0.1;
Ubar = Jbar^(1/3) * expm( t*( direction(1)*X{1} + direction(2)*X{2} + direction(3)*X{3} + direction(4)*X{4} + direction(5)*X{5} ) )
%% set initial values
xi = zeros(N,1);                                                           % initial guess of coefficients
W = zeros(N_qp,1);                                                         % energy
P = zeros(N_qp*9,1);                                                       % first Piola-Kirchhoff stress
C = zeros(N_qp*81,1);                                                      % stiffness w.r.t. deformation gradient
Jacobi = zeros(N,N);                                                       % Jacobi matrix
ct_newt = 0;                                                               % counts the number of Newton iterations
%% Newton loop
while ct_newt<20
    ct_newt = ct_newt + 1;                                                 % increment and display Newton iteration
    F = repmat(reshape(Ubar',9,1),N_qp,1) + B*xi;                          % localize F (cf. equation (18))
    MaterialLaw(N_qp_phase1);                                              % evaluate W, P, C at every quadrature point (cf. equations (24), (25), (26))
    r = B' * (P.*w9);                                                      % compute residual (cf. equation (28))
    if norm(r,2)<1e-7, break;                                              % norm of residual small enough? then stop loop.
    else,              fprintf(1,'Newton iteration #%2i: residual norm = %4.2e\n', ct_newt, norm(r,2)); end                                          % check norm of residual
    Jacobi = AssembleJacobian();                                           % build the Jacobi matrix, cf. equation (29)
    xi = xi - (Jacobi\r);                                                  % update the coefficient vector, cf. equation (30)
end
%% effective stress
Pbar = reshape(sum( reshape(P,9,N_qp) .* reshape(w9,9,N_qp) , 2 ),3,3)'
