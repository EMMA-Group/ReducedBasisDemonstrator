%% Contains a single matrix with helpful indices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% CreateIndexMatrix.m
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



% Let M be a 9x9-matrix stored as a row-major 81-vector. Then its
% components can be accessed as M(a) with 1<=a<=81. If it was in 9x9-format
% then that would be M(b,c) with a=9*(b-1)+c and 1<=b,c<=9. If it was in
% 3x3x3x3-format, it would be M(i,j,k,l) with [i,j,k,l]=indexmatrix(a,:)

    indexmatrix = ...
       [1, 1, 1, 1;
        1, 1, 1, 2;
        1, 1, 1, 3;
        1, 1, 2, 1;
        1, 1, 2, 2;
        1, 1, 2, 3;
        1, 1, 3, 1;
        1, 1, 3, 2;
        1, 1, 3, 3;%
        1, 2, 1, 1;
        1, 2, 1, 2;
        1, 2, 1, 3;
        1, 2, 2, 1;
        1, 2, 2, 2;
        1, 2, 2, 3;
        1, 2, 3, 1;
        1, 2, 3, 2;
        1, 2, 3, 3;%
        1, 3, 1, 1;
        1, 3, 1, 2;
        1, 3, 1, 3;
        1, 3, 2, 1;
        1, 3, 2, 2;
        1, 3, 2, 3;
        1, 3, 3, 1;
        1, 3, 3, 2;
        1, 3, 3, 3;%
        2, 1, 1, 1;
        2, 1, 1, 2;
        2, 1, 1, 3;
        2, 1, 2, 1;
        2, 1, 2, 2;
        2, 1, 2, 3;
        2, 1, 3, 1;
        2, 1, 3, 2;
        2, 1, 3, 3;%
        2, 2, 1, 1;
        2, 2, 1, 2;
        2, 2, 1, 3;
        2, 2, 2, 1;
        2, 2, 2, 2;
        2, 2, 2, 3;
        2, 2, 3, 1;
        2, 2, 3, 2;
        2, 2, 3, 3;%
        2, 3, 1, 1;
        2, 3, 1, 2;
        2, 3, 1, 3;
        2, 3, 2, 1;
        2, 3, 2, 2;
        2, 3, 2, 3;
        2, 3, 3, 1;
        2, 3, 3, 2;
        2, 3, 3, 3;%
        3, 1, 1, 1;
        3, 1, 1, 2;
        3, 1, 1, 3;
        3, 1, 2, 1;
        3, 1, 2, 2;
        3, 1, 2, 3;
        3, 1, 3, 1;
        3, 1, 3, 2;
        3, 1, 3, 3;%
        3, 2, 1, 1;
        3, 2, 1, 2;
        3, 2, 1, 3;
        3, 2, 2, 1;
        3, 2, 2, 2;
        3, 2, 2, 3;
        3, 2, 3, 1;
        3, 2, 3, 2;
        3, 2, 3, 3;%
        3, 3, 1, 1;
        3, 3, 1, 2;
        3, 3, 1, 3;
        3, 3, 2, 1;
        3, 3, 2, 2;
        3, 3, 2, 3;
        3, 3, 3, 1;
        3, 3, 3, 2;
        3, 3, 3, 3];%
    