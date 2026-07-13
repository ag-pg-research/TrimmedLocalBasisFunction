function [f, A] = construct_polynorm_basis(m, k)
% Function prepares normalized polynomial basis functions of a form 
% f_l(j) = (j / k)^(l-1)
%
% Inputs:
% m     number of polynomial basis functions
% k     half-width of the analysis window (admissible delay in number of
%       samples)
%
% Outputs:
% f     matrix (m x 2k+1) of system parameters estimates
% A     matrix (m x m) allowing for recursive calculations
%
% Author:   Artur Gancza, Gdansk University of Technology, Poland
% email:    artgancz@pg.edu.pl
    t = (-k:1:k) / k;
    f = ones(m, 2*k+1);
    A = eye(m);
    for ii = 2:1:m
        f(ii, :) = f(ii-1, :) .* t;
        for jj = 1:1:ii-1
            A(ii, jj) = nchoosek(ii-1, ii-jj) / (-k)^(ii-jj);
        end
    end
end