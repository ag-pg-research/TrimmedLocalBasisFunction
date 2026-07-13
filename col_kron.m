function C = col_kron(A, B)
% The function performs a Kronecker product between two column vectors.
% This function is faster than the Matlab "kron" function used in the
% general case.
%
% Inputs:
% A     first column vector;
% B     second column vector;
%
% Outputs:
% C     Kronecker product of two column vectors
%
% Author:   Artur Gancza, Gdansk University of Technology, Poland
% email:    artgancz@pg.edu.pl
    n = length(A);
    m = length(B);
    C = ones(n*m, 1);
    for ii = 1:1:n
        C((ii-1)*m+1:ii*m) = A(ii)*B;
    end
end