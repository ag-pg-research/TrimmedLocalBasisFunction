function B = eye_kron(n, A)
% The function performs a Kronecker product between an identity matrix and
% any other matrix. This function is faster than the Matlab "kron" function 
% used in the general case.
%
% Inputs:
% n     size of the identity matrix;
% A     matrix (m1 x m2)
%
% Outputs:
% B     Kronecker product of identity matrix and any other matrix
%
% Author:   Artur Gancza, Gdansk University of Technology, Poland
% email:    artgancz@pg.edu.pl
    [m1, m2] = size(A);
    B = zeros(n*m1, n*m2);
    for ii = 1:1:n
        B((ii-1)*m1+1:ii*m1, (ii-1)*m2+1:ii*m2) = A;
    end
end