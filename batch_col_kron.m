function C = batch_col_kron(A, B)
% The function performs a Kronecker product between two batches of column 
% vectors. The result is C = [A(1,:).*B; ...;
%                             A(n, :).*B]            
%
% Inputs:
% A     matrix of first batch of column vectors;
% B     matrix of second batch of column vectors;
%
% Outputs:
% C     Kronecker product of two batches of column vectors
%
% Author:   Artur Gancza, Gdansk University of Technology, Poland
% email:    artgancz@pg.edu.pl
    [n, T] = size(A);
    m = size(B, 1);
    C = ones(n*m, T);
    for ii = 1:1:n
        C((ii-1)*m+1:ii*m, :) = A(ii, :).*B;
    end
end