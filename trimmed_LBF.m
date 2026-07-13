function theta = trimmed_LBF(y, u, n, k, m, delta)
% Function estimates least squares of trimmed residuals combined with basis
% expansion method (trimmed local basis function) for the FIR model.
%
% Inputs:
% y     column vector containing output signal (T x 1)
% u     column vector of input signal (T x 1)
% n     FIR model order
% k     half-width of the analysis window (admissible delay in number of
%       samples)
% m     number of polynomial basis functions
% delta number of samples for trimming
%
% Outputs:
% theta matrix (n x T) of system parameters estimates
%
% Author:   Artur Gancza, Gdansk University of Technology, Poland
% email:    artgancz@pg.edu.pl
    T = length(y);
    K = 2*k+1;
    nm = n*m;
    K_tilde = K - delta;
    [f, A] = construct_polynorm_basis(m, k);
    fp = A\f(:, K);

    % Regression vectors
    phi = zeros(n, T);
    for ni = 1:1:n
        phi(ni, ni:T) = u(1:T-ni+1).';
    end
    
    theta = zeros(n, T);
    F0 = eye_kron(n, f(:,k+1)');
    
    psi = batch_col_kron(phi(:, 1:K), f);
    P = psi*psi';
    p = psi*conj(y(1:K));
    beta_0 = P \ p;
    epsilon = y(1:K).' - beta_0'*psi;
    [~, id] = sort(abs(epsilon));
    id_tr = id(1:K_tilde);
    P = psi(:, id_tr)*psi(:, id_tr)';
    p = psi(:, id_tr)*conj(y(id_tr));
    beta_0 = P \ p;
    epsilon = y(1:K).' - beta_0'*psi;
    [~, id] = sort(abs(epsilon));
    id_tr = id(1:K_tilde);
    for t = k+1:1:T-k
        local_y = y(t-k:t+k);
        P = psi(:, id_tr)*psi(:, id_tr)';
        p = psi(:, id_tr)*conj(local_y(id_tr));
        % Vector of hyperparameters
        beta = P\p;
        theta(:, t) = F0*beta;

        e = y(t-k+1:t+k).' - beta'*psi(:, 2:K);

        if t < T-k
            psi(m+1:nm, :) = psi(1:nm-m, :);
            psi(1:m, :) = phi(1, t-k+1:t+k+1) .* f;
            ep = y(t+k+1) - beta'*col_kron(phi(:, t+k+1), fp);
            epsilon = [e, ep];
            [~, id] = sort(abs(epsilon));
            id_tr = id(1:K_tilde);
        end
    end     
end