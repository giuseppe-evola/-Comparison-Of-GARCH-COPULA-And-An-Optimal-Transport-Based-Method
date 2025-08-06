function [pi_opt, pi_matrix] = sinkhorn_mot_4(mu1, mu2, mu3, mu4, C, epsilon, maxIter, tol)

% SINKHORN_MOT_4_LOG Solves the MOT problem for 4 marginals in log domain
% Same inputs/outputs as original, but computations done in log domain

[n1, n2, n3, n4] = size(C);
assert(length(mu1) == n1 && length(mu2) == n2 && length(mu3) == n3 && length(mu4) == n4, ...
    'Marginal dimensions must match the cost tensor.');
assert(abs(sum(mu1) - 1) < 1e-8 && abs(sum(mu2) - 1) < 1e-8 && ...
       abs(sum(mu3) - 1) < 1e-8 && abs(sum(mu4) - 1) < 1e-8, ...
    'Marginals must sum to 1.');

mu1 = mu1(:); mu2 = mu2(:); mu3 = mu3(:); mu4 = mu4(:);

% Initialize in log domain
log_u1 = zeros(n1, 1); 
log_u2 = zeros(n2, 1); 
log_u3 = zeros(n3, 1); 
log_u4 = zeros(n4, 1);

% Scale and convert cost to log kernel
C_scaled = C / mean(C(:));
log_K = -C_scaled / epsilon;

% Convert marginals to log domain
log_mu1 = log(mu1);
log_mu2 = log(mu2);
log_mu3 = log(mu3);
log_mu4 = log(mu4);

for iter = 1:maxIter
    log_u1_old = log_u1; 
    log_u2_old = log_u2; 
    log_u3_old = log_u3; 
    log_u4_old = log_u4;

    % Update log_u1
    for i = 1:n1
        temp = squeeze(log_K(i,:,:,:));
        log_scaled = temp + ...
            reshape(log_u2, [1, n2, 1, 1]) + ...
            reshape(log_u3, [1, 1, n3, 1]) + ...
            reshape(log_u4, [1, 1, 1, n4]);
        % Use log-sum-exp trick for numerical stability
        max_val = max(log_scaled(:));
        log_marginal = max_val + log(sum(exp(log_scaled(:) - max_val)));
        log_u1(i) = log_mu1(i) - log_marginal;
    end

    % Update log_u2
    for j = 1:n2
        temp = squeeze(log_K(:,j,:,:));
        log_scaled = temp + ...
            reshape(log_u1, [n1, 1, 1, 1]) + ...
            reshape(log_u3, [1, 1, n3, 1]) + ...
            reshape(log_u4, [1, 1, 1, n4]);
        max_val = max(log_scaled(:));
        log_marginal = max_val + log(sum(exp(log_scaled(:) - max_val)));
        log_u2(j) = log_mu2(j) - log_marginal;
    end

    % Update log_u3
    for k = 1:n3
        temp = squeeze(log_K(:,:,k,:));
        log_scaled = temp + ...
            reshape(log_u1, [n1, 1, 1, 1]) + ...
            reshape(log_u2, [1, n2, 1, 1]) + ...
            reshape(log_u4, [1, 1, 1, n4]);
        max_val = max(log_scaled(:));
        log_marginal = max_val + log(sum(exp(log_scaled(:) - max_val)));
        log_u3(k) = log_mu3(k) - log_marginal;
    end

    % Update log_u4
    for l = 1:n4
        temp = squeeze(log_K(:,:,:,l));
        log_scaled = temp + ...
            reshape(log_u1, [n1, 1, 1]) + ...
            reshape(log_u2, [1, n2, 1]) + ...
            reshape(log_u3, [1, 1, n3]);
        max_val = max(log_scaled(:));
        log_marginal = max_val + log(sum(exp(log_scaled(:) - max_val)));
        log_u4(l) = log_mu4(l) - log_marginal;
    end

    % Check convergence in log domain
    err = max(abs([
        log_u1 - log_u1_old;
        log_u2 - log_u2_old;
        log_u3 - log_u3_old;
        log_u4 - log_u4_old
    ]));
    disp(err)
    if err < tol
        break;
    end
end

% Final transport plan in log domain
log_pi_matrix = zeros(size(log_K));
for i = 1:n1
    for j = 1:n2
        for k = 1:n3
            for l = 1:n4
                log_pi_matrix(i,j,k,l) = log_K(i,j,k,l) + ...
                    log_u1(i) + log_u2(j) + log_u3(k) + log_u4(l);
            end
        end
    end
end

% Convert back from log domain and normalize
pi_matrix = exp(log_pi_matrix);
pi_matrix = pi_matrix / sum(pi_matrix(:));
pi_opt = pi_matrix(:);

end