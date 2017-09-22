function [pVec] =  iterate_p_e(S, delta, k, r, freq)
% delta = 10^-3; % convergence threshold
% k = 3;
% s - vector of raw s_i signals, concatenated to a row vector?
% J - all the training people
% output - pVec - vector of final p_i's

N = size(S,2);%length(S); %??

% initialize 
t = 1;
p_prev = sqrt(N)/N;
p_diff  = 100000;
while abs(p_diff) >= delta
     for j = 1:J
    p_squared = p_prev.^2;
    P = diag(p_squared);
    sigma  = S'*P*S;

    % PCA on sigma
    [U,SS,~] = svd(sigma);
    phi = U(:,k);  %  ?
    % use U or V to compute the eigenvectors?
    E = (((S))*phi)*phi'; % what is the size of E?
    e_sum = 0;
    m = length(E/N);
    for n_e = 1:m
        e_hat = E((N*(n_e-1)+1):(n_e*N));
        e_sum = e_sum + e_hat ;
    end

    e_next = 1/N * e_sum;
    i = 1; % i or j??
    s_i = S(:,i);
    % p_next = % update with e_next ???????????????????????
     [goodness] = getGoodness(s_i, e_next, freq, r);
     
    argmax_p = dot(p_prev, goodness(s,e));



    p_diff = p_next - p_prev;
    p_prev = p_next;
    end
end