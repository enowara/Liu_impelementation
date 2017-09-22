function p = iterate_p_e(S, delta, k, J, N)
% k = 3;
% delta - convergence threshold
% S = t X N*J
% J; % number of subjects
% N = size(Sj,2);%size(S,2)/J;% number of grid points
% initialize p_prev (p(0)) - p(t-1)
% N = 10;
p_previ = sqrt(N)/N;
p_prevVec = p_previ*ones(1,N);
e_nextJ = []; %?
gVecj = [];
t = 1; 
counter = 1;
pConvex = 10000;
dif =sum(abs(p_prevVec - pConvex'));
while dif >= delta
    for j = 1:J 
        
        Sj = S(:,((j-1)*N+1):(j*N));  % T x N   %S(:,1:N);
        if length(find(isnan(Sj))) > 0 
            continue
        else
        T = size(Sj,1);
        P = diag(p_prevVec.^2);  % size T x T ( or N X N)

        Sigma = Sj*P*Sj'; % N x N b/c (NxT) (TxT) (TxN)
        [U, ~, ~] = svd(Sigma);
        phi = U(:,1:k);    % N x 3

        E_hat = Sj' * phi * phi';    % I x N
        % E_hat = [e_hati....e_hatN];

        %update e_next
        e_next = (mean(E_hat,1)');   % length I  %sum(e_hati)/N

        e_nextJ = [e_nextJ e_next];  % concatenate e_next for all people J?
    end
    end

 for j = 1:J
    g = getGoodness(Sj, e_next);   % g should be a vector   1 x N
    gVecj = [gVecj; g];      % matrix J x N
 end

    % update p_next
    % [~, pIdx] = max(g_whole*p_prevVec');  % sum among all people, update p for all people J
    % p_next = p_prevVec(pIdx);
    % diff = abs(p_next -p_prev ) ;
    % p_prev = p_next;
% for j = 1:J
gVecjAve = sum(gVecj,1);    % sum or average?

    n = N;
    clear pConvex
    cvx_begin quiet
        variable pConvex(n)
        
        minimize(-sum(gVecj*pConvex));
        subject to
        norm(pConvex) <= 1;
    cvx_end
% end 
    dif = sum(abs(p_prevVec - pConvex'));
    p_prevVec = pConvex';
    counter = counter +1;
end   
p = p_prevVec;
end

