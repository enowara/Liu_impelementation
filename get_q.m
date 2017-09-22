function q = get_q(pVec, N)

% define i and jk - to do C(N, 2)
q=[];
for i = 1:N
    p_i = pVec(i);
    q_i = [];
    for j = 1:N
        p_j = pVec(j);
        q_ij = p_i*p_j;
%         x_ij = xTemp(i,j);    % each element of NxN x
%         q_ij = x_ij*p_i*p_j;
        q_i = [q_i q_ij];
    end
    q = [q; q_i];
end

% append back into NxN matrix q
end
