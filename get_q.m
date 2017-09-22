function q = get_q(pVec, N)

% q should be size ignore repeating pairs 
q=[];
for i = 1:N % for each ROI
    p_i = pVec(i); 
    q_i = [];
    for j = i:N  % instead of 1:N, then don't keep repeating pairs 2,3 = 3,2
                 % instead N^2 features, get  C(N,2) +N
        p_j = pVec(j);
        q_ij = p_i*p_j;  % becuase assume ROIs are independent of each other, here i,j are ROIs
%         x_ij = xTemp(i,j);    % each element of NxN x
%         q_ij = x_ij*p_i*p_j;
        q_i = [q_i q_ij];
    end
    q = [q; q_i];
end

% append back into NxN matrix q
end
