function RBF_q = getRBF(U,V, Q) 
gamma = 1;

RBF_q = exp(-gamma*(sqrt((x1-x2)' * Q * (x1-x2)))^.2);
end







Dq = sqrt(x1-x2)' * Q * (x1-x2);
RBF_qxi_xj = exp(-gamma*Dq(x1, x2)^2)  % a number for each i and j ?

% for i 
%     for j 
% % generate i, j 
% % xi = x(i,:); % row vector
% % xj = x(j,:);
% 
% qi = q(i,:); % row vec of q
% Q = diag(qi);
% Dq = sqrt(xi-xj)' * Q * (xi-xj);
% RBF_qxi_xj = exp(-gamma*Dq(xi, xj)^2)  % a number for each i and j ?
% q_i = [q_i q_ij];
% RBF_qxi_xj
%     end
% end
end

