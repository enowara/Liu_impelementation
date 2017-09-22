function numEigenVec = find_90_perc_eigen(s, perc_var)

% input
% s - matrix of eigenvalues
% perc_var - percent variance to keep, if 90% , should be inputted as 0.9

% output
% numEigenVec - number of eienvectors to keep

i = 1;
lambdas = diag(s); % eigenvalues
sum_temp = 0;

var_90 = perc_var*sum(lambdas(:));
if lambdas(1) > var_90% check if the first eigenvector already contributes to 90% of var
    numEigenVec = 1;
else
    while sum_temp < var_90
        sum_temp1 = lambdas(i) + lambdas(i+1);
        sum_temp = sum_temp + sum_temp1;

        i = i+1;
    end

    numEigenVec = i+1; % number of eigenvectors to keep, add one b/c last sum had i + (i+1)
    
end

end