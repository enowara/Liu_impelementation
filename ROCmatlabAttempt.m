
Train an SVM classifier on the same sample data. Standardize the data.

mdlSVM = fitcsvm(pred,resp,'Standardize',true);
Compute the posterior probabilities (scores).

mdlSVM = fitPosterior(mdlSVM);
[~,score_svm] = resubPredict(mdlSVM);
The second column of score_svm contains the posterior probabilities of bad radar returns.

Compute the standard ROC curve using the scores from the SVM model.

% [Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(resp,score_svm(:,mdlSVM.ClassNames),'true');                                                                                                                                                                                                                                               
[Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(Yts,scores_SVM(:,2), 1);                                                                                                                                                                                                                                               
figure, plot(Xsvm, Ysvm)
title('test data')
AUCsvm

[Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(Ytr,scores_SVM_post(:,2), 1);                                                                                                                                                                                                                                               
figure, plot(Xsvm, Ysvm)
title('test data')
AUCsvm