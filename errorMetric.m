% error metrics
% 
% AUC, ROC
matNames = dir([['Nov30SVM_MasksShuffled/'] ['*.mat']]);

for r = 1:length(matNames)
    load(['Nov30SVM_MasksShuffled/' matNames(r).name])

    % Compute the posterior probabilities (scores).

    % mdlSVM = fitPosterior(SVMModel);
    % [~,score_svm] = resubPredict(mdlSVM);

    [ScoreCVSVMModel,ScoreParameters]  = fitPosterior(SVMModel,Xtr,Ytr);
    [labelSVMPost,scorePost] = predict(ScoreCVSVMModel,Xts);

    [X,Y,T,AUC] = perfcurve(Yts,score(:,2),1) ;
    AUC
    figure
    plot(X,Y)
end

% By default, X is false
%     positive rate, FPR and Y is
%     true positive rate, TPR

%% false fake rate and false liveness rate
% false positives : 

% concatenate all fake and live labels and get a single FLR, etc per person
% or get a number for each person, keep and plot AUC?

% for each test person
LabelFakeIDX = find(Yts == 0);    % observations that are fake
FalsePos = length(find(labelSVM(LabelFakeIDX) == 1));  % indices
TrueNeg = length(find(labelSVM(LabelFakeIDX) == 0));

FFR = FalsePos / (FalsePos + TrueNeg) ; % # between 0 and 1

%%%%%%
LabelLiveIDX = find(Yts == 1);  
FalseNeg = length(find(labelSVM(LabelLiveIDX) == 0)); 
TruePos = length(find(labelSVM(LabelLiveIDX) == 1));

FLR = FalseNeg / (FalseNeg + TruePos);

HTER = (FFR+FLR)/2;

FFR_all = [];
FLR_all = [];
EER_all = [];
 % load each person's p prediction with SVM and RDF
 for p = 1:15 
    SVMpred = labelsSVM{p} ;
    label_p = Ytss{p};
    
    LabelFakeIDX = find(label_p == 0); 
    FalsePos = length(find(SVMpred(LabelFakeIDX) == 1));  % indices
    TrueNeg = length(find(SVMpred(LabelFakeIDX) == 0));
    FFR = FalsePos / (FalsePos + TrueNeg) ;
    
    LabelLiveIDX = find(label_p == 1);  
    FalseNeg = length(find(SVMpred(LabelLiveIDX) == 0)); 
    TruePos = length(find(SVMpred(LabelLiveIDX) == 1));
    FLR = FalseNeg / (FalseNeg + TruePos);

    HTER = (FFR+FLR)/2;
    
    FFR_all = [FFR_all; FFR];
    FLR_all = [FLR_all; FLR];
    EER_all = [EER_all; HTER];

 end
 
%  ROC
plot(FFR_all, FLR_all, '*')
% AUC
  
    EE - equal error rate, when acceptance and rejection errors are equal
    
    