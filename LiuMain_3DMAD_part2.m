% after part 1 has already been run, which computed the S matrices with raw
% 
datte = '10-02';
saveLiuFolder = ['LiuMasks/' datte '/'];
mkdir(saveLiuFolder)

N = 16;
tic
%% split for LOOV 
    % implement the LOOV strategy. Loop over all people in the dataset
    % except one when training? 
startTestPerson1 = [1:5:85];
endTestPerson1 = [1:5:85] + 4;

startTestPerson2 = [];
startTestPerson3 = [];
endTestPerson2 = [];
endTestPerson3 = [];
    
allPeople = startTestPerson1(1):endTestPerson1(end); 
%     allPeople = startTestPerson1(1):endTestPerson3(end);
%     pEnd = 15;
pEnd = 17; 

% LOAD S matrices
load([saveLiuFolder 'Masks-LiuData-1.mat'])
S1 = S;
Mlist1 = Mlist;
load([saveLiuFolder 'Masks-LiuData-2.mat'])
S2 = S;
Mlist2 = Mlist;
load([saveLiuFolder 'Masks-LiuData-3.mat'])
Sfake_p = S;
Mlist3 = Mlist;
Mlist_fake_p = Mlist3;

% SliveAll = [S1 S2];

predictionAllSVM = [];
predictionAllSVMLive = [];
predictionAllSVMFake = [];
testPeople = [];
labelsSVM = [];
predtests = [];
Ytss = [];

for p = 1:pEnd
    
    % if Replay testPersonInit(1:4); change indicing a little
    
    testPerson1 = startTestPerson1(p):endTestPerson1(p);
    if isempty(startTestPerson2) ~= 1
        testPerson2 = startTestPerson2(p):endTestPerson2(p);
    else
        testPerson2 = [];
    end
     if isempty(startTestPerson2) ~= 1
        testPerson3 = startTestPerson3(p):endTestPerson3(p);
     else
         testPerson3 = [];
     end
%     testPerson = [testPerson1 testPerson2 testPerson3];
    testPersonInit = [testPerson1 testPerson2 testPerson3];
    
    testPerson = testPersonInit;%(1:4); % videos to leave out for testing, if live, 
                                          % those will not be considered
                                          % for learning p and q
    trainPeople = setdiff(allPeople, testPerson); % keep all videos for training 
                                                % and learning p and q
                                            % except j_leave ones
% get e and p iteratively
    
% TODO: fix indicing to loop over correct columns and keep correct Mlist

% find training people indices for N ROIs% Str_pStart = (trainPeople-1)*N+1;
Str_pStart = (trainPeople-1)*N+1;
Str_pEnd = (trainPeople)*N;
Str_idx = [];
for ll = 1:length(Str_pStart)
    Str_i = Str_pStart(ll):Str_pEnd(ll);
    Str_idx = [Str_idx Str_i];
end
      
    Slive_p1 = S1(:,Str_idx);% include each tr persons 16 ROIs
    Slive_p2 = S2(:,Str_idx);
    Slive_p_tr = [Slive_p1 Slive_p2];
    
    Sfake_p_tr = Sfake_p(:,Str_idx);
    
% find testing people indices     
Sts_pStart = (testPerson-1)*N+1;
Sts_pEnd = (testPerson)*N;
Sts_idx = [];
for ll = 1:length(Sts_pStart)
    Sts_i = Sts_pStart(ll):Sts_pEnd(ll);
    Sts_idx = [Sts_idx Sts_i];
end    

    Slive_p1_ts = S1(:,Sts_idx);% include each tr persons 16 ROIs
    Slive_p2_ts = S2(:,Sts_idx);
    Slive_p_ts = [Slive_p1_ts Slive_p2_ts];
    
    Sfake_p_ts = Sfake_p(:,Sts_idx);
    
    Mlist_live_p1_tr = Mlist1([trainPeople(1):trainPeople(end)], :);
    Mlist_live_p2_tr = Mlist2([trainPeople(1):trainPeople(end)], :);
    Mlist_live_p_tr = [Mlist_live_p1_tr; Mlist_live_p2_tr];
    
    Mlist_fake_p_tr = Mlist3([trainPeople(1):trainPeople(end)], :);
    
    Mlist_live_p1_ts = Mlist1([testPerson(1):testPerson(end)], :);
    Mlist_live_p2_ts = Mlist2([testPerson(1):testPerson(end)], :);
    Mlist_live_p_ts = [Mlist_live_p1_ts; Mlist_live_p2_ts];
    
    Mlist_fake_p_ts = Mlist3([testPerson(1):testPerson(end)], :);
    
    delta = 10^-3; % convergence threshold
        % r = 3; % bpm error toleration
    %     TODO: fix this part to determine indices defined as m list to leave out
%     k = 3; % keep top 3 eigenvectors, convert to 90 % variance instead
    alpha = 0.9; % in percent of variance

    warning off
    [pVec] =  iterate_p_e(Slive_p_tr, delta, N, alpha); % only computed for live 
%     % previously when keeping repeating pairs, q was 16 x 16, but it should
%     % be 1 x C(N,2)+N and then diagonalized, then q 16 x 16 was vectorized
%     % and diagonalized anyway to 16^2 x 16^2. Now, it'll be lower
%     % dimensional
    q = get_q(pVec,N); % only computed for live 
    Q =diag(q(:));   % Q = R'*R;
%     % R = sqrtm(Q);
% save([saveLiuFolder 'Masks-LiuData-Live-p-' num2str(p) '.mat'], 'pVec', 'q')
  
%% SVM

liveFolders = 1:2;
fakeFolders = 3;


testPersonLiv = testPerson; % change if Replay
testPersonAt = testPerson;
 [score_posterior, score, Yts, Ytr, labelSVM, predictionSVM, predictionSVMLive, ...
     predictionSVMFake] = SVM_Liu(saveLiuFolder, N, Q, liveFolders, ...
     fakeFolders, testPerson, testPersonLiv, testPersonAt, trainPeople, ...
     Slive_p_tr, Sfake_p_tr, Slive_p_ts, Sfake_p_ts, ...
     Mlist_live_p_tr, Mlist_fake_p_tr, Mlist_live_p_ts, Mlist_fake_p_ts);
     
     
   % append results from each LOOV run
scores_SVM_postcell{p} = score_posterior;
scores_SVMcell{p} = num2cell(score);
Ytsscell{p} = Yts;
Ytrscell{p} = Ytr;
testPeople = [testPeople; testPerson];
labelsSVMcell{p} = labelSVM;
predictionAllSVM = [predictionAllSVM; predictionSVM];
predictionAllSVMLive = [predictionAllSVMLive; predictionSVMLive];
predictionAllSVMFake = [predictionAllSVMFake; predictionSVMFake];
end % end p, for each LOOV person

% sum up accuracy over all test people    
        predictionAverageSVM = sum(predictionAllSVM)/length(predictionAllSVM);
        disp([num2str(predictionAverageSVM) '% Average SVM accuracy']);

        predictionAverageSVMLive = sum(predictionAllSVMLive)/length(predictionAllSVMLive);
        disp([num2str(predictionAverageSVMLive) '% Average Live SVM accuracy']);
        predictionAverageSVMFake = sum(predictionAllSVMFake)/length(predictionAllSVMFake);
        disp([num2str(predictionAverageSVMFake) '% Average Fake SVM accuracy']);

save([saveLiuFolder 'SVMresults_all_p.mat'], 'scores_SVM_postcell', 'scores_SVMcell', ...
    'Ytsscell', 'Ytrscell', 'testPeople', 'labelsSVMcell', 'predictionAllSVM', ...
    'predictionAllSVMLive', 'predictionAllSVMFake', 'predictionAverageSVM', ...
    'predictionAverageSVMLive', 'predictionAverageSVMFake')


toc