% after part 1 has already been run, which computed the S matrices with raw
% 
datte = '10-03';
saveLiuFolder = ['LiuReplay/' datte '/'];
mkdir(saveLiuFolder)

N = 16;
% tic
%% split for LOOV 
    % implement the LOOV strategy. Loop over all people in the dataset
    % except one when training? 
    startTestPerson1 = (1:4:(60));
    startTestPerson2 = (1:4:(60))+60;
    startTestPerson3 = (1:2:(60))+120 ;
    
    endTestPerson1 =   ((1:4:(60)) + 3 );
    endTestPerson2 =     (((1:4:(60))+60) +3);
    endTestPerson3 =     (((1:2:(60))+120) +1); 
    
    allPeople = startTestPerson1(1):endTestPerson3(end);

%     allPeople = startTestPerson1(1):endTestPerson3(end);
%     pEnd = 15;
pEnd = 15; 

% LOAD S matrices
load([saveLiuFolder 'Replay-LiuData-1.mat'])
S1 = S;
Mlist1 = Mlist;
NameList1 = NameListAll;
clear NameListAll

load([saveLiuFolder 'Replay-LiuData-2.mat'])
S2 = S;
Mlist2 = Mlist;
NameList2 = NameListAll;
clear NameListAll

load([saveLiuFolder 'Replay-LiuData-3.mat'])
S3 = S;
Mlist3 = Mlist;
NameList3 = NameListAll;
clear NameListAll

NameList_Live = [NameList1];
NameList_Attack = [NameList2; NameList3];

% Mlist_fake_p = Mlist3;

% SliveAll = [S1 S2];

predictionAllSVM = [];
predictionAllSVMLive = [];
predictionAllSVMFake = [];
testPeople = [];
labelsSVM = [];
predtests = [];
Ytss = [];


% choose which light and motion set up to train and test on
train_cases_All = {'all', 'specific'}; 
light_condition_All = {'adverse', 'controlled', 'all_light'};
attack_All = {'photo', 'video'};

for t1 = 1:length(train_cases_All)
    train_cases = train_cases_All{t1};
    for t2 = 1:length(light_condition_All)
        light_condition = light_condition_All{t2};
        for t3 = 1:length(attack_All)
            attack = attack_All{t3};

for p = 1%:pEnd
    
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

% choose the subset for testing -  live

testPeople_Liv = [];
    for tt = testPersonInit 
        TsNameList = NameList_Live(tt).name;
        if  strcmp(light_condition,'all_light')
            idx_match = TsNameList; 
       elseif strcmp(light_condition,'adverse')
           idx_match = strfind(TsNameList, 'adverse');
       elseif strcmp(light_condition,'controlled')
           idx_match = strfind(TsNameList, 'controlled');
        end
        
       if isempty(idx_match)
           idx_match = 0;
       end
       testPeople_Liv = [testPeople_Liv idx_match];
    end

    testPeople_idx_keep1 = find((testPeople_Liv) ~=0); % choose non zero indices
    testPersonLiv = intersect(testPersonInit,testPeople_idx_keep1);
% choose the subset for testing -  face attacks

    testPeople_At = [];
    for tt = testPersonInit 
        TsNameList = NameList_Attack(tt).name;
        if strcmp(attack,'photo') & strcmp(light_condition,'all_light')
            idx_match = strfind(TsNameList, 'photo'); % match all the name strings containing photo
       elseif strcmp(attack,'video') & strcmp(light_condition,'all_light')
           idx_match = strfind(TsNameList, 'video');
       elseif strcmp(attack,'photo') & strcmp(light_condition,'adverse')
           idx_match = strfind(TsNameList, 'photo_adverse');
       elseif strcmp(attack,'video') & strcmp(light_condition,'adverse')
           idx_match = strfind(TsNameList, 'video_adverse');
       elseif strcmp(attack,'photo') & strcmp(light_condition,'controlled')
           idx_match = strfind(TsNameList, 'photo_controlled');
       elseif strcmp(attack,'video') & strcmp(light_condition,'controlled')
           idx_match = strfind(TsNameList, 'video_controlled');
        end
       if isempty(idx_match)
           idx_match = 0;
       end
       testPeople_At = [testPeople_At idx_match];
    end

    testPeople_idx_keep2 = find((testPeople_At) ~=0); % choose non zero indices
    testPersonAt = intersect(testPersonInit,testPeople_idx_keep2);
    
    testPerson = [testPersonLiv testPersonAt]; % combine correct indices from live and attack for testing
    testPerson = sort(testPerson, 'ascend');

if strcmp(train_cases,'all')
% if keeping all situations for training
    trainPeople = setdiff(allPeople, testPerson);
    trainPeopleLiv = trainPeople;
    trainPeopleAt = trainPeople;
elseif strcmp(train_cases,'specific')
    trainPeopleInit = setdiff(allPeople, testPerson); % m values

    
    % if keeping only specific scenarios for training -  live
    trainPeople_Match1 = [];
    for tt = trainPeopleInit
        TrNameList = NameList_Live(tt).name;
        if  strcmp(light_condition,'all_light')
           idx_match = trainPeopleInit; % match all the name strings containing photo
       elseif strcmp(light_condition,'adverse')
           idx_match = strfind(TrNameList, 'adverse');
       elseif strcmp(light_condition,'controlled')
           idx_match = strfind(TrNameList, 'controlled');
       end
       if isempty(idx_match)
           idx_match = 0;
       end
       trainPeople_Match1 = [trainPeople_Match1 idx_match];
    end

    trainPeople_idx_keep1 = find((trainPeople_Match1) ~=0); % choose non zero indices
    trainPeopleLiv = setdiff(trainPeople_idx_keep1, testPerson);
    
    % if keeping only specific scenarios for training -  attacks
    
    trainPeople_Match2 = [];
    for tt = trainPeopleInit
        TrNameList = NameList_Attack(tt).name;
        if strcmp(attack,'photo') & strcmp(light_condition,'all_light')
            idx_match = strfind(TrNameList, 'photo'); % match all the name strings containing photo
       elseif strcmp(attack,'video') & strcmp(light_condition,'all_light')
           idx_match = strfind(TrNameList, 'video');
       elseif strcmp(attack,'photo') & strcmp(light_condition,'adverse')
           idx_match = strfind(TrNameList, 'photo_adverse');
       elseif strcmp(attack,'video') & strcmp(light_condition,'adverse')
           idx_match = strfind(TrNameList, 'video_adverse');
       elseif strcmp(attack,'photo') & strcmp(light_condition,'controlled')
           idx_match = strfind(TrNameList, 'photo_controlled');
       elseif strcmp(attack,'video') & strcmp(light_condition,'controlled')
           idx_match = strfind(TrNameList, 'video_controlled');
        end
       if isempty(idx_match)
           idx_match = 0;
       end
       trainPeople_Match2 = [trainPeople_Match2 idx_match];
    end

    trainPeople_idx_keep2 = find((trainPeople_Match2) ~=0); % choose non zero indices
    trainPeopleAt = setdiff(trainPeople_idx_keep2, testPerson);
    
    trainPeople = [trainPeopleLiv trainPeopleAt];
end

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
%     Slive_p2 = S2(:,Str_idx);
%     Slive_p_tr = [Slive_p1 Slive_p2];
    Slive_p_tr = Slive_p1;
    
    Sfake_p_tr1 = S2(:,Str_idx);
    Sfake_p_tr2 = S3(:,Str_idx);
    
    Sfake_p_tr = [Sfake_p_tr1 Sfake_p_tr2];
    
    
% find training people indices     
Sts_pStart = (testPerson-1)*N+1;
Sts_pEnd = (testPerson)*N;
Sts_idx = [];
for ll = 1:length(Sts_pStart)
    Sts_i = Sts_pStart(ll):Sts_pEnd(ll);
    Sts_idx = [Sts_idx Sts_i];
end    

    Slive_p1_ts = S1(:,Sts_idx);% include each tr persons 16 ROIs
    Slive_p_ts = Slive_p1_ts;
%     Slive_p2_ts = S2(:,Sts_idx);
%     Slive_p_ts = [Slive_p1_ts Slive_p2_ts];

    Sfake_p_ts1 = S2(:,Sts_idx);
    Sfake_p_ts2 = S3(:,Sts_idx);
    
    Sfake_p_ts = [Sfake_p_ts1 Sfake_p_ts2];
    
    Mlist_live_p1_tr = Mlist1([trainPeople(1):trainPeople(end)], :);
%     Mlist_live_p2_tr = Mlist2([trainPeople(1):trainPeople(end)], :);
%     Mlist_live_p_tr = [Mlist_live_p1_tr; Mlist_live_p2_tr];
    
    Mlist_live_p_tr = Mlist_live_p1_tr;
    
    Mlist_fake_p2_tr = Mlist2([trainPeople(1):trainPeople(end)], :);
    Mlist_fake_p3_tr = Mlist3([trainPeople(1):trainPeople(end)], :);
    
    Mlist_fake_p_tr = [Mlist_fake_p2_tr; Mlist_fake_p3_tr];
    
    Mlist_live_p1_ts = Mlist1([testPerson(1):testPerson(end)], :);
%     Mlist_live_p2_ts = Mlist2([testPerson(1):testPerson(end)], :);
%     Mlist_live_p_ts = [Mlist_live_p1_ts; Mlist_live_p2_ts];
    
    Mlist_live_p_ts = Mlist_live_p1_ts;

    Mlist_fake_p2_ts = Mlist2([testPerson(1):testPerson(end)], :);    
    Mlist_fake_p3_ts = Mlist3([testPerson(1):testPerson(end)], :);
    
    Mlist_fake_p_ts = [Mlist_fake_p2_ts; Mlist_fake_p3_ts];
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

liveFolders = 1; 
fakeFolders = 2;%:3;  % if only 2 - fixed, 3 - handheld

if fakeFolders == 2
    motion_set_up = 'fixed';
elseif fakeFolders == 3
    motion_set_up = 'handheld';
end
    
labelsSVM = [];
predtests = [];
Ytss = [];
orderTrAll = [];
orderTsAll = []; 

% testPersonLiv = testPerson; % change if Replay
% testPersonAt = testPerson;
 [score_posterior, score, Yts, Ytr, labelSVM, predictionSVM, predictionSVMLive, ...
     predictionSVMFake] = SVM_Liu(saveLiuFolder, N, Q, liveFolders, ...
     fakeFolders, testPersonLiv, testPersonAt, trainPeopleLiv, trainPeopleAt, ...
     Slive_p_tr, Sfake_p_tr, Slive_p_ts, Sfake_p_ts, ...
     Mlist_live_p_tr, Mlist_fake_p_tr, Mlist_live_p_ts, Mlist_fake_p_ts);
     
     
   % append results from each LOOV run
scores_SVM_postcell{p} = score_posterior;
scores_SVMcell{p} = num2cell(score);
Ytsscell{p} = Yts;
Ytrscell{p} = Ytr;
testPeople = [testPeople testPerson];
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

save([saveLiuFolder '_' motion_set_up '_' train_cases '_' light_condition '_' attack '_' 'SVMresults_all_p.mat'], 'scores_SVM_postcell', 'scores_SVMcell', ...
    'Ytsscell', 'Ytrscell', 'testPeople', 'labelsSVMcell', 'predictionAllSVM', ...
    'predictionAllSVMLive', 'predictionAllSVMFake', 'predictionAverageSVM', ...
    'predictionAverageSVMLive', 'predictionAverageSVMFake')

        end % end t3
    end % end t2
end % end t1