    function [scores_SVM_postcell, scores_SVMcell, Ytsscell, Ytrscell, testPeople, labelsSVMcell, ...
        orderTscell]  = SVM_LiuTogether(saveLiuFolder, N, q, liveFolders, fakeFolders)
 % initialize 
    predictionAllSVM = [];
    predictionAllSVMLive = [];
    predictionAllSVMFake = [];
    testPeople = [];
    labelsSVM = [];
    predtests = [];
    Ytss = [];
    orderTrAll = [];
    orderTsAll = []; 
        
%     startTestPerson1 = (1:4:(60));
%     startTestPerson2 = (1:4:(60))+60;
%     startTestPerson3 = (1:2:(60))+120 ;
%     
%     endTestPerson1 =   ((1:4:(60)) + 3 );
%     endTestPerson2 =     (((1:4:(60))+60) +3);
%     endTestPerson3 =     (((1:2:(60))+120) +1); 
%     
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
   
   for p = 1:pEnd;  % p = 1:17 if 3DMAD, 15 per folder for Replay - test, train, devel - different people
    
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

%     trainPeople = setdiff(allPeople, testPerson);
    % get training live %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    xZero = 0;
    for ii = 1:3
               PdataLtr = [];
          for f = liveFolders
              testPersonLiv = testPersonInit;%(1:4);
              trainPeople = setdiff(allPeople, testPersonLiv);
              if ii == 1
                load(['Train-LiuData-' num2str(f) '.mat'])
              elseif ii == 2
                load(['Test-LiuData-' num2str(f) '.mat'])
              elseif ii == 3
               load(['Devel-LiuData-' num2str(f) '.mat'])
                 end     
                Q = diag(q(:));
                R = sqrtm(Q);
            for ff = 1:length(trainPeople)
            m = trainPeople(ff);
            try
                % get indices from Mlist
                mrow = [f m];
                idx = find(ismember(Mlist,mrow,'rows') == 1);
                % load Sj from S
                Sj = S(:,((idx-1)*N+1):(idx*N)); 
                xj = getRho(Sj, N); % without goodness metric, should be NxN or vector N^2
                %% SVM with RBF kernel
                xjkernel = xj*R;
                 if sum(isnan(xjkernel)) > 0
                    continue
%                 if sum(xjkernel) <= 0.00001
%                     xZero = xZero+1;
%                 end
                % make Xts and Xtr matrices 
                 else
                   PdataLtr = double([PdataLtr; xjkernel]);
                end
                if sum(sum(isnan(PdataLtr))) > 0
                    disp(['PdataLtr is NaN for' num2str(p) '-' num2str(ii)])
                
                 end
            catch
            continue
            end
            end
          end
%           sum(PdataLtr)
          PdataLtrCell{ii} = PdataLtr;
    end  
%         xZero  
          % get testing live %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii = 1:3
                PdataLts = [];
          for f = liveFolders
              testPersonLiv = testPersonInit;%(1:4);
%               trainPeople = setdiff(allPeople, testPersonLiv);
              if ii == 1
                load(['Train-LiuData-' num2str(f) '.mat'])
              elseif ii == 2
                load(['Test-LiuData-' num2str(f) '.mat'])
              elseif ii == 3
               load(['Devel-LiuData-' num2str(f) '.mat'])
                 end       
                Q = diag(q(:));
                R = sqrtm(Q);
            for ff = 1:length(testPersonLiv)
            m = testPersonLiv(ff);
             try
                % get indices from Mlist
                mrow = [f m];
                idx = find(ismember(Mlist,mrow,'rows') == 1);
                % load Sj from S
                Sj = S(:,((idx-1)*N+1):(idx*N));
                xj = getRho(Sj, N); % without goodness metric, should be NxN or vector N^2
                %% SVM with RBF kernel
                xjkernel = xj*R;
                if sum(isnan(xjkernel)) > 0
                    continue
                % make Xts and Xtr matrices 
                else
                    PdataLts = [PdataLts; xjkernel];
                end
                if sum(sum(isnan(PdataLts))) > 0
                    disp(['PdataLts is NaN for' num2str(p) '-' num2str(ii)])
                 end
            catch
            continue
            end
            end
          end
%           sum(PdataLts)
       PdataLtsCell{ii} = PdataLts;
    end
            % get training fake %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       for ii = 1:3
             
               PdataFtr = [];
          for f = fakeFolders
%                  if strcmp(attack,'photo')
                    testPersonAt = testPersonInit;%([1,2,5,6,9,10]);
                    trainPeople = setdiff(allPeople, testPersonInit);
%                 elseif strcmp(attack,'video')
%                     testPersonAt = testPersonInit([3,4,7,8]);
%                     trainPeople = setdiff(allPeople, testPersonInit);
%                  end
                if ii == 1
                load(['Train-LiuData-' num2str(f) '.mat'])
              elseif ii == 2
                load(['Test-LiuData-' num2str(f) '.mat'])
              elseif ii == 3
               load(['Devel-LiuData-' num2str(f) '.mat'])
                 end     
              
                Q = diag(q(:));
                R = sqrtm(Q);
            for ff = 1:length(trainPeople)
            m = trainPeople(ff);
            try
                % get indices from Mlist
                mrow = [f m];
                idx = find(ismember(Mlist,mrow,'rows') == 1);
                % load Sj from S
                Sj = S(:,((idx-1)*N+1):(idx*N));
                xj = getRho(Sj, N); % without goodness metric, should be NxN or vector N^2
                %% SVM with RBF kernel
                xjkernel = xj*R;
                if sum(isnan(xjkernel)) > 0
                    continue
                % make Xts and Xtr matrices 
                else
                    PdataFtr = [PdataFtr; xjkernel];
                end
                if sum(sum(isnan(PdataFtr))) > 0
                    disp(['PdataFtr is NaN for' num2str(p) '-' num2str(ii)])
                 end
            catch
            continue
            end
            end
          end
%           sum(PdataFtr)
        PdataFtrCell{ii} = PdataFtr;
          end
       
% get testing live %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
          

 for ii = 1:3
    PdataFts = [];
          for f = fakeFolders
%               if strcmp(attack,'photo')
                    testPersonAt = testPersonInit;%([1,2,5,6,9,10]);
%                     trainPeople = setdiff(allPeople, testPersonAt);
%                 elseif strcmp(attack,'video')
%                     testPersonAt = testPersonInit([3,4,7,8]);
%                     trainPeople = setdiff(allPeople, testPersonAt);
%               end
                 if ii == 1
                load(['Train-LiuData-' num2str(f) '.mat'])
              elseif ii == 2
                load(['Test-LiuData-' num2str(f) '.mat'])
              elseif ii == 3
               load(['Devel-LiuData-' num2str(f) '.mat'])
                 end     
                              Q = diag(q(:));
                R = sqrtm(Q);
            for ff = 1:length(testPersonAt)
            m = testPersonAt(ff);
            try
                % get indices from Mlist
                mrow = [f m];
                idx = find(ismember(Mlist,mrow,'rows') == 1);
                % load Sj from S
                Sj = S(:,((idx-1)*N+1):(idx*N));
                xj = getRho(Sj, N); % without goodness metric, should be NxN or vector N^2
                %% SVM with RBF kernel
                xjkernel = xj*R;
                if sum(isnan(xjkernel)) > 0
                    continue
                % make Xts and Xtr matrices 
                else
                    PdataFts = [PdataFts; xjkernel];
                end
                  if sum(sum(isnan(PdataFts))) > 0
                    disp(['PdataFts is NaN for' num2str(p) '-' num2str(ii)])
                 end
              catch
            continue
            end
          end
          end
%           sum(PdataFts)
          PdataFtsCell{ii} = PdataFts;
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  if not empty

%% set labels
    PdataLtrTotalTemp = [PdataLtrCell{1} ; PdataLtrCell{2} ; PdataLtrCell{3}];
    PdataFtrTotalTemp = [PdataFtrCell{1} ; PdataFtrCell{2} ; PdataFtrCell{3}];
%     size(PdataLtrTotalTemp)
%     size(PdataFtrTotalTemp)
%     size(PdataLtrCell{1})
%     size(PdataFtrCell{1})
        scores_SVM_postI = [];
        scores_SVMI = [];
        YtssI = [];
        YtrsI = [];
        labelsSVMI = [];
        orderTsI = [];
         
        for ii = 1:3
            if ii == 1
                PdataLtrTotal = [PdataLtrTotalTemp; PdataLtsCell{2}; PdataLtsCell{3}];
                PdataFtrTotal = [PdataFtrTotalTemp; PdataFtsCell{2}; PdataFtsCell{3}];
                PdataLtsTotal = [PdataLtsCell{1}];
                PdataFtsTotal = [PdataFtsCell{1}];
                
            elseif ii == 2                
                PdataLtrTotal = [PdataLtrTotalTemp; PdataLtsCell{1}; PdataLtsCell{3}];
                PdataFtrTotal = [PdataFtrTotalTemp; PdataFtsCell{1}; PdataFtsCell{3}];
                PdataLtsTotal = [PdataLtsCell{2}];
                PdataFtsTotal = [PdataFtsCell{2}];
            elseif ii == 3      
                PdataLtrTotal = [PdataLtrTotalTemp; PdataLtsCell{1}; PdataLtsCell{2}];
                PdataFtrTotal = [PdataFtrTotalTemp; PdataFtsCell{1}; PdataFtsCell{2}];
                PdataLtsTotal = [PdataLtsCell{3}];
                PdataFtsTotal = [PdataFtsCell{3}];
            end
            
        YtrL = ones(size(PdataLtrTotal,1), 1);
        YtrF = zeros(size(PdataFtrTotal,1), 1);
        YtsL = ones(size(PdataLtsTotal,1), 1);
        YtsF = zeros(size(PdataFtsTotal,1), 1);
        
        XtrTemp = [PdataLtrTotal; PdataFtrTotal];
        XtsTemp = [PdataLtsTotal; PdataFtsTotal];

        % combine. No shuffling if LOOV?
        YtrTemp = [YtrL; YtrF];
        YtsTemp = [YtsL; YtsF];

        XYtrTemp = [XtrTemp YtrTemp];
        s = RandStream('mt19937ar','Seed',sum(100*clock));
        orderTri = randperm(s, size(XYtrTemp,1));
        XYtr = XYtrTemp(orderTri,:);

        Xtr = XYtr(:,1:(end-1));  % is 582, should be 1200?
        Ytr = XYtr(:,end);

        XYtsTemp = [XtsTemp YtsTemp];
        s = RandStream('mt19937ar','Seed',sum(100*clock));
        orderTsi = randperm(s, size(XYtsTemp,1));
        XYts = XYtsTemp(orderTsi,:);

        Xtr = XYtr(:,1:(end-1));
        Ytr = XYtr(:,end);      % should be longer than 6?

        Xts = XYts(:,1:(end-1));
        Yts = XYts(:,end);
    %% SVM
        SVMModel = fitcsvm(Xtr,Ytr,'KernelFunction','linear','Standardize',true);
        [labelSVM,score] = predict(SVMModel,Xts);
        predictionSVM = (length(find(labelSVM==Yts))/length(Yts))*100;
        predictionAllSVM = [predictionAllSVM; predictionSVM];    

        % prediction for live and fake separately 
        LiveIdx = find(Yts == 1);
        FakeIdx = find(Yts == 0);
        labelLive = labelSVM(LiveIdx); %label(1:end/2);
        labelFake = labelSVM(FakeIdx); %label((end/2+1):end);
        YtsLive = Yts(LiveIdx); %Yts(1:end/2);
        YtsFake = Yts(FakeIdx); %Yts((end/2+1):end);
        
        SVMModel2 = fitPosterior(SVMModel);
        [~,score_posterior] = resubPredict(SVMModel2);
        if length(YtsLive) == 0 || length(YtsFake) ==0
%             continue 
        break
        end

        predictionSVMLive = (length(find(labelLive==YtsLive))/length(YtsLive))*100;
        predictionAllSVMLive = [predictionAllSVMLive; predictionSVMLive]; 
        predictionSVMFake = (length(find(labelFake==YtsFake))/length(YtsFake))*100;
        predictionAllSVMFake = [predictionAllSVMFake; predictionSVMFake];

        scores_SVM_postI = [scores_SVM_postI; score_posterior];
        scores_SVMI = [scores_SVMI; score];
        YtssI = [YtssI; Yts];
        YtrsI = [YtrsI; Ytr];
%         testPeopleI = [testPeople; testPersonInit];
        labelsSVMI = [labelsSVMI; labelSVM];
%         predtestsI = [predtestsI; predtest];
        orderTsI = [orderTsI; orderTsi'];
        
  end
        scores_SVM_postcell{p} = scores_SVM_postI;
        scores_SVMcell{p} = num2cell(scores_SVMI);
%         scores_RDFcell{p} = scores_RDFI;
        Ytsscell{p} = YtssI;
        Ytrscell{p} = YtrsI;
        testPeople = [testPeople; testPersonInit];
        labelsSVMcell{p} = labelsSVMI;
%         predtestscell{p} = predtestsI;
        orderTscell{p} = orderTsI;
   end
    
      % sum up accuracy over all test people    
        predictionAverageSVM = sum(predictionAllSVM)/length(predictionAllSVM);
        disp([num2str(predictionAverageSVM) '% Average SVM accuracy']);

        predictionAverageSVMLive = sum(predictionAllSVMLive)/length(predictionAllSVMLive);
        disp([num2str(predictionAverageSVMLive) '% Average Live SVM accuracy']);
        predictionAverageSVMFake = sum(predictionAllSVMFake)/length(predictionAllSVMFake);
        disp([num2str(predictionAverageSVMFake) '% Average Fake SVM accuracy']);
 save([saveLiuFolder 'LiuMasks.mat'])

        
   end
          