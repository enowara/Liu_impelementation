    function [scores_SVM_postcell, scores_SVMcell, Ytsscell, Ytrscell, testPeople, labelsSVMcell, ...
        orderTscell]  = SVM_Liu(saveLiuFolder, N, liveFolders, fakeFolders, attack)
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
        
    startTestPerson1 = (1:4:(60));
    startTestPerson2 = (1:4:(60))+60;
    startTestPerson3 = (1:2:(60))+120 ;
    
    endTestPerson1 =   ((1:4:(60)) + 3 );
    endTestPerson2 =     (((1:4:(60))+60) +3);
    endTestPerson3 =     (((1:2:(60))+120) +1); 
    
    allPeople = startTestPerson1(1):endTestPerson3(end);
    pEnd = 15;
   
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
   
              
    PdataLtr = [];
          for f = liveFolders
              testPersonLiv = testPersonInit(1:4);
              trainPeople = setdiff(allPeople, testPersonLiv);
              load(['Train-LiuData-' num2str(f) '.mat'])
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
                % make Xts and Xtr matrices 
                PdataLtr = [PdataLtr; xjkernel];
            catch
            continue
            end
            end
          end
            
          
          % get testing live %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
    PdataLts = [];
          for f = liveFolders
              testPersonLiv = testPersonInit(1:4);
              load(['Train-LiuData-' num2str(f) '.mat'])
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
                % make Xts and Xtr matrices 
                PdataLts = [PdataLts; xjkernel];
            catch
            continue
            end
            end
          end
          
              % get training fake %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
               PdataFtr = [];
          for f = fakeFolders
                 if strcmp(attack,'photo')
                    testPersonAt = testPersonInit([1,2,5,6,9,10]);
                    trainPeople = setdiff(allPeople, testPersonInit);
                elseif strcmp(attack,'video')
                    testPersonAt = testPersonInit([3,4,7,8]);
                    trainPeople = setdiff(allPeople, testPersonInit);
                 end
         
              load(['Train-LiuData-' num2str(f) '.mat'])
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
                % make Xts and Xtr matrices 
                PdataFtr = [PdataFtr; xjkernel];
            catch
            continue
            end
            end
          end
          
          % get testing live %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
          
    PdataFts = [];
          for f = fakeFolders
              if strcmp(attack,'photo')
                    testPersonAt = testPersonInit([1,2,5,6,9,10]);
%                     trainPeople = setdiff(allPeople, testPersonAt);
                elseif strcmp(attack,'video')
                    testPersonAt = testPersonInit([3,4,7,8]);
%                     trainPeople = setdiff(allPeople, testPersonAt);
                 end
              load(['Train-LiuData-' num2str(f) '.mat'])
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
                % make Xts and Xtr matrices 
                PdataFts = [PdataFts; xjkernel];
              catch
            continue
            end
          end
          end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  if not empty

%% set labels
        if isempty(PdataLtr) || isempty(PdataFtr) || isempty(PdataLts) || isempty(PdataFts)
            break % skip if there is data missing! otherwise - overfits
        end
        YtrL = ones(size(PdataLtr,1), 1);
        YtrF = zeros(size(PdataFtr,1), 1);
        YtsL = ones(size(PdataLts,1), 1);
        YtsF = zeros(size(PdataFts,1), 1);

        % combine. No shuffling if LOOV?
        YtrTemp = [YtrL; YtrF];
        YtsTemp = [YtsL; YtsF];

        XtrTemp = [PdataLtr; PdataFtr];
        XtsTemp = [PdataLts; PdataFts];

        XYtrTemp = [XtrTemp YtrTemp];
        s = RandStream('mt19937ar','Seed',sum(100*clock));
        orderTri = randperm(s, size(XYtrTemp,1));
        XYtr = XYtrTemp(orderTri,:);

        Xtr = XYtr(:,1:(end-1));
        Ytr = XYtr(:,end);

        XYtsTemp = [XtsTemp YtsTemp];
        s = RandStream('mt19937ar','Seed',sum(100*clock));
        orderTsi = randperm(s, size(XYtsTemp,1));
        XYts = XYtsTemp(orderTsi,:);

        Xtr = XYtr(:,1:(end-1));
        Ytr = XYtr(:,end);

        Xts = XYts(:,1:(end-1));
        Yts = XYts(:,end);
        
        %% SVM
        SVMModel = fitcsvm(Xtr,Ytr,'KernelFunction','rbf','Standardize',true);
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

        predictionSVMLive = (length(find(labelLive==YtsLive))/length(YtsLive))*100;
        predictionAllSVMLive = [predictionAllSVMLive; predictionSVMLive]; 
        predictionSVMFake = (length(find(labelFake==YtsFake))/length(YtsFake))*100;
        predictionAllSVMFake = [predictionAllSVMFake; predictionSVMFake]; 
        
        scores_SVM_postcell{p} = score_posterior;
        scores_SVMcell{p} = num2cell(score);
        Ytsscell{p} = Yts;
        Ytrscell{p} = Ytr;
%         testPeople = [testPeople; testPerson];
        testPeople = [testPeople; testPersonInit];

        labelsSVMcell{p} = labelSVM;
        orderTscell{p} = orderTsi;
   end
        % sum up accuracy over all test people    
        predictionAverageSVM = sum(predictionAllSVM)/length(predictionAllSVM);
        disp([num2str(predictionAverageSVM) '% Average SVM accuracy']);

        predictionAverageSVMLive = sum(predictionAllSVMLive)/length(predictionAllSVMLive);
        disp([num2str(predictionAverageSVMLive) '% Average Live SVM accuracy']);
        predictionAverageSVMFake = sum(predictionAllSVMFake)/length(predictionAllSVMFake);
        disp([num2str(predictionAverageSVMFake) '% Average Fake SVM accuracy']);
save([saveLiuFolder attack 'LiuReplay.mat'])
