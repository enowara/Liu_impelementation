function  [score_posterior, score, Yts, Ytr, labelSVM, predictionSVM, predictionSVMLive, ...
     predictionSVMFake] = SVM_LiuTogether(saveLiuFolder, N, Q, liveFolders, ...
     fakeFolders, testPerson, testPersonLiv, testPersonAt, trainPeople, Slive_p_tr, ...
     Sfake_p_tr, Slive_p_ts, Sfake_p_ts, Mlist_live_p_tr, Mlist_fake_p_tr, ...
     Mlist_live_p_ts, Mlist_fake_p_ts)
    
    PdataLtr = [];
    R = sqrtm(Q); % computed once for each p, 
%                     q only depends on Slive which only changes with p
    for f = liveFolders
        for ff = 1:length(trainPeople)
            m = trainPeople(ff);
            try
                % get indices from Mlist
                mrow = [f m];
                idx = find(ismember(Mlist_live_p_tr,mrow,'rows') == 1);
                % load Sj from S
                Sj = Slive_p_tr(:,((idx-1)*N+1):(idx*N)); 
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
        for ff = 1:length(testPersonLiv)
            m = testPersonLiv(ff);
            try
            % get indices from Mlist
            mrow = [f m];
            idx = find(ismember(Mlist_live_p_ts,mrow,'rows') == 1);
            % load Sj from S
            Sj = Slive_p_ts(:,((idx-1)*N+1):(idx*N));
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
        % if replay
%                  if strcmp(attack,'photo')
%                     testPersonAt = testPersonInit([1,2,5,6,9,10]);
%                     trainPeople = setdiff(allPeople, testPersonInit);
%                 elseif strcmp(attack,'video')
%                     testPersonAt = testPersonInit([3,4,7,8]);
%                     trainPeople = setdiff(allPeople, testPersonInit);
%                  end            
        for ff = 1:length(trainPeople)
            m = trainPeople(ff);
            try
                % get indices from Mlist
                mrow = [f m];
                idx = find(ismember(Mlist_fake_p_tr,mrow,'rows') == 1);
                % load Sj from S
                Sj = Sfake_p_tr(:,((idx-1)*N+1):(idx*N));
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
%               if strcmp(attack,'photo')
%                     testPersonAt = testPersonInit([1,2,5,6,9,10]);
% %                     trainPeople = setdiff(allPeople, testPersonAt);
%                 elseif strcmp(attack,'video')
%                     testPersonAt = testPersonInit([3,4,7,8]);
% %                     trainPeople = setdiff(allPeople, testPersonAt);
%               end
              for ff = 1:length(testPersonAt)
                  m = testPersonAt(ff);
                  try
                    % get indices from Mlist
                    mrow = [f m];
                    idx = find(ismember(Mlist_fake_p_ts,mrow,'rows') == 1);
                    % load Sj from S
                    Sj = Sfake_p_ts(:,((idx-1)*N+1):(idx*N));
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
%             break % skip if there is data missing! otherwise - overfits
            disp('data missing')
        end
        YtrL = ones(size(PdataLtr,1), 1);
        YtrF = zeros(size(PdataFtr,1), 1);
        YtsL = ones(size(PdataLts,1), 1);
        YtsF = zeros(size(PdataFts,1), 1);

        % combine live and fake
        
        Ytr = [YtrL; YtrF];
        Yts = [YtsL; YtsF];

        Xtr = [PdataLtr; PdataFtr];
        Xts = [PdataLts; PdataFts];
        
        % No shuffling if LOOV
%         YtrTemp = [YtrL; YtrF];
%         YtsTemp = [YtsL; YtsF];
% 
%         XtrTemp = [PdataLtr; PdataFtr];
%         XtsTemp = [PdataLts; PdataFts];

%         XYtrTemp = [XtrTemp YtrTemp];
%         s = RandStream('mt19937ar','Seed',sum(100*clock));
%         orderTri = randperm(s, size(XYtrTemp,1));
%         XYtr = XYtrTemp(orderTri,:);

%         Xtr = XYtr(:,1:(end-1));
%         Ytr = XYtr(:,end);

%         XYtsTemp = [XtsTemp YtsTemp];
%         s = RandStream('mt19937ar','Seed',sum(100*clock));
%         orderTsi = randperm(s, size(XYtsTemp,1));
%         XYts = XYtsTemp(orderTsi,:);

%         Xtr = XYtr(:,1:(end-1));
%         Ytr = XYtr(:,end);

%         Xts = XYts(:,1:(end-1));
%         Yts = XYts(:,end);
        
        %% SVM
        SVMModel = fitcsvm(Xtr,Ytr,'KernelFunction','rbf','Standardize',true);
        [labelSVM,score] = predict(SVMModel,Xts);
        predictionSVM = (length(find(labelSVM==Yts))/length(Yts))*100;

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
        predictionSVMFake = (length(find(labelFake==YtsFake))/length(YtsFake))*100;
   end
        