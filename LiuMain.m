% Liu et al. 2016 ECCV implementation
%% load  videos

folderMain = '/media/ewa/SH/DatasetsAntiSpoof/3dmadDirectories/';

% fileNameList = dir([[folderMain folderEnd] [ '*.avi']]); 
% folderMain r= '/media/ewa/SH/DatasetsAntiSpoof/ReplayAttackDirectories/';
% folderReal = 'devel/real/';
% folderAttackf = 'devel/attack/fixed/';  % both contain photo, vid in adverse and controlled
% folderAttackh = 'devel/attack/hand/';

 for f = 1:3
     S = [];
     Mlist = [];
     
     folderEnd = ['Data0' num2str(f) 'Keep/'];

%     if f == 1
%         folderEnd = folderReal;
%     end
%     if f == 2
%         folderEnd = folderAttackf;
% 
%     end
%     if f == 3
%         folderEnd = folderAttackh;
%     end
    fileNameList = dir([[folderMain folderEnd] ['*' '.avi']]); %

    for i =1:length(fileNameList)
        imgCells{i} = fileNameList(i).name;  
    end
    
    [cs,index] = sort_nat(imgCells,'ascend');
    img_names = cs;
   for  m = 1:length(fileNameList)
    try
            vidName = img_names{m};       
            % read in the videos
            v = VideoReader([folderMain folderEnd vidName]);
            videoLength = v.Duration;
            videoRate = v.FrameRate;
            numFrame = videoLength*videoRate;
            width = v.Width;
            height = v.Height; 
            Fps = v.FrameRate;

            frames = read(v);

            vg = frames; %(:,:,2,:);
%             vg = permute(vg,[1,2,4,3]);   
t=1; % first frame
    firstFrame = vg(:,:,t);
    load(['/media/ewa/Data/PreliminaryResultsToCleanUp/FromBPADCameraVitals2016/Data/3DMAD1stFrame/Data0' num2str(f) 'Dlib/dLib-' vidName(1:end-4) '_C.avi' '.png.mat'])
%      load(['../Data/IDAP1stFrameNewDevel/' num2str(f) 'Dlib/dLib-' vidName '.png.mat'])

    firstPoints = pointsResized;
    
    
%% split into small regions
[facialRegions] = generateFacialROIs(firstPoints);   % each cell is each ROI

%% track
 whole = [firstPoints(1:68,1), firstPoints(1:68,2)];

    % whole face - lots of missing parts
    % between eyebrows and chin
    eyebrows1 = [firstPoints(18:22,1), firstPoints(18:22,2)];
    eyebrows2 = [firstPoints(23:27,1), firstPoints(23:27,2)];
    eyebrows = [eyebrows1; eyebrows2];
    
    chin = [firstPoints(1:17,1), firstPoints(1:17,2)];
    
        eyebrowsWidth = max(eyebrows(:,1)) - min(eyebrows(:,1));
    forehead1 = [eyebrows(:,1), (eyebrows(:,2)-round(eyebrowsWidth/2))];
    forehead2 = [eyebrows(:,1), (eyebrows(:,2))];

% plot([whole_face(:,1)], [whole_face(:,2)], 'r.');

whole_face1 = [forehead1(:,1); chin(end:-1:1,1)];
whole_face2 = [forehead1(:,2); chin(end:-1:1,2)];
whole_face = [(whole_face1(:)) whole_face2];

regionMask_whole_face = roipoly(firstFrame, whole_face(:,1), whole_face(:,2));

smallmask =  regionMask_whole_face; % regionMask_whole_face; %foreheadMasksubsampled;% faceROIsubsampled; %foreheadMask;
rr =  4;
pointsList = KLTtrackerMASK(vg,smallmask, rr);%%  rPPG detection using Haan and Jeanne

%% select facial ROIs
pulseXY = [];
N = length(facialRegions);
for n = 1:N
    %ismember(round(facialRegions{n}), round(pointsList(:,:,:)));
%     gridROI = pointsList(:,round(facialRegions{n}),:);
    tempReg = round(facialRegions{n});
    temp = inpolygon(pointsList(1,:,1), pointsList(1,:,2), tempReg(:,1), tempReg(:,2));
    tempIDX = find(temp ==1);
    gridROI = pointsList(:,tempIDX,:); % keep tracked points that are inside a given ROI
% raw PPG in R, G, B
    chR = 1;
    chG = 2;
    chB = 3;
    nn = 10;
    % CHROM PPG extraction
    [PPGrawR] = getPPGperChannel(chR, vg, gridROI, nn);
    [PPGrawG] = getPPGperChannel(chG, vg, gridROI, nn);
    [PPGrawB] = getPPGperChannel(chB, vg, gridROI, nn);

    % Chrom
    % raw PPG for a grid region 
    pulseXYii = getPulseXY(PPGrawR, PPGrawG, PPGrawB, Fps); % per grid ROI
    pulseXYi = mean(pulseXYii,2);  % average spatially for each small facial ROI
    pulseXY=  [pulseXY pulseXYi]; % per person
end

if  length(find(isnan(pulseXY))) > 0 
    continue
else 
    S = [S pulseXY];   % for all people together
    Mlist = [Mlist; f m]; % keep track of all the videos and from which folder, to split for LOOV tr and ts
end 
    catch 
        continue
    end
   end % end m
save(['Masks-LiuData-' num2str(f) '.mat'], 'S', 'Mlist')
 end % end f
   % iterate through videos, each time leaving out different ones for
   % validation 

%    TODO: Implement the same splitting scheme as in my LOOV_SVM code to
%    leave out all videos of the same person, instead of just one video at
%    a time but we can try one video at a time too

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
load('Masks-LiuData-1.mat')
S1 = S;
Mlist1 = Mlist;
load('Masks-LiuData-2.mat')
S2 = S;
Mlist2 = Mlist;
load('Masks-LiuData-3.mat')
S3 = S;
Mlist3 = Mlist;

% SliveAll = [S1 S2];

predictionAllSVM = [];
predictionAllSVMLive = [];
predictionAllSVMFake = [];
testPeople = [];
labelsSVM = [];
predtests = [];
Ytss = [];
orderTrAll = [];
orderTsAll = []; 
    
for p = 1:pEnd
    
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
    

    Slive_p1 = S1(:,[trainPeople(1)*N:trainPeople(end)*N]);  % include each tr persons 16 ROIs
    Slive_p2 = S2(:,[trainPeople(1)*N:trainPeople(end)*N]);
    Slive_p = [Slive_p1 Slive_p1];
    
    delta = 10^-3; % convergence threshold
        % r = 3; % bpm error toleration
    %     TODO: fix this part to determine indices defined as m list to leave out
%     k = 3; % keep top 3 eigenvectors, convert to 90 % variance instead
    alpha = 0.9; % in percent of variance
    N = 16;
%     warning o ff
    [pVec] =  iterate_p_e(Slive_p, delta, N, alpha); % only computed for live 
    % previously when keeping repeating pairs, q was 16 x 16, but it should
    % be 1 x C(N,2)+N and then diagonalized, then q 16 x 16 was vectorized
    % and diagonalized anyway to 16^2 x 16^2. Now, it'll be lower
    % dimensional
    q = get_q(pVec,N); % only computed for live 
    % Q =diag(q(:));   % Q = R'*R;
    % R = sqrtm(Q);
save(['Masks-LiuData-Live-p-' num2str(p) '.mat'], 'pVec', 'q')
  
   
%% SVM

saveLiuFolder = 'LiuMasks/';
mkdir(saveLiuFolder)
liveFolders = 1:2;
fakeFolders = 3;
N = 16;
% attack = 'photo';
 [scores_SVM_postcell, scores_SVMcell, Ytsscell, Ytrscell, testPeople, labelsSVMcell, ...
       orderTscell]  = SVM_LiuTogether(saveLiuFolder, N, q, liveFolders, fakeFolders, testPerson, trainPeople);
end % end p, for each LOOV person
