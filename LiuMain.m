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
end
%% get e and p iteratively
% load('Test-LiuData-1.mat')
% Sts = S;
% load('Train-LiuData-1.mat')
% Str = S;
% load('Devel-LiuData-1.mat')
% Sdev = S;
% 
% Slive = [Str Sts Sdev];

if f == 1 || f ==2 % live videos
    Slive = S;
    delta = 10^-3; % convergence threshold
    % alpha = 0.6; %(% of variance preserved)
    % r = 3; % bpm error toleration
    J = size(Slive,2); %1;% size(Slive,2)/3;
    
    % implement the LOOV strategy. Loop over all people in the dataset
    % except one when training? 
    
    k = 3; % keep top 3 eigenvectors, convert to 90 % variance instead
    alpha = 0.9; % in percent of variance
    N = 16;
    warning off
    [pVec] =  iterate_p_e(Slive, delta, k, J, N, alpha); % only computed for live 
    q = get_q(pVec,N); % only computed for live 
    % Q =diag(q(:));   % Q = R'*R;
    % R = sqrtm(Q);
    save(['Masks-LiuData-' num2str(f) '.mat'], 'S', 'Mlist', 'pVec', 'q')
else
    save(['Masks-LiuData-' num2str(f) '.mat'], 'S', 'Mlist', 'pVec', 'q')
end



end
%% SVM

saveLiuFolder = 'LiuMasks/';
mkdir(saveLiuFolder)
liveFolders = 1:2;
fakeFolders = 3;
N = 16;
% attack = 'photo';
 [scores_SVM_postcell, scores_SVMcell, Ytsscell, Ytrscell, testPeople, labelsSVMcell, ...
       orderTscell]  = SVM_LiuTogether(saveLiuFolder, N, q, liveFolders, fakeFolders);
   
