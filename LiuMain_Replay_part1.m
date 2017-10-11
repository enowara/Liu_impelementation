% Liu et al. 2016 ECCV implementation
%% load  videos

folderMain = '/media/ewa/SH/DatasetsAntiSpoof/ReplayAttackDirectories/';
folderReal_Alll = {'train/real/', 'test/real/', 'devel/real/'};
folderAttackf_All = {'train/attack/fixed/', 'test/attack/fixed/', 'devel/attack/fixed/'};  
folderAttackh_All = {'train/attack/hand/', 'test/attack/hand/', 'devel/attack/hand/'};

DlibFolder = {'IDAP1stFrameNew', 'IDAP1stFrameNewTest', 'IDAP1stFrameNewDevel'};

datte = '10-03';
saveLiuFolder = ['LiuReplay/' datte '/'];
mkdir(saveLiuFolder)

for f = 1:3
    NameListAll = [];
    S = [];
    Mlist = []; % combine S and Mlist matrices for each train, test, devel folder
% combine from train, test, devel and split by types of attacks 
    for FF = 1:3
        folderReal = folderReal_Alll{FF};
        folderAttackf = folderAttackf_All{FF};
        folderAttackh = folderAttackh_All{FF};

    %  for f = 1:3
         S_init = [];
         Mlist_init = [];

        if f == 1
            folderEnd = folderReal;
        end
        if f == 2
            folderEnd = folderAttackf;

        end
        if f == 3
            folderEnd = folderAttackh;
        end
        fileNameList = dir([[folderMain folderEnd] ['*' '.mov']]); %

        clear imgCells
        
        for i =1:length(fileNameList)
            imgCells{i} = fileNameList(i).name;  
        end
    
        [cs,index] = sort_nat(imgCells,'ascend');
        img_names = cs;
        
        % save all video names to later choose tr and ts sets with specific
        % set ups
        if f == 1
            NameListAll = [NameListAll; fileNameList];
            
        elseif f == 2 
            NameListAll = [NameListAll; fileNameList];
            
        elseif f == 3
            NameListAll = [NameListAll; fileNameList];
            
        end
            
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
        load(['/media/ewa/Data/PreliminaryResultsToCleanUp/FromBPADCameraVitals2016/Data/' DlibFolder{FF} '/' num2str(f) 'Dlib/dLib-' vidName '.png.mat'])

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
        N = size(facialRegions,2);
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
            
            % only keep the length corresponding to shortest video in dataet
            PPGrawR = PPGrawR(1:229,:);
            PPGrawG = PPGrawG(1:229,:);
            PPGrawB = PPGrawB(1:229,:);

            % Chrom
            % raw PPG for a grid region 
            pulseXYii = getPulseXY(PPGrawR, PPGrawG, PPGrawB, Fps); % per grid ROI
            pulseXYi = mean(pulseXYii,2);  % average spatially for each small facial ROI
            pulseXY=  [pulseXY pulseXYi]; % per person
        end

        if  length(find(isnan(pulseXY))) > 0 
            continue
        else 
            
            S_init = [S_init pulseXY];   % for all people together
            Mlist_init = [Mlist_init; f m]; % keep track of all the videos and from which folder, to split for LOOV tr and ts
        end 
    catch 
        continue
        end
    disp([num2str(f) '-' num2str(FF) '-' num2str(m)])
   end % end m
    S = [S S_init]; % are the videos going to be the same length to put them together here
    Mlist = [Mlist; Mlist_init];
    end % end FF
    
save([saveLiuFolder 'Replay-LiuData-' num2str(f) '.mat'], 'S', 'Mlist', 'NameListAll')
 end % end f
 
   % iterate through videos, each time leaving out different ones for
   % validation 

%    TODO: Implement the same splitting scheme as in my LOOV_SVM code to
%    leave out all videos of the same person, instead of just one video at
%    a time but we can try one video at a time too

