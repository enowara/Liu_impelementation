% Liu et al. 2016 ECCV implementation
%% load  videos

filename =''; %'01_01_01.avi';
vid_phot = 'phot'; %'vid';
fileFormat = '.png';% '';


mainFolderVids = '/media/ewa/SH/DatasetsAntiSpoof/EwaPPGdataCollection/NewerHandvsFixedNotxt/';
mainFolderPnts = '/media/ewa/SH/DatasetsAntiSpoof/EwaPPGdataCollection/NewerHandvsFixed/';

datte = '10_10';
saveLiuFolder = ['/home/ewa/Dropbox (Rice Scalable Health)/DocumentsUbuntu/' ...
        'Liveness_Detection_Security/Results/BPADCameraVitals2016/LiuEwaData/' datte '/'];
mkdir(saveLiuFolder)

vidList = dir(mainFolderVids);

datte = '10-11';
saveLiuFolder = ['LiuEwaData/' datte '/'];
mkdir(saveLiuFolder)

S = [];
Mlist = [];

for vidNum = 3%:length(vidList)% [4:6, 8:length(vidList)]
    
    vidName = vidList(vidNum).name; % EwaLiveCamFixed/cam1/';
    fullPath2File = [mainFolderVids vidName '/cam1/'];
    pointsFile = [mainFolderPnts vidName '.txt'];


%     addpath('/home/ewa/Dropbox (Rice Scalable Health)/DocumentsUbuntu/Liveness_Detection_Security/Code/temp_helper_functions/altmany-export_fig-5be2ca4/');
%     addpath('/home/ewa/Dropbox (Rice Scalable Health)/DocumentsUbuntu/Liveness_Detection_Security/Code/BPAD2017_Fall/')

    %% split videos into shorter time windows
    
    % first check how long each video is
    vid_length = length(dir(fullPath2File));
    tim_win = 10*30;
    counter = 0;
    
    S_init = [];
    Mlist_init = [];
    
    for tt = 1%:tim_win:vid_length-tim_win
        counter = counter+1;
        vid_start = tt;
        vid_finish = tt+tim_win - 1;%[]; % remove this many frames from the end of the video
        vid_end = vid_length - tim_win - vid_start - 1;
        %% read video 
        [~, vidSin_out, Fps] = read_video(filename, fullPath2File, vid_phot, fileFormat, vid_start, vid_end);
cd '/home/ewa/Dropbox (Rice Scalable Health)/DocumentsUbuntu/Liveness_Detection_Security/Code/BPAD2017_Fall'
        % %% load landmarks
        img_names = 'frame_00000.png';

        [pointsResizedAll] = dlib2FacialPoints(fullPath2File,img_names, pointsFile);
        
        firstPoints = pointsResizedAll;
        firstFrame = vidSin_out(:,:,:,1);
cd '/home/ewa/Dropbox (Rice Scalable Health)/DocumentsUbuntu/Liveness_Detection_Security/Code/BPADCameraVitals2016/Liu_Implementation'   
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
        pointsList = KLTtrackerMASK(vidSin_out,smallmask, rr);%%  rPPG detection using Haan and Jeanne

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
            [PPGrawR] = getPPGperChannel(chR, vidSin_out, gridROI, nn);
            [PPGrawG] = getPPGperChannel(chG, vidSin_out, gridROI, nn);
            [PPGrawB] = getPPGperChannel(chB, vidSin_out, gridROI, nn);

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

    disp([num2str(f) '-' num2str(counter)])
   end % end counter
    S = [S S_init]; % are the videos going to be the same length to put them together here
    Mlist = [Mlist; Mlist_init];
    
save([saveLiuFolder 'Replay-LiuData-' num2str(f) '.mat'], 'S', 'Mlist')
 end % end f
 
   % iterate through videos, each time leaving out different ones for
   % validation 

%    TODO: Implement the same splitting scheme as in my LOOV_SVM code to
%    leave out all videos of the same person, instead of just one video at
%    a time but we can try one video at a time too

