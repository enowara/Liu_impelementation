function [vidSin_out1, vidCol_out1, Fps1] = read_video(filename1, fullPath2File1, vid_phot1, fileFormat1, vid_start1, vid_end1)

% given filename and file location, read in the video and keep only the
% channel of interest

% input:
% filename1 - name of the video to be read
% folderMain1 - main folder where all videos are stored
% fullPath2File1 - specific folder where this video of interest is stored
% vid_start1 - frame from which the video should start being read
% vid_end1 - frame until which the video should end being read - subtract
% it from the total video length

% output:
% vidCol_out1 - original (RGB or IR) video - some applications may require
% RGB
% vidSin_out1 - signle channel read video, green or IR 

% if video has three channels - RGB, select green channel

% elseif video has only one channel - IR or G, just keep that channel


%% debugging parameters
% filename1 =  filename;
% fullPath2File1 =  fullPath2File;
% vid_phot1 =  vid_phot;
% fileFormat1 =  fileFormat;
% vid_start1 =  vid_start;
% vid_end1 = vid_end;
% Fps1 =  Fps;

addpath('/home/ewa/Dropbox (Rice Scalable Health)/DocumentsUbuntu/Liveness_Detection_Security/Code/temp_helper_functions');
%% if the input is a video
if strcmp(vid_phot1, 'vid')
%     disp('videoMode')
    if nargin == 4 % if video start time is not provided, start at frame 1
        vid_start1 = 0; % 1 - first second? 1/30 - first frame
    end

%     if strcmp(filename, '') ~= 1  % if filename is not empty, that is the case for images
        
    v = VideoReader([fullPath2File1 filename1], 'CurrentTime', vid_start1);
    Fps1 = v.FrameRate;
    video_Length = Fps1*v.Duration;  % the whole video duration is frame rate times number of seconds
    Height = v.Height;
    Width = v.Width; 
    
    if nargin == 4 || 5 % if video end time is not provided, or neither start and end times, set the end time to be the last frame
        
%         if ndims(video) == 4 
%             vid_end1 = size(video,4);
%         elseif ndims(video) == 3
%             vid_end1 = size(video,3);
%         end

    vid_end1 = video_Length;
            
    end
    
%     if vid_end1 == 0   % 0 can be an argument passed to make the video full length
%         vid_end1 = video_Length;    
%     end
    
        
    counter = 0;
% initialize video object
video2Save_init = zeros(Height, Width, 3, (vid_start1:floor(video_Length - vid_end1)));
    while counter <  size(video2Save_init,3) %hasFrame(v)
     counter = counter+1;
%     for counter = vid_start1:v.Duration-1
        temp_video = readFrame(v);
        video2Save(:,:,:,counter) = temp_video;  % otherwise the temp_video gets overwritten and ends up being a single frame
    end

% if counter < floor(vid_end1) % it means that there are less frames that were read than what it says the duration of the video is
%     video2Save = video2Save_init(:,:,:,1:counter);
% end
    % if video is 4-D - w,h, rgb, time, select only the green channel
    if ndims(video2Save) == 4 % if number of dimensions is 4, it means the video is RGB
        videoGreen = permute(video2Save(:,:,2,:), [1 2 4 3]);
    else
        videoGreen = video2Save; % if IR, just one channel, keep
    end

    vidCol_out1 = video2Save;
    vidSin_out1 = videoGreen;

%% if the input is a stream of photos
elseif strcmp(vid_phot1, 'phot')
%         disp('photoMode')

    % first check if the images were already read in and concatenated into
    % an array
% if exist([folderMain1 folderVid1 'images.mat'])
%     ('exists')
%     load([folderMain1 folderVid1 'images.mat'])
% else % if not then follow these steps to read in each frame and concatenate them together
        
    if nargin <= 4 % if video start time is not provided, start at frame 1
        vid_start1 = 1; % 1 - first frame
    end

    % filename = strcat(filePrefix,num2str(vid_start1-1,fileFormat),filetype);
    % filepath = fullfile(dataDir,filename);

    % read in all images in the directory of the given format, ex. .png
    fileList_init = dir([[[fullPath2File1]] ['*' fileFormat1]]);
    % sort the image names to make sure they are in the right order
    for i =1:length(fileList_init)
        imgCells{i} = fileList_init(i).name;  
    end
     
    [cs,index] = sort_nat(imgCells,'ascend');
    fileList = cs;

    % read in the first frame (1 or vid_start1) 
    tempImg = imread([fullPath2File1 fileList{vid_start1}]);    % first image
    Height = size(tempImg,1);
    Width = size(tempImg,2);
    
    
    if nargin <= 5 % if video end time is not provided, or neither start and end times, set the end time to be the last frame
        vid_end1 = length(fileList);
    end
    
    numFrame = length(vid_start1:length(fileList_init) - vid_end1); % figure out how many frames are to be read = length of fileList.names?
    
    % initialize an image array with the size of the first image
        
    if isa(tempImg,'uint8')
        rawImage = zeros(Height,Width, 3, numFrame,'uint8');
        rawImageG = zeros(Height,Width, numFrame,'uint8');
    elseif isa(tempImg,'uint8')
        rawImage = zeros(Height,Width, 3, numFrame,'uint16');
        rawImageG = zeros(Height,Width, numFrame,'uint8');
%     else
%         ?
    end

    %rawImage = zeros(Height,Width, numFrame);
    counter = 0;
    for t = vid_start1:(vid_start1+numFrame)
        counter = counter+1;
        Img = imread([fullPath2File1 fileList{t}]);    % first image

        if ndims(Img) == 3
            img_green = Img(:,:,2); % Green channel 
            img_RGB = Img;
            % append each image to the initialized image array
            rawImageG(:,:,counter) = img_green; 
            rawImage(:,:,:,counter) = img_RGB; 
        else
            img_green = Img;
            img_RGB = Img;
            % append each image to the initialized image array
            rawImageG(:,:,counter) = img_green; 
            rawImage(:,:,counter) = img_RGB;
        end

        
        
    end 

    % save
%     save([folderMain1 folderVid1 'images.mat'], 'rawImage');
% end % if file does not exist

vidSin_out1 = rawImageG;
vidCol_out1 = rawImage; % for now keep color image output to be just the green channel but later it can be changed to concatenate raw RGB imgs 

Fps1 = 30;
disp('Fps is automatically set to 30 by default')
% TODO - save color image stream
% TODO - check Fps for images, usually 30 
end
end
