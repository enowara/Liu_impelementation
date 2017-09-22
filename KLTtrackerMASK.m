function [pointsList] = KLTtrackerMASK(images, regionMask, rr)
% KLTtracker
% % debugging:
% images = vg;
% regionMask = faceROI;
ROISizeGrid =rr;
images = images(:,:,2,:);
images = permute(images,[1,2,4,3]); 

first_frame = images(:,:,1);
% first_frame = imresize(images(:,:,1), [240 320]);

numFrame = size(images, 3);


widthImg = size(first_frame, 1);
heightImg = size(first_frame, 2);

widthImg = size(first_frame, 1);
heightImg = size(first_frame, 2);
% set ROI
width = floor(widthImg/(ROISizeGrid/2));
height = floor(heightImg/(ROISizeGrid/2));

xmin = 1;
ymin =1;

% figure, imshow(regionMask)
% hold on
% initialize grid points, get each pixel within the ROI, counter should = w*h
counter = 0;
pointsList_init = [];
for i = 1:width-1
    for j = 1:height-1
        pixelX = xmin + (i-1)*ROISizeGrid/2 + floor(ROISizeGrid/2);  %+ i;  
        pixelY = ymin + (j-1)*ROISizeGrid/2 + floor(ROISizeGrid/2); %+ j;
        if regionMask(pixelX,pixelY)
%             regionMask(pixelX,pixelY)
           counter = counter+1;

%             plot(pixelY, pixelX, 'r *')
% [pixelX, pixelY]
%             pause(2)

            pointsList_init(counter, 1) = pixelY;
            pointsList_init(counter, 2) = pixelX;
        end
    end
end

ROIboxTr = [min(pointsList_init(:,1)), min(pointsList_init(:,2)), ... 
    (max(pointsList_init(:,1)) - min(pointsList_init(:,1))), ... 
    (max(pointsList_init(:,2)) - min(pointsList_init(:,2))) ];

% % check where the grid points are
% testImg = first_frame;
% testImg = insertMarker(testImg, [pointsList_init], '+', 'Color', 'blue');
% testImg = insertShape(testImg, 'Rectangle', ROIboxTr, 'Color', 'red');
% 
% figure, imshow(testImg)

% detect feature poitns
pointsStr = detectMinEigenFeatures(first_frame, 'MinQuality', 0.001, 'ROI', ROIboxTr);

% select the strongest K points 
numKLT = round(0.85*size(pointsStr,1));
pointsStr = pointsStr.selectStrongest(numKLT); 

% reject points outside the polygon 
in_sel = inpolygon(pointsStr.Location(:,1), pointsStr.Location(:,2),pointsList_init(:,1),pointsList_init(:,2)); 
points = pointsStr.Location(in_sel,:); 

% displayImage = first_frame;
% displayImage = insertShape(displayImage, 'Rectangle', ROIboxTr, 'Color', 'red');
% displayImage = insertMarker(displayImage, points, '+', 'Color', 'white');
% figure, imshow(displayImage);


%%

% create a tracker object
tracker = vision.PointTracker('MaxBidirectionalError', 1);
% initialize the tracker
initialize(tracker, points, first_frame);
% psize = size(pointsList_init);
% pointsList = zeros(psize(1), 2, numFrame);
pointsList = zeros(numFrame, size(pointsList_init, 1), 2);
pointsListpoints = pointsList_init;
oldPoints = points;
% track the points:
for i = 1:numFrame

%             if (rem(i,100)<1)
%                 disp(['percentage complete = ',num2str(i*100/numFrame)])
%             end
      next_frame = images(:,:,i);  % imresize(images(:,:,i), [240 320]);
      [points, isFound] = step(tracker, images(:,:,i)); % (imresize(images(:,:,i), [240 320])));
      visiblePoints = points(isFound, :);
      oldInliers = oldPoints(isFound, :);
      
      if size(visiblePoints, 1) >= 20 % need at least 2 points
        % Estimate the geometric transformation between the old points
        % and the new points and eliminate outliers
           [xform] = estimateGeometricTransform(...
           oldInliers, visiblePoints, 'similarity');%, 'MaxDistance', 2);
        % Apply the transformation to the bounding box points
        pointsListpoints = transformPointsForward(xform, pointsListpoints);
        xform.T;
        % Display tracked points
        next_frame = insertMarker(next_frame, visiblePoints, '+', ...
            'Color', 'white');
%         imshow(next_frame);
        
        pointsList(i,:,:) = pointsListpoints;
        % Reset the points
        oldPoints = visiblePoints;
        setPoints(tracker, oldPoints);
        i;
       
      else
          disp('not tracking')
      end
      
      if round(permute(pointsList(i,:,:), [2 3 1])) == round(pointsList_init(:,:))
%           disp('All same...')
      end
end
        % Clean up
        release(tracker);
        
        % check how the ROI was transformed
%         testROI = images(:,:,numFrame); %imresize(images(:,:,numFrame), [240 320]);%  % last frame
%         testROI = insertShape(testROI, 'Rectangle', ROIboxTr, 'Color', 'red');
%         testROI = insertMarker(testROI, pointsList_init, '+', 'Color', 'yellow');
%         testROI = insertMarker(testROI, permute(pointsList(numFrame, :,:), [2 3 1]), '+', 'Color', 'blue');

%        figure, imshow(testROI);

% end 
