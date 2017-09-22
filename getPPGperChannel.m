function [PPGraw] = getPPGperChannel(ch, filteredImages, gridPoints, n)
% get a raw estimate of the PPG in R, , B channels

filteredImages = filteredImages(:,:,ch,:);     % choose channel ch - rgb
filteredImages = permute(filteredImages,[1,2,4,3]);

numFrame = size(filteredImages,3);
imageHeight = size(filteredImages,1);
imageWidth = size(filteredImages,2);
% initialize PPG
PPG = zeros(floor(numFrame), size(gridPoints,2));

%     PPGraw = getPPGperChannel(ch)
%     function getPPGperChannel(ch)

for k = 1:floor(numFrame)
      rawPerfussion = (filteredImages(:,:,k));             
         for c = 1:size(gridPoints,2) 
            % Get the tracked gridPoints 
            cxStart = round(gridPoints(k, c,1));
            cyStart = round(gridPoints(k, c,2)); 
            if (cyStart-n)<1
                yStart = 1;
            else
                yStart = (cyStart-n);
            end 

            if (cxStart-n)<1
                xStart = 1;
            else
                xStart = (cxStart-n);
            end 

            if (cyStart+(n-1)) > imageHeight
                yEnd = imageHeight;
            else
                yEnd = (cyStart+(n-1));
            end

            if (cxStart+(n-1)) > imageWidth 
                xEnd = imageWidth ;
            else 
                xEnd = (cxStart+(n-1)) ;
            end           

%         pause(2)
if yStart > yEnd || xStart > xEnd
    disp('ISSUE PPG')
end
            rawPPG = rawPerfussion(yStart:yEnd, xStart:xEnd);            
            rawPPG = mean(rawPPG,1); 
            rawPPG = mean(rawPPG,2);
            
            PPGraw(k,c) = rawPPG(:); 

         end
end

    

end
        
%          % plot
%         figure
%         for c = 1:size(gridPoints,2)
%             hold on
%         plot((1:floor(numFrame)), PPG(:, c), 'b')
%         end
