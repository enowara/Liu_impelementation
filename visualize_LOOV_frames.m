% [scores_SVM_postcell, scores_SVMcell, Ytsscell, Ytrscell, testPeople, labelsSVMcell, ...
%        orderTscell]  = SVM_LiuTogether(saveLiuFolder, N, q, liveFolders, fakeFolders, testPerson, trainPeople);
%    
   
folderMain = '/media/ewa/SH/DatasetsAntiSpoof/3dmadDirectories/';


figure,
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
    trainPeople = setdiff(allPeople, testPerson); 
    
   % visualize video frames of LOOV, tr and ts subjects
   for f = 1:3
       folderEnd = ['Data0' num2str(f) 'Keep/'];

        fileNameList = dir([[folderMain folderEnd] ['*' '.avi']]); %

        for i =1:length(fileNameList)
            imgCells{i} = fileNameList(i).name;  
        end

        [cs,index] = sort_nat(imgCells,'ascend');
        img_names = cs;

       for m_ts = trainPeople  % trainPeople
            vidName = img_names{m_ts};       
            % read in the videos
            v = VideoReader([folderMain folderEnd vidName]);
            frames = read(v);
            firstFrame = frames(:,:,:,1);
            imshow(firstFrame)
            title([num2str(f) ' - ' num2str(m_ts) 'test person p = ' num2str(p)])
%             pause(1)
            drawnow
       end
   end
end
   
