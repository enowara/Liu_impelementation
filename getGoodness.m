function [goodness] = getGoodness(Sj, e_i)

goodness = [];
% debugging
% s=1;
% s_i = Sj(:,s);
% e_i = e_next;
    L = size(Sj,1);   % time duration of PPG signals 
    Fs= 30;%v.FrameRate;
    freq =Fs*(0:(L/2))/L; 

for s = 1:size(Sj,2)
    s_i = Sj(:,s);

    r  = 3/60; % Hz % 3 bpm 
    % find max freq of ground truth signal e 
    s_hat = abs(fft(s_i));

    Ye_i = fft(e_i);
    [~,fIDX] = max(abs(Ye_i));
    fHR = freq(fIDX);

    numerator = sum(s_hat((fIDX-r) : (fIDX+r)));
    denominator1 = sum(s_hat) -  numerator;
    goodnessi = numerator/denominator1;
    goodness = [goodness, goodnessi];

end