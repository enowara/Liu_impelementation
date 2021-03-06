function PulseXY = getpulseXY(PPGr, PPGg, PPGb, Fps)
% Computes a motion robust pulse signal based on chrominance-based method from
%     Gerard's 2013 paper

% debugging
% PPGr = PPGrawR;
% PPGg = PPGrawG;
% PPGb = PPGrawB;

L = size(PPGr,1);

% PPG_red: non-filtered raw red channel recordings
for c = 1:size(PPGr,2)

%     Rn(:,c) = PPGr(:,c)./mean(PPGr,2);   % mean over a window
%     Gn(:,c) = PPGg(:,c)./mean(PPGg,2);
%     Bn(:,c) = PPGb(:,c)./mean(PPGb,2);
% for tt = 1:floor(L/(2.5*Fps))  % take 2.5 seconds in frames
%     meanShortR = mean(PPGr(((tt-1)*75+1) : tt*75, 2 )); 
%     Rn(((tt-1)*75+1 : tt*75),c) = PPGr(((tt-1)*75+1 : tt*75),c)./meanShortR;
%     meanShortG = mean(PPGg(((tt-1)*75+1) : tt*75, 2 )); 
%     Gn(((tt-1)*75+1 : tt*75),c) = PPGg(((tt-1)*75+1 : tt*75),c)./meanShortG;
%     meanShortB = mean(PPGb(((tt-1)*75+1) : tt*75, 2 )); 
%     Bn(((tt-1)*75+1 : tt*75),c) = PPGb(((tt-1)*75+1 : tt*75),c)./meanShortB;
% end
timeWin = floor(2.5*Fps);
for tt = 1:timeWin:L-timeWin % take 2.5 seconds in frames
    meanShortR = mean(PPGr(tt:tt+timeWin)); 
    Rn((tt:tt+timeWin),c) = PPGr((tt:tt+timeWin),c)./meanShortR;
    meanShortG = mean(PPGg(tt:tt+timeWin)); 
    Gn((tt:tt+timeWin),c) = PPGg((tt:tt+timeWin),c)./meanShortG;
    meanShortB = mean(PPGb(tt:tt+timeWin)); 
    Bn((tt:tt+timeWin),c) = PPGb((tt:tt+timeWin),c)./meanShortB;
end

% if the difference between the current averaged signal length is much
% smaller than the actual length
if (L-size(Rn,1)) > 1*Fps
    % add the last end part to the PPG signal
    meanShortR = mean(PPGr(size(Rn,1):L)); 
    Rn((size(Rn,1):L),c) = PPGr((size(Rn,1):L),c)./meanShortR;
    meanShortG = mean(PPGg(size(Gn,1):L)); 
    Gn((size(Gn,1):L),c) = PPGg((size(Gn,1):L),c)./meanShortG;
    meanShortB = mean(PPGb(size(Bn,1):L)); 
    Bn((size(Bn,1):L),c) = PPGb((size(Bn,1):L),c)./meanShortB;
end

end % end c

X_s = 3*Rn -2*Gn;
Y_s = 1.5*Rn +Gn - 1.5*Bn;

% bandpass filter the signals
if Fps == 30
    load('../highpass_05_30.mat') % replace with bandpass filters for 25 fps
    load('../lowpass_5_30.mat')
    for c = 1:size(PPGr,2)
        PPG2filt0 =  X_s(:,c);
        PPG2filt1 = filtfilt(lowpass_5,1,PPG2filt0);  %low pass
        PPG2filt2 = filtfilt(highpass_05,IIR_part,PPG2filt1); % high pass
        X_f(:,c) = PPG2filt2;

        PPG2filt0 =  Y_s(:,c);
        PPG2filt1 = filtfilt(lowpass_5,1,PPG2filt0);  %low pass
        PPG2filt2 = filtfilt(highpass_05,IIR_part,PPG2filt1); % high pass
        Y_f(:,c) = PPG2filt2;
    end 

elseif Fps == 25
    load('../EwaHighPass.mat')
    load('../EwaLowPass.mat')
    for c = 1:size(PPGr,2)
        PPG2filt0 =  X_s(:,c);
        PPG2filt1 = filtfilt(b, 1, PPG2filt0);  %low pass
        PPG2filt2 = filtfilt(b_i, a_i,PPG2filt1); % high pass
        X_f(:,c) = PPG2filt2;

        PPG2filt0 =  Y_s(:,c);
        PPG2filt1 = filtfilt(b, 1, PPG2filt0);  %low pass
        PPG2filt2 = filtfilt(b_i, a_i, PPG2filt1); % high pass
        Y_f(:,c) = PPG2filt2;
    end 
end


alpha = (X_f)./(Y_f);
S = X_f - alpha.*Y_f;
PulseXY = S;% X_f ;
end