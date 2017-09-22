function PulseXY = pulseXY(PPGr, PPGg, PPGb)
% debugging
% PPGr = PPGrawR;
% PPGg = PPGrawG;
% PPGb = PPGrawB;

L = size(PPGr,1);

% Computes a motion robust pulse signal based on chrominance-based method from
%     Gerard's 2013 paper

% PPG_red: non-filtered raw red channel recordings
for c = 1:size(PPGr,2)

%     Rn(:,c) = PPGr(:,c)./mean(PPGr,2);   % mean over a window
%     Gn(:,c) = PPGg(:,c)./mean(PPGg,2);
%     Bn(:,c) = PPGb(:,c)./mean(PPGb,2);
for tt = 1:L/75
    meanShortR = mean(PPGr(((tt-1)*75+1) : tt*75, 2 )); 
    Rn(((tt-1)*75+1 : tt*75),c) = PPGr(((tt-1)*75+1 : tt*75),c)./meanShortR;
    meanShortG = mean(PPGg(((tt-1)*75+1) : tt*75, 2 )); 
    Gn(((tt-1)*75+1 : tt*75),c) = PPGg(((tt-1)*75+1 : tt*75),c)./meanShortG;
    meanShortB = mean(PPGb(((tt-1)*75+1) : tt*75, 2 )); 
    Bn(((tt-1)*75+1 : tt*75),c) = PPGb(((tt-1)*75+1 : tt*75),c)./meanShortB;
end


end
X_s = 3*Rn -2*Gn;
Y_s = 1.5*Rn +Gn - 1.5*Bn;

% filter the signals
load('../highpass_05_30.mat')
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


alpha = (X_f)./(Y_f);
S = X_f - alpha.*Y_f;
PulseXY = S;% X_f ;
end