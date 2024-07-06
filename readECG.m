clc
clear
close all
%%

fileName = 'sample_data';
data = edfread([fileName,'.edf']);
info = edfinfo([fileName,'.edf']);

%%
EKG_entry = 11;
fs = info.NumSamples/seconds(info.DataRecordDuration);
fs = fs(EKG_entry);
ECG=[];
for i = 1:height(data)
    ECG = [ECG; data.(EKG_entry){i}];
end
%%
t = (1:length(ECG))/fs;
plot(t,ECG);
xlabel('time(s)')
%% interpolate bad period
ECG(4120*fs+1:4960*fs)=0;
ECG(25890*fs+1:26200*fs)=0;
plot(t,ECG)
%% filp if reversed
ECG = -ECG;
%% apply some filter to clean up the data
fECG = freq_filt(ECG,fs,[2,80],[1,10],'bandpass');
fECG = freq_filt(fECG,fs,60,2,'notch');
%% wavelet to amplify R peak
sig = cwt(fECG, 0.08*fs, 'db3')';
plot(t,sig); hold on
plot(t,fECG)
%% find peaks
findpeaks(sig,fs,'MinPeakHeight',200,'MinPeakDistance',0.4,'MinPeakProminence',0.01);
%% 

[~,peaks]=findpeaks(sig,fs,'MinPeakHeight',250,'MinPeakDistance',0.5,'MinPeakProminence',0.01);
RR = diff(peaks);
plot(RR)
%%
csvwrite([fileName,'_HR.csv'], RR);