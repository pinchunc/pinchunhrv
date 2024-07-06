clear all
clc
close all
SID = 'sample_data'
%load([SID '_RR.mat'])
filename = ([SID '_RR.txt'])
%filenameid = strrep(filename,'.txt','')
load([SID '_PC.mat'])
rr = readtable(filename);
Res.CNT = struct;
%Res.CNT.rate.ECG1 = 250;
Res.HRV = struct;
Res.HRV.Data.T_RR = table2array(rr(:,3));
Res.CNT.Offset = table2array(rr(1,3))
RES=Res;
fs=250;

RR_tot_ind=round((RES.HRV.Data.T_RR-RES.CNT.Offset));
RR_tot=(RR_tot_ind(2:end)-RR_tot_ind(1:end-1))./fs;

RR_tot_ind(1)=[];
rj=find(RR_tot>1.7 | RR_tot<0.5);
RR_tot(rj)=[];
RR_tot_ind(rj)=[];

RR_tot_time=RR_tot_ind./fs;
RRts=spline(RR_tot_ind./fs,RR_tot,1/fs:0.25:RR_tot_ind(end)/fs);
RRts_time=1/fs:0.25:RR_tot_ind(end)/fs;

uw=1.5; %undisturbed minutes
aw=5; %window size
epl=30 % epoch length in sec
mrkr=stageData.stages;
mdiv=floor(length(mrkr));
t=find((mrkr(2:end)-mrkr(1:end-1))~=0); %transition of stage
t=cat(1,0,t);
t=cat(1,t,length(mrkr));
bnd=zeros(length(t)-1,2);
for i=1:length(t)-1
    bnd(i,:)=[t(i)+1 t(i+1)];
end
bmn=(bnd(:,2)-bnd(:,1)+1).*epl./60;
bndstg=mrkr(bnd(:,1));
mintof=find(bmn>(uw+aw));
analepochs_mrk=bndstg(mintof);
analepochs_dur=bmn(mintof);
ContEpoch=floor((analepochs_dur-uw)./aw); %column 5
analepochs_bnd=bnd(mintof,:);
ll=length(analepochs_dur);
for j=1:length(analepochs_dur)
    t1=(analepochs_bnd(j,1)-1)*epl+uw*60;
    start_time{j,1}=[num2str(floor((t1)/60)) ':' num2str((t1)-60*floor((t1)/60))];
    t3=analepochs_bnd(j,2)*epl;
    t2=t1+60*aw*floor((t3-t1)/(aw*60));
    x=RR_tot(find(RR_tot_time>t1 & RR_tot_time<=t2));
    x=reshape(x,1,length(x));
    xts=RRts(find(RRts_time>t1 & RRts_time<=t2)); 
    dxx=1000*abs(x(2:end)-x(1:end-1));
    if analepochs_mrk(j)==1
        stage{j,1}='N1'; % col4
    elseif analepochs_mrk(j)==2
        stage{j,1}='N2';
    elseif analepochs_mrk(j)==3
        stage{j,1}='N3';
    elseif analepochs_mrk(j)==5
        stage{j,1}='REM';
    elseif analepochs_mrk(j)==0
        stage{j,1}='Wake';
    elseif analepochs_mrk(j)==-1
        stage{j,1}='NoStage';
    elseif analepochs_mrk(j)==7
        stage{j,1}='NoStage';
    end
    x(find(x>2.5))=[]; x(find(x<0.35))=[];
    HRmean(j,1)=mean(60./x);
    RRmean(j,1)=mean(x*1000);
    SDNN(j,1)=std(x*1000);
    SDSD(j,1)=std(dxx);
    RMSSD(j,1)=sqrt(sum(dxx.^2)./(length(dxx)-1));
    NN50(j,1)=length(find(dxx>50));
    pNN50(j,1)=length(find(dxx>50))/length(dxx)*100;
    [px,f]=pyulear(detrend(xts*1000),16,0:1/1000:1,4);
    VLF(j,1)=sum(px(find(f>=0.0033 & f<0.04)))/1000;
    LF(j,1)=sum(px(find(f>=0.04 & f<0.15)))/1000;
    HF(j,1)=sum(px(find(f>=0.15 & f<=0.4)))/1000;
    Total(j,1)=sum(px)/1000;
    LFovHF(j,1)=LF(j,1)/HF(j,1);
    LFnu(j,1)=LF(j,1)/(LF(j,1)+HF(j,1));
    HFnu(j,1)=HF(j,1)/(LF(j,1)+HF(j,1));
    [~,kk]=max(px(find(f>=0.04 & f<0.15)));
    LFpk(j,1)=0.04+f(kk);
    [~,kk]=max(px(find(f>=0.15 & f<=0.4)));
    HFpk(j,1)=0.15+f(kk);
    MPF(j,1)=sum(f.*px')/sum(px);
    Epoch(j,1)=j;
    ContEpochs(j,1)=ContEpoch(j);
    EpochDuration(j,1)=analepochs_dur(j);
end
output = table(start_time,Epoch,stage,ContEpochs,EpochDuration,RRmean,HRmean,SDNN,SDSD,...
    RMSSD,NN50,pNN50,VLF,LF,HF,Total,LFovHF,LFnu,HFnu,LFpk,...
    HFpk,MPF);
writetable(output,[SID '_hrv.csv'],'Delimiter',',','WriteVariableNames',true);