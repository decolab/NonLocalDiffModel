clear all;
path2=[ '../../Turbulence/Basics'];
addpath(path2);
path3=[ '../../Nonequilibrium/'];
addpath(genpath(path3));
path4=[ '../../Tenet/TENET/'];
addpath(genpath(path4));

NSUB=1003; %% 1003
N=1000;

load sc_schaefer_MK.mat
C=sc_schaefer;
C=C/max(max(C));

Isubdiag = find(tril(ones(N),-1));
indexregion=1:1000;

load schaefercog.mat;
load SCFClongrange.mat;

%%% Define Training and test data
% Parameters of the data
TR=0.72;  % Repetition Time (seconds)
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

load (['hcp_REST_LR_schaefer1000.mat']);

FCemp2=zeros(NSUB,N,N);

for sub=1:NSUB
    sub
    clear events;
    ts=subject{sub}.schaeferts ;
    ts=ts(indexregion,:);

    %% Transition matrix empirical

    for seed=1:N
        ts(seed,:)=demean(detrend(ts(seed,:)));
        if sum(isnan(ts(seed,:)))==0
            signal_filt2=(filtfilt(bfilt,afilt,ts(seed,:)));
            tss=demean(signal_filt2(50:end-50));
            tfilt(seed,:)=tss;
        else
            tfilt(seed,:)=0*tfilt(seed-1,:);
        end
    end
    FCemp2(sub,:,:)=corrcoef(tfilt','rows','complete');
end

FCemp=squeeze(nanmean(FCemp2));

%% FC model

NTRIALS=1000;
NREP=200;
epsilon=1400;

%% HARM

Kmatrix=zeros(N,N);
for i=1:N
    for j=1:N
        dij2=sum((SchaeferCOG(i,:)-SchaeferCOG(j,:)).^2);
        Kmatrix(i,j)=exp(-dij2/epsilon);
    end
end
Dmatrix=diag(sum(Kmatrix,2));
Pmatrix=inv(Dmatrix)*Kmatrix;

%% CHARM SC
Thorizont=2;
Kmatrix=zeros(N,N);

for i=1:N
    for j=1:N
        dij2=sum((SchaeferCOG(i,:)-SchaeferCOG(j,:)).^2);
        Kmatrix(i,j)=exp(complex(0,1)*dij2/epsilon);
    end
end

Ktr_t=Kmatrix^Thorizont;
Ptr_t=abs(Ktr_t).^2;

Dmatrix=diag(sum(Ptr_t,2));
PmatrixC=inv(Dmatrix)*Ptr_t;

FCsim3=zeros(NTRIALS,N,N);
FCsim2=zeros(NREP,N,N);
FCsimC3=zeros(NTRIALS,N,N);
FCsimC2=zeros(NREP,N,N);

for tr=1:NREP
    tr
    for sub=1:NTRIALS
        tssim(:,1)=zeros(N,1);
        tssimC(:,1)=zeros(N,1);
        ini=ceil(rand*N);
        while ini==555 || ini==908
            ini=ceil(rand*N);
        end
        tssim(ini,1)=1;
        tssimC(ini,1)=1;
        for tt=2:1000
            tssim(:,tt)=zeros(N,1);
            while sum(tssim(:,tt))==0
                for i=1:N
                    if tssim(i,tt-1)==1
                        for j=1:N
                            if rand<Pmatrix(i,j)
                                tssim(j,tt)=1;
                            end
                        end
                    end
                end
            end

            tssimC(:,tt)=zeros(N,1);
            while sum(tssimC(:,tt))==0
                for i=1:N
                    if tssimC(i,tt-1)==1
                        for j=1:N
                            if rand<PmatrixC(i,j)
                                tssimC(j,tt)=1;
                            end
                        end
                    end
                end
            end
        end
        FCsim3(sub,:,:)=corrcoef(tssim');
        FCsimC3(sub,:,:)=corrcoef(tssimC');
    end
    FCsim2(tr,:,:)=squeeze(nanmean(FCsim3));
    FCsimC2(tr,:,:)=squeeze(nanmean(FCsimC3));
end
FCsim=squeeze(nanmean(FCsim2));
FCsimC=squeeze(nanmean(FCsimC2));

cc=corrcoef(FCemp(Isubdiag),FCsim(Isubdiag),'rows','complete');
cc(2)

cc=corrcoef(FCemp(Isubdiag),FCsimC(Isubdiag),'rows','complete');
cc(2)

nn=1;
NW=10;
for i=1:NW:200
    FCs=squeeze(nanmean(FCsim2(i:i+NW-1,:,:)));
    FCsC=squeeze(nanmean(FCsimC2(i:i+NW-1,:,:)));
    cc=corrcoef(FCemp(Isubdiag),FCs(Isubdiag),'rows','complete');
    corrFCempHARM(nn)=cc(2);
    cc=corrcoef(FCemp(Isubdiag),FCsC(Isubdiag),'rows','complete');
    corrFCempCHARM(nn)=cc(2);
    fittFCempHARM(nn)=(nanmean((FCemp(Isubdiag)-FCs(Isubdiag)).^2));
    fittFCempCHARM(nn)=(nanmean((FCemp(Isubdiag)-FCsC(Isubdiag)).^2));

    gbcharm=squeeze(nanmean(FCs));
    gbccharm=squeeze(nanmean(FCsC));
    gbcemp=squeeze(nanmean(FCemp));
    cc=corrcoef(gbcemp,gbcharm,'rows','complete');
    corrGBCempHARM(nn)=cc(2);    
    fittGBCempHARM(nn)=(nanmean((gbcemp-gbcharm).^2));
    cc=corrcoef(gbcemp,gbccharm,'rows','complete');
    corrGBCempCHARM(nn)=cc(2);    
    fittGBCempCHARM(nn)=(nanmean((gbcemp-gbccharm).^2));
    %%events
    % cc=corrcoef(FCempeve(Isubdiag),FCs(Isubdiag),'rows','complete');
    % corrFCempeveHARM(nn)=cc(2);
    % cc=corrcoef(FCempeve(Isubdiag),FCsC(Isubdiag),'rows','complete');
    % corrFCempeveCHARM(nn)=cc(2);
    % fittFCempeveHARM(nn)=(nanmean((FCempeve(Isubdiag)-FCs(Isubdiag)).^2));
    % fittFCempeveCHARM(nn)=(nanmean((FCempeve(Isubdiag)-FCsC(Isubdiag)).^2));
    %%
    nn=nn+1;
end

figure(1)
boxplot([corrFCempHARM' corrFCempCHARM'])
ranksum(corrFCempHARM,corrFCempCHARM)

figure(2)
boxplot([fittFCempHARM'  fittFCempCHARM'])
ranksum(fittFCempHARM,fittFCempCHARM)

figure(3)
boxplot([corrGBCempHARM' corrGBCempCHARM'])
ranksum(corrGBCempHARM,corrGBCempCHARM)

figure(4)
boxplot([fittGBCempHARM'  fittGBCempCHARM'])
ranksum(fittGBCempHARM,fittGBCempCHARM)
% figure(3)
% boxplot([corrFCempeveHARM' corrFCempeveCHARM'])
% ranksum(corrFCempeveHARM,corrFCempeveCHARM)
% 
% figure(4)
% boxplot([fittFCempeveHARM' fittFCempeveCHARM'])
% ranksum(fittFCempeveHARM,fittFCempeveCHARM)

save results_FCmodel.mat FCemp FCsim2 FCsimC2;
