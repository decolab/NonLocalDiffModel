clear all;
path2=[ '../../Turbulence/Basics'];
addpath(path2);
path3=[ '../../Nonequilibrium/'];
addpath(genpath(path3));
path4=[ '../../Tenet/TENET/'];
addpath(genpath(path4));

NSUB=1000;
NGroup=20;
N=1000;
epsilon=1400;

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

load (['hcp_RELATIONAL_LR_schaefer1000.mat']);

Pstatesemp2=zeros(NSUB,N);
nsub=1;

for sub=1:NGroup:NSUB-NGroup
    sub
    Pm2=zeros(N,N);

    for sub2=sub:sub+NGroup-1
        clear events;
        ts=subject{sub2}.schaeferts ;
        ts=ts(indexregion,:);

        %% Transition matrix empirical

        for seed=1:N
            ts(seed,:)=demean(detrend(ts(seed,:)));
            if sum(isnan(ts(seed,:)))==0
                signal_filt2=(filtfilt(bfilt,afilt,ts(seed,:)));
                tss=demean(signal_filt2(50:end-50));
                ev1=tss>std(tss)+mean(tss);
                ev2=[0 ev1(1:end-1)];
                events(seed,:)=(ev1-ev2)>0;
            else
                events(seed,:)=0*events(seed-1,:);
            end
        end

        for seed=1:N
            lisev=find(events(seed,:)==1);
            for t=lisev
                for t2=t+1:size(events,2)
                    lista=find(events(:,t2)==1);
                    if isempty(lista)
                        Pm2(seed,seed)=Pm2(seed,seed)+1;
                    else
                        Pm2(seed,lista)=Pm2(seed,lista)+1;
                        break;
                    end
                end
            end
        end
    end
    Dmatrix=diag(sum(Pm2,2));
    Dmatrix(555,555)=Dmatrix(554,554);
    Dmatrix(908,908)=Dmatrix(907,907);
    Pmatrixemp=inv(Dmatrix)*Pm2;
    Pmatrixemp100=Pmatrixemp^50;
    Pstatesemp=Pmatrixemp100(1,:);
    Pstatesemp([555 908])=[];

    %% HARM SC
    Kmatrix=zeros(N,N);
    for i=1:N
        for j=1:N
            dij2=sum((SchaeferCOG(i,:)-SchaeferCOG(j,:)).^2);
            Kmatrix(i,j)=exp(-dij2/epsilon);
        end
    end
    Dmatrix=diag(sum(Kmatrix,2));
    Pmatrix=inv(Dmatrix)*Kmatrix;
    Pmatrix100=Pmatrix^50;
    Pstates=Pmatrix100(1,:);

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
    PmatrixC100=PmatrixC^50;
    PstatesC=PmatrixC100(1,:);

    Pstates([555 908])=[];
    PstatesC([555 908])=[];

    KLfitt(nsub)=-log(nansum(sqrt(Pstatesemp.*Pstates)))
    cc=corrcoef(Pstates,Pstatesemp,'rows','complete');
    Corrfitt(nsub)=cc(2)

    KLfittC(nsub)=-log(nansum(sqrt(Pstatesemp.*PstatesC)))
    cc=corrcoef(PstatesC,Pstatesemp,'rows','complete');
    CorrfittC(nsub)=cc(2)
    Pstatesempsub(nsub,:)=Pstatesemp;
    Pstatessub(nsub,:)=Pstates;
    PstatesCsub(nsub,:)=PstatesC;
    nsub=nsub+1;
end

Pstatesempm=nanmean(Pstatesempsub);
Pstatesm=nanmean(Pstatessub);
PstatesCm=nanmean(PstatesCsub);

figure(1);
boxplot([KLfitt' KLfittC']);
ranksum(KLfitt,KLfittC)

figure(2);
boxplot([Corrfitt' CorrfittC']);
ranksum(Corrfitt,CorrfittC)


save results_model_CHARM_subjects_RELATIONAL.mat KLfitt Corrfitt KLfittC CorrfittC Pstatesempm Pstatesm PstatesCm;
