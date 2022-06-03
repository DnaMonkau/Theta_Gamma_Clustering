%%%%%%%This code (includes Step1, Step2, Step3)show how cluster individual theta cycle of LFP recorded in hippocampal CA1 pyrmidal layer
%%%%%%%into four clusters:slow gamma; median gamma; early high gamma and late high gamma through k-means clustering
%%%%%%%author: Lu Zhang;math2437@hotmail.com; tested with Matlab 2017b;

%%%%toolbox needed:
%%%%Toolbox 1:buzsaki Lab code for processing hc-11 data; not necessary if use your own
%%%%data:https://github.com/buzsakilab/buzcode

%%%%Toolbox 2:Comunity Clustering;not necessary if only try k-means
%%%%clustering;http://netwiki.amath.unc.edu/GenLouvain/GenLouvain
% 

%%%%%%%Step 1 includes loading hc-11 data;please download the data first if not use your own data;
%% https://crcns.org/data-sets/hc/hc-11/about-hc-11

%%%%%%%Step 2 prepare data samples for clustering; 
%%%%%%%LFP, filtered LFP in theta band,LFP theta phase is needed;



%%%%Step 1; prepare the data from hc-11 data set; This is not necesary if
%%%%you have your own data instead.
%%%%samprate=1250;%%sampling rate
%%%%%%%%%Load LFP Chan 30 from example data in hc-11 data set
clear all
data = load('Theta_Values.mat');
% for data in file
samprate=1000;%%%%%sampling rate of LFP in .eeg filestructs
% s=[];
for j =1:16
%     disp(j)
    lfp = transpose(data.structs.REM_Data{j});
    [len, extra]= size(lfp);
    for l =1:extra
    lfp(:,l)=cast(lfp(:,l), 'double');
    end
    %%%GET FILTERED THETA
    
    amp = data.structs.Amplitude{j};
    Thetaphase = data.structs.Phase{j};%%%%%%%%%instanous theta phase;
    lfpTheta = data.structs.Filtered{j};
%     disp(lfpTheta)
    
    %%%%Step 2; creat frequency-phase spectrum for each individual theta cycle
    WaveParam.ntw=11;      %11 points smoothing time-window 
    WaveParam.nsw=3;        %3 points smoothing fre-window 
    WaveParam.DownSample=4; %DownSampling parameter;
    WaveParam.wname='morl'; %Morlet Wavelet
    WaveParam.Samplingrate=samprate;
    WaveParam.Freq=[20:2:120]; %Frequeny band interested in Hz
    WaveParam.Zscore=1;  %Normalize waveletPower by zscore in time domain for a given freqency point;
    timerange=ones([1, extra]);
    t_list=[];
    for l = 1:extra
    
        t=length(lfp(:,l));   %%%%%Using only the Linear Track Trial ;
        t_list=[t_list,t];
    %%%%Analysis period from 0 to 40000 )seconds;Using large end (15000) could ensure using the whole data.
    end
    timerange=[timerange;t_list];
    % timerange=timerange{:,:};
%     disp(timerange)
    PhaseBinNum=20; %%%%%%Divid each theta cycle by 20 phase bin;
    PhaseStep=2*pi/PhaseBinNum;
    PhaseBin=-pi:PhaseStep:pi;
    PhaseBin=(PhaseBin(1:end-1)+PhaseBin(2:end))/2;
%     disp(size(lfp));
    Thetaphase=transpose(Thetaphase);
    lfpTheta = transpose(lfpTheta);
    ParamOut = SpecThetaExtract(lfp,lfpTheta,Thetaphase,samprate,timerange,WaveParam,PhaseStep);
    Div=ParamOut.Div;
    Sample = ParamOut.sample';
%     disp(size(Sample));
    s{j}=Sample;
%%%%cluster the first 5000 samples;
SampleCorr=corr(ParamOut.sample);
AdjMat=SampleCorr+1;  %%%%%%%%%%Non-negative weighted correlation.
clear SampleCorr;
for itt=1:size(AdjMat,1)
    AdjMat(itt,itt)=0;%%%%%%diagonal zeros for Adjcent matrix
end
k = full(sum(AdjMat));
twom = sum(k);
%         B = full(Adj - gamma*k'*k/twom);
tic
k = full(sum(AdjMat));
twom = sum(k);
gamma = 1;  %%Community Clustering Paramter%%
limit =10000; %%memory consideration for community clustering 

B = @(i) AdjMat(:,i) - gamma*k'*k(i)/twom;
[b{j},c{j}]=community_louvain(AdjMat);
end
save('Matrix120.mat','s', '-mat')
save('Cluster120.mat', 'b', '-mat')

