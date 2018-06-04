function aE_UP_DOWN_detect_field(valtozok,dirs,xlsdata)    
overwrite=valtozok.overwrite;

% Parameters
% LFP traces pre-processing
sfLFP = 1000;   % LFP decimated to 1KHz
lowPassCutOff = 500; % Hz
% UP & DOWN states detection parameters
stateDurTh_ms = 100;    % 100 ms
stateIntervTh_ms = 50;  % 50 ms
nSTD_UP = 2;
nSTD_DOWN = 2;
stateDurTh = stateDurTh_ms*1e-3*sfLFP;
stateIntervTh = stateIntervTh_ms*1e-3*sfLFP;
% Frequency bands for LFP analysis
fBands = [.2 2;2 4];
fHigh = [20 100];
% Likelihood computation
theta = [193.33289; 189.12116]./360*2*pi; % for hilbert: come from analysis of intracellular traces --> 0-1, 1-3 Hz
nBins = 100;
phaseMethod = 'hilbert';
% phaseMethod = 'interpol';
% Combined
combFlag = 1;   % this parameter must be set to 1 if you want to include info about high frequencies
fAlphaBeta = [10 40];
winRmsFilt = 5*1e-3;
winSmoothFilt = 50*1e-3;
win_samples = winRmsFilt*sfLFP;
winSmooth_samples = winSmoothFilt*sfLFP;
if mod(winSmooth_samples,2)==0
    winSmooth_samples = winSmooth_samples+1;
end
% Transition-triggered analysis of LFP
halfTransTriggWin = 0.5*sfLFP; % samples


segmentlength=valtozok.segmentlength;
for filei=1:length({xlsdata.ID})
        a=dir([dirs.statedir,xlsdata(filei).ID,'.mat']);
        if xlsdata(filei).field==1 & (isempty(a) | overwrite==1)
        disp(['exporting state transitions from  ', xlsdata(filei).ID])
        %      [Selection,ok] = listdlg('ListString',{xlsdata.ID},'ListSize',[300 600]); % az XLS file alapján kiválasztjuk, hogy melyik file összes mérésén szeretnénk végigmenni
        temp=load([dirs.bridgeddir,xlsdata(filei).ID]);
        bridgeddata=temp.bridgeddata;
        statedata.UP=struct;
        statedata.DOWN=struct;
        for sweep=1:length(bridgeddata)
            progressbar(sweep/length(bridgeddata))
            si=bridgeddata(sweep).si;
            time=[1:length(bridgeddata(sweep).y)]*si-si;
            duration=time(end)-time(1);
            if duration>segmentlength
                Y=bridgeddata(sweep).y;
                sf=1/si;
                lowPassCutOffNorm = lowPassCutOff/sf;
                [bLow,aLow] = ellip(2,0.1,40,lowPassCutOffNorm);
                % Pre-processing of LFP trace
                % Low-pass filtering of LFP
                LFPseries_mV = Y.*1e+3;    % 3rd column
                LFPdataFilt = filtfilt(bLow,aLow,LFPseries_mV); % LP filter, cut-off 500 Hz
                % Downsampling to 1 KHz
                nSamples = length(LFPdataFilt);
                nSamplesDec = ceil(nSamples/(sf/sfLFP));
                timeDec = (1/sfLFP).*(1:nSamplesDec);
                LFPdata_dec = decimate(LFPdataFilt,round(sf/sfLFP));
                % Filter LFP in the delta frequency range
                [b_deltaFilt,a_deltaFilt] = ellip(2,0.1,40,[0.1 4]./sfLFP);
                deltaFiltLFP = filtfilt(b_deltaFilt,a_deltaFilt,LFPdata_dec);
                if valtozok.plotthestuff==1
                    %% Plot raw signals and synchronization index
                    figure(3)
                    clf
                    h(1)=subplot(4,1,1);
                    plot(timeDec,deltaFiltLFP)
                    h(2)=subplot(4,1,2);
                    plot(timeDec,LFPdata_dec,'k');
                    xlim([0 max(timeDec)])
                    title('LFP')
                    xlabel('Time [s]')
                    ylabel('Voltage [mV]')
                    hold all
                end
                %% Filter LFP in different frequency bands
                LFPdataDecFilt = zeros(length(LFPdata_dec),size(fBands,1));
                for bb = 1:size(fBands,1)
                    if fBands(bb,1) == 0 % lowpass
                        [b_LFPBandFilt,a_LFPBandFilt] = ellip(2,0.1,40,fBands(bb,2)./sfLFP,'low');
                    else    % bandpass
                        [b_LFPBandFilt,a_LFPBandFilt] = ellip(2,0.1,40,fBands(bb,:)./sfLFP);
                    end
                    LFPdataDecFilt(:,bb) = filtfilt(b_LFPBandFilt,a_LFPBandFilt,LFPdata_dec);
                end
                %% For each frequency band, approximate L_t
                k = zeros(size(LFPdataDecFilt,1),size(fBands,1));
                phi = zeros(size(LFPdataDecFilt,1),size(fBands,1));
                L_t = zeros(size(LFPdataDecFilt,1),size(fBands,1));
                for bb = 1:size(fBands,1)
                    [k(:,bb),phi(:,bb)] = UP_DOWN_DET_compInstPhaseAmpl(LFPdataDecFilt(:,bb),phaseMethod);
                    phi(phi(:,bb)<=0,bb) = phi(phi(:,bb)<=0,bb)+2*pi; % rescale from 0 to 2*pi
                    phi(:,bb) = round(phi(:,bb)./(2*pi/360)).*(2*pi/360); % round to 1 degree resolution
                    L_t(:,bb) = cos(phi(:,bb)-theta(bb));
                end
                %% Compute k_high
                [b_high,a_high] = ellip(2,0.1,40,fHigh./sfLFP);
                LFPdataHigh = filtfilt(b_high,a_high,LFPdata_dec);
                [k_high,phi_high] = UP_DOWN_DET_compInstPhaseAmpl(LFPdataHigh,'hilbert');
                %% Derive S_LFP
                K = k./repmat(k_high'+sum(k,2),1,size(fBands,1));
                S_LFP = (sum(K.*L_t,2)+1)./2;
                %% Compute k_alphaBeta
                if combFlag==1
                    [b_alphaBeta,a_alphaBeta] = ellip(2,0.1,40,fAlphaBeta./sfLFP);
                    LFPdataAlphaBeta = filtfilt(b_alphaBeta,a_alphaBeta,LFPdata_dec);
                    rmsLFPdata = rmsFilt(LFPdataAlphaBeta,win_samples);
                    smoothLFPdata = smooth(rmsLFPdata,winSmooth_samples);
                    S_alphaBeta = smoothLFPdata./max(smoothLFPdata);
                    % Compute S_comb
                    S_comb = (S_LFP+S_alphaBeta)./2;
                else
                    S_comb = S_LFP;
                end
                %% Fit distribution of S_comb by GM
                options = statset('Display','final','MaxIter',500);
                try
                    GMmodels1 = gmdistribution.fit(S_comb,3,'Options',options);
                    GMmodels2= gmdistribution.fit(S_comb,3,'Options',options);
                    GMmodels3 = gmdistribution.fit(S_comb,3,'Options',options);
                    [~,selectedidx]=min([GMmodels1.NlogL,GMmodels2.NlogL,GMmodels3.NlogL]);
                    if selectedidx==1
                        GMmodel=GMmodels1;
                    elseif selectedidx==2
                        GMmodel=GMmodels2;
                        elseif selectedidx==3
                            GMmodel=GMmodels3;
                    end
                    
                catch err
                    keyboard;
%                     continue;
                end
                tonan=[GMmodel.PComponents]<.05;
                mu = GMmodel.mu;
                mu(tonan)=NaN;
                sigma = squeeze(GMmodel.Sigma);
                sigma(tonan)=NaN;
                if valtozok.plotthestuff==1
                    %%
                    figure(2)
                    clf
                    [hAmpl,bins] = hist(S_comb,nBins);
                    bar(bins,hAmpl./length(S_comb).*100);
                    ylabel('# occurrences')
                    xlabel('S_{LFP}')
                    hold all
                    yFit = pdf(GMmodel,bins');
                    plot(bins,yFit,'r');
                end
                %% Compute thresholds for UP & DOWN state detection
                [~,idxUP] = nanmax(mu);
                [~,idxDOWN] = nanmin(mu);
                thUP = mu(idxUP)-nSTD_UP*sqrt(sigma(idxUP));
                thDOWN = mu(idxDOWN)+nSTD_DOWN*sqrt(sigma(idxDOWN));
                if valtozok.plotthestuff==1
                    figure(2)
                    line([thUP thUP],[0 max(yFit)],'LineStyle','--','LineWidth',1,'Color','c')
                    line([thDOWN thDOWN],[0 max(yFit)],'LineStyle','--','LineWidth',1,'Color','g')
                end
                %% Detect UP and DOWN states from LFP trace and save results
                [UP_states_DET, DOWN_states_DET] = UP_DOWN_DET_detectStates(S_comb, thUP, thDOWN, stateIntervTh, stateDurTh);
                UP_state_freq = length(UP_states_DET)./(length(LFPdata_dec)/sfLFP);
                DOWN_state_freq = length(DOWN_states_DET)./(length(LFPdata_dec)/sfLFP);
                mUPstateDur = mean((UP_states_DET(:,2)-UP_states_DET(:,1))./sfLFP);
                stdUPstateDur = std((UP_states_DET(:,2)-UP_states_DET(:,1))./sfLFP);
                mDOWNstateDur = mean((DOWN_states_DET(:,2)-DOWN_states_DET(:,1))./sfLFP);
                stdDOWNstateDur = std((DOWN_states_DET(:,2)-DOWN_states_DET(:,1))./sfLFP);
                UPDOWNstates_stat = [UP_state_freq DOWN_state_freq mUPstateDur stdUPstateDur mDOWNstateDur stdDOWNstateDur];
                if valtozok.plotthestuff==1
                    figure(3)
                    subplot(4,1,2)
                    hold on
                    for uu = 1:size(UP_states_DET,1)
                        plot(timeDec(UP_states_DET(uu,1):UP_states_DET(uu,2)),LFPdata_dec(UP_states_DET(uu,1):UP_states_DET(uu,2)),'Color','r','LineWidth',2)
                    end
                    for dd = 1:size(DOWN_states_DET,1)
                        plot(timeDec(DOWN_states_DET(dd,1):DOWN_states_DET(dd,2)),LFPdata_dec(DOWN_states_DET(dd,1):DOWN_states_DET(dd,2)),'Color','b','LineWidth',2)
                    end
                    %% Plot S_comb & thresholds
                    figure(3)
                    h(3)=subplot(4,1,3);
                    plot(timeDec,S_comb)
                    hold all
                    line(xlim,[thUP thUP],'Color','r','LineWidth',2,'LineStyle','--');
                    line(xlim,[thDOWN thDOWN],'Color','b','LineWidth',2,'LineStyle','-.');
                    h(4)=subplot(4,1,4);
                    plot(timeDec,radtodeg(phi))
                    
                    linkaxes(h,'x')
                    pause
                end
                

                downsampletimes=round(sf/sfLFP);
                for upi=1:size(UP_states_DET,1)
                    if isempty(fieldnames(statedata.UP))
                        NEXT=1;
                    else
                        NEXT=length(statedata.UP)+1;
                    end
                    idx=UP_states_DET(upi,1):UP_states_DET(upi,2);
                    statedata.UP(NEXT).type='UP';
                    statedata.UP(NEXT).onseth=UP_states_DET(upi,1)*downsampletimes;
                    statedata.UP(NEXT).endh=UP_states_DET(upi,2)*downsampletimes;
                    statedata.UP(NEXT).onsett=time(idx(1)*downsampletimes)+bridgeddata(sweep).realtime;
                    statedata.UP(NEXT).endt=time((idx(end)-1)*downsampletimes)+bridgeddata(sweep).realtime;
                    statedata.UP(NEXT).midt=median(timeDec(idx));
                    statedata.UP(NEXT).duration=length(idx)*si*downsampletimes;
                    statedata.UP(NEXT).value=mean(LFPdataDecFilt(idx));
                    statedata.UP(NEXT).sweepnum=sweep;
                end
                for downi=1:size(DOWN_states_DET,1)
                    if isempty(fieldnames(statedata.DOWN))
                        NEXT=1;
                    else
                        NEXT=length(statedata.DOWN)+1;
                    end
                    idx=DOWN_states_DET(downi,1):DOWN_states_DET(downi,2);
                    statedata.DOWN(NEXT).type='DOWN';
                    statedata.DOWN(NEXT).onseth=DOWN_states_DET(downi,1)*downsampletimes;
                    statedata.DOWN(NEXT).endh=DOWN_states_DET(downi,2)*downsampletimes;
                    statedata.DOWN(NEXT).onsett=time(idx(1)*downsampletimes)+bridgeddata(sweep).realtime;
                    statedata.DOWN(NEXT).endt=time((idx(end)-1)*downsampletimes)+bridgeddata(sweep).realtime;
                    statedata.DOWN(NEXT).midt=median(timeDec(idx));
                    statedata.DOWN(NEXT).duration=length(idx)*si*downsampletimes;
                    statedata.DOWN(NEXT).value=mean(LFPdataDecFilt(idx));
                    statedata.DOWN(NEXT).sweepnum=sweep;
                end
            end
        end
        save([dirs.statedir,xlsdata(filei).ID],'statedata');
        end
    end
   