%%

close all
clear all
reloadxls=0;
projectnames={'CB1elphys','InVivo','Persistent-ChRstim','persistent firing','bleb recording'};
% projectnum=3;

% [projectnum,ok] = listdlg('ListString',projectnames,'ListSize',[300 600]); % az XLS file alapján kiválasztjuk, hogy melyik file összes mérésén szeretnénk végigmenni
projectdata.owbridge=0;
projectdata.owbridge=0;
projectdata.owevent=0;
h = aE_projectselector(projectnames);
uiwait(h);
projectnum=projectdata.projectnum;


% alapadatok


if projectnum==1;
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/Human_rosehip/CB1elphys/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'Figures/'];
    amplifier='HEKA';
    xlsdata=aE_readxls([dirs.basedir,'cb1elphys.xls']);
elseif projectnum==2;
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/_persistent/_InVivo/'];
    %     dirs.rawdir=[locations.tgtardir,'AXONdata/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    dirs.statedir=[dirs.basedir,'State/'];
    dirs.taxonomydir=[locations.tgtardir,'ANALYSISdata/marci/_Taxonomy/persistent_invivo/'];
    dirs.PSDdir=[dirs.basedir,'PSD/'];
    dirs.PSDdir_high=[dirs.basedir,'PSD_high/'];
    dirs.PSDdir_log=[dirs.basedir,'PSD_log/'];
    dirs.videodir=[dirs.basedir,'Videodata/'];
    dirs.brainstatedir=[dirs.basedir,'BrainState/'];
    % dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    % dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    % dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'Figures/'];
    dirs.breathingdir=[dirs.basedir,'Breathing/'];
    dirs.breathingdir_psd=[dirs.basedir,'Breathing_PSD/'];
    dirs.offsetdir=[dirs.basedir,'Offset/'];
    dirs.cross_spectrum_breathing=[dirs.basedir,'Cross_spectrum_with_breathing/'];
    if reloadxls==1
        xlsdata=aE_readxls([dirs.basedir,'invivodata.xls']);
    else
        load([dirs.basedir,'xlsdata'])
    end
    amplifier='HEKA';
elseif projectnum==3;
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/_persistent/_ChRstim/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    %     dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    %     dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    %     dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'Figures/'];
    amplifier='HEKA';
    xlsdata=aE_readxls([dirs.basedir,'ChRstimdata_windows.xls']);
elseif projectnum==4
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.EMdir,'ANALYSISdata/marci/_persistent/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.taxonomydir=[locations.tgtardir,'ANALYSISdata/marci/_Taxonomy/persistent_slice/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    dirs.PSDdir=[dirs.basedir,'PSD/'];
    dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'figures/'];
    amplifier='HEKA';
    xlsdata=aE_readxls([dirs.basedir,'persistentdata_windows.xls']);
    dirs.v0distdir=[dirs.basedir,'v0_dist/'];
elseif projectnum==5
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/_persistent/_BlebRecording/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    %     dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    %     dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    %     dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'Figures/'];
    amplifier='HEKA';
    xlsdata=aE_readxls([dirs.basedir,'blebdata_windows.xls']);
end
if projectdata.doaAPstatistics==1
    xlsdata = persistent_aAPstatistics(dirs,xlsdata);
end
if projectdata.inspecttraces==1
    aE_InspectTraces(dirs,xlsdata);
    return
elseif projectdata.selectvideoROIs==1
    aE_videoanalyzer_ROIselector(dirs, xlsdata);
    return
end
%% create neurontaxonomy xls file
if isfield(dirs,'taxonomydir')
    path=[locations.matlabstuffdir,'NotMine/20130227_xlwrite/'];
    javaaddpath([path 'poi_library/poi-3.8-20120326.jar']);
    javaaddpath([path 'poi_library/poi-ooxml-3.8-20120326.jar']);
    javaaddpath([path 'poi_library/poi-ooxml-schemas-3.8-20120326.jar']);
    javaaddpath([path 'poi_library/xmlbeans-2.3.0.jar']);
    javaaddpath([path 'poi_library/dom4j-1.6.1.jar']);
    path=[locations.matlabstuffdir,'20130227_xlwrite/'];
    javaaddpath([path 'poi_library/poi-3.8-20120326.jar']);
    javaaddpath([path 'poi_library/poi-ooxml-3.8-20120326.jar']);
    javaaddpath([path 'poi_library/poi-ooxml-schemas-3.8-20120326.jar']);
    javaaddpath([path 'poi_library/xmlbeans-2.3.0.jar']);
    javaaddpath([path 'poi_library/dom4j-1.6.1.jar']);
    taxdata=struct;
    for xlsi=1:length(xlsdata)
        if xlsdata(xlsi).field==0 && xlsdata(xlsi).juxta==0
            if isempty(fieldnames(taxdata))
                NEXT=1;
            else
                NEXT=length(taxdata)+1;
            end
            gsc=xlsdata(xlsi).G_S_C;
            hyps=strfind(gsc,'_');
            g=num2str(gsc(1:hyps(1)-1));
            s=num2str(gsc(hyps(1)+1:hyps(2)-1));
            c=num2str(gsc(hyps(2)+1:end));
            
            taxdata(NEXT).anatgroup='0';
            taxdata(NEXT).g=g;
            taxdata(NEXT).s=s;
            taxdata(NEXT).c=c;
            taxdata(NEXT).ID=num2str(xlsi);
            taxdata(NEXT).fname=xlsdata(xlsi).HEKAfname;
            if any(strfind(taxdata(NEXT).fname,','))
                eddig=strfind(taxdata(NEXT).fname,',');
                eddig=eddig(1)-1;
                taxdata(NEXT).fname=taxdata(NEXT).fname(1:eddig);
            end
        end
    end
    outmatrix=[{taxdata.anatgroup}',{taxdata.ID}',{taxdata.fname}',{taxdata.g}',{taxdata.s}',{taxdata.c}'];
    delete([dirs.taxonomydir,'/taxonomydata.xls']);
    xlwrite([dirs.taxonomydir,'/taxonomydata.xls'],outmatrix,'Sheet 1',['A1']);
end

%%
if strcmp(amplifier,'AXON')
    aE_exportAXONdata(dirs,xlsdata,projectdata.owexport)
else
    %% export Raw data from HEKA
    aE_exportrawHEKAdata(dirs,xlsdata,projectdata.owexport)
    %% generate PGF data, bridge balancing
    aE_generatePGF_bridge_HEKA(dirs,xlsdata,projectdata.owbridge)
end
%% exporting breathing
if isfield(dirs,'breathingdir') & isfield(xlsdata,'Thermosensor_channel')
    for xlsidx=1:length(xlsdata)
        a=dir([dirs.breathingdir,xlsdata(xlsidx).ID,'.mat']);
        if xlsdata(xlsidx).field == 1 & ~isnan(xlsdata(xlsidx).Thermosensor_channel) & (isempty(a) | overwrite==1)
            rawdata=HEKAexportbytime_main(xlsdata(xlsidx).HEKAfname,xlsdata(xlsidx).setup,xlsdata(xlsidx).Thermosensor_channel,xlsdata(xlsidx).startT,xlsdata(xlsidx).endT);
            save([dirs.breathingdir,xlsdata(xlsidx).ID],'rawdata','xlsdata','xlsidx','-v7.3')
            disp([xlsdata(xlsidx).ID,'breathing exported'])
        end
    end
end

%% finding events
valtozok.overwrite=projectdata.owevent;
valtozok.overwritebefore=datenum(2017,12,29);%;datenum(datetime('today'));%;
valtozok.plotit=0;
valtozok.threshholdaveragetime=15;%s
% valtozok.mindvpdt=1;
valtozok.minampl=.0001;
valtozok.apampl=.005;
% valtozok.maxfilteredhw=.01;
valtozok.maxdecayriseratio=5;
valtozok.maxrisetime=.01;
valtozok.minaphw=.0001;
valtozok.maxapwidth=.005;
valtozok.filtermovingtime=.0025;
valtozok.diffmovingt=.0005;
valtozok.steptime=.0005; %s
valtozok.eventminsdval=3;
valtozok.apthreshval=10;
parallelcount=4;
aE_findevents(valtozok,dirs,parallelcount,xlsdata)
paralleldata.count=NaN;
while isnan(paralleldata.count) | ~isempty(paralleldata.files)
    paralleldata.files=dir(dirs.eventparaleldir);
    paralleldata.files([paralleldata.files.isdir])=[];
    pause(3)
    paralleldata.prevcount=paralleldata.count;
    paralleldata.count=length(paralleldata.files);
    if paralleldata.prevcount~=paralleldata.count
        disp(['waiting for ',num2str(paralleldata.count),' eventfinding scripts to finish'])
    end
end
% %% Breathing-field coherence
% destdir=[dirs.basedir,'Coherence/breathing-field/'];
% for xlsi=length(xlsdata):-1:1
%     if xlsdata(xlsi).field==1 
%         fname=[xlsdata(xlsi).HEKAfname,'_',xlsdata(xlsi).G_S_C];
%         a=dir([destdir,fname,'.mat']);
%         if (isempty(a) | valtozok.overwrite==1) & ~isnan(xlsdata(xlsi).Thermosensor_channel)
%             sourcedir = dirs.breathingdir;
%             variablename = 'rawdata';
%             temp=load([sourcedir,fname],variablename);
%             breathing=temp.(variablename);
%             sourcedir =dirs.bridgeddir;
%             variablename = 'bridgeddata';
%             temp=load([sourcedir,fname],variablename);
%             field=temp.(variablename);
%             clear temp
%             %%
%             coherencedata=struct;
%             for fieldsweep=1:length(field)
%                 %%
%                realtime=field(fieldsweep).realtime;
%                si=field(fieldsweep).si;
%                breathingsweep=find([breathing.realtime]==realtime);
%                sweeplength=length(field(fieldsweep).y)*si;
%                if ~isempty(breathingsweep) & sweeplength > 3
%                    time=[1:length(field(fieldsweep).y)]*si;
%                    y_field=field(fieldsweep).y;
%                    y_breath=breathing(breathingsweep).y;
%                    %%
%                    fs=1/si;
%                    f=[.5:.1:10.5];
%                    window=hanning(round(4/si));
%                    noverlap=(round(3.6/si));
%                    [Cxy,F] = mscohere(y_field,y_breath,window,noverlap,f,fs);
%                    pause
%                end
%             end
%         end
%     end
% end


%% extracting PSD

valtozok.overwrite=0;
valtozok.analyseonlyfield=0;
valtozok.downsamplenum='auto';

parameters=struct;
parameters.min=.5; %minimal frequency for decomposition
parameters.max=10;% - maximal frequency for decomposition
parameters.step=.05; %- frequency stepsize
parameters.scale=1;% - scale type, can be linear (1) or logarithmic (2)
parameters.wavenumber=9;% - number of waves in wavelet
parameters.waveletlength=30;
parameters.addtaper=1;
parameters.taperlength=30;

valtozok.parameters=parameters;

aE_PSD_export(dirs,xlsdata,valtozok);
 


%% extracting PSD - high frequency
if isfield(dirs,'PSDdir_high')
    valtozok=struct;
    valtozok.overwrite=0;
    valtozok.analyseonlyfield=1;
    valtozok.downsamplenum='auto';
    valtozok.high=1;
    valtozok.onlyawake=1;
    parameters=struct;
    parameters.min=.5; %minimal frequency for decomposition
    parameters.max=220;% - maximal frequency for decomposition
    parameters.step=1; %- frequency stepsize
    parameters.scale=1;% - scale type, can be linear (1) or logarithmic (2)
    parameters.wavenumber=9;% - number of waves in wavelet
    parameters.waveletlength=30;
    parameters.addtaper=1;
    parameters.taperlength=30;
    
    valtozok.parameters=parameters;
    aE_PSD_export(dirs,xlsdata,valtozok);
end

%% extracting PSD - log
if isfield(dirs,'PSDdir_log')
    valtozok=struct;
    valtozok.overwrite=0;
    valtozok.analyseonlyfield=1;
    valtozok.downsamplenum='auto';
    valtozok.log=1;
    valtozok.onlyawake=1;
    parameters=struct;
    parameters.min=.5; %minimal frequency for decomposition
    parameters.max=200;% - maximal frequency for decomposition
    parameters.step=1; %- frequency stepsize
    parameters.scale=2;% - scale type, can be linear (1) or logarithmic (2)
    parameters.wavenumber=9;% - number of waves in wavelet
    parameters.waveletlength=30;
    parameters.addtaper=1;
    parameters.taperlength=30;
    
    valtozok.parameters=parameters;
    aE_PSD_export(dirs,xlsdata,valtozok);
end
%% extracting PSD for breathing
if isfield(dirs,'PSDdir_breathing')
    valtozok=struct;
    valtozok.breathing=1;
    valtozok.overwrite=0;
    valtozok.analyseonlyfield=1;
    valtozok.downsamplenum='auto';
    
    parameters=struct;
    parameters.min=.5; %minimal frequency for decomposition
    parameters.max=10;% - maximal frequency for decomposition
    parameters.step=.05; %- frequency stepsize
    parameters.scale=1;% - scale type, can be linear (1) or logarithmic (2)
    parameters.wavenumber=9;% - number of waves in wavelet
    parameters.waveletlength=30;
    parameters.addtaper=1;
    parameters.taperlength=30;
    
    valtozok.parameters=parameters;
    
    aE_PSD_export(dirs,xlsdata,valtozok);
    
    %% extracting wavelet cross spectrum between breathing and field
    valtozok=struct;
    valtozok.analyseonlyfield=1;
    valtozok.downsamplenum='auto';
    valtozok.overwrite=1;
    parameters=struct;
    parameters.min=.5; %minimal frequency for decomposition
    parameters.max=10;% - maximal frequency for decomposition
    parameters.step=.05; %- frequency stepsize
    parameters.scale=1;% - scale type, can be linear (1) or logarithmic (2)
    parameters.wavenumber=9;% - number of waves in wavelet
    parameters.waveletlength=30;
    parameters.addtaper=1;
    parameters.taperlength=30;
    
    valtozok.parameters=parameters;
    
    aE_cross_spectrum_export(dirs,xlsdata,valtozok);
end

%% extracting movement info from video

valtozok.sampleinterval=.01;
valtozok.overwrite=0;
for xlsi=length(xlsdata):-1:1
    if any(strfind(xlsdata(xlsi).anaesthesia,'awake'))
        valtozok.setupname=xlsdata(xlsi).setup;
        valtozok.filename=xlsdata(xlsi).HEKAfname;
        aE_analyzevideo_main(valtozok,dirs);
    end
end

%% PCA on ROIs
overwrite=0;
ROIdir=[dirs.videodir,'ROIs/'];
movementdir=[dirs.videodir,'movement/'];
files=dir(movementdir);
files([files.isdir])=[];
for filei=1:length(files)
    a=dir([ROIdir,files(filei).name]);
    if ~isempty(a)
    load([ROIdir,files(filei).name]);
    %     if true%~isfield(ROIdata,'PCA')
    load([movementdir,files(filei).name]);
    for videoi=1:length(rawvideo)
        %         if videoi==1
        %             bigvideo=rawvideo(videoi).vid;
        %         else
        %             bigvideo=cat(3,bigvideo,rawvideo(videoi).vid);
        %         end
        bigvideo=rawvideo(videoi).vid;
        %%
        for ROIi=1:length(ROIdata)
            if isempty(ROIdata(ROIi).mask)
                disp(['no ROI selected for ROI:',ROIdata(ROIi).ROIname,' in file ',files(filei).name]);
            elseif overwrite==1 | ~isfield(ROIdata,'PCA') | ~isfield(ROIdata(ROIi).PCA,'mask') | length(ROIdata(ROIi).PCA)<videoi | isempty(ROIdata(ROIi).PCA(videoi).mask) | ROIdata(ROIi).PCA(videoi).mask~=ROIdata(ROIi).mask | ~isfield(ROIdata(ROIi).PCA(videoi),'uniquePC')
                ROIdata(ROIi).PCA(videoi).mask=ROIdata(ROIi).mask;
                bigmask=repmat(ROIdata(ROIi).mask,1,1,size(bigvideo,3));
                ROIvideo=bigvideo;
                ROIvideo(~bigmask)=0;
                %                 tic
                linearROIvideo=reshape(ROIvideo,size(ROIvideo,1)*size(ROIvideo,2),size(ROIvideo,3));
                %%
                borders=[];
                whattoprincomp=((double(linearROIvideo')));
                whattoprincomp=diff(zscore(whattoprincomp));
                %%
                % %                 borders(1,:) = prctile(whattoprincomp,2);
                % %                 borders(2,:) = prctile(whattoprincomp,98);
                % %                 borders(3,:) = diff(borders);
                % % %                 borders=borders;
                % %                 whattoprincomp=whattoprincomp-repmat(borders(1,:),size(whattoprincomp,1),1);
                % %                 whattoprincomp=whattoprincomp./repmat(borders(3,:),size(whattoprincomp,1),1);
                % %                 whattoprincomp(whattoprincomp>1)=1;
                % %                 whattoprincomp(whattoprincomp<0)=0;
                % %                 whattoprincomp=diff(whattoprincomp);
                % %                 whattoprincompold=whattoprincomp;
                % %                 rowstodel=find(isnan(whattoprincomp(end,:)));
                
                rowstodel=find(sum(whattoprincomp)==0);
                whattoprincomp(:,rowstodel)=[];
                
                %%
                [coeff,score,latent]=princomp(whattoprincomp);
                ROIdata(ROIi).PCA(videoi).PC=zeros(size(score,1),10);
                eddig=min(10,size(score,2));
                ROIdata(ROIi).PCA(videoi).PC(:,1:eddig)=score(:,1:eddig);
                ROIdata(ROIi).PCA(videoi).latent=latent;
                ROIdata(ROIi).PCA(videoi).cumweigths=cumsum(latent)/sum(latent);
                ROIdata(ROIi).PCA(videoi).weigths=(latent)/sum(latent);
                %                 toc
                %         %%
                %         scorelayout=[];
                %         close all
                %         princompstoshow=600;
                %         figure(10)
                %         clf
                %         for i=1:princompstoshow
                %             %                 subplot(princompstoshow,2,(i-1)*2+1)
                %             subplot(1,2,1)
                %             plot((score(:,i)))
                %             %                 subplot(princompstoshow,2,(i-1)*2+2)
                %             subplot(1,2,2)
                % %             scorelayout(:,:,i)=reshape(coeff(:,i),size(bigvideo,1),size(bigvideo,2));
                % %             imagesc(scorelayout(:,:,i))
                %             %                 caxis([-.04 .04])
                %             pause
                %         end
                %%
                if strcmp(ROIdata(ROIi).ROIname,'Body')
                    powervals=[];
                    for pci=1:10
                        parameters=struct;
                        parameters.min=1; %minimal frequency for decomposition
                        parameters.max=5;% - maximal frequency for decomposition
                        parameters.step=.1; %- frequency stepsize
                        parameters.scale=1;% - scale type, can be linear (1) or logarithmic (2)
                        parameters.wavenumber=7;% - number of waves in wavelet
                        parameters.waveletlength=10;
                        parameters.addtaper=1;
                        parameters.taperlength=5;
                        parameters.interval=mode(diff(videodata(videoi).time));
                        [powerMatrix, frequencyVector]  = calculateTFDPower(ROIdata(ROIi).PCA(videoi).PC(:,pci), parameters);
                        powervals(pci)=median(max(powerMatrix));
%                         figure(1)
%                         clf
%                         subplot(3,1,1)
%                         imagesc(1:length(ROIdata(ROIi).PCA(videoi).PC(:,pci)),frequencyVector,powerMatrix);
%                         set(gca,'YDir','normal');
%                         subplot(3,1,2)
%                         plot(max(powerMatrix));
%                         hold on
%                         plot(max(powerMatrix)*0+median(max(powerMatrix)),'r-')
%                         subplot(3,1,3)
%                         plot(ROIdata(ROIi).PCA(videoi).PC(:,pci))
%                         pause
                    end
                    [~,breathindex]=max(powervals);
                    ROIdata(ROIi).PCA(videoi).breath=ROIdata(ROIi).PCA(videoi).PC(:,breathindex);
                end
            end
        end
        %% PCA based PC selection
        ennyiPCkell=3;
        PCs=struct;
        for ROIi=1:length(ROIdata)
            PCs(ROIi).pc=ROIdata(ROIi).PCA(videoi).PC(:,1:ennyiPCkell);
        end
        whattoprincomp=zscore([PCs.pc]);
        [coeff,score,latent]=princomp(whattoprincomp);
        osszidx=1:length(ROIdata)*ennyiPCkell;
        for ROIi=4:length(ROIdata)
            idxmost=(ROIi-1)*ennyiPCkell+1:ennyiPCkell+(ROIi-1)*ennyiPCkell;
            idxrest=osszidx;
            idxrest(idxmost)=[];
            idxbckground=[-ennyiPCkell+1:0]+ennyiPCkell*length(ROIdata);
            [bestcompcoeff,bestcompidx]=max(abs(coeff(idxmost,:)),[],2);
            
            [correspondingbackgroundcoeff,corrbcgidx]=max(abs(coeff(idxbckground,bestcompidx)));
            [correspondingrestcoeff,~]=max(abs(coeff(idxrest,bestcompidx)));
            
            coeffbackgroundratio=bestcompcoeff./correspondingbackgroundcoeff';
            coeffrestratio=bestcompcoeff./correspondingrestcoeff';
            
            [~,idx]=max(coeffbackgroundratio);
            ROIdata(ROIi).PCA(videoi).bestPC=ROIdata(ROIi).PCA(videoi).PC(:,idx);
            [~,idx]=max(coeffrestratio);
            ROIdata(ROIi).PCA(videoi).uniquePC=ROIdata(ROIi).PCA(videoi).PC(:,idx);
            
%             figure(1)
%             clf
%             for princompi=1:ennyiPCkell
%                 subplot(ennyiPCkell,1,princompi)
%                 plot(whattoprincomp(:,idxmost(princompi)))
%                 hold on
%                 plot(whattoprincomp(:,idxbckground(corrbcgidx(princompi))),'r-')
%             end
            
            %%% itt kiválasztom azt a principális komponenst, ami egy adott
            %%% ROI-ra jellemző, de a háttérre nem
            %             return
        end
        
        %%
        %         for i=1:100
        %
        %                 subplot(2,2,1)
        %                 plot((score(:,i)))
        %                xlabel('Frames')
        %                 ylabel('master PCA score')
        %
        %                 subplot(2,2,2)
        %                 plot([1:ennyiPCkell*length(ROIdata)]/ennyiPCkell+1,(coeff(:,i)))
        %                 xlabel('ROIs')
        %                 ylabel('PCA coefficient')
        %
        %                 [~,idx]=max(abs(coeff(:,i)));
        %                 subplot(2,2,3)
        %                 plot(whattoprincomp(:,idx));
        %
        %                 xlabel('Frames')
        %                 ylabel('original PCA score')
        %                 subplot(2,2,4)
        %                 imagesc(abs(coeff'))
        %                 set(gca,'YDir','normal')
        %                 pause
        %             end
    end
    disp(['PCA for ', files(filei).name, 'is done'])
    save([ROIdir,files(filei).name],'ROIdata');
    else
        disp(['no ROIS selected for ', files(filei).name])
    end 
end
    %% processing Gaspar's pupil diameter data

savedir=[dirs.videodir,'eye/'];
for xlsi=1:length(xlsdata)
    setupname=xlsdata(xlsi).setup;
    filename=xlsdata(xlsi).HEKAfname(1:6);
    pupilsizedir=[locations.tgtardir,'VIDEOdata/',setupname,'/pupilsize/'];
    files=dir(pupilsizedir);
    files([files.isdir])=[];
    a=dir([savedir,xlsdata(xlsi).HEKAfname,'.mat']);
    if isempty(a)
        if length(files)>0
            for i=1:length(files)
                hyps=strfind(files(i).name,'_');
                files(i).filename=[files(i).name(hyps(1)-2:hyps(1)-1),files(i).name(hyps(1)+1:hyps(2)-1),files(i).name(hyps(2)+1:hyps(3)-1)];
                files(i).timestamp=[files(i).name(1:hyps(end)+2)];
                files(i).ext=files(i).name(end-2:end);
            end
            neededfiles=find(strcmp({files.filename},filename));
            if length(neededfiles)>0
                time=[];
                diameter=[];
                %             flash=[];
                for i=1:length(neededfiles)
                    %                 %%old method
                    %                 txt=load([pupilsizedir,files(neededfiles(i)).name]);
                    %                  time=[time;txt(:,3)];
                    %                 diameter=[diameter;txt(:,4)];
                    %                 flash=[flash;txt(:,5)];
                    %% new method - string read
                    fid=fopen([pupilsizedir,files(neededfiles(i)).name],'r');
                    txt=textscan(fid,'%s');
                    txt=txt{1};
                    fclose(fid);
                    dottedtimeidx = find(~cellfun(@isempty,regexp(txt,':..:')));
                    timenow=txt(dottedtimeidx+1);
                    diameternow=txt(dottedtimeidx+2);
                    timenow_num=zeros(size(timenow));
                    diameternow_num=zeros(size(diameternow));
                    for idxi=1:length(timenow)
                        if any(strfind(timenow{idxi},','))
                            timenow{idxi}(strfind(timenow{idxi},','))='.';
                        end
                        timenow_num(idxi)=str2num(timenow{idxi});
                        if any(strfind(diameternow{idxi},','))
                            diameternow{idxi}(strfind(diameternow{idxi},','))='.';
                        end
                        diameternow_num(idxi)=str2num(diameternow{idxi});
                    end
                    
                    time=[time;timenow_num];
                    diameter=[diameter;diameternow_num];
                    %%
                    
                end
                [time,ix]=sort(time);
                diameter=diameter(ix);
                %             flash=flash(ix);
                pupildata.diameter=diameter;
                pupildata.time=time;
                %             pupildata.flash=flash;

                save([savedir,xlsdata(xlsi).HEKAfname],'pupildata');
                disp([xlsdata(xlsi).HEKAfname,' - pupil diameter extraction is done'])
            end
        end
    end
end

%% normalizing pupil size and movement data
videopercentiles=struct;
alldata=struct;
files=dir([dirs.videodir,'eye/']);
files([files.isdir])=[];
for filei=1:length(files)
    
    load([dirs.videodir,'ROIs/',files(filei).name]);
    neededroi=find(strcmp({ROIdata.ROIname},'Eye'));
    stats=regionprops(ROIdata(neededroi).mask,'MajorAxisLength','MinorAxisLength');
    alldata(filei).MajorAxisLength=stats.MajorAxisLength;
    alldata(filei).MinorAxisLength=stats.MinorAxisLength;
    
    load([dirs.videodir,'eye/',files(filei).name]);
    load([dirs.videodir,'movement/',files(filei).name],'videodata');
    xlsidx=find(strcmp({xlsdata.HEKAfname},files(filei).name(1:end-4)),1,'first');
    
    moviename=[videodata(end).timestamp,'.avi'];
    setupname=xlsdata(xlsidx).setup;
    locations=marcicucca_locations;
    moviefiletoplay=[locations.tgtardir,'VIDEOdata/',setupname,'/',moviename];
    movieobj=VideoReader([moviefiletoplay]);
    szorzo=movieobj.Height/size(videodata(1).originalpic_all,1);
    %%
    
    alldata(filei).pupildiameter=pupildata.diameter'/alldata(filei).MinorAxisLength/szorzo*2;
    pupildiametersnow=sort(alldata(filei).pupildiameter);
    percentile95=pupildiametersnow(round(.95*length(pupildiametersnow)));
    if percentile95>1
        alldata(filei).pupildiameter=alldata(filei).pupildiameter/percentile95;
    end
    alldata(filei).pupildiameter_raw=pupildata.diameter';
%     alldata(filei).pupildiameter=pupildata.diameter';
    alldata(filei).pupildiameter=moving(alldata(filei).pupildiameter,2)';
    alldata(filei).pupildiameter_raw=moving(alldata(filei).pupildiameter_raw,2)';
    
    %%
    
    
end
sorteddiameter=sort([alldata.pupildiameter]);
pupilpercentiles=zeros(100,1);
for i=1:100
    pupilpercentiles(i)=sorteddiameter(round(length(sorteddiameter)*(i/100)));
end
videopercentiles.pupilpercentiles=pupilpercentiles;

allmovdata=struct;
files=dir([dirs.videodir,'movement/']);
files([files.isdir])=[];
for filei=1:length(files)
    load([dirs.videodir,'movement/',files(filei).name],'videodata');
    for fieldi=1:length(videodata);
        if isempty(fieldnames(allmovdata))
            NEXT=1;
        else
            NEXT=length(allmovdata)+1;
        end
        if length(videodata(fieldi).movement_all)>2
            allmovdata(NEXT).movementall=moving(videodata(fieldi).movement_all,2)';
        end
    end
end
sortedmovdata=sort([allmovdata.movementall]);
movementpercentiles=zeros(100,1);
for i=1:100
    movementpercentiles(i)=sortedmovdata(round(length(sortedmovdata)*(i/100)));
end
videopercentiles.movementpercentiles=movementpercentiles;
save([dirs.videodir,'percentiles'],'videopercentiles');
return

%% looking for slow oscillations on the field
window=5;%s
windowstep=2.5;
for xlsidx=length(xlsdata):-1:1
    ID=xlsdata(xlsidx).ID;
    load([dirs.PSDdir,ID])
    peakdata=struct;
    for sweep=1:length(PSDdata)
        sweephossz=length(PSDdata(sweep).y)*PSDdata(sweep).si_powerMatrix;
        if sweephossz>=window
            time=[1:length(PSDdata(sweep).y)]*PSDdata(sweep).si_powerMatrix-PSDdata(sweep).si_powerMatrix;
            PSDdata(sweep).powerMatrix=(double(PSDdata(sweep).powerMatrix)-PSDdata(sweep).compress_offset)*PSDdata(sweep).compress_multiplier;
            PSDmedian=zeros(size(PSDdata(sweep).powerMatrix,1),round(sweephossz/windowstep)-1);
            PSDmediantime=1:size(PSDmedian,2)*windowstep;
            
            for wini=1:round(sweephossz/windowstep)-1
                [~,ettol]=min(abs(time-((wini)*windowstep-window/2)));
                [~,eddig]=min(abs(time-((wini)*windowstep+window/2)));
                PSDmedian(:,wini)=nanmedian(PSDdata(sweep).powerMatrix(:,ettol:eddig),2);
                [pks,locs,w,p]=findpeaks(PSDmedian(:,wini));
                [p,idx]=sort(p,'descend');
                pks=pks(idx);
                locs=locs(idx);
                w=w(idx);
                idxnow=locs(1:3);
                
                
                if isempty(fieldnames(peakdata))
                    next=1;
                else
                    next=length(peakdata)+1;
                end
                peakdata(next).peakval=pks(1:3)';
                peakdata(next).peakwidth=w(1:3)';
                peakdata(next).prominence=p(1:3)';
                peakdata(next).freq=PSDdata(sweep).frequencyVector(locs(1:3));
                
%                 figure(2)
%                 clf
%                 plot(PSDmedian(:,wini));
%                 hold on
%                 plot(idxnow,PSDmedian(idxnow,wini),'ro')
%                 pause
            end
            figure(1)
            clf
            subplot(3,2,1)
            plot(time,PSDdata(sweep).y);
            subplot(3,2,2)
            plot(PSDdata(sweep).frequencyVector,PSDmedian);
            subplot(3,2,3)
            imagesc(time,PSDdata(sweep).frequencyVector,PSDdata(sweep).powerMatrix);
            set(gca,'YDir','normal');
            colormap linspecer
            ylabel('Frequency (Hz)')
            xlabel('Time (s)')
            subplot(3,2,4)
            imagesc(PSDmediantime,PSDdata(sweep).frequencyVector,PSDmedian);
            set(gca,'YDir','normal');
            colormap linspecer
            ylabel('Frequency (Hz)')
            xlabel('Time (s)')
            
            subplot(3,2,5)
            semilogy([peakdata.freq],[peakdata.prominence],'ko')
            
            pause
        end
    end
    
end








%%

valtozok=struct;
valtozok.movingvindowsize=3;
valtozok.movingvindowstep=.5; %seconds for downsampling and median filtering
valtozok.timeborders=[0 inf];
valtozok.frequencyrange=[1 4];
valtozok.PSDonfield=true;
valtozok.minsweeplength=valtozok.movingvindowsize/2;
xlsnum=find(strcmp({xlsdata.ID},'1705301rm_1_2_1'));

aE_V0_vs_PSD(dirs,xlsdata,xlsnum,valtozok)

%
%% UP-DOWN transitions
if projectnum==2
    valtozok=struct;
    valtozok.overwrite=0;
    valtozok.segmentlength=5;
    valtozok.plotthestuff=0;
    aE_UP_DOWN_detect(valtozok,dirs,xlsdata);
    
    valtozok=struct;
    valtozok.overwrite=1;
    valtozok.segmentlength=5;
    valtozok.plotthestuff=1;
    aE_UP_DOWN_detect_field(valtozok,dirs,xlsdata);
end

%% plot state transitions
timeborders=[69070, 69510];
timeborders=[0, 76500];
timeborders=[0, inf];
minimumstateduration=0;
[Selection,ok] = listdlg('ListString',{xlsdata.ID},'ListSize',[300 600]);
ID=xlsdata(Selection).ID;
load([dirs.eventdir,ID],'eventdata');
load([dirs.bridgeddir,ID],'bridgeddata','stimdata');
load([dirs.statedir,ID],'statedata');
statedata.UP=statedata.UP([statedata.UP.onsett]>timeborders(1)&[statedata.UP.onsett]<timeborders(2));
statedata.DOWN=statedata.DOWN([statedata.DOWN.onsett]>timeborders(1)&[statedata.DOWN.onsett]<timeborders(2));
% get transition data
timebefore=1;
timeafter=1;
apdata=eventdata(strcmp({eventdata.type},'AP'));
prevsweepnum=0;
for api=1:length(apdata)
    sweepnum=apdata(api).sweepnum;
    if sweepnum~=prevsweepnum
        y=bridgeddata(sweepnum).y;
        si=bridgeddata(sweepnum).si;
        yfiltered=moving(y,5);
        dyfiltered=diff(yfiltered)/si;
        prevsweepnum=sweepnum;
        stepback=round(.0002/si);
    end
    onseth=apdata(api).onseth;
    maxh=apdata(api).maxh;
    threshh=maxh-stepback;
    while dyfiltered(threshh)>5 & threshh>2
        threshh=threshh-1;
    end
    threshv=yfiltered(threshh);
    apdata(api).threshv=threshv;
    apdata(api).APamplitude=apdata(api).maxval-apdata(api).threshv;
    %     if threshv>-.04% & threshv<-.04
    %         figure(1)
    %         clf
    %         hold on
    %         plot(yfiltered,'r-')
    %         plot(y,'k-')
    %
    %         plot(maxh,yfiltered(maxh),'ro')
    %         plot(threshh,yfiltered(threshh),'ko')
    %         plot(onseth,yfiltered(onseth),'kx')
    %         pause
    %     end
end
figure(1)
clf
% hist([apdata.threshv],100)
plot([apdata.threshv],[apdata.APamplitude],'ko')
[x,~]=ginput(1);
aapdata=apdata([apdata.threshv]<x);
sapdata=apdata([apdata.threshv]>x);
epdata=eventdata(strcmp({eventdata.type},'ep'));
ipdata=eventdata(strcmp({eventdata.type},'ip'));
statestocheck={'UP','DOWN'};
%
for statei=1:length(statestocheck)
    statedatanow=statedata.(statestocheck{statei});
    statedatanow=statedatanow([statedatanow.duration]>minimumstateduration);
    transitiondata=struct;
    for i=1:length(statedatanow)
        sweepnum=statedatanow(i).sweepnum;
        transitiont=statedatanow(i).onsett;
        transitionh=statedatanow(i).onseth;
        si=bridgeddata(sweepnum).si;
        stepback=round(timebefore/si);
        stepforward=round(timeafter/si);
        y=bridgeddata(sweepnum).y;
        if transitionh>stepback & length(y)>transitionh+stepforward
            if isempty(fieldnames(transitiondata))
                NEXT=1;%% recording statistics
                baselineSDbinsize=5;
                baselineSDbinnum=round(30*60/baselineSDbinsize);
                recstats=struct;
                files={xlsdata.ID};
                for filei=1:length(files)
                    if xlsdata(filei).field==0
                        if isempty(fieldnames(recstats))
                            NEXT=1;
                        else
                            NEXT=length(recstats)+1;
                        end
                        ID=files{filei};
                        load([dirs.bridgeddir,ID],'bridgeddata','lightdata','stimdata');
                        recstats(NEXT).ID=ID;
                        recstats(NEXT).RS=nanmedian([lightdata.RS]);
                        recstats(NEXT).anaesth=xlsdata(filei).anaesthesia;
                        if lightdata(end).realtime>lightdata(1).realtime;
                            recstats(NEXT).recordinglength=lightdata(end).realtime-lightdata(1).realtime;
                        else
                            recstats(NEXT).recordinglength=lightdata(end).realtime-lightdata(1).realtime+24*3600;
                        end
                        recstats(NEXT).recordinglength= recstats(NEXT).recordinglength+bridgeddata(end).si*length(bridgeddata(end).y);
                        %baselineSD statistics
                        baselineSD=NaN(baselineSDbinnum,1);
                        if recstats(NEXT).recordinglength>=20*60
                            for bini=1:baselineSDbinnum
                                startt=(bini-1)*baselineSDbinsize+lightdata(1).realtime;
                                endt=(bini)*baselineSDbinsize+lightdata(1).realtime;
                                neededsweepnum=find([lightdata.realtime]>=startt & [lightdata.realtime]<=endt);
                                sweepSD=NaN;
                                if ~isempty(neededsweepnum)
                                    if neededsweepnum(1)>1
                                        neededsweepnum=[neededsweepnum(1)-1,neededsweepnum];
                                    end
                                    sweepSD=nan(size(neededsweepnum));
                                    for sweepi=1:length(neededsweepnum)
                                        sweepnum=neededsweepnum(sweepi);
                                        if strcmp(stimdata(sweepnum).Amplifiermode,'C-Clamp') &  any(strfind(bridgeddata(sweepnum).channellabel,'Vmon')) & std(stimdata(sweepnum).y)==0
                                            [b,a]=butter(1,500/(1/bridgeddata(sweepnum).si)/2,'low');
                                            y=bridgeddata(sweepnum).y;
                                            y=filtfilt(b,a,y);
                                            x=[1:length(bridgeddata(sweepnum).y)]*bridgeddata(sweepnum).si+bridgeddata(sweepnum).realtime-bridgeddata(sweepnum).si;
                                            neededidx=find(x>=startt & x<=endt);
                                            sweepSD(sweepi)=nanstd(y(neededidx));
                                        end
                                    end
                                end
                                baselineSD(bini)=nanmedian(sweepSD);
                            end
                            %             figure(2)
                            %             clf
                            %             plot([1:baselineSDbinnum]*baselineSDbinsize/60,baselineSD,'ko-')
                            %             pause
                        end
                        recstats(NEXT).baselineSD=baselineSD;
                    end
                end
                
                baselineSDtimevector=[1:baselineSDbinnum]*baselineSDbinsize/60;
                newbaselineSDtimevector=[1:1:30];
                newbaselineSDtimevectorcenters=newbaselineSDtimevector-mode(diff(newbaselineSDtimevector));
                newbaselineSDtimevectorsteps=mode(diff(newbaselineSDtimevector));
                for i=1:length(recstats)
                    recstats(i).newbaselineSD=[];
                    for stepi=1:length(newbaselineSDtimevectorcenters)
                        idx=find(baselineSDtimevector>=newbaselineSDtimevectorcenters(stepi)-newbaselineSDtimevectorsteps/2 & baselineSDtimevector<=newbaselineSDtimevectorcenters(stepi)+newbaselineSDtimevectorsteps/2);
                        recstats(i).newbaselineSD(stepi)=nanmean(recstats(i).baselineSD(idx));
                    end
                    recstats(i).newbaselineSD=(recstats(i).newbaselineSD/nanmedian(recstats(i).newbaselineSD(find(newbaselineSDtimevector>15))))';
                end
                figure(1)
                clf
                subplot(3,1,1)
                [nall,xbin]=hist([recstats.recordinglength]/60,[2.5:5:62.5]);
                [nketamine,~]=hist([recstats(strcmp({recstats.anaesth},'ketamine xylazine')).recordinglength]/60,[2.5:5:62.5]);
                [nchloral,~]=hist([recstats(strcmp({recstats.anaesth},'chloral hydrate')).recordinglength]/60,[2.5:5:62.5]);
                [hAxes,hBar,hLine] = plotyy(xbin,[nketamine;nchloral]',xbin,cumsum([nketamine;nchloral]'),'bar','plot');
                set(hBar,'BarLayout','stacked')
                colororder=get(hAxes,'colororder');
                set(hBar(2),'FaceColor',colororder{1}(2,:));
                set(hLine,'LineWidth',2)
                legend({'ketamine xylazine','chloral hydrate'})
                set(gca,'Xtick',[0:5:60])
                xlabel('min')
                ylabel('# of cells')
                title('recording length')
                subplot(3,1,2)
                [nall,xbin]=hist([recstats.RS]/10^6,[5:10:95]);
                [nketamine,~]=hist([recstats(strcmp({recstats.anaesth},'ketamine xylazine')).RS]/10^6,[5:10:95]);
                [nchloral,~]=hist([recstats(strcmp({recstats.anaesth},'chloral hydrate')).RS]/10^6,[5:10:95]);
                [hAxes,hBar,hLine] = plotyy(xbin,[nketamine;nchloral]',xbin,cumsum([nketamine;nchloral]'),'bar','plot');
                set(hBar,'BarLayout','stacked')
                colororder=get(hAxes,'colororder');
                set(hBar(2),'FaceColor',colororder{1}(2,:));
                set(hLine,'LineWidth',2)
                legend({'ketamine xylazine','chloral hydrate'})
                set(gca,'Xtick',[0:10:90])
                xlabel('MOhm')
                ylabel('# of cells')
                title('median RS')
                subplot(3,1,3)
                hold on
                for i=1:length(recstats)
                    needed=find(~isnan(recstats(i).newbaselineSD));
                    if strcmp(recstats(i).anaesth,'ketamine xylazine')
                        colorka=colororder{1}(1,:);
                    elseif strcmp(recstats(i).anaesth,'chloral hydrate')
                        colorka=colororder{1}(2,:);
                    end
                    if ~isempty(needed)
                        plot(newbaselineSDtimevector(needed),recstats(i).newbaselineSD(needed),'Color',colorka,'LineWidth',2);
                    end
                end
                ylim([0 1.2])
                xlabel('time (min)')
                ylabel('relative baseline SD')
            else
                NEXT=length(transitiondata)+1;
            end
            neededAPs=find([apdata.maxtime]>transitiont-timebefore &[apdata.maxtime]<transitiont+timeafter);
            neededaAPs=find([aapdata.maxtime]>transitiont-timebefore &[aapdata.maxtime]<transitiont+timeafter);
            neededsAPs=find([sapdata.maxtime]>transitiont-timebefore &[sapdata.maxtime]<transitiont+timeafter);
            neededeps=find([epdata.maxtime]>transitiont-timebefore &[epdata.maxtime]<transitiont+timeafter);
            neededips=find([ipdata.maxtime]>transitiont-timebefore &[ipdata.maxtime]<transitiont+timeafter);
            transitiondata(NEXT).transitiont=transitiont;
            transitiondata(NEXT).time=[-stepback:stepforward]'*si;
            transitiondata(NEXT).y=y(transitionh-stepback:transitionh+stepforward)';
            
            if ~isempty(neededAPs)
                transitiondata(NEXT).APtimes=[apdata(neededAPs).maxtime]-transitiont;
                transitiondata(NEXT).APnum=length(transitiondata(NEXT).APtimes);
            else
                transitiondata(NEXT).APtimes=[];
                transitiondata(NEXT).APnum=0;
            end
            if ~isempty(neededaAPs)
                transitiondata(NEXT).aAPtimes=[aapdata(neededaAPs).maxtime]-transitiont;
                transitiondata(NEXT).aAPnum=length(transitiondata(NEXT).aAPtimes);
            else
                transitiondata(NEXT).aAPtimes=[];
                transitiondata(NEXT).aAPnum=0;
            end
            if ~isempty(neededsAPs)
                transitiondata(NEXT).sAPtimes=[sapdata(neededsAPs).maxtime]-transitiont;
                transitiondata(NEXT).sAPnum=length(transitiondata(NEXT).sAPtimes);
            else
                transitiondata(NEXT).sAPtimes=[];
                transitiondata(NEXT).sAPnum=0;
            end
            if ~isempty(neededeps)
                transitiondata(NEXT).eptimes=[epdata(neededeps).maxtime]-transitiont;
                transitiondata(NEXT).epnum=length(transitiondata(NEXT).eptimes);
            else
                transitiondata(NEXT).eptimes=[];
                transitiondata(NEXT).epnum=0;
            end
            if ~isempty(neededips)
                transitiondata(NEXT).iptimes=[ipdata(neededips).maxtime]-transitiont;
                transitiondata(NEXT).ipnum=length(transitiondata(NEXT).iptimes);
            else
                transitiondata(NEXT).iptimes=[];
                transitiondata(NEXT).ipnum=0;
            end
        end
    end
    Transitiondata.(statestocheck{statei})=transitiondata;
end
% plot
figure(1)
clf
subplot(4,2,1)
needed=find([Transitiondata.UP.aAPnum]>0);
plot([Transitiondata.UP(needed).time],[Transitiondata.UP(needed).y]);
xlim([-timebefore,timeafter])
title('DOWN to UP transition')
subplot(4,2,2)
needed=find([Transitiondata.DOWN.aAPnum]>0);
plot([Transitiondata.DOWN(needed).time],[Transitiondata.DOWN(needed).y]);
xlim([-timebefore,timeafter])
title('UP to DOWN transition')
subplot(4,2,3)
hist([Transitiondata.UP.sAPtimes],[-timebefore:0.05:timeafter])
xlim([-timebefore,timeafter])
ylabel('somatic APnum');
subplot(4,2,4)
hist([Transitiondata.DOWN.sAPtimes],[-timebefore:0.05:timeafter])
xlim([-timebefore,timeafter])
ylabel('somatic AP num');
subplot(4,2,5)
hist([Transitiondata.UP.aAPtimes],[-timebefore:0.05:timeafter])
xlim([-timebefore,timeafter])
ylabel('axonal APnum');
subplot(4,2,6)
hist([Transitiondata.DOWN.aAPtimes],[-timebefore:0.05:timeafter])
xlim([-timebefore,timeafter])
ylabel('axonal AP num');
subplot(4,2,7)
hist([Transitiondata.UP.eptimes],[-timebefore:0.05:timeafter])
xlim([-timebefore,timeafter])
ylabel('EPnum');
subplot(4,2,8)
hist([Transitiondata.DOWN.eptimes],[-timebefore:0.05:timeafter])
xlim([-timebefore,timeafter])
ylabel('EPnum');


%%
apersratio=[APdata.aAPnum]./[APdata.sAPnum];
[~,idx]=sort(apersratio,'descend');
APdata=APdata(idx);
% for filei=1:length(files)
%     fname=files(filei).name(1:end-4);
%     load([sortedeventdir,fname],'eventdata');
%     xlsidx=find(strcmp({xlsdata.ID},fname));
%     
%     for xlsfieldi=1:length(fieldek)
%         APdata(filei).(fieldek{xlsfieldi})=xlsdata(xlsidx).(fieldek{xlsfieldi});
%     end
% %     a=dir()
%     APdata(filei).aAPnum=sum([eventdata.axonalAP]);
%     APdata(filei).sAPnum=sum([eventdata.somaticAP]);
% end



%%
if projectnum==4
    %% defining stimepochs and spike clusters
    valtozok_stimepochs.overwrite=projectdata.owstimepoch;
    valtozok_stimepochs.histbins=[0:.001:.5];
    valtozok_stimepochs.sdtimesval_persistentgroup=4;
    valtozok_stimepochs.maxskewness_persistentgroup=2;
    valtozok_stimepochs.maxisiwithingroup=.1;
    valtozok_stimepochs.minapnum_persistentgroup=1;
    valtozok_stimepochs.minAPtimediff_stimepoch=10;
    valtozok_stimepochs.minAPnum_stimepoch=20;
    valtozok_stimepochs.plotpersistentgroups=0;
    persistent_definestimepoch(dirs,valtozok_stimepochs)
    
    %% determining TAU and V0
    valtozok_tau.timeconstanthossz=.100; %s
    valtozok_tau.stepstocheck=20;%*2
    valtozok_states.overwrite=0;
    valtozok_states.ordfiltorder=.2;
    valtozok_states.ordfiltlength=.02; %s
    files=dir(dirs.bridgeddir);
    files([files.isdir])=[];
    for fnum=1:length(files)
        progressbar(fnum/length(files),[],[])
        fname=files(fnum).name;
        a=dir([dirs.v0distdir,fname]);
        if isempty(a) | valtozok_states.overwrite==1
            taudata=persistent_gettau(valtozok_tau,dirs,fname);
            persistent_getv0dist(dirs,valtozok_states,valtozok_tau,taudata,fname);
        end
    end
end
if  projectnum==5
    %% checking somatically triggered axonal events
    timebefore=.003;
    timeafter=.005;
    cutofffreq=[3000];
    filterorder=1;
    for i=1:length(xlsdata)
        if xlsdata(i).startT>xlsdata(i).endT
            xlsdata(i).endTmodified=xlsdata(i).endT+86400;
        else
            xlsdata(i).endTmodified=xlsdata(i).endT;
        end
    end
    blebidxes=find([xlsdata.axonal]);
    fieldidxes=find([xlsdata.field]);
    icdxes=find(~[xlsdata.field]&~[xlsdata.axonal]);
    [icSelection,ok] = listdlg('ListString',{xlsdata(icdxes).ID},'ListSize',[300 600]); % az XLS file alapján kiválasztjuk, hogy melyik file összes mérésén szeretnénk végigmenni
    for ici=1:length(icSelection)
        icxlsidx=icdxes(icSelection(ici));
        blebSelection=zeros(size(blebidxes));
        for blebi=1:length(blebidxes)
            blebidx=blebidxes(blebi);
            if strcmp(xlsdata(icxlsidx).HEKAfname,xlsdata(blebidx).HEKAfname) & xlsdata(blebidx).startT<xlsdata(icxlsidx).endTmodified & xlsdata(icxlsidx).startT<xlsdata(blebidx).endTmodified
                blebSelection(blebi)=1;
            end
        end
        blebSelection=blebidxes(find(blebSelection));
        for blebi=1:length(blebSelection)
            blebxlsidx=blebSelection(blebi);
            disp(xlsdata(blebxlsidx).ID)
            ic=struct;
            bleb=struct;
            ic=load([dirs.bridgeddir,xlsdata(icxlsidx).ID]);
            load([dirs.eventdir,xlsdata(icxlsidx).ID],'eventdata');
            ic.eventdata=eventdata;
            bleb=load([dirs.bridgeddir,xlsdata(blebxlsidx).ID]);
            for blebsweepnum=1:length(bleb.bridgeddata) %filtering
                if strcmp(xlsdata(blebxlsidx).ID,'1702083rm_2_20_1')
                    bleb.bridgeddata(blebsweepnum).realtime=bleb.bridgeddata(blebsweepnum).realtime+24*3600;
                end
                si=bleb.bridgeddata(blebsweepnum).si;
                ninq=.5/si;
                if length(cutofffreq)==1
                    [b,a] = butter(filterorder,cutofffreq/ninq,'low');
                else
                    [b,a] = butter(filterorder,cutofffreq/ninq,'bandpass');
                end
                bleb.bridgeddata(blebsweepnum).y_filt=filter(b,a,bleb.bridgeddata(blebsweepnum).y);
            end
            apidxes=find(strcmp({ic.eventdata.type},'AP'));
            APdata=struct;
            NEXT=0;
            for api=1:length(apidxes)
                eventidx=apidxes(api);
                si=round(ic.eventdata(eventidx).si*10^6)/10^6;
                icsweepnum=ic.eventdata(eventidx).sweepnum;
                maxh=ic.eventdata(eventidx).maxh;
                maxtime=ic.eventdata(eventidx).maxtime;
                stimulated=ic.eventdata(eventidx).stimulated;
                baselineval=ic.eventdata(eventidx).baselineval;
                amplitude=ic.eventdata(eventidx).amplitude;
                halfwidth=ic.eventdata(eventidx).halfwidth;
                blebsweepnum=find([bleb.bridgeddata.realtime]==ic.bridgeddata(icsweepnum).realtime);
                stepback=round(timebefore/si);
                stepforward=round(timeafter/si);
                
                if ~isempty(blebsweepnum) && maxh>stepback && stepforward+maxh<length(ic.bridgeddata(icsweepnum).y)
                    NEXT=NEXT+1;
                    APdata(NEXT).maxtime=maxtime;
                    APdata(NEXT).si=si;
                    APdata(NEXT).time=[-stepback:stepforward]'*si*1000;
                    APdata(NEXT).ic_y=ic.bridgeddata(icsweepnum).y(maxh-stepback:maxh+stepforward)';
                    APdata(NEXT).bleb_y=bleb.bridgeddata(blebsweepnum).y(maxh-stepback:maxh+stepforward)';
                    APdata(NEXT).bleb_y_filt=bleb.bridgeddata(blebsweepnum).y_filt(maxh-stepback:maxh+stepforward)';
                    APdata(NEXT).bleb_y_baseline=mean(APdata(NEXT).bleb_y(1:stepback));
                    APdata(NEXT).bleb_channellabel=bleb.lightdata(blebsweepnum).channellabel;
                    APdata(NEXT).bleb_preamplnum=bleb.lightdata(blebsweepnum).preamplnum;
                    APdata(NEXT).bleb_Amplifiermode=bleb.lightdata(blebsweepnum).Amplifiermode;
                    APdata(NEXT).stimulated=stimulated;
                    APdata(NEXT).baselineval=median(APdata(NEXT).ic_y(1:stepback));
                    APdata(NEXT).amplitude=amplitude;
                    APdata(NEXT).halfwidth=halfwidth;
                end
            end
            si=max([APdata.si]);
            for NEXT=1:length(APdata) %resampling
                if APdata(NEXT).si<si
                    APdata(NEXT).time=resample(APdata(NEXT).time,APdata(NEXT).si*10^6,si*10^6);
                    APdata(NEXT).ic_y=resample(APdata(NEXT).ic_y,APdata(NEXT).si*10^6,si*10^6);
                    APdata(NEXT).bleb_y=resample(APdata(NEXT).bleb_y,APdata(NEXT).si*10^6,si*10^6);
                    APdata(NEXT).bleb_y_filt=resample(APdata(NEXT).bleb_y_filt,APdata(NEXT).si*10^6,si*10^6);
                    APdata(NEXT).si=si;
                end
            end
        end
    end
    APdataOriginal=APdata;
    %% visualize bleb data
    APdata=APdataOriginal;
    recordingmode='V-Clamp';
    baselinepercentilerange=[.01 .99];
    baselinestdmin=3;
    % detect amplitudes
    %      figure(1)
    %     subplot(4,4,i*4+2);
    %     [x,~]=ginput(2);
    x=[-2,2];
    baselinex=[-2, 0];
    for api=1:length(APdata)
        bleb_y_baseline=APdata(api).bleb_y_baseline;
        time=APdata(api).time;
        bleb_y=APdata(api).bleb_y;
        bleb_y_filt=APdata(api).bleb_y_filt;
        ic_y=APdata(api).ic_y;
        idxes=find(time>min(x)&time<max(x));
        baselineidxes=find(time>min(baselinex)&time<max(baselinex));
        if strcmp(recordingmode,'V-Clamp')
            [peakv,peakh]=max(bleb_y_filt(idxes));
            peakv=(peakv-bleb_y_baseline);
            peakh=peakh+idxes(1);
            [baselinepeakv,~]=min(bleb_y_filt(baselineidxes));
            baselinepeakv=-(baselinepeakv-bleb_y_baseline);
        elseif strcmp(recordingmode,'C-Clamp')
            [peakv,peakh]=min(bleb_y_filt(idxes));
            peakv=-(peakv-bleb_y_baseline);
            peakh=peakh+idxes(1);
            [baselinepeakv,~]=max(bleb_y_filt(baselineidxes));
            baselinepeakv=(baselinepeakv-bleb_y_baseline);
        end
        APdata(api).bleb_amplitude=peakv;
        APdata(api).bleb_amplitude_t=time(peakh);
        APdata(api).bleb_baselineamplitude=baselinepeakv;
    end
    
    baselineamplitudes=sort([APdata(strcmp({APdata.bleb_Amplifiermode},recordingmode)).bleb_baselineamplitude]);
    
    [nbase,xbincenters]=hist(baselineamplitudes,100);
    baselineamplitudes=baselineamplitudes(round(baselinepercentilerange(1)*length(baselineamplitudes)):round(baselinepercentilerange(2)*length(baselineamplitudes)));
    [nbase2,~]=hist(baselineamplitudes,xbincenters);
    minamplitude=mean(baselineamplitudes)+std(baselineamplitudes)*baselinestdmin;
    
    todel=find([APdata.bleb_amplitude]<minamplitude&strcmp({APdata.bleb_Amplifiermode},recordingmode) );
    
    APdata(todel)=[];
    
    close all
    figure(33)
    clf
    bar(xbincenters,nbase,'r')
    hold on
    bar(xbincenters,nbase2,'b')
    plot([minamplitude,minamplitude],[0, max(nbase)],'r-','LineWidth',3)
    figure(1)
    clf
    hold on
    figure(2)
    clf
    hold on
    for i =0:3
        figure(1)
        needed=(strcmp({APdata.bleb_Amplifiermode},recordingmode));
        if i==0 % persistent spikelet
            cim='axonal spikelet';
            needed=find(needed&[APdata.stimulated]==0 & [APdata.amplitude]<.06 & [APdata.halfwidth]>300*10^-6); %
        elseif i==1 % persistent spike
            cim='axonal spike';
            needed=find(needed&[APdata.stimulated]==0 & [APdata.amplitude]>.07 &  [APdata.baselineval]<-.00); %& [APdata.halfwidth]>300*10^-6
        elseif i==2 % somatic spike
            cim='somatic spike - hyperpolarized v0';
            needed=find(needed&[APdata.stimulated]==1 & [APdata.amplitude]>.07&  [APdata.baselineval]<-.07);
        elseif i==3 % somatic spike
            cim='somatic spike - depolarized v0';
            needed=find(needed&[APdata.stimulated]==1 & [APdata.amplitude]>.07&  [APdata.baselineval]>-.05);
        end
        figure(1)
        hsom(i+1)=subplot(4,4,i*4+1);
        hold on
        plot([APdata(needed).time],[APdata(needed).ic_y])
        plot(mean([APdata(needed).time],2),mean([APdata(needed).ic_y],2),'k-','LineWidth',3)
        title(cim)
        hax(i+1)=subplot(4,4,i*4+2);
        hold on
        plot([APdata(needed).time],bsxfun(@(a,b)a-b,[APdata(needed).bleb_y_filt],[APdata(needed).bleb_y_baseline]))
        plot(mean([APdata(needed).time],2),mean(bsxfun(@(a,b)a-b,[APdata(needed).bleb_y_filt],[APdata(needed).bleb_y_baseline]),2),'k-','LineWidth',3)
        if isfield(APdata,'bleb_amplitude')
            [~,amplbins]=hist([APdata(strcmp({APdata.bleb_Amplifiermode},recordingmode)).bleb_amplitude],100);
            [~,ampltbins]=hist([APdata(strcmp({APdata.bleb_Amplifiermode},recordingmode)).bleb_amplitude_t],100);
            subplot(4,4,i*4+3);
            hist([APdata(needed).bleb_amplitude],amplbins)
            xlabel('peak amplitude')
            subplot(4,4,i*4+4);
            hist([APdata(needed).bleb_amplitude_t],ampltbins)
            xlabel('peak latency')
        end
        linkaxes(hsom,'xy')
        linkaxes(hax,'xy')
        subplot(4,4,i*4+2);
        axis tight
        subplot(4,4,i*4+1);
        axis tight
        if isfield(APdata,'bleb_amplitude')
            figure(2)
            if i<2
                colorka='rx';
            else
                colorka='ko';
            end
            subplot(3,1,1)
            hold on
            plot([APdata(needed).maxtime],[APdata(needed).bleb_amplitude],colorka)
            ylabel('peak amplitude')
            subplot(3,1,2)
            hold on
            plot([APdata(needed).maxtime],[APdata(needed).bleb_amplitude_t],colorka)
            ylabel('peak latency')
            xlabel('time')
            subplot(3,1,3)
            hold on
            plot([APdata(needed).bleb_amplitude_t],[APdata(needed).bleb_amplitude],colorka)
            xlabel('peak latency')
            ylabel('peak amplitude')
        end
    end
    
end
% return
% %% uj cucc
%
% files=dir(dirs.bridgeddir);
% files([files.isdir])=[];
% sweepdata=struct;
% progressbar('preloading files')
% for fi=1:length(files)
%     load([dirs.bridgeddir,files(fi).name],'lightdata');
%     if isempty(fieldnames(sweepdata))
%         sweepdata=lightdata;
%     else
%         sweepdata=[sweepdata,lightdata];
%     end
%     progressbar(fi/length(files))
% end
%
%
% %% puffnlightstim
% aE_persistent_puffnlightstim %this script plots puffing and light stim experiments.. ap waveforms and onsets are analysed - base values are needed from this main script
%% check electrotonic and chemical connectivity
valtozok.plot.dpi=150;
valtozok.plot.xcm=20;
valtozok.plot.ycm=14;
valtozok.plot.betumeret=14;
valtozok.plot.betutipus='Arial';
valtozok.plot.axesvastagsag=2;
valtozok.plot.xinch=valtozok.plot.xcm/2.54;
valtozok.plot.yinch=valtozok.plot.ycm/2.54;
valtozok.plot.xsize=valtozok.plot.dpi*valtozok.plot.xinch;
valtozok.plot.ysize=valtozok.plot.dpi*valtozok.plot.yinch;


valtozok.gj_baselinelength=.010;
valtozok.gj_baselinelengthend=.08;
valtozok.gj_minlinelength=.05;
valtozok.gj_mincurrampl=-10*10^-12;

valtozok.noAPbeforetheevent=1; %s
valtozok.noAPaftertheevent=.05; %s
valtozok.pairedpulseneeded=0; %boolean
valtozok.pairedpulsedelay=NaN;%.06; %s
valtozok.pairedpulsejitter=.05; %s
valtozok.baselinelength=0.025; %s
valtozok.psplength=.15; %s
valtozok.filterorder=3;
valtozok.cutofffreq=1000;
valtozok.drugwashintime=120;
valtozok.maxy0baselinedifference=.0005;
valtozok.discardpostsweepswithap=1;
valtozok.postrecordingmode='C-Clamp';%'C-Clamp' or 'V-Clamp' or 'any'
valtozok.prerecordingmode='C-Clamp';%'C-Clamp' or 'V-Clamp' or 'any'

aE_checkGJandChemicalSynapse(valtozok,xlsdata,dirs)
%% plotting IV
sweepbordersrelativetorheobase=[-1,2];
scalebarx=[.1];
scalebary=[.05];
valtozok.plot.betumeret=8;
valtozok.plot.axesvastagsag=2;
valtozok.plot.xcm=20;
valtozok.plot.ycm=14;
valtozok.plot.dpi=150;
xinch=valtozok.plot.xcm/2.54;
yinch=valtozok.plot.ycm/2.54;
valtozok.plot.xsize=valtozok.plot.dpi*xinch;
valtozok.plot.ysize=valtozok.plot.dpi*yinch;


dothesecond=zeros(size(xlsdata));
files=dir([dirs.taxonomydir,'IVs']);
files([files.isdir])=[];
datafiles=dir([dirs.taxonomydir,'datafiles']);
datafiles([datafiles.isdir])=[];
for filenum=length(files):-1:1
    %     pause
    fname=files(filenum).name;
    datafileidx=[];
    for datai=1:length(datafiles)
        if any(strfind(datafiles(datai).name,fname))
            datafileidx=datai;
        end
    end
    if ~isempty(datafileidx)
        load([dirs.taxonomydir,'IVs/',files(filenum).name]);
        load([dirs.taxonomydir,'datafiles/',datafiles(datafileidx).name]);
        rheobasesweep=find(cellStruct.apNums>0,1,'first');
        si=mode(diff(iv.time));
        [b,a]=butter(3,15000/(1/mode(diff(iv.time)))/2,'low');
        [bb,aa]=butter(3,1000/(1/mode(diff(iv.time)))/2,'low');
        sweepborders=sweepbordersrelativetorheobase+rheobasesweep;
        sweepnums=[1,[sweepborders(1):sweepborders(2)]];
        vs=[filtfilt(bb,aa,iv.v1)];
        for sweepi=2:length(sweepnums)
            if sweepnums(sweepi)<=iv.sweepnum
                vs(:,sweepi)=[filtfilt(b,a,iv.(['v',num2str(sweepnums(sweepi))]))];
                difi=diff(vs(:,sweepi-1:sweepi)');
                vs(:,sweepi)=vs(:,sweepi)-min(difi)+.01;
            end
        end
        
        figure(33)
        clf
        hold on
        plot(iv.time,vs,'k-','LineWidth',2)
        medv=median(vs(:));
        maxt=max(iv.time);
        plot([0 0]+maxt,[0 scalebary]+medv,'k-','LineWidth',5)
        plot([0 scalebarx]+maxt,[0 0]+medv,'k-','LineWidth',5)
        
        text(maxt+scalebarx/2,medv+scalebary/2,{[num2str(scalebary*1000), 'mV'],[num2str(scalebarx*1000), 'ms']})
        axis tight
        set(gca,'LineWidth',valtozok.plot.axesvastagsag,'FontSize',valtozok.plot.betumeret,'Position',[1/valtozok.plot.xcm 1/valtozok.plot.ycm 1-2/valtozok.plot.xcm 1-2/valtozok.plot.ycm])
        axis off
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 valtozok.plot.xsize/valtozok.plot.dpi valtozok.plot.ysize/valtozok.plot.dpi])
        
        %     print(gcf,[dirs.figuresdir,'/IVs/IV_',xlsdata(prenum).ID,'.jpg'],'-djpeg',['-r',num2str(valtozok.plot.dpi)])
        
        fname([strfind(fname,'.mat'),strfind(fname,'.mat')+1,strfind(fname,'.mat')+2,strfind(fname,'.mat')+3])=[];
        
        saveas(gcf,[dirs.figuresdir,'IV/',fname],'jpg')
        saveas(gcf,[dirs.figuresdir,'IV/',fname],'pdf')
    end
end
