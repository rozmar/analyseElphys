%%
close all
% clear all
projectnames={'CB1elphys','InVivo','Persistent-ChRstim','persistent firing'};
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
    dirs.rawdir=[locations.tgtardir,'AXONdata/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
    % dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    % dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    % dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    % dirs.figuresdir=[dirs.basedir,'Figures/'];
    xlsdata=aE_readxls([dirs.basedir,'invivodata.xls']);
    amplifier='AXON';
elseif projectnum==3;
    overwrite=0;
    locations=marcicucca_locations;
    dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/_persistent/_ChRstim/'];
    dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
    dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
    dirs.eventdir=[dirs.basedir,'Events/'];
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
    dirs.eventdir=[dirs.basedir,'Events/'];
    dirs.eventparaleldir=[dirs.basedir,'Events/paralel/'];
    dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
    dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
    dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];
    dirs.figuresdir=[dirs.basedir,'figures/'];
    amplifier='HEKA';
    xlsdata=aE_readxls([dirs.basedir,'persistentdata_windows.xls']);
end
%%
if strcmp(amplifier,'AXON')
    aE_exportAXONdata(dirs,xlsdata,projectdata.owexport)
else
    %% export Raw data from HEKA
    aE_exportrawHEKAdata(dirs,xlsdata,projectdata.owexport)
    %% generate PGF data, bridge balancing
    overwrite=1;
    aE_generatePGF_bridge_HEKA(dirs,xlsdata,projectdata.owbridge)
end

%% finding events
valtozok.overwrite=projectdata.owevent;
valtozok.plotit=0;
valtozok.threshholdaveragetime=30;%s
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
aE_findevents(valtozok,dirs,parallelcount)

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
    
end
return
%% uj cucc

files=dir(dirs.bridgeddir);
files([files.isdir])=[];
sweepdata=struct;
progressbar('preloading files')
for fi=1:length(files)
    load([dirs.bridgeddir,files(fi).name],'lightdata');
    if isempty(fieldnames(sweepdata))
        sweepdata=lightdata;
    else
        sweepdata=[sweepdata,lightdata];
    end
    progressbar(fi/length(files))
end


%% puffnlightstim
aE_persistent_puffnlightstim %this script plots puffing and light stim experiments.. ap waveforms and onsets are analysed - base values are needed from this main script
%% check electrotonic and chemical connectivity
valtozok.plot.dpi=600;
valtozok.plot.xcm=20;
valtozok.plot.ycm=14;
valtozok.plot.betumeret=8;
valtozok.plot.axesvastagsag=2;

valtozok.plot.xsize=valtozok.plot.dpi*xinch;
valtozok.plot.ysize=valtozok.plot.dpi*yinch;


valtozok.gj_baselinelength=.010;
valtozok.gj_baselinelengthend=.08;
valtozok.gj_minlinelength=.05;
valtozok.gj_mincurrampl=-10*10^-12;

valtozok.noAPbeforetheevent=1; %s
valtozok.noAPaftertheevent=.05; %s
valtozok.pairedpulseneeded=0; %boolean
valtozok.pairedpulsedelay=.06; %s
valtozok.pairedpulsejitter=.01; %s
valtozok.baselinelength=0.025; %s
valtozok.psplength=.15; %s
valtozok.filterorder=3;
valtozok.cutofffreq=1500;
valtozok.drugwashintime=120;
valtozok.maxy0baselinedifference=.0005;
valtozok.discardpostsweepswithap=1;


aE_checkGJandChemicalSynapse(valtozok,xlsdata,dirs)
%% plotting IV

valtozok.plot.betumeret=8;
valtozok.plot.axesvastagsag=2;


xinch=valtozok.plot.xcm/2.54;
yinch=valtozok.plot.ycm/2.54;
valtozok.plot.xsize=valtozok.plot.dpi*xinch;
valtozok.plot.ysize=valtozok.plot.dpi*yinch;


dothesecond=zeros(size(xlsdata));

for prenum=47:length(xlsdata)
    %     pause
    fname=[xlsdata(prenum).ID,'.mat'];
    HEKAfname=xlsdata(prenum).HEKAfname;
    if any(strfind(HEKAfname,','))
        comma=strfind(HEKAfname,',');
        HEKAfname=HEKAfname(1:comma(1)-1);
    end
    load([locations.tgtardir,'MATLABdata/IV/',xlsdata(prenum).setup,'/',HEKAfname,'.mat']);
    
    gsc=xlsdata(prenum).G_S_C;
    hyps=strfind(gsc,'_');
    commas=strfind(gsc,',');
    g=num2str(gsc(1:hyps(1)-1));
    if isempty(commas)
        s=num2str(gsc(hyps(1)+1:hyps(2)-1));
    elseif dothesecond(prenum)==1 & length(commas)==1
        s=num2str(gsc(commas(1)+1:hyps(2)-1));
    elseif dothesecond(prenum)==1 & length(commas)>1
        s=num2str(gsc(commas(1)+1:commas(2)-1));
    else
        s=num2str(gsc(hyps(1)+1:commas(1)-1));
    end
    c=num2str(gsc(hyps(2)+1:end));
    iv=iv.(['g',g,'_s',s,'_c',c]);
    si=mode(diff(iv.time));
    [b,a]=butter(3,15000/(1/mode(diff(iv.time)))/2,'low');
    [bb,aa]=butter(3,1000/(1/mode(diff(iv.time)))/2,'low');
    %     for ii=1:iv.sweepnum
    %         if ii<5
    %             iv.(['v',num2str(ii)])=filtfilt(bb,aa,iv.(['v',num2str(ii)]));
    %         else
    %             iv.(['v',num2str(ii)])=filtfilt(b,a,iv.(['v',num2str(ii)]));
    %         end
    %     end
    sweeplevonas=xlsdata(prenum).IV_sweeplevonas;
    if strcmp(sweeplevonas,'NaN')
        sweeplevonas=NaN;
    else
        sweeplevonas=str2num(sweeplevonas);
    end
    if isnan(sweeplevonas)
        figure(3)
        clf
        subplot(5,1,1)
        hold on;
        plot(iv.time,iv.v1,'k-','LineWidth',2)
        plot(iv.time,iv.(['v',num2str(iv.sweepnum-4)]),'k-','LineWidth',2);
        axis tight
        title(xlsdata(prenum).ID)
        subplot(5,1,2)
        hold on;
        plot(iv.time,iv.v1,'k-','LineWidth',2)
        plot(iv.time,iv.(['v',num2str(iv.sweepnum-3)]),'k-','LineWidth',2);
        axis tight
        title(['Ca: ',num2str(xlsdata(prenum).Ca),' mM    Mg:',num2str(xlsdata(prenum).Mg),' mM'])
        subplot(5,1,3)
        hold on;
        plot(iv.time,iv.v1,'k-','LineWidth',2)
        plot(iv.time,iv.(['v',num2str(iv.sweepnum-2)]),'k-','LineWidth',2);
        axis tight
        subplot(5,1,4)
        hold on;
        plot(iv.time,iv.v1,'k-','LineWidth',2)
        plot(iv.time,iv.(['v',num2str(iv.sweepnum-1)]),'k-','LineWidth',2);
        axis tight
        subplot(5,1,5)
        hold on;
        plot(iv.time,iv.v1,'k-','LineWidth',2)
        plot(iv.time,iv.(['v',num2str(iv.sweepnum)]),'k-','LineWidth',2);
        axis tight
        pause
        sweeplevonas=0;
    end
    
    %
    figure(33)
    clf
    hold on
    plot(iv.time,filtfilt(bb,aa,iv.v1),'k-','LineWidth',2)
    plot(iv.time,filtfilt(b,a,iv.(['v',num2str(iv.sweepnum-(sweeplevonas))])),'k-','LineWidth',2);
    axis tight
    ylim([-.13 .050])
    xlim([0 1])
    %     ylimitek(:,i)=get(gca,'Ylim');
    set(gca,'LineWidth',valtozok.plot.axesvastagsag,'FontSize',valtozok.plot.betumeret,'Position',[1/valtozok.plot.xcm 1/valtozok.plot.ycm 1-2/valtozok.plot.xcm 1-2/valtozok.plot.ycm])
    axis off
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 valtozok.plot.xsize/valtozok.plot.dpi valtozok.plot.ysize/valtozok.plot.dpi])
    print(gcf,[dirs.figuresdir,'/IVs/IV_',xlsdata(prenum).ID,'.jpg'],'-djpeg',['-r',num2str(valtozok.plot.dpi)])
end
figure(1)
clf
hold on

plot([.5 .5],[-.050 -.03],'k-','LineWidth',5)
plot([.5 .6],[-.050 -.050],'k-','LineWidth',5)
ylim([-.13 .050])
xlim([0 1])
axis off

set(gca,'Position',[0 0 1 1])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 valtozok.plot.xsize/valtozok.plot.dpi valtozok.plot.ysize/valtozok.plot.dpi])
print(gcf,[dirs.figuresdir,'/IVs/IVscalebar_100ms_20mV.jpg'],'-djpeg',['-r',num2str(valtozok.plot.dpi)])
