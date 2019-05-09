%ANALYSELPHYS_MAIN.m file is an ordered collection of scripts. Their purpouse
%is to export and analyze electrophysiological data obtained during single
%channel extracellular or intracellular recordings.
%   -(1) aE_projectselector_main.m where you can choose a project to analyze (GUI)
%   this script also generates the directory structure (dirs) and the metadata
%   of the selected project (xlsdata)
%   based on the user selection the following scripts are run:
%   -(2) persistent_statistics_main.m - collecting and saving statistics in a .mat 
%   file (slow oscillation power, aAP number and frequency) and adds these to
%   xlsdata variable as fields
%   -(3) persistent_statistics_loader.m loads the statistics that were earlier generated
%   (alternative to persistent_statistics_main,)
%   -(4) starts GUIs for inspection and analysis of exported data: 
%                   -aE_InspectTraces.m: elphys, video, breathing, PSD..
%                   etc. data can be loaded and inspected, zoomed,
%                   filtered. From this GUI multiple subGUIs can be started
%                   for detailed analysis.
%                   -aE_videoanalyzer_ROIselector.m: select ROIs on mouse face 
%                   videos done during in vivo recordings
%   -(5) creates an xls file for the neurontaxonomy script based on the input
%   xls file (this makes the taxonomy_routines.m script able to analyze the IVs 
%   of the cells)
%   -(6) aE_exportrawHEKAdata.m exports raw HEKA data from .dat files. This 
%   script uses batch communication to control a Fitmaster program and export 
%   the data into .mat files. The exported .mat files are read in and arranged
%   in a struct array. The exported data is then saved in the
%   dirs.rawexporteddir folder (defined in aE_projectselector_main.m).
%   -(7) aE_generatePGF_bridge_HEKA.m further processes the data obtained
%   by aE_exportrawHEKAdata.m:  - generates a stimulus struct where the
%                               indexes correspond to the indexes in the 
%                               raw exported data.
%                               - calculates series resistance in each
%                               sweep and applies offline bridge balance
%                               compensation. The data is saved in a new
%                               sruct similarly to the raw data.
%   The two newly generated structs (bridgeddata and stimdata) are saved in
%   the dirs.bridgeddata folder (defined in aE_projectselector_main.m).
%   -(8) if a thermo sensor was placed in front of the mouse during an in 
%   vivo experiment, the voltage fluctuation generated by the sensor and
%   recorded by the HEKA amplifier is also exported similarly to aE_exportrawHEKAdata.m
%   The exported raw data is saved in the dirs.breathingdir folder (defined 
%   in aE_projectselector_main.m).
%   -(9) aE_findevents.m tries to extract PSPs and APs. The metadata of the 
%   events are stored in chronological order in a struct (eventdata). With
%   the help of the eventdata struct the events can be located in the
%   rawdata or bridgeddata structs.
%   -(10) aE_PSD_export.m runs a wavelet analysis on the recording and
%   saves the resulting power spectrum density matrices as a struct array with
%   the same indexing as rawdata, bridgeddata or stimdata structs. The
%   resulting struct is saved in the dirs.PSDdir directory.
%   -(11) the mean, min, max, median of PSDs is extracted for easier
%   display in the aE_InspectTraces.m script.
%   -(12) wavelet anaysis is done for breathing as well but with different
%   parameters. The results are saved in the dirs.breathingdir_psd folder.
%   -(13) aE_analyzevideo_main.m loads the videos recorded during the
%   experiment. The videos are downsampled and saved. The average changes in pixel
%   intensity is also measured and saved. The results are saved in the
%   dirs.videodir folder.
%   -(14) aE_analyzevideo_PCA_on_ROIs.m does PCA on the ROIs (selected in 
%   aE_videoanalyzer_ROIselector.m) on the videos (exported in aE_analyzevideo_main.m)
%   This script is currently not working properly, the results are not
%   helping a lot..
%   -(15) aE_analyzevideo_pupilsize_load.m loads the pupil sizes generated
%   by a Labview scrip (Gáspár Oláh). The data is stored in the
%   dirs.videodir folder, in the 'eye' subfolder.
%   -(16) aE_analyzevideo_normalize_pupil_movement.m normalizes pupil
%   diameters and movements accross experiments. This makes comparing
%   experiments possible. The data is saved in the dirs.videodir folder, in
%   the percentiles subfolder.
%   -(17) all the IVs are plotted and saved in the dirs.figuresdir folder, 
%   in the IV subfolder.
%   -(18) aE_checkGJandChemicalSynapse.m averages presynaptic action
%   potentials, hyper- and depolarizing square pulses and averages them so
%   PSPs and GJ currents will be apparent. The generated figures are saved
%   in the dirs.figuresdir folder
%
% See also aE_projectselector_main, aE_InspectTraces, aE_exportrawHEKAdata,
% aE_generatePGF_bridge_HEKA, aE_findevents, aE_PSD_export,
% aE_analyzevideo_main, aE_checkGJandChemicalSynapse
close all
clear all
locations=marcicucca_locations;
%% (1) select project
[dirs,xlsdata,projectdata,amplifier]=aE_projectselector_main;
%% (2) perform statistics (SO, aAP)
if projectdata.dostatistics==1
    persistent_statistics_main(dirs,xlsdata)
end
%% (3) load statistics
xlsdata=persistent_statistics_loader(dirs,xlsdata,projectdata);
%% (4) start inspecttraces or videoanalyzer
if projectdata.inspecttraces==1
    aE_InspectTraces(dirs,xlsdata);
    return
elseif projectdata.selectvideoROIs==1
    aE_videoanalyzer_ROIselector(dirs, xlsdata);
    return
end
%% (5) create neurontaxonomy xls file
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
            if projectdata.projectnum==4 %1-2:rodent; 3-4:human; 1,3: no rhythmic; 2,4:rhythmic  
                taxdata(NEXT).anatgroup=num2str(strcmpi(xlsdata(xlsi).species,'human')*2+any(strfind(xlsdata(xlsi).persistentfiring,'rhythmic'))+1);
            else
                taxdata(NEXT).anatgroup='0';
            end
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
%% (6) and (7)
if strcmp(amplifier,'AXON') % (6) this is a simple script that uses abfload.m to load data generated with axon amplifiers 
    aE_exportAXONdata(dirs,xlsdata,projectdata.owexport)
else
    % (6) export Raw data from HEKA
    aE_exportrawHEKAdata(dirs,xlsdata,projectdata.owexport)
    % (7) generate PGF data, bridge balancing
    aE_generatePGF_bridge_HEKA(dirs,xlsdata,projectdata.owbridge)
end
%% (8) exporting breathing
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

%% (9) finding events
valtozok.overwrite=projectdata.owevent;
valtozok.overwritebefore=datenum(datetime('today'));%;datenum(2018,12,29);%;
valtozok.plotit=0;
valtozok.sdval=3;
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
parallelcount=2;
aE_findevents(valtozok,dirs,parallelcount,xlsdata)
paralleldata.count=NaN;
while isnan(paralleldata.count) | ~isempty(paralleldata.files) % this loop only tacks the progress of the aE_findevents.m script
    paralleldata.files=dir(dirs.eventparaleldir);
    paralleldata.files([paralleldata.files.isdir])=[];
    pause(3)
    paralleldata.prevcount=paralleldata.count;
    paralleldata.count=length(paralleldata.files);
    if paralleldata.prevcount~=paralleldata.count
        disp(['waiting for ',num2str(paralleldata.count),' eventfinding scripts to finish'])
    end
end

%% (10) extracting PSD
% lower frequency band - linear
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
% extracting PSD - high frequency - linear
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
    parameters.step=2; %- frequency stepsize
    parameters.scale=1;% - scale type, can be linear (1) or logarithmic (2)
    parameters.wavenumber=9;% - number of waves in wavelet
    parameters.waveletlength=30;
    parameters.addtaper=1;
    parameters.taperlength=30;
    
    valtozok.parameters=parameters;
    aE_PSD_export(dirs,xlsdata,valtozok);
end

% extracting PSD - logarithmic
if isfield(dirs,'PSDdir_log')
    valtozok=struct;
    valtozok.overwrite=0;
    valtozok.analyseonlyfield=0;
    valtozok.downsamplenum='auto';
    valtozok.log=1;
    valtozok.onlyawake=0;
    parameters=struct;
    parameters.min=.5; %minimal frequency for decomposition
    parameters.max=200;% - maximal frequency for decomposition
    parameters.step=4; %- frequency stepsize
    parameters.scale=2;% - scale type, can be linear (1) or logarithmic (2)
    parameters.wavenumber=9;% - number of waves in wavelet
    parameters.waveletlength=30;
    parameters.addtaper=1;
    parameters.taperlength=30;
    
    valtozok.parameters=parameters;
    aE_PSD_export(dirs,xlsdata,valtozok);
end

%% (11) PSD stats
dirfields=fieldnames(dirs);
dirstodo=dirfields(~cellfun(@isempty,regexp(dirfields,'PSDdir')));
for diri=1:length(dirstodo)
    dirnow=dirs.(dirstodo{diri});
    destdir=[dirnow,'stats/'];
    files=dir(dirnow);
    files([files.isdir])=[];
    for filei=1:length(files)
        fname=files(filei).name;
        a=dir([dirnow,'stats/',fname]);
        if isempty(a)
            load([dirnow,fname]);
            PSDdata_stats=rmfield(PSDdata,{'powerMatrix','y','si_powerMatrix'});
            if isfield(PSDdata_stats,'compress_offset')
                PSDdata_stats=rmfield(PSDdata_stats,{'compress_offset','compress_multiplier'});
            end
            for sweepi=1:length(PSDdata);
                if isfield(PSDdata,'compress_offset')
                    PSDdata(sweepi).powerMatrix=double(PSDdata(sweepi).powerMatrix)*PSDdata(sweepi).compress_multiplier+PSDdata(sweepi).compress_offset;
                end
                PSDdata_stats(sweepi).mu=mean(PSDdata(sweepi).powerMatrix,2);
                PSDdata_stats(sweepi).sigma=std(PSDdata(sweepi).powerMatrix,[],2);
                PSDdata_stats(sweepi).min=min(PSDdata(sweepi).powerMatrix,[],2);
                PSDdata_stats(sweepi).max=max(PSDdata(sweepi).powerMatrix,[],2);
                PSDdata_stats(sweepi).median=median(PSDdata(sweepi).powerMatrix,2);
                PSDdata_stats(sweepi).length=length(PSDdata(sweepi).y)*PSDdata(sweepi).si_powerMatrix;
            end
            save([destdir,fname],'PSDdata_stats','-v7.3');
            disp(['PSD stat export from ',fname, ' is done'])
        end
    end
end

%% (12) extracting PSD for breathing
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

%% (13) extracting movement info from video
if isfield(dirs,'videodir')
    valtozok.sampleinterval=.01;
    valtozok.overwrite=0;
    for xlsi=1:length(xlsdata)%length(xlsdata):-1:1
        if any(strfind(xlsdata(xlsi).anaesthesia,'awake'))
            valtozok.setupname=xlsdata(xlsi).setup;
            valtozok.filename=xlsdata(xlsi).HEKAfname;
            aE_analyzevideo_main(valtozok,dirs);
        end
    end
end
%% (14) PCA on ROIs
if isfield(dirs,'videodir')
    overwrite=0;
    aE_analyzevideo_PCA_on_ROIs(dirs,xlsdata,overwrite)
end
%% (15) processing Gaspar's pupil diameter data
if isfield(dirs,'videodir')
    aE_analyzevideo_pupilsize_load(dirs,xlsdata)
end
%% (16) normalizing pupil size and movement data
if isfield(dirs,'videodir')
    aE_analyzevideo_normalize_pupil_movement(dirs,xlsdata)
end
%% (17) plotting IV
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

%% (18) check electrotonic and chemical connectivity
valtozok.plot.dpi=150;
valtozok.plot.xcm=20;
valtozok.plot.ycm=14;
valtozok.plot.betumeret=8;
valtozok.plot.betutipus='Arial';
valtozok.plot.axesvastagsag=2;
valtozok.plot.xinch=valtozok.plot.xcm/2.54;
valtozok.plot.yinch=valtozok.plot.ycm/2.54;
valtozok.plot.xsize=valtozok.plot.dpi*valtozok.plot.xinch;
valtozok.plot.ysize=valtozok.plot.dpi*valtozok.plot.yinch;
valtozok.plot.savejpg=1;
valtozok.plot.savefig=0;

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
[Selection,ok] = listdlg('ListString',{xlsdata.ID},'ListSize',[300 600]); % az XLS file alapján kiválasztjuk, hogy melyik file összes mérésén szeretnénk végigmenni

for xlsi=length(Selection):-1:1 %xlsnum=1:length(Selection)%going throught potential presynaptic cells
    xlsnum=Selection(xlsi);
    close all
    prenum=xlsnum;%Selection(xlsnum);
    aE_checkGJandChemicalSynapse(valtozok,xlsdata,dirs,prenum);
end

