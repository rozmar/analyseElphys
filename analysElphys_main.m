%% alapadatok
close all
clear all
overwrite=0;
locations=marcicucca_locations;
dirs.basedir=[locations.tgtardir,'ANALYSISdata/marci/Human_rosehip/CB1elphys/'];
dirs.rawexporteddir=[dirs.basedir,'Exported_raw/'];
dirs.bridgeddir=[dirs.basedir,'Bridged_stim/'];
dirs.eventdir=[dirs.basedir,'Events/'];
dirs.onlyAPeventdir=[dirs.basedir,'Events_onlyAP/'];
dirs.grpupedeventdir=[dirs.basedir,'Events_grouped/'];
dirs.stimepochdir=[dirs.basedir,'Stimepochs/'];

xlsdata=aE_readxls([dirs.basedir,'cb1elphys.xls']);

%% export Raw data from HEKA
for xlsidx=1:length(xlsdata)
    a=dir([dirs.rawexporteddir,xlsdata(xlsidx).ID,'.mat']);
    if isempty(a) | overwrite==1
        rawdata=HEKAexportbytime_main(xlsdata(xlsidx).HEKAfname,xlsdata(xlsidx).setup,xlsdata(xlsidx).Channel,xlsdata(xlsidx).startT,xlsdata(xlsidx).endT);
        save([dirs.rawexporteddir,xlsdata(xlsidx).ID],'rawdata','xlsdata','xlsidx')
        disp([xlsdata(xlsidx).ID,' done'])
    else
        disp([xlsdata(xlsidx).ID,' already done.. skipped'])
    end
end
xlsdataold=xlsdata;

%% generate PGF data, bridge balancing
overwrite=0;
plotRSvalues=0;
plotbridgedsweeps=0;
RSbaselinelength=.00005; %ennyi időt átlagol össze a feszültség megmérésekor
RSrisetime=.00005; %ennyi időt hagy ki az áram injekcióját követően
poolRStime=30; %ebben az idoablakban atlagolja ossze az RS-t a sweep-ek kozott
files=dir(dirs.rawexporteddir);
files([files.isdir])=[];
for fnum=length(files):-1:1%1:length(files)
    fname=files(fnum).name(1:end-4);
%     load([dirs.rawexporteddir,files(fnum).name],'xlsidx');
    xlsidx=find(strcmp({xlsdata.ID},fname));
    if isempty(xlsidx)
        disp(['xls file és filenevek közti összetűzés'])
        pause
    end
    a=dir([dirs.bridgeddir,xlsdata(xlsidx).ID,'.mat']);
    if isempty(a) | overwrite==1
        load([dirs.rawexporteddir,files(fnum).name]);
        xlsdata=xlsdataold;

        %     rawdata=temp.rawdata;
        %     xlsdata=temp.xlsdata;
        %     xlsidx=temp.xlsidx;
        stimdata=aE_generatePGF_calculateRS(rawdata,plotRSvalues,RSbaselinelength,RSrisetime,poolRStime);  % generating PGF data, calculating RS
        %bridge balancing
        bridgeddata=struct;
        for sweepnum=1:length(rawdata)
            time=[0:rawdata(sweepnum).si:rawdata(sweepnum).si*(length(stimdata(sweepnum).y)-1)];
            bridgeddata(sweepnum).y=rawdata(sweepnum).y-stimdata(sweepnum).yforbridge*stimdata(sweepnum).RS;
            bridgeddata(sweepnum).si=rawdata(sweepnum).si;
            bridgeddata(sweepnum).realtime=rawdata(sweepnum).realtime;
            bridgeddata(sweepnum).timertime=rawdata(sweepnum).timertime;
            bridgeddata(sweepnum).channellabel=rawdata(sweepnum).channellabel;
            if plotbridgedsweeps==1 & max(diff(stimdata(sweepnum).y))>100*10^-12
                figure(12211)
                clf
                subplot(3,1,1)
                plot(time,rawdata(sweepnum).y)
                hold on
                plot(time(stimdata(sweepnum).segmenths),rawdata(sweepnum).y(stimdata(sweepnum).segmenths),'ro')
                title('raw data')
                subplot(3,1,2)
                plot(time,stimdata(sweepnum).y)
                ylim([-950*10^-12 950*10^-12])
                title('stimulus')
                subplot(3,1,3)
                plot(time,bridgeddata(sweepnum).y)
                title('bridged data')
                pause
            end
        end
        %bridge balancing
        save([dirs.bridgeddir,xlsdata(xlsidx).ID],'stimdata','bridgeddata','RSbaselinelength','RSrisetime','poolRStime','xlsdata','xlsidx','-v7.3')
        disp([xlsdata(xlsidx).ID,' done (bridge_stim)'])
    else
        disp([xlsdata(xlsidx).ID,' already done .. skipped (bridge_stim)'])
    end
end
% return
% % % %% search for NAN data ezek jottek ki eddig: 1410293rm_4_1_3.mat 1211291rm_4_3_4.mat
% % % startdir=dirs.bridgeddir;
% % % files=dir(startdir);
% % % files([files.isdir])=[];
% % % for i=119: length(files)
% % %     fname=files(i).name;
% % %     load([startdir,fname],'bridgeddata');
% % %     nanperc(i)=sum(isnan([bridgeddata.y]))/length([bridgeddata.y]);
% % %     disp([fname,'  -  ',num2str(nanperc(i))]);
% % %     
% % % end
%% finding events
valtozok.overwrite=0;
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

aE_findevents(valtozok,dirs)

% return