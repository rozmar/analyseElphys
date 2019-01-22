function aE_generatePGF_bridge_HEKA(dirs,xlsdata,overwrite)
xlsdataold=xlsdata;
plotRSvalues=0;
plotbridgedsweeps=0;
RSbaselinelength=.00005; %ennyi időt átlagol össze a feszültség megmérésekor
RSrisetime=.00005;%.00005; %ennyi időt hagy ki az áram injekcióját követően
poolRStime=30; %ebben az idoablakban atlagolja ossze az RS-t a sweep-ek kozott
files=dir(dirs.rawexporteddir);
files([files.isdir])=[];
%progressbar('generating PGF and bridge balancing')
for fnum=1:length(files)
    fname=files(fnum).name(1:end-4);
    %     load([dirs.rawexporteddir,files(fnum).name],'xlsidx');
    if isfield(xlsdataold,'ID')
        xlsidx=find(strcmp({xlsdataold.ID},fname));
        ID=xlsdataold(xlsidx).ID;
    else
        ID=fname;
        xlsidx=NaN;
    end
    
    %     if isempty(xlsidx)
    %         disp(['xls file és filenevek közti összetűzés .. v'])
    %     end
    a=dir([dirs.bridgeddir,ID,'.mat']);
    if isempty(a) | overwrite==1
        temp=load([dirs.rawexporteddir,files(fnum).name]);
        rawdata=temp.rawdata;
        xlsdata=xlsdataold;
        % realtime past midnight error
        idxtocorrect=find([rawdata.realtime]<rawdata(1).realtime);
        for idxi=1:length(idxtocorrect)
            rawdata(idxtocorrect(idxi)).realtime=rawdata(idxtocorrect(idxi)).realtime+24*3600;
        end
        %     rawdata=temp.rawdata;
        %     xlsdata=temp.xlsdata;
        %     xlsidx=temp.xlsidx;
        stimdata=aE_generatePGF_calculateRS(rawdata,plotRSvalues,RSbaselinelength,RSrisetime,poolRStime);  % generating PGF data, calculating RS
        
        %bridge balancing
        bridgeddata=struct;
        lightdata=struct;
        for sweepnum=1:length(rawdata)
            time=[0:rawdata(sweepnum).si:rawdata(sweepnum).si*(length(stimdata(sweepnum).y)-1)];
            bridgeddata(sweepnum).y=rawdata(sweepnum).y-stimdata(sweepnum).yforbridge*stimdata(sweepnum).RS;
            bridgeddata(sweepnum).si=rawdata(sweepnum).si;
            bridgeddata(sweepnum).realtime=rawdata(sweepnum).realtime;
            bridgeddata(sweepnum).timertime=rawdata(sweepnum).timertime;
            bridgeddata(sweepnum).channellabel=rawdata(sweepnum).channellabel;
            if isfield(rawdata, 'stimulation')
                bridgeddata(sweepnum).stimulation=rawdata(sweepnum).stimulation;
            end
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
            lightdata(sweepnum).si=rawdata(sweepnum).si;
            lightdata(sweepnum).realtime=rawdata(sweepnum).realtime;
            lightdata(sweepnum).timertime=rawdata(sweepnum).timertime;
            lightdata(sweepnum).channellabel=rawdata(sweepnum).channellabel;
            lightdata(sweepnum).preamplnum=stimdata(sweepnum).preamplnum;
            lightdata(sweepnum).Amplifiermode=stimdata(sweepnum).Amplifiermode;
            lightdata(sweepnum).segmenths=stimdata(sweepnum).segmenths;
            lightdata(sweepnum).RS=stimdata(sweepnum).RS;
            lightdata(sweepnum).ID=ID;
        end
        %bridge balancing
        xlsdata=xlsdataold;
        save([dirs.bridgeddir,ID],'lightdata','stimdata','bridgeddata','RSbaselinelength','RSrisetime','poolRStime','xlsdata','xlsidx','-v7.3')
        disp([ID,' done (bridge_stim)'])
    else
        disp([ID,' already done .. skipped (bridge_stim)'])
    end
    progressbar(fnum/length(files))
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