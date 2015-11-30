function aE_exportAXONdata(dirs,xlsdata,overwrite)
for xlsidx=1:length(xlsdata)
    expname=xlsdata(xlsidx).AXONfname;
    cellname=xlsdata(xlsidx).Cellname;
    setupname=xlsdata(xlsidx).setup;
    
    % expname='proba';
    % cellname='a';
    rawdirnow=[dirs.rawdir,setupname,'/',expname,'/',cellname,'/'];
    files=dir(rawdirnow);
    files([files.isdir])=[];
    bridgeddata=struct;
    stimdata=struct;
    for filei=1:length(files)
%         [~,~,h]=abfload([rawdirnow,files(filei).name],'info') ;
        [temp,si,h]=abfload([rawdirnow,files(filei).name],'channels',{xlsdata(xlsidx).ChannelRec,xlsdata(xlsidx).ChannelStim});%,'channels','a');%channelstodo(chi)
        for sweepi=1:size(temp,3)
            if isempty(fieldnames(bridgeddata))
                NEXT=1;
            else
                NEXT=length(bridgeddata)+1;
            end
            bridgeddata(NEXT).y=temp(:,1,sweepi)'/1000;
            bridgeddata(NEXT).si=si/1000000;
            bridgeddata(NEXT).realtime=h.recTime(1)+(sweepi)*h.dataPtsPerChan*h.si/1000000;
            bridgeddata(NEXT).endtime=h.recTime(2)+(sweepi)*h.dataPtsPerChan*h.si/1000000;
            bridgeddata(NEXT).channellabel=xlsdata(xlsidx).ChannelRec;
            
            stimdata(NEXT).y=round(temp(:,2,sweepi)'/20)/1000000000000*20;
            stimdata(NEXT).segmenths=find(abs(diff(stimdata(NEXT).y))>0);
%             stimdata(NEXT).yforbridge=temp(:,2,sweepi)'/1000000000000;
            stimdata(NEXT).channellabel=xlsdata(xlsidx).ChannelStim;
            stimdata(NEXT).Amplifiermode='C-Clamp';
        end
    end
    save([dirs.bridgeddir,xlsdata(xlsidx).ID],'stimdata','bridgeddata','xlsdata','xlsidx','-v7.3')
end