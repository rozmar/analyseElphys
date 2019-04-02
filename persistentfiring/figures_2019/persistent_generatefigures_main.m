%% figure_generation_2019_main
%% in vitro human aAP stats
% close all
species=unique({xlsdata.species});
for xlsi=1:length(xlsdata)
    if any(strfind(xlsdata(xlsi).persistentfiring,'tonic'))
        xlsdata(xlsi).persistent_tonic=1;
    else
        xlsdata(xlsi).persistent_tonic=0;
    end
    if any(strfind(xlsdata(xlsi).persistentfiring,'sporadic'))
        xlsdata(xlsi).persistent_sporadic=1;
    else
        xlsdata(xlsi).persistent_sporadic=0;
    end
    if any(strfind(xlsdata(xlsi).persistentfiring,'rhythmic'))
        xlsdata(xlsi).persistent_rhythmic=1;
    else
        xlsdata(xlsi).persistent_rhythmic=0;
    end
    if any(strfind(xlsdata(xlsi).persistentfiring,'not_tested'))
        xlsdata(xlsi).persistent_tested=0;
    else
        xlsdata(xlsi).persistent_tested=1;
    end
    if any(strfind(xlsdata(xlsi).persistentfiring,'dunno'))
        xlsdata(xlsi).persistent_dunno=1;
    else
        xlsdata(xlsi).persistent_dunno=0;
    end
    if any(strfind(xlsdata(xlsi).persistentfiring,'none'))
        xlsdata(xlsi).persistent_nopersistent=1;
    else
        xlsdata(xlsi).persistent_nopersistent=0;
    end
end
cellnums=struct;
ishuman=strcmp({xlsdata.species},'Human')& ~[xlsdata.field];
islayer1=~cellfun(@any,regexpi({xlsdata.celltype},'L2')) & ~cellfun(@any,regexpi({xlsdata.celltype},'pyr'));
neededcells=ishuman & islayer1;
cellnums.all=sum(neededcells);
cellnums.persistent_dunno=sum(neededcells&[xlsdata.persistent_dunno]==1);
cellnums.persistent_tested=sum(neededcells&[xlsdata.persistent_dunno]==0 & [xlsdata.persistent_tested]==1);
cellnums.persistent_rhythmic=sum(neededcells&[xlsdata.persistent_dunno]==0 & [xlsdata.persistent_tested]==1 & [xlsdata.persistent_rhythmic]==1);
cellnums.persistent_tonic=sum(neededcells&[xlsdata.persistent_dunno]==0 & [xlsdata.persistent_tested]==1 & [xlsdata.persistent_tonic]==1);
cellnums.persistent_sporadic=sum(neededcells&[xlsdata.persistent_dunno]==0 & [xlsdata.persistent_tested]==1 & [xlsdata.persistent_sporadic]==1);
cellnums.persistent_sporadic_tonic=sum(neededcells&[xlsdata.persistent_dunno]==0 & [xlsdata.persistent_tested]==1 & [xlsdata.persistent_sporadic]==1 & [xlsdata.persistent_tonic]==1);
cellnums.persistent_sporadic_rhythmic=sum(neededcells&[xlsdata.persistent_dunno]==0 & [xlsdata.persistent_tested]==1 & [xlsdata.persistent_sporadic]==1 & [xlsdata.persistent_rhythmic]==1);
cellnums.persistent_sporadic_only=cellnums.persistent_sporadic-cellnums.persistent_sporadic_rhythmic-cellnums.persistent_sporadic_tonic;
cellnums.persistent_no_aAP=sum(neededcells&[xlsdata.persistent_dunno]==0 & [xlsdata.persistent_tested]==1 & [xlsdata.persistent_nopersistent]==1);
% tested=cellnums.persistent_rhythmic+cellnums.persistent_tonic+cellnums.persistent_nopersistent+cellnums.persistent_sporadic-cellnums.persistent_sporadic_rhythmic-cellnums.persistent_sporadic_tonic;
piechartorder={'persistent_sporadic_only','persistent_no_aAP','persistent_sporadic_tonic','persistent_tonic','persistent_rhythmic','persistent_sporadic_rhythmic'};
piechartstring={};
cellnumvec=zeros(length(piechartorder),1);
for i=1:length(piechartorder)
    cellnumvec(i)=cellnums.(piechartorder{i});
    piechartstring{i} = [strrep(strrep(piechartorder{i},'_',' '),'persistent ',''),' (n = ',num2str(cellnumvec(i)),')'];
end
figure(9)
clf
p=pie(cellnumvec);%,piechartstring
% title(['persistent firing in n=',num2str(cellnums.persistent_tested),' human interneurons'])
% színezés és sráfozás
colors=struct;
colors.no_aAP=[1 1 1];%[255,237,160]/256;
colors.tonic=[254,178,76]/256;
colors.rhythmic=[240,59,32]/256;
colors.colororder=[];
for i=1:length(piechartorder)
    if any(regexpi(piechartorder{i},'rhythmic'))
        colors.colororder(:,i)=colors.rhythmic;
    elseif any(regexpi(piechartorder{i},'tonic'))
        colors.colororder(:,i)=colors.tonic;
    else
        colors.colororder(:,i)=colors.no_aAP;
    end
end
patches=p.findobj('type','Patch');
patchestostripe=find(cellfun(@any,strfind(piechartorder,'sporadic')));
for i = 1:length(patches)
    if any(i==patchestostripe)
        patches(i)=hatchfill(patches(i), 'single', 45, 10,colors.colororder(:,i));
    else
        patches(i).FaceColor=colors.colororder(:,i);
    end
end
[~,hlegend]=legend(['no persistent firing (n = ',num2str(cellnums.persistent_no_aAP),')'],['sporadic aAPs (n = ',num2str(cellnums.persistent_sporadic),')'],['tonic persistent firing (n = ',num2str(cellnums.persistent_tonic),')'],['rhythmic persistent firing (n = ',num2str(cellnums.persistent_rhythmic),')'],'Location','southoutside');
PatchInLegend = hlegend.findobj('type', 'Patch');
PatchInLegend(1).FaceColor=colors.no_aAP;
PatchInLegend(2).FaceColor=[1 1 1];
PatchInLegend(3).FaceColor=colors.tonic;
PatchInLegend(4).FaceColor=colors.rhythmic;

legend boxoff 
% hatchfill(PatchInLegend(4), 'single', 45, 5,colors.no_aAP); % this
% doesn't work for some reason...fuck

valtozok.xcm=5;
valtozok.ycm=5;
valtozok.ycm_current=1;
valtozok.fontsize=8;
valtozok.fonttype='Helvetica';
valtozok.axesvastagsag=1;
valtozok.dpi=900;
set(gca,'LineWidth',valtozok.axesvastagsag,'FontSize',valtozok.fontsize,'Fontname',valtozok.fonttype,'Units','normalized','Position',[.25 .25 .5 .5])
set(gcf,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[valtozok.xcm/.5 valtozok.ycm/.5]+2,'PaperPosition',[2 2 valtozok.xcm/.5 valtozok.ycm/.5])
saveas(gcf,[dirs.figuresdir,'Human_aAP_stats.pdf'])
print(gcf,[dirs.figuresdir,'Human_aAP_stats.jpg'],'-djpeg',['-r',num2str(valtozok.dpi)])

%% human persistent anatomy
PFtypes={'rhythmic','sporadic','tonic','nopersistent'};
cellanatomy=struct;
for pftypei=1:length(PFtypes)
    pftype=PFtypes{pftypei};
    ishuman=strcmp({xlsdata.species},'Human')& ~[xlsdata.field];
    islayer1=~cellfun(@any,regexpi({xlsdata.celltype},'L2')) & ~cellfun(@any,regexpi({xlsdata.celltype},'pyr'));
    neededcells=ishuman & islayer1;
    rhythmic_celltypes={xlsdata(neededcells&[xlsdata.persistent_dunno]==0 & [xlsdata.persistent_tested]==1 & [xlsdata.(['persistent_',pftype])]==1).celltype};
    IDs={xlsdata(neededcells&[xlsdata.persistent_dunno]==0 & [xlsdata.persistent_tested]==1 & [xlsdata.(['persistent_',pftype])]==1).ID};
    unique(rhythmic_celltypes);
    cellanatomy.(pftype).unidentified=0;
    cellanatomy.(pftype).ngf=0;
    cellanatomy.(pftype).nonngf=0;
    cellanatomy.(pftype).no_recovery=0;
    for i=1:length(rhythmic_celltypes)
        if any(regexpi(rhythmic_celltypes{i},'nonngf'))
            cellanatomy.(pftype).nonngf=cellanatomy.(pftype).nonngf+1;
        elseif any(regexpi(rhythmic_celltypes{i},'ngf'))
            cellanatomy.(pftype).ngf=cellanatomy.(pftype).ngf+1;
        elseif any(regexpi(rhythmic_celltypes{i},'none'))
            cellanatomy.(pftype).no_recovery=cellanatomy.(pftype).no_recovery+1;
        else
                    disp([IDs{i},': ',rhythmic_celltypes{i}])
            cellanatomy.(pftype).unidentified=cellanatomy.(pftype).unidentified+1;
        end
    end
    
    
    piechartorder={'no_recovery','nonngf','unidentified','ngf'};%unique(rhythmic_celltypes);%{'ngf','L1interneuron','dunno','none'};
    todel=zeros(size(piechartorder));
    for i=1:length(piechartorder)
        if cellanatomy.(pftype).(piechartorder{i})<1
           todel(i)=1;
        end
    end
    piechartorder(find(todel))=[];
    piechartstring={};
    cellnumvec=zeros(length(piechartorder),1);
    for i=1:length(piechartorder)
        cellnumvec(i)=cellanatomy.(pftype).(piechartorder{i});
        piechartstring{i} = [strrep(piechartorder{i},'_',' '),' (n = ',num2str(cellnumvec(i)),')'];
    end
    figure(10)
    clf
    p=pie(cellnumvec);
    
    title(['anatomy of ',pftype,' cells'])
    
    colors=struct;
    colors.no_recovery=[1 1 1];%[255,237,160]/256;
    colors.non_ngf=[254,178,76]/256;
    colors.ngf=[240,59,32]/256;
    colors.unidentified=[255,237,160]/256;
    colors.colororder=[];
    for i=1:length(piechartorder)
        if any(regexpi(piechartorder{i},'nonngf'))
            colors.colororder(:,i)=colors.non_ngf;
        elseif any(regexpi(piechartorder{i},'ngf'))
            colors.colororder(:,i)=colors.ngf;
        elseif any(regexpi(piechartorder{i},'unidentified'))
            colors.colororder(:,i)=colors.unidentified;
        else
            colors.colororder(:,i)=colors.no_recovery;
        end
    end
    patches=p.findobj('type','Patch');
    patchestostripe=find(cellfun(@any,strfind(piechartorder,'sporadic')));
    for i = 1:length(patches)
        if any(i==patchestostripe)
            patches(i)=hatchfill(patches(i), 'single', 45, 10,colors.colororder(:,i));
        else
            patches(i).FaceColor=colors.colororder(:,i);
        end
    end
    offset=0;
    [~,hlegend]=legend(['no recovery (n = ',num2str(cellanatomy.(pftype).no_recovery),')'],['non NGF (n = ',num2str(cellanatomy.(pftype).nonngf),')'],['NGF (n = ',num2str(cellanatomy.(pftype).ngf),')'],['unidentified (n = ',num2str(cellanatomy.(pftype).unidentified),')'],'Location','southoutside');
    PatchInLegend = hlegend.findobj('type', 'Patch');
    if cellanatomy.(pftype).no_recovery>0
        PatchInLegend(1+offset).FaceColor=colors.no_recovery;
    else
        offset=offset-1;
    end
    if cellanatomy.(pftype).nonngf>0
        PatchInLegend(2+offset).FaceColor=colors.non_ngf;
    else
        offset=offset-1;
    end
    if cellanatomy.(pftype).ngf>0
        PatchInLegend(3+offset).FaceColor=colors.ngf;
    else
        offset=offset-1;
    end
    if cellanatomy.(pftype).unidentified>0
        PatchInLegend(4+offset).FaceColor=colors.unidentified;
    else
        offset=offset-1;
    end
    legend boxoff
    
    valtozok.xcm=5;
    valtozok.ycm=5;
    valtozok.ycm_current=1;
    valtozok.fontsize=8;
    valtozok.fonttype='Helvetica';
    valtozok.axesvastagsag=1;
    valtozok.dpi=900;
    set(gca,'LineWidth',valtozok.axesvastagsag,'FontSize',valtozok.fontsize,'Fontname',valtozok.fonttype,'Units','normalized','Position',[.25 .25 .5 .5])
    set(gcf,'PaperUnits','centimeters','PaperPositionMode','manual','PaperSize',[valtozok.xcm/.5 valtozok.ycm/.5]+2,'PaperPosition',[2 2 valtozok.xcm/.5 valtozok.ycm/.5])
    saveas(gcf,[dirs.figuresdir,'Human_aAP_',pftype,'_anatomy_stats.pdf'])
    print(gcf,[dirs.figuresdir,'Human_aAP_',pftype,'_anatomy_stats.jpg'],'-djpeg',['-r',num2str(valtozok.dpi)])
end
%% in vivo recording statistics 
persistent_generatefigures_recording_statistics % takes quite a long time.. recording lengths, RS values etc
%% in vivo aAP statistics
close all
cellnums=struct;
minimumrectime=180;
isawake=strcmp({xlsdata.anaesthesia},'awake-ketamine');% strcmp({xlsdata.anaesthesia},'awake-isoflurane');% | 
neededcells=isawake & [xlsdata.field]==0 & [xlsdata.recordinglength]>=minimumrectime;
cellnums.all=sum(neededcells);
cellnums.with_aAP=sum(neededcells&[xlsdata.axonalAPnum]>0);
cellnums.without_aAP=sum(neededcells&[xlsdata.axonalAPnum]==0);
cellnums.with_persistent_aAP=sum(neededcells&[xlsdata.axonalAPnum_persistent]>0 & [xlsdata.axonalAPnum_sporadic]==0);
cellnums.with_sporadic_aAP=sum(neededcells&[xlsdata.axonalAPnum_sporadic]>0 & [xlsdata.axonalAPnum_persistent]==0);
cellnums.with_sporadic_and_persistent_aAP=sum(neededcells&[xlsdata.axonalAPnum_sporadic]>0 & [xlsdata.axonalAPnum_persistent]>0);
piechartorder={'without_aAP','with_persistent_aAP','with_sporadic_and_persistent_aAP','with_sporadic_aAP'};
piechartstring={};
cellnumvec=zeros(length(piechartorder),1);
for i=1:length(piechartorder)
    cellnumvec(i)=cellnums.(piechartorder{i});
    piechartstring{i} = [strrep(piechartorder{i},'_',' '),' (n = ',num2str(cellnumvec(i)),')'];
end
if sum(cellnumvec)~=cellnums.all
    disp('ERROR!! - the cell numbers doesn''t add up in the piechart!!')
end
figure(8)
clf
hist([xlsdata(neededcells).recordinglength]/60,[.5:1:80])
xlim([0 60])
title('recording length')
xlabel('time (min)')
ylabel('count')
saveas(gcf,[dirs.figuresdir,'invivo_awake_recording_length'],'pdf')
saveas(gcf,[dirs.figuresdir,'invivo_awake_recording_length'],'jpg')
figure(9)
clf
pie(cellnumvec,piechartstring);
title(['aAP occurence in n=',num2str(cellnums.all),' cells'])
saveas(gcf,[dirs.figuresdir,'invivo_awake_aAP_occurence'],'pdf')
saveas(gcf,[dirs.figuresdir,'invivo_awake_aAP_occurence'],'jpg')
cellswithsporadicaAP=neededcells&[xlsdata.axonalAPnum_sporadic]>0;
[aAPcounts_sporadic,aAPbincenters]=hist([xlsdata(cellswithsporadicaAP).axonalAPnum_sporadic],.5:1:30.5);
figure(10)
bar(aAPbincenters,aAPcounts_sporadic)
xlabel('aAP count')
ylabel('number of cells')
title('aAP numbers in experiments with sporadic aAP')
saveas(gcf,[dirs.figuresdir,'invivo_awake_aAP_num_sporadic'],'pdf')
saveas(gcf,[dirs.figuresdir,'invivo_awake_aAP_num_sporadic'],'jpg')

figure(12)
plot([xlsdata(cellswithsporadicaAP).recordinglength]/60,[xlsdata(cellswithsporadicaAP).axonalAPnum_sporadic],'ko')
xlabel('recording length (min)')
ylabel('sporadic aAP num')
saveas(gcf,[dirs.figuresdir,'invivo_awake_rec_length_vs_aAPnum'],'pdf')
saveas(gcf,[dirs.figuresdir,'invivo_awake_rec_length_vs_aAPnum'],'jpg')
cellswithspersistentaAP=neededcells&[xlsdata.axonalAPnum_persistent]>0;
[aAPcounts_persistent,aAPbincenters]=hist([xlsdata(cellswithspersistentaAP).axonalAPnum_persistent],5:10:1000);
figure(11)
bar(aAPbincenters,aAPcounts_persistent)
xlabel('aAP count')
ylabel('number of cells')
title('aAP numbers in experiments with persistent aAP')

%% ISI, packet, etc
% neededIDs={''}
minimumrectime=180;
isawake=strcmp({xlsdata.anaesthesia},'awake-ketamine');% strcmp({xlsdata.anaesthesia},'awake-isoflurane');% | 
neededcells=find(isawake & [xlsdata.field]==0 & [xlsdata.recordinglength]>=minimumrectime);
aAPisis=[];
sAPisis=[];
for xlsi=1:length(xlsdata)
    if ~isempty(xlsdata(xlsi).aAPtimes) & length(xlsdata(xlsi).aAPtimes)>1 & any(xlsi==neededcells)
        aAPisis=[aAPisis,diff(xlsdata(xlsi).aAPtimes)];
        sAPisis=[sAPisis,diff(xlsdata(xlsi).sAPtimes)];
    end
end
figure(7)
clf
subplot(2,1,1)
[n,bincenters]=hist(aAPisis,[0:1:600]);
bar(bincenters,n)
% hold on
% plot(bincenters,sum(n)-cumsum(n))
set(gca,'yscale','log')
% set(gca,'xscale','log')
% ylim([0 100])
xlim([0 120])
xlabel('ISI (s)')
title('aAP ISI distribution')
subplot(2,1,2)
[n,bincenters]=hist(sAPisis,[0:1:600]);
bar(bincenters,n)
% hold on
% plot(bincenters,sum(n)-cumsum(n))
set(gca,'yscale','log')
% set(gca,'xscale','log')
% ylim([0 100])
xlim([0 120])
xlabel('ISI (s)')
title('sAP ISI distribution')
%
figure(6)
subplot(2,2,1)
packetnums=[xlsdata(neededcells).aAP_packet_num];
packetnums(packetnums==0)=[];
hist(packetnums,[1:20])
xlabel('packet number in each recording')
ylabel('cell count')
title('Packet statistics')
subplot(2,2,2)
packetaAPnum=[xlsdata(neededcells).aAP_packet_APnum];
hist(packetaAPnum,[1:2:50]);
axis tight
xlim([0 50])
xlabel('aAP number in each packet')
ylabel('packet count')
subplot(2,2,3)
packetlength=[xlsdata(neededcells).aAP_packet_length];
packetlength(packetlength==0)=[];
hist(packetlength,[0:120]);
axis tight
xlim([0 120])
xlabel('packet length (s)')
ylabel('packet count')
subplot(2,2,4)
packetlength=[xlsdata(neededcells).aAP_packet_length];
packetaAPnum=[xlsdata(neededcells).aAP_packet_APnum];
packetfreq=packetaAPnum./packetlength;
plot(packetlength,packetaAPnum,'ko')
xlabel('Packet length (s)')
ylabel('aAP number in each packet')
set(gca,'xscale','log')
set(gca,'yscale','log')

%%
packetapmaxtimes={xlsdata(neededcells).aAP_packet_APmaxtimes};
packetapnums=cellfun(@length,packetapmaxtimes);
uniquepacketapnums=unique(packetapnums);
figure(5)
clf
packetaAPnum=[xlsdata(neededcells).aAP_packet_APnum];


%%
cellcoord_lat=[xlsdata(neededcells).Cranio_Lat]+(([xlsdata(neededcells).locationX]-[xlsdata(neededcells).Cranio_center_X])/1000).*[xlsdata(neededcells).Lateral_dir_X];
cellcoord_AP=[xlsdata(neededcells).Cranio_AP]+(([xlsdata(neededcells).locationy]-[xlsdata(neededcells).Cranio_center_Y])/1000).*[xlsdata(neededcells).Rostral_dir_Y];
cranio_lat=[xlsdata(neededcells).Cranio_Lat];
cranio_AP=[xlsdata(neededcells).Cranio_AP];
figure(13)
plot(cellcoord_lat,cellcoord_AP,'ko')
figure(14)
hist3([cranio_lat',cranio_AP'],[5,5]);
xlabel('lateral (mm)')
xlim([0 2])
ylabel('anterior (mm)')
% hist(cranio_AP)
% plot(cranio_lat,cranio_AP,'ko')

%%
figure(11)
clf
plot([xlsdata(cellswithsporadicaAP).recordinglength]/60,[xlsdata(cellswithsporadicaAP).axonalAPnum_sporadic],'ko')
%%





[aAPcounts,aAPbincenters]=hist([xlsdata(neededcells).axonalAPnum],0:1:2000);
[aAPcounts_sporadic,aAPbincenters]=hist([xlsdata(neededcells).axonalAPnum_sporadic],0:1:2000);
[aAPcounts_persistent,aAPbincenters]=hist([xlsdata(neededcells).axonalAPnum_persistent],0:1:2000);
aAPcounts(1)=[];
aAPcounts_sporadic(1)=[];
aAPcounts_persistent(1)=[];
aAPbincenters(1)=[];
cumaAPcounts=cumsum(aAPcounts);%/sum(aAPcounts);
cumaAPcounts_sporadic=cumsum(aAPcounts_sporadic);%/sum(aAPcounts_sporadic);
cumaAPcounts_persistent=cumsum(aAPcounts_persistent);%/sum(aAPcounts_persistent);
figure(10)
clf
plot(aAPbincenters,cumaAPcounts,'k-','LineWidth',3)
hold on
plot(aAPbincenters,cumaAPcounts_sporadic,'g-','LineWidth',3)
plot(aAPbincenters,cumaAPcounts_persistent,'r-','LineWidth',3)
ylim([0 max(cumaAPcounts)])
xlabel('aAP num')
ylabel('number of cells')





%% plotting the location of the craniotomies
cellcoord_lat=[xlsdata.Cranio_Lat]+(([xlsdata.locationX]-[xlsdata.Cranio_center_X])/1000).*[xlsdata.Lateral_dir_X];
cellcoord_AP=[xlsdata.Cranio_AP]+(([xlsdata.locationy]-[xlsdata.Cranio_center_Y])/1000).*[xlsdata.Rostral_dir_Y];
for i=1:length(xlsdata)
    xlsdata(i).cellcoord_lat=cellcoord_lat(i);
    xlsdata(i).cellcoord_AP=cellcoord_AP(i);
end
searchfields={'anaesthesia'};
searchstrs={'awake'};
ezaz=true(size(xlsdata));
for j=1:length(searchfields)
    searchfield=searchfields{j};
    searchstr=searchstrs{j};
    for i=1:length(xlsdata)
        if isnumeric(xlsdata(i).(searchfield))
            if ~xlsdata(i).(searchfield)==searchstr
                ezaz(i)=false;
            end
        else
            if ~any(regexp(xlsdata(i).(searchfield),searchstr))
                ezaz(i)=false;
            end
        end
    end
end
needed=find(ezaz);
for xlsi=1:length(needed)
    xlsidx=needed(xlsi);
    if ~isnan(xlsdata(xlsidx).Cranio_AP)
        
    end
end

figure(3)
clf
% plot([xlsdata(needed).Cranio_Lat],[xlsdata(needed).Cranio_AP],'ro')
hold on
needed=find(ezaz & [xlsdata.field]==1);
plot3([xlsdata(needed).cellcoord_lat],[xlsdata(needed).cellcoord_AP],-[xlsdata(needed).locationz],'bo')
needed=find(ezaz & [xlsdata.field]==0);
plot3([xlsdata(needed).cellcoord_lat],[xlsdata(needed).cellcoord_AP],-[xlsdata(needed).locationz],'rx')
xlabel('Lateral')
ylabel('Anterior')
zlabel('depth')

figure(4)
clf
needed=find(ezaz & [xlsdata.field]==0); 
plot3([xlsdata(needed).cellcoord_lat],[xlsdata(needed).cellcoord_AP],[xlsdata(needed).axonalAPnum],'bo')
xlabel('Lateral')
ylabel('Anterior')
zlabel('aAPnum')

figure(5)
beteg=[xlsdata.axonalAPfreq]>.09;
jo=[xlsdata.axonalAPfreq]<.09 ;
clf
plot3([xlsdata(beteg).cellcoord_lat],[xlsdata(beteg).cellcoord_AP],[xlsdata(beteg).axonalAPnum],'ro')
hold on
plot3([xlsdata(jo).cellcoord_lat],[xlsdata(jo).cellcoord_AP],[xlsdata(jo).axonalAPnum],'bo')
xlabel('Lateral')
ylabel('Anterior')
zlabel('aAPnum')
title('Persistent and spoarid APs')

figure(6)
clf
semilogy([xlsdata(beteg).Cranio_AP],[xlsdata(beteg).axonalAPnum],'ro')
hold on
semilogy([xlsdata(jo).Cranio_AP],[xlsdata(jo).axonalAPnum],'bo')
xlabel('Cranio location (anterior)')
ylabel('aAPnum')

%%
%% plotting aAP freqencies


pfinductionidx = find([xlsdata.stimulatedsAPnum]>400);
etalonidx=find(strcmp({xlsdata.ID},'1705301rm_1_2_1'));
xvals=1./[xlsdata.aAPisi_med];
yvals=[xlsdata.axonalAPnum];
figure(11)
clf
loglog(xvals,yvals,'ko')
hold on
loglog(xvals(pfinductionidx),yvals(pfinductionidx),'ro')
loglog(xvals(etalonidx),yvals(etalonidx),'bo')
xlabel('median instantenous aAP frequency')
ylabel('aAP number')
[x,y]=ginput(1);
distances=pdist([x,xvals;y,yvals]');
[~,idx]=nanmin(distances(1:length(xlsdata)));
disp([xlsdata(idx).ID,'    anterior:',num2str(xlsdata(idx).Cranio_AP),'    lateral:',num2str(xlsdata(idx).Cranio_Lat)])
loglog(xvals(idx),yvals(idx),'rx')
%%
isawake=strcmp({xlsdata.anaesthesia},'awake-isoflurane') | strcmp({xlsdata.anaesthesia},'awake-ketamine');
needed = ~isnan([xlsdata.sAPisi_med]) & ~isnan([xlsdata.aAPisi_med]) & [xlsdata.axonalAPnum]>5 & isawake;
figure(12)
clf
loglog(1./[xlsdata(needed).sAPisi_med],1./[xlsdata(needed).aAPisi_med],'ko')
hold on
betegidx=([xlsdata.axonalAPnum]>90);
loglog(1./[xlsdata(betegidx&needed).sAPisi_med],1./[xlsdata(betegidx&needed).aAPisi_med],'ro')
xlabel('median instantenous sAP frequency')
ylabel('median instantenous aAP frequency')
%% aAP timing 
isawake=strcmp({xlsdata.anaesthesia},'awake-isoflurane') | strcmp({xlsdata.anaesthesia},'awake-ketamine');
timestep=1;
hosszak=[3,5,10,15,20];

figure(33)
clf

for i=1:length(hosszak)
    needed= [xlsdata.axonalAPnum]>0 & [xlsdata.axonalAPnum]<90 & isawake &[xlsdata.recordinglength]>hosszak(i)*60;
    neededidx=find(needed);
    ns=[];
    
    subplot(length(hosszak),2,i*2-1)
    hist([xlsdata(needed).aAPtimes]/60,[timestep/2:timestep:40])
    axis tight
    xlim([0 hosszak(i)])
    
    for ii=1:length(neededidx)
        [ns(ii,:),~]=hist([xlsdata(neededidx(ii)).aAPtimes]/60,[timestep/2:timestep:40]);
    end
    subplot(length(hosszak),2,i*2)
    bar([timestep/2:timestep:40],mean(ns));
%     hold on
%     errorbar([timestep/2:timestep:40],mean(ns),std(ns)/sqrt(length(neededidx)),'k.','LineWidth',2);
    xlim([0 hosszak(i)])
    ylim([0 0.6])
    annotation('textbox',[.9 1-(i/length(hosszak)*.85) .1 .05],'String',['n = ',num2str(length(neededidx))],'EdgeColor','none')
end
subplot(length(hosszak),2,i*2)
xlabel('Time from obtaining whole cell (min)')
subplot(length(hosszak),2,i*2-1)
xlabel('Time from obtaining whole cell (min)')
subplot(length(hosszak),2,1)
title('registered aAP number in each minute')
subplot(length(hosszak),2,2)
title('average aAP number per minute')
%% aAP brainstate dependence
aAPnumbins=[0:1:60];
timespentbins=[-600:60:1200];
figure(34)
clf
figure(33)
clf

brainstatenames_varname={'Slow_wave_sleep', 'REM_sleep', 'Quiet_wakefulness','Active_wakefulness'};
brainstatenames={'Slow wave sleep', 'REM sleep', 'Quiet wakefulness','Active wakefulness'};
needed= [xlsdata.axonalAPnum]>0 & [xlsdata.axonalAPnum]<90 & isawake;
for brainstatei=1:length(brainstatenames_varname)
    brainstatename_var=brainstatenames_varname{brainstatei};
    brainstatename=brainstatenames{brainstatei};
    figure(34)
    subplot(length(brainstatenames_varname)+1,1,brainstatei);
    hist([xlsdata(needed).([brainstatename_var,'_aAPnum'])],aAPnumbins);
    title(brainstatename)
    figure(33)
    subplot(length(brainstatenames_varname)+1,1,brainstatei);
    hist([xlsdata(needed).([brainstatename_var,'_timespent'])],timespentbins);
    title(brainstatename)
end
figure(34)
subplot(length(brainstatenames_varname)+1,1,brainstatei+1);
hist([xlsdata(needed).brainstateless_aAPnum],aAPnumbins);
xlabel('aAP num per recording in each brain state')
title('brainstateless')

figure(33)
subplot(length(brainstatenames_varname)+1,1,brainstatei+1);
hist([xlsdata(needed).brainstateless_timespent],timespentbins);
xlabel('time spent in each brain state per recording')
title('brainstateless')

     
     
