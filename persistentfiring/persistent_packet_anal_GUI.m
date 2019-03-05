function persistent_packet_anal_GUI(dirs,xlsdata)
data.dirs=dirs;
data.xlsdata=xlsdata;
data.xlsidx=0;
[~,order]=sort([xlsdata.aAP_packet_num],'descend');%axonalAPnum
aAPnum=[xlsdata(order).aAP_packet_num];%axonalAPnum
data.IDs={xlsdata(order).ID};
data.IDs=[data.IDs(~isnan(aAPnum)),data.IDs(isnan(aAPnum))];
aAPnum=[aAPnum(~isnan(aAPnum)),aAPnum(isnan(aAPnum))];
IDstring=data.IDs;
for i=1:length(data.IDs)
    IDstring{i}=[IDstring{i},' - ',num2str(aAPnum(i)),' packets'];
end
IDstring=['no file selected',IDstring];

data.handles.master=figure('Visible','off','Position',[0,0,1000,600]);
data.handles.axes1=axes('Position',[.05,.68, .6 ,.3]);
data.handles.axes2=axes('Position',[.05,.38, .6 ,.3]);
data.handles.axes3=axes('Position',[.05,.08, .6 ,.3]);

data.handles.fileselector=uicontrol('Style','popupmenu','String',IDstring,'Units','normalized','Position',[.7,.9,.3,.1],'Callback',@fileselector_Callback);
data.handles.sourceselector=uicontrol('Style','popupmenu','String',{'packets','single aAPs'},'Units','normalized','Position',[.7,.85,.3,.1],'Callback',@fileselector_Callback);
data.handles.packetselector=uicontrol('Style','listbox','String','no packets yet','max',10,'min',1,'Units','normalized','Position',[.7,.6,.3,.2],'Callback',@packetselector_Callback);


data.handles.caxis_containter=uipanel('Title','Caxis','Units','normalized','Position',[.7 .05 .3 .15]);
data.handles.caxis_lower=uicontrol(data.handles.caxis_containter,'Style','slider','SliderStep', [.01, 0.1],'Min', 0, 'Max', 100, 'Value', 0,'Units','normalized','Position',[0,.1,1,.2],'Callback', @setcaxis);
data.handles.caxis_upper=uicontrol(data.handles.caxis_containter,'Style','slider','SliderStep', [.01, 0.1],'Min', 0, 'Max', 100, 'Value', 100,'Units','normalized','Position',[0,.4,1,.2],'Callback', @setcaxis);
data.handles.caxis_lower_str=uicontrol(data.handles.caxis_containter,'Style','text','String','nan','Units','normalized','Position',[0,.6,.4,.3]);
data.handles.caxis_upper_str=uicontrol(data.handles.caxis_containter,'Style','text','String','nan','Units','normalized','Position',[0.5,.6,.4,.3]);
data.handles.PSD_baseline_containter=uipanel('Title','PSD','Units','normalized','Position',[.7 .2 .3 .2]);
data.handles.PSD_baseline_switch=uicontrol(data.handles.PSD_baseline_containter,'Style','checkbox','String','Calculate baseline','Units','normalized','Position',[.05 .3 .95 .15],'Callback', @plotthedata_Callback);
data.handles.PSD_baseline_start=uicontrol(data.handles.PSD_baseline_containter,'Style','edit','String',-20,'Units','normalized','Position',[.05 0.1 0.3 0.2],'Callback', @plotthedata_Callback);
data.handles.PSD_baseline_end=uicontrol(data.handles.PSD_baseline_containter,'Style','edit','String',-10,'Units','normalized','Position',[.4 0.1 0.3 0.2],'Callback', @plotthedata_Callback);
data.handles.PSD_averagingmethod=uicontrol(data.handles.PSD_baseline_containter,'Style','popupmenu','String',{'mean','median'},'Units','normalized','Position',[.05 .8 .35 .15],'Callback', @plotthedata_Callback);
data.handles.elfiz_baseline_containter=uipanel('Title','Baseline','Units','normalized','Position',[.7 .4 .3 .2]);
data.handles.elfiz_baseline_switch=uicontrol(data.handles.elfiz_baseline_containter,'Style','checkbox','String','Calculate baseline','Units','normalized','Position',[.05 .3 .95 .15],'Callback', @plotthedata_Callback);
data.handles.elfiz_baseline_start=uicontrol(data.handles.elfiz_baseline_containter,'Style','edit','String',-0.5,'Units','normalized','Position',[.05 0.1 0.3 0.2],'Callback', @plotthedata_Callback);
data.handles.elfiz_baseline_end=uicontrol(data.handles.elfiz_baseline_containter,'Style','edit','String',0,'Units','normalized','Position',[.4 0.1 0.3 0.2],'Callback', @plotthedata_Callback);
data.handles.plotproperties_containter=uipanel('Title','Plot properties','Units','normalized','Position',[.7 .8 .3 .10]);
data.handles.plotproperties_start=uicontrol(data.handles.plotproperties_containter,'Style','edit','String',-20,'Units','normalized','Position',[.05 0.5 0.3 0.5]); %,'Callback', @plotthedata_Callback
data.handles.plotproperties_end=uicontrol(data.handles.plotproperties_containter,'Style','edit','String',20,'Units','normalized','Position',[.4 0.5 0.3 0.5]);%,'Callback', @plotthedata_Callback

% uicontrol('Style','text','String','Caxis','Units','normalized','Position',[.7 .13 .1 .1]);
% align([handles.axes1,handles.fileselector],'Center','None')
data.handles.master.Visible='on';
guidata(data.handles.master,data);
end

function updateGUI(data)
if data.xlsidx>0
    if isfield(data.aAPvicinity,'aAPsinpacket')
        packetstring={};
        for i=1:length(data.aAPvicinity)
            aAPn=length(data.aAPvicinity(i).aAPsinpacket);
            duration=round((data.aAPvicinity(i).aAPsinpacket(end).maxtime-data.aAPvicinity(i).aAPsinpacket(1).maxtime)*1000)/1000;
            timebefore=round(data.aAPvicinity(i).timeback);
            timeafter=round(data.aAPvicinity(i).timeforward);
            si=round(data.aAPvicinity(i).si*1e6);
            packetstring{i}=['aAPs: ',num2str(aAPn),' dur: ', num2str(duration),'s ','before: ',num2str(timebefore),' after: ',num2str(timeafter)];
        end
        data.handles.packetselector.Value=1;
        data.handles.packetselector.String=packetstring;
    else
        packetstring={};
        for i=1:length(data.aAPvicinity)
            realtime=round((data.aAPvicinity(i).realtime*1000))/1000;
%             duration=round((data.aAPvicinity(i).aAPsinpacket(end).maxtime-data.aAPvicinity(i).aAPsinpacket(1).maxtime)*1000)/1000;
            timebefore=round(data.aAPvicinity(i).timeback);
            timeafter=round(data.aAPvicinity(i).timeforward);
            si=round(data.aAPvicinity(i).si*1e6);
            packetstring{i}=['time: ',num2str(realtime),'before: ',num2str(timebefore),' after: ',num2str(timeafter)];
        end
        data.handles.packetselector.Value=1;
        data.handles.packetselector.String=packetstring;
%         data.handles.packetselector.String=1:length(data.aAPvicinity);
    end
    disp('loller')
else
    data.handles.packetselector.Value=1;
    data.handles.packetselector.String='no packets';
end
end

function plotthedata_Callback(hObject,eventdata)
data=guidata(hObject);
data=plotthedata(data);
guidata(hObject,data);
end

function data=plotthedata(data)
selectedpackets=data.handles.packetselector.Value;
axes(data.handles.axes1)
cla
hold all
% for i=1:length(selectedpackets)
%     plot(data.aAPvicinity(selectedpackets(i)).time_elphys,data.aAPvicinity(selectedpackets(i)).ic)
% end

x=reshape([data.aAPvicinity(selectedpackets).time_elphys],length(data.aAPvicinity(selectedpackets(1)).time_elphys),length(selectedpackets));
y=reshape([data.aAPvicinity(selectedpackets).ic],length([data.aAPvicinity(selectedpackets(1)).ic]),length(selectedpackets));

baselineedges=sort([str2double(data.handles.elfiz_baseline_start.String),str2double(data.handles.elfiz_baseline_end.String)]);
if data.handles.elfiz_baseline_switch.Value & ~any(isnan(baselineedges));
    baselineidx=find(data.aAPvicinity(1).time_elphys>=baselineedges(1) & data.aAPvicinity(1).time_elphys<=baselineedges(2)); 
    mu=nanmean(y(baselineidx,:),1);
    y=y-repmat(mu,size(y,1),1);
end

plotedges=sort([str2double(data.handles.plotproperties_start.String),str2double(data.handles.plotproperties_end.String)]);
if ~any(isnan(plotedges));
    neededidx=find(data.aAPvicinity(1).time_elphys>=plotedges(1) & data.aAPvicinity(1).time_elphys<=plotedges(2)); 
    x=x(neededidx,:);
    y=y(neededidx,:);
end

if length(selectedpackets)>1
    shadedErrorBar(nanmean(x,2),nanmean(y,2),nanstd(y,0,2))
end
plot(nanmean(x,2),nanmean(y,2),'k-','LineWidth',1)

axis tight

axes(data.handles.axes2)
cla
hold all
% for i=1:length(selectedpackets)
%     plot(data.aAPvicinity(selectedpackets(i)).time_elphys,data.aAPvicinity(selectedpackets(i)).field)
% end
x=reshape([data.aAPvicinity(selectedpackets).time_elphys],length(data.aAPvicinity(selectedpackets(1)).time_elphys),length(selectedpackets));
y=reshape([data.aAPvicinity(selectedpackets).field],length([data.aAPvicinity(selectedpackets(1)).field]),length(selectedpackets));
baselineedges=sort([str2double(data.handles.elfiz_baseline_start.String),str2double(data.handles.elfiz_baseline_end.String)]);
if data.handles.elfiz_baseline_switch.Value & ~any(isnan(baselineedges));
    baselineidx=find(data.aAPvicinity(1).time_elphys>=baselineedges(1) & data.aAPvicinity(1).time_elphys<=baselineedges(2)); 
    mu=nanmean(y(baselineidx,:),1);
    y=y-repmat(mu,size(y,1),1);
end
plotedges=sort([str2double(data.handles.plotproperties_start.String),str2double(data.handles.plotproperties_end.String)]);
if ~any(isnan(plotedges));
    neededidx=find(data.aAPvicinity(1).time_elphys>=plotedges(1) & data.aAPvicinity(1).time_elphys<=plotedges(2)); 
    x=x(neededidx,:);
    y=y(neededidx,:);
end
if length(selectedpackets)>1
    shadedErrorBar(nanmean(x,2),nanmean(y,2),nanstd(y,0,2))
end
plot(nanmean(x,2),nanmean(y,2),'k-','LineWidth',1)

axis tight

axes(data.handles.axes3)
cla
time=data.aAPvicinity(1).time_PSD_field;
firstidx=find(~cellfun(@isempty,{data.aAPvicinity(selectedpackets).PSD_field}),1,'first');
temp=reshape([data.aAPvicinity(selectedpackets).PSD_field],[size([data.aAPvicinity(selectedpackets(firstidx)).PSD_field],1),size([data.aAPvicinity(selectedpackets(firstidx)).PSD_field],2),sum(~cellfun(@isempty,{data.aAPvicinity(selectedpackets).PSD_field}))]);%length(selectedpackets)
%%
if ~isempty(temp)
    tempp=temp;
    tempp(tempp<1e-10)=NaN;
    PSDbaselineedges=sort([str2double(data.handles.PSD_baseline_start.String),str2double(data.handles.PSD_baseline_end.String)]);
    if data.handles.PSD_baseline_switch.Value & ~any(isnan(PSDbaselineedges));
        baselineidx=find(data.aAPvicinity(1).time_PSD_field>=PSDbaselineedges(1) & data.aAPvicinity(1).time_PSD_field<=PSDbaselineedges(2));
        mu=nanmean(tempp(:,baselineidx,:),2);
        sigma=nanstd(tempp(:,baselineidx,:),[],2);
        tempp=tempp-repmat(mu,1,size(tempp,2),1);
        tempp=tempp./repmat(sigma,1,size(tempp,2),1);
        temp=tempp;
    end
    
    plotedges=sort([str2double(data.handles.plotproperties_start.String),str2double(data.handles.plotproperties_end.String)]);
    if ~any(isnan(plotedges));
        neededidx=find(time>=plotedges(1) & time<=plotedges(2));
        time=time(neededidx);
        temp=temp(:,neededidx,:);
    end
    %%
    averagingmethod=data.handles.PSD_averagingmethod.String;
    averagingmethod=averagingmethod{data.handles.PSD_averagingmethod.Value};
    switch averagingmethod
        case 'mean'
            toplot=nanmean(temp,3);
        case 'median'
            toplot=nanmedian(temp,3);
        case 'std'
            toplot=nanstd(temp,[],3);
    end
    
    imagesc(time,data.aAPvicinity(1).frequencyVector_field,toplot);
    set(gca,'YDir','normal');
    colormap linspecer
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    
    
    set(gca, 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
    ytickidxs=[];
    yticknow=get(gca,'YTick');
    ylimnow=get(gca,'Ylim');
    
    linearyticks = linspace(ylimnow(1),ylimnow(end),length(data.aAPvicinity(1).frequencyVector_field));
    for ticki=1:length(yticknow)
        [~,ytickidxs(ticki)]=min(abs(linearyticks-yticknow(ticki)));
    end
    yticklabelnow=data.aAPvicinity(1).frequencyVector_field(ytickidxs);
    set(gca,'Yticklabel',round(yticklabelnow*100)/100)
    
    linkaxes([data.handles.axes1,data.handles.axes2,data.handles.axes3],'x');
    
    % caxisline=sort(toplot(:));
    % caxisline(isnan(caxisline))=[];
    % data.caxisvals=[caxisline(round(length(caxisline)*.01)),caxisline(round(length(caxisline)*.99))];
    data.caxisvals=[nanmin(toplot(:)),nanmax(toplot(:))];
    
    
    caxismultipliers=[data.handles.caxis_lower.Value,data.handles.caxis_upper.Value]/100;
    caxis(sort([data.caxisvals(1)+caxismultipliers(1)*diff(data.caxisvals),data.caxisvals(1)+caxismultipliers(2)*diff(data.caxisvals) ]));
    data.handles.caxis_lower_str.String=min([data.caxisvals(1)+caxismultipliers(1)*diff(data.caxisvals),data.caxisvals(1)+caxismultipliers(2)*diff(data.caxisvals) ]);
    data.handles.caxis_upper_str.String=max([data.caxisvals(1)+caxismultipliers(1)*diff(data.caxisvals),data.caxisvals(1)+caxismultipliers(2)*diff(data.caxisvals) ]);
end
end



function fileselector_Callback(hObject,eventdata)
data = guidata(hObject);
selectedid = data.handles.fileselector.Value-1;%get(hObject,'Value')-1;
if data.handles.sourceselector.Value==1
    directory=[data.dirs.basedir,'Vicinity_packet/'];
elseif data.handles.sourceselector.Value==2
    directory=[data.dirs.basedir,'Vicinity_aAP/'];
end
xlsidx=find(strcmp({data.xlsdata.ID},data.IDs{selectedid}));
a=dir([directory,data.xlsdata(xlsidx).ID,'.mat']);
if selectedid>0 & ~isempty(a)
    
    data.xlsidx=xlsidx;
    disp(['loading ',data.xlsdata(xlsidx).ID]);
    hObject.Enable='off';
    pause(.3)
    load([directory,data.xlsdata(xlsidx).ID]);
    
    %%
    if isfield(aAPvicinity,'PSD_field_compress_multiplier')
        for i=1:length(aAPvicinity)
            aAPvicinity(i).PSD_field=double(aAPvicinity(i).PSD_field)*aAPvicinity(i).PSD_field_compress_multiplier+aAPvicinity(i).PSD_field_compress_offset;
        end
    end
    %%
    data.aAPvicinity=aAPvicinity;
    hObject.Enable='on';
    
     data.handles.plotproperties_start.String=(max([aAPvicinity.aAP_interval_before]))*-1;
     data.handles.plotproperties_end.String=(max([aAPvicinity.aAP_interval_after]));
else
    data.xlsidx=0;
end
updateGUI(data);
guidata(hObject,data);
end
function packetselector_Callback(hObject,eventdata)
data = guidata(hObject);
data=plotthedata(data);
guidata(hObject,data);
end
function setcaxis(hObject,eventdata)
data = guidata(hObject);
caxismultipliers=[data.handles.caxis_lower.Value,data.handles.caxis_upper.Value]/100;
axes(data.handles.axes3)
caxis(sort([data.caxisvals(1)+caxismultipliers(1)*diff(data.caxisvals),data.caxisvals(1)+caxismultipliers(2)*diff(data.caxisvals) ]));
data.handles.caxis_lower_str.String=min([data.caxisvals(1)+caxismultipliers(1)*diff(data.caxisvals),data.caxisvals(1)+caxismultipliers(2)*diff(data.caxisvals) ]);
data.handles.caxis_upper_str.String=max([data.caxisvals(1)+caxismultipliers(1)*diff(data.caxisvals),data.caxisvals(1)+caxismultipliers(2)*diff(data.caxisvals) ]);
% disp('lol')
end
