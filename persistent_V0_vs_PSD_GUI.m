function persistent_V0_vs_PSD_GUI(dirs,xlsdata,xlsnum,valtozok)
a=dir([dirs.v0_vs_PSDdir,xlsdata(xlsnum).ID,'.mat']);
if isempty(a)
    if nargin<4
        valtozok.export=struct;
        valtozok.export.movingvindowsize=3;
        valtozok.export.movingvindowstep=.5; %seconds for downsampling and median filtering
        valtozok.export.PSDonfield=false;
        valtozok.export.minsweeplength=valtozok.export.movingvindowsize/2;

        valtozok.plot.v0_start=-90;
        valtozok.plot.v0_step=0.01;
        valtozok.plot.v0_end=-40;
        valtozok.plot.v0_window=1;
        valtozok.plot.time_start=0;
        valtozok.plot.time_end=999999;
        valtozok.plot.freq_start=0;
        valtozok.plot.freq_end=10;
    end
    [v0_PSD_data,valtozok.export] = aE_V0_vs_PSD(dirs,xlsdata,xlsnum,valtozok.export);
else
    load([dirs.v0_vs_PSDdir,xlsdata(xlsnum).ID]);
    disp('loadin v0 PSD data!')
    if ~isfield(valtozok,'plot')
        valtozok.plot.v0_start=-90;
        valtozok.plot.v0_step=0.01;
        valtozok.plot.v0_end=-40;
        valtozok.plot.v0_window=1;
        valtozok.plot.time_start=0;
        valtozok.plot.time_end=999999;
        valtozok.plot.freq_start=0;
        valtozok.plot.freq_end=10;
    end
end
data.dirs=dirs;
data.xlsdata=xlsdata;
data.xlsnum=xlsnum;
data.valtozok=valtozok;
data.v0_PSD_data=v0_PSD_data;
data.valtozokplot.time_start=round(min([v0_PSD_data.time])*10)/10;
data.valtozokplot.time_end=round(max([v0_PSD_data.time])*10)/10;
V0percentilesneeded=v0_PSD_data(1).v0_percentiles;
data.handles.master=figure('Visible','off','Position',[0,0,1000,600]);
data.handles.axes1=axes('Position',[.05,.55, .35 ,.4]);
data.handles.axes2=axes('Position',[.45,.55, .35 ,.4]);
data.handles.axes3=axes('Position',[.05,.05, .35 ,.4]);
data.handles.axes4=axes('Position',[.45,.05, .35 ,.4]);
spacing=.06;
toploc=.95;
data.handles.controls_containter=uipanel('Title','Plot details','Units','normalized','Position',[.85 .5 .15 .45]);
uicontrol(data.handles.controls_containter,'Style','text','String','V0 percentile','Units','normalized','Position',[.05,toploc,.9,spacing]);
data.handles.percentileselector=uicontrol(data.handles.controls_containter,'Style','popupmenu','String',V0percentilesneeded,'Units','normalized','Position',[.05,toploc-spacing,.9,spacing],'Callback',@analandplot_Callback);
uicontrol(data.handles.controls_containter,'Style','text','String','function done on PSD','Units','normalized','Position',[.05,toploc-spacing*2,.9,spacing]);
data.handles.PSD_function_selector=uicontrol(data.handles.controls_containter,'Style','popupmenu','String',{'max','mean','median'},'Units','normalized','Position',[.05,toploc-spacing*3,.9,spacing],'Callback',@analandplot_Callback);
data.handles.Zscoreselector=uicontrol(data.handles.controls_containter,'Style','checkbox','String','Zscore','Units','normalized','Position',[.05,toploc-spacing*4,.9,spacing],'Callback',@analandplot_Callback);
uicontrol(data.handles.controls_containter,'Style','text','String','brainstate','Units','normalized','Position',[.05,toploc-spacing*5,.9,spacing]);
if isempty(v0_PSD_data(1).statenames)
    brainstates='All';
else
    brainstates=v0_PSD_data(1).statenames;
end
data.handles.PSD_brainstate_selector=uicontrol(data.handles.controls_containter,'Style','popupmenu','String',brainstates,'Units','normalized','Position',[.05,toploc-spacing*6,.9,spacing],'Callback',@analandplot_Callback);
uicontrol(data.handles.controls_containter,'Style','text','String','V0 start','Units','normalized','Position',[.0,toploc-spacing*7,.25,spacing]);
uicontrol(data.handles.controls_containter,'Style','text','String','V0 step','Units','normalized','Position',[.25,toploc-spacing*7,.25,spacing]);
uicontrol(data.handles.controls_containter,'Style','text','String','V0 end','Units','normalized','Position',[.5,toploc-spacing*7,.25,spacing]);
uicontrol(data.handles.controls_containter,'Style','text','String','V0 window','Units','normalized','Position',[.75,toploc-spacing*7,.25,spacing]);
data.handles.V0_start=uicontrol(data.handles.controls_containter,'Style','edit','String','-90','Units','normalized','Position',[.00,toploc-spacing*8,.25,spacing],'Callback',@analandplot_Callback);
data.handles.V0_step=uicontrol(data.handles.controls_containter,'Style','edit','String','0.01','Units','normalized','Position',[.25,toploc-spacing*8,.25,spacing],'Callback',@analandplot_Callback);
data.handles.V0_end=uicontrol(data.handles.controls_containter,'Style','edit','String','-60','Units','normalized','Position',[.50,toploc-spacing*8,.25,spacing],'Callback',@analandplot_Callback);
data.handles.V0_window=uicontrol(data.handles.controls_containter,'Style','edit','String','1','Units','normalized','Position',[.75,toploc-spacing*8,.25,spacing],'Callback',@analandplot_Callback);
uicontrol(data.handles.controls_containter,'Style','text','String','Start time','Units','normalized','Position',[.0,toploc-spacing*9,.5,spacing]);
uicontrol(data.handles.controls_containter,'Style','text','String','End time','Units','normalized','Position',[.5,toploc-spacing*9,.5,spacing]);
data.handles.time_start=uicontrol(data.handles.controls_containter,'Style','edit','String',num2str(round(min([v0_PSD_data.time])*10)/10),'Units','normalized','Position',[0,toploc-spacing*10,.5,spacing],'Callback',@analandplot_Callback);
data.handles.time_end=uicontrol(data.handles.controls_containter,'Style','edit','String',num2str(round(max([v0_PSD_data.time])*10)/10),'Units','normalized','Position',[.5,toploc-spacing*10,.5,spacing],'Callback',@analandplot_Callback);
uicontrol(data.handles.controls_containter,'Style','text','String','Start freq','Units','normalized','Position',[.0,toploc-spacing*11,.5,spacing]);
uicontrol(data.handles.controls_containter,'Style','text','String','End freq','Units','normalized','Position',[.5,toploc-spacing*11,.5,spacing]);
data.handles.freq_start=uicontrol(data.handles.controls_containter,'Style','edit','String',num2str(min(v0_PSD_data(1).frequencyVector)),'Units','normalized','Position',[0,toploc-spacing*12,.5,spacing],'Callback',@analandplot_Callback);
data.handles.freq_end=uicontrol(data.handles.controls_containter,'Style','edit','String',num2str(max(v0_PSD_data(1).frequencyVector)),'Units','normalized','Position',[.5,toploc-spacing*12,.5,spacing],'Callback',@analandplot_Callback);
data.handles.plotbutton=uicontrol(data.handles.controls_containter,'Style','pushbutton','String','Plot!','Units','normalized','Position',[.05,toploc-spacing*13,.9,spacing],'Callback',@analandplot_Callback_pushbotton);

data.handles.export_containter=uipanel('Title','Exporting','Units','normalized','Position',[.85 0 .15 .45]);
uicontrol(data.handles.export_containter,'Style','text','String','window step','Units','normalized','Position',[.0,toploc-spacing,.33,spacing]);
uicontrol(data.handles.export_containter,'Style','text','String','window size','Units','normalized','Position',[.33,toploc-spacing,.33,spacing]);
uicontrol(data.handles.export_containter,'Style','text','String','min sweep length','Units','normalized','Position',[.66,toploc-spacing,.33,spacing]);
data.handles.export_movingvindowstep=uicontrol(data.handles.export_containter,'Style','edit','String','0.5','Units','normalized','Position',[0,toploc-spacing*2,.33,spacing]);
data.handles.export_movingvindowsize=uicontrol(data.handles.export_containter,'Style','edit','String','2','Units','normalized','Position',[.33,toploc-spacing*2,.33,spacing]);
data.handles.export_minsweeplength=uicontrol(data.handles.export_containter,'Style','edit','String','1','Units','normalized','Position',[.66,toploc-spacing*2,.33,spacing]);
data.handles.export_PSDonfield=uicontrol(data.handles.export_containter,'Style','checkbox','String','field PSD','Units','normalized','Position',[.05,toploc-spacing*3,.9,spacing]);
data.handles.plotbutton=uicontrol(data.handles.export_containter,'Style','pushbutton','String','Reexport','Units','normalized','Position',[.05,toploc-spacing*5,.9,spacing],'Callback',@export_Callback);

data.handles.plotbutton=uicontrol('Style','pushbutton','String','Save','Units','normalized','Position',[.85,.02,.15,.02],'Callback',@savedata_Callback);
data.handles.master.Visible='on';
guidata(data.handles.master,data);
updategui(data);

end

function updategui(data)
valtozok=data.valtozok;
data.handles.export_movingvindowsize.String=valtozok.export.movingvindowsize;
data.handles.export_movingvindowstep.String=valtozok.export.movingvindowstep;
data.handles.export_minsweeplength.String=valtozok.export.minsweeplength;
data.handles.export_PSDonfield.Value=valtozok.export.PSDonfield;
data.handles.V0_start.String=valtozok.plot.v0_start;
data.handles.V0_end.String=valtozok.plot.v0_end;
data.handles.V0_step.String=valtozok.plot.v0_step;
data.handles.V0_window.String=valtozok.plot.v0_window;
data.handles.time_start.String=valtozok.plot.time_start;
data.handles.time_end.String=valtozok.plot.time_end;
data.handles.freq_start.String=valtozok.plot.freq_start;
data.handles.freq_end.String=valtozok.plot.freq_end;
end

function analandplot(data)
%%
percentileidx=data.handles.percentileselector.Value;

v0_PSD_data=data.v0_PSD_data;

icVs=[v0_PSD_data.icVs];
icV=icVs(percentileidx,:);
time=[v0_PSD_data.time];
stateidx=[v0_PSD_data.stateidx];

functionneeded=data.handles.PSD_function_selector.String;
functionneeded=functionneeded{data.handles.PSD_function_selector.Value};
brainstateneeded=data.handles.PSD_brainstate_selector.String;
if ~ischar(brainstateneeded)
    brainstateneeded=brainstateneeded{data.handles.PSD_brainstate_selector.Value};
end
switch functionneeded
    case 'max'
        PSD_to_use=[v0_PSD_data.PSD_max];
    case 'mean'
        PSD_to_use=[v0_PSD_data.PSD_mean];
    case 'median'
        PSD_to_use=[v0_PSD_data.PSD_median];
        
end
%%
if data.handles.Zscoreselector.Value% & ~isempty(data.PSDdata_stats)
    %     PSDdata_stats=data.PSDdata_stats;
    %     neededwaves=[PSDdata_stats.length]>1;
    %         mus=[PSDdata_stats(neededwaves).mu];%
    % %         mus=mus(neededfrequencies,:);%(neededwaves)
    %         mu=median(mus,2);
    %         sigmas=[PSDdata_stats(neededwaves).sigma];%(neededwaves)
    % %         sigmas=sigmas(neededfrequencies,:);
    %         sigma=median(sigmas,2);
    %         bigmatrix_Z=PSD_to_use-repmat(mu,1,size(PSD_to_use,2));
    %         bigmatrix_Z=bigmatrix_Z./repmat(sigma,1,size(PSD_to_use,2));
    %         bigmatrix_Z(PSD_to_use==0)=0;
    %         PSD_to_use=bigmatrix_Z;
    PSD_to_use=zscore(PSD_to_use,0,2);
    
end
%%

v0bins=[str2double(data.handles.V0_start.String):str2double(data.handles.V0_step.String):str2double(data.handles.V0_end.String)]/1000;
binsize=str2double(data.handles.V0_window.String)/1000;

% PSD_to_use=zscore(PSD_median,0,2);


V0_PSD_bindata=struct;
frequencyVector=v0_PSD_data(1).frequencyVector;
statenames=v0_PSD_data(1).statenames;
valtozok.timeborders=[str2num(data.handles.time_start.String),str2num(data.handles.time_end.String)];

if isempty(statenames)
    eddig=0;
else
    eddig=length(statenames);
end
for statei=0:eddig
    if statei==0
        needed=stateidx<inf;
    else
        needed=stateidx==statei;
    end
    if ~isempty(valtozok.timeborders) & diff(valtozok.timeborders)>0
        needed=needed&time>min(valtozok.timeborders) & time<max(valtozok.timeborders);
    end
    powervals=nan(length(frequencyVector),length(v0bins));
    powervalsSD=nan(length(frequencyVector),length(v0bins));
    n=nan(size(v0bins));
    for bini=1:length(v0bins)
        idx=icV>=v0bins(bini)-binsize/2 & icV<v0bins(bini)+binsize/2;
        powervals(:,bini)=nanmean(PSD_to_use(:,idx&needed),2);
        powervalsSD(:,bini)=nanstd(PSD_to_use(:,idx&needed),1,2);
        n(bini)=sum(idx&needed);
    end
    if statei==0
        V0_PSD_bindata(statei+1).statename='All';
    else
        V0_PSD_bindata(statei+1).statename=statenames{statei};
    end
    V0_PSD_bindata(statei+1).v0bins=v0bins;
    V0_PSD_bindata(statei+1).binsize=binsize;
    V0_PSD_bindata(statei+1).powervals=powervals;
    V0_PSD_bindata(statei+1).powervalsSD=powervalsSD;
end

stateidx=find(strcmp({V0_PSD_bindata.statename},brainstateneeded));
v0bins=V0_PSD_bindata(stateidx).v0bins;
binsize=V0_PSD_bindata(stateidx).binsize;
powervals=V0_PSD_bindata(stateidx).powervals;
powervalsSD=V0_PSD_bindata(stateidx).powervalsSD;
startfreq=str2num(data.handles.freq_start.String);
endfreq=str2num(data.handles.freq_end.String);


axes(data.handles.axes1)
cla
imagesc(v0bins*1000,frequencyVector,powervals)
set(gca,'YDir','normal');
colormap linspecer
xlabel('Voltage (mV)')
ylabel('Frequency (Hz)')
axis tight
xlimits=get(gca,'Xlim');
hold on
plot(xlimits,[startfreq,startfreq],'k-','LineWidth',2)
plot(xlimits,[endfreq,endfreq],'k-','LineWidth',2)
title('mean')
axes(data.handles.axes2)
cla
imagesc(v0bins*1000,frequencyVector,powervalsSD)
set(gca,'YDir','normal');
colormap linspecer
xlabel('Voltage (mV)')
ylabel('Frequency (Hz)')
axis tight
xlimits=get(gca,'Xlim');
hold on
plot(xlimits,[startfreq,startfreq],'k-','LineWidth',2)
plot(xlimits,[endfreq,endfreq],'k-','LineWidth',2)
title('sd')

axes(data.handles.axes3)
bar(v0bins*1000,n);
xlim(xlimits);
xlabel('Voltage (mV)')
ylabel('# of points averaged')
title(V0_PSD_bindata(stateidx).statename)% data = guidata(hObject);
% analandplot(data);


axes(data.handles.axes4)
neededfreqidx=frequencyVector>=startfreq &frequencyVector<=endfreq;
shadedErrorBar(v0bins*1000,nanmean(powervals(neededfreqidx,:),1),nanmean(powervalsSD(neededfreqidx,:),1));
end

function analandplot_Callback(hObject,eventdata)
%not plotting at every click
% data = guidata(hObject);
% analandplot(data);
end
function analandplot_Callback_pushbotton(hObject,eventdata)
data = guidata(hObject);
analandplot(data);
end

function savedata_Callback(hObject,eventdata)
data = guidata(hObject);
data=updatevaltozok(data);
v0_PSD_data=data.v0_PSD_data;
valtozok=data.valtozok;
save([data.dirs.v0_vs_PSDdir,data.xlsdata(data.xlsnum).ID],'v0_PSD_data','valtozok');
disp('V0-PSD data saved')
end

function data=updatevaltozok(data)
valtozok=data.valtozok;
valtozok.export.movingvindowsize=str2num(data.handles.export_movingvindowsize.String);
valtozok.export.movingvindowstep=str2num(data.handles.export_movingvindowstep.String);
valtozok.export.minsweeplength=str2num(data.handles.export_minsweeplength.String);
valtozok.export.PSDonfield=data.handles.export_PSDonfield.Value;
valtozok.plot.v0_start=str2num(data.handles.V0_start.String);
valtozok.plot.v0_end=str2num(data.handles.V0_end.String);
valtozok.plot.v0_step=str2num(data.handles.V0_step.String);
valtozok.plot.v0_window=str2num(data.handles.V0_window.String);
valtozok.plot.time_start=str2num(data.handles.time_start.String);
valtozok.plot.time_end=str2num(data.handles.time_end.String);
valtozok.plot.freq_start=str2num(data.handles.freq_start.String);
valtozok.plot.freq_end=str2num(data.handles.freq_end.String);
data.valtozok=valtozok;
end

function export_Callback(hObject,eventdata)
data = guidata(hObject);
data=updatevaltozok(data);
valtozok=data.valtozok;
[v0_PSD_data,valtozok] = aE_V0_vs_PSD(data.dirs,data.xlsdata,data.xlsnum,data.valtozok.export);
data.valtozok.plot.time_start=round(min([v0_PSD_data.time])*10)/10;
data.valtozok.plot.time_end=round(max([v0_PSD_data.time])*10)/10;
data.valtozok.export=valtozok;
data.v0_PSD_data=v0_PSD_data;
updategui(data)
guidata(hObject,data);
end