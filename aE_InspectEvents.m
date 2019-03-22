function varargout = aE_InspectEvents(varargin)
% AE_INSPECTEVENTS MATLAB code for aE_InspectEvents.fig
%      AE_INSPECTEVENTS, by itself, creates a new AE_INSPECTEVENTS or raises the existing
%      singleton*.
%
%      H = AE_INSPECTEVENTS returns the handle to a new AE_INSPECTEVENTS or the handle to
%      the existing singleton*.
%
%      AE_INSPECTEVENTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AE_INSPECTEVENTS.M with the given input arguments.
%
%      AE_INSPECTEVENTS('Property','Value',...) creates a new AE_INSPECTEVENTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before aE_InspectEvents_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to aE_InspectEvents_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help aE_InspectEvents

% Last Modified by GUIDE v2.5 05-Mar-2019 14:00:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @aE_InspectEvents_OpeningFcn, ...
                   'gui_OutputFcn',  @aE_InspectEvents_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before aE_InspectEvents is made visible.
function aE_InspectEvents_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to aE_InspectEvents (see VARARGIN)

% Choose default command line output for aE_InspectEvents
handles.output = hObject;
handles.data.bridgeddata=varargin{1};
handles.data.stimdata=varargin{2};
eventdata_orig=varargin{3};
fieldek=fieldnames(eventdata_orig);
for fieldi=1:length(fieldek)
    fieldnev=fieldek{fieldi};
    if length([eventdata_orig.(fieldnev)])<length(eventdata_orig)
        [eventdata_orig((cellfun(@isempty,{eventdata_orig.(fieldnev)}))).(fieldnev)]=deal(NaN);
    end
end
%%
if ~isfield(eventdata_orig,'RS')
    progressbar('Assigning RS to each AP.. please wait')
    for sweepi=1:length(handles.data.stimdata)
        RS=handles.data.stimdata(sweepi).RS;
        [eventdata_orig([eventdata_orig.sweepnum]==sweepi).RS]=deal(RS);
        progressbar(sweepi/length(handles.data.stimdata))
    end
end
%%
handles.data.eventdata_orig=eventdata_orig;
handles.data.eventdata=eventdata_orig;
handles.data.eventadata_selectedidxes=1:length(handles.data.eventdata_orig);
handles.data.rules=struct;
handles.data.plotdata(1).X='t';
handles.data.plotdata(1).Y='v';
handles.data.plotdata(2).X='t';
handles.data.plotdata(2).Y='dv';
handles.data.plotdata(3).X='v';
handles.data.plotdata(3).Y='dv';
handles.data.plotdata(4).X='t';
handles.data.plotdata(4).Y='ddv';

set(handles.edit4,'String','1')
set(handles.edit3,'String','1')
set(handles.edit2,'String','1')
set(handles.edit1,'String','');
set(handles.listbox2,'String','');

set(handles.popupmenu11,'String',{'','abs'});
set(handles.popupmenu10,'String',{'gaussian sigma','boxcar'});
set(handles.popupmenu9,'String',{'aAP','sAP','not AP'});
set(handles.popupmenu8,'String',{'histogram','plot VS'});
set(handles.popupmenu7,'String',fieldnames(handles.data.eventdata));
set(handles.popupmenu6,'String',{'=','~=','<','>'});
set(handles.popupmenu5,'String',fieldnames(handles.data.eventdata));
set(handles.popupmenu4,'String',[1:4]);
set(handles.popupmenu3,'String',fieldnames(handles.data.eventdata));
set(handles.listbox1,'Max',inf,'Min',0);
set(handles.popupmenu1,'String',{'t','v'});
set(handles.popupmenu2,'String',{'v','dv','ddv','dddv'});
handles = updateGUI(handles);
guidata(hObject, handles);

function handles = updateGUI(handles)
%%
selectedeventdatafield=get(handles.popupmenu3,'String');
selectedeventdatafield=selectedeventdatafield{get(handles.popupmenu3,'Value')};
if length(handles.data.eventdata)==0
    uniquevals='';
else
    if ~ischar(handles.data.eventdata(1).(selectedeventdatafield))
        uniquevals=unique([handles.data.eventdata.(selectedeventdatafield)]);
        uniquevals=mat2cell(uniquevals,1,length(uniquevals));
    else
        uniquevals=unique({handles.data.eventdata.(selectedeventdatafield)});
    end
end
set(handles.listbox1,'Value',1);
set(handles.listbox1,'String',uniquevals);

if ~ isempty(fieldnames(handles.data.rules)) & length(handles.data.rules)>0
    
    for i=1:length(handles.data.rules)
        if isnumeric(handles.data.rules(i).value)
            stringtoshow{i}=[handles.data.rules(i).selectedeventdatafield,' ',handles.data.rules(i).comparefunc,' ',num2str(handles.data.rules(i).value)];
        else
            stringtoshow{i}=[handles.data.rules(i).selectedeventdatafield,' ',handles.data.rules(i).comparefunc,' ',handles.data.rules(i).value];
        end
    end
else
    stringtoshow='no rules';
    
end
set(handles.listbox2,'String',stringtoshow);

function handles = applyrules(handles)
if ~ isempty(fieldnames(handles.data.rules))
    kell=true(size(handles.data.eventdata_orig));
    for i=1:length(handles.data.rules)
        if isnumeric(handles.data.rules(i).value)
            if strcmp(handles.data.rules(i).modifyfunc,'abs')
                valvector=abs([handles.data.eventdata_orig.(handles.data.rules(i).selectedeventdatafield)]);
            else
                valvector=[handles.data.eventdata_orig.(handles.data.rules(i).selectedeventdatafield)];
            end
            if strcmp(handles.data.rules(i).comparefunc,'=')
                kell=kell&  valvector== handles.data.rules(i).value;
            elseif strcmp(handles.data.rules(i).comparefunc,'~=')
                kell=kell& valvector ~= handles.data.rules(i).value;
            elseif strcmp(handles.data.rules(i).comparefunc,'>')
                kell=kell& valvector > handles.data.rules(i).value;
            elseif strcmp(handles.data.rules(i).comparefunc,'<')
                kell=kell& valvector < handles.data.rules(i).value;
            end
        else
            if strcmp(handles.data.rules(i).comparefunc,'=')
                kell=kell& strcmpi({handles.data.eventdata_orig.(handles.data.rules(i).selectedeventdatafield)}, handles.data.rules(i).value);
            else
                kell=kell& ~strcmpi({handles.data.eventdata_orig.(handles.data.rules(i).selectedeventdatafield)}, handles.data.rules(i).value);
            end
        end
    end
    handles.data.eventdata=handles.data.eventdata_orig(kell);
    handles.data.eventadata_selectedidxes=find(kell);
else
    handles.data.eventdata=handles.data.eventdata_orig;
    handles.data.eventadata_selectedidxes=1:length(handles.data.eventdata_orig);
end

function handles=updateeventdata(handles)
 eventdata_orig=handles.data.eventdata_orig;
%  if length(fieldnames(eventdata_orig)) ~= length(fieldnames(handles.data.eventdata))
     fieldnevek=fieldnames(handles.data.eventdata);
     for fieldi=1:length(fieldnevek)
         fieldnev=fieldnevek{fieldi};
         if ~any(strcmp(fieldnev,fieldnames(eventdata_orig)))
             [eventdata_orig.(fieldnev)] = deal(NaN);
             disp(['new field added to eventdata: ',fieldnev])
         end
     end
     set(handles.popupmenu7,'String',fieldnames(handles.data.eventdata));
     set(handles.popupmenu5,'String',fieldnames(handles.data.eventdata));
     set(handles.popupmenu3,'String',fieldnames(handles.data.eventdata));
%  end
 eventdata_orig(handles.data.eventadata_selectedidxes)=handles.data.eventdata;
 for fieldi=1:length(fieldnevek)
     fieldnev=fieldnevek{fieldi};
     if length([eventdata_orig.(fieldnev)])<length(eventdata_orig)
         [eventdata_orig((cellfun(@isempty,{eventdata_orig.(fieldnev)}))).(fieldnev)]=deal(NaN);
         disp(['field filled with NANs: ',fieldnev])
     end
 end
 handles.data.eventdata_orig=eventdata_orig;

function handles = updateDATA(handles,doitall)
%%

% valtozok.timeback=.001;
% valtozok.timeforward=.001;
% valtozok.movingn=3;
valtozok.timebackforthreshold=.002;
valtozok.timeback=str2num(get(handles.edit2,'String'))/1000;
valtozok.timeforward=str2num(get(handles.edit3,'String'))/1000;
valtozok.movingt=round(str2num(get(handles.edit4,'String')));

filtstring=get(handles.popupmenu10,'String');
filtstring=filtstring{get(handles.popupmenu10,'Value')};
if any(strfind(filtstring,'gauss'))
    valtozok.filter='gauss';
elseif any(strfind(filtstring,'boxcar'))
    valtozok.filter='boxcar';
end

if valtozok.movingt<1
    valtozok.movingt=1;
end
% if strcmp(valtozok.filter,'gauss')
%     hgauss=fspecial('gaussian',[valtozok.movingn*3,1],valtozok.movingn/6);
% end
valtozok.depolratewindow=.001;
valtozok.baselineVwindow=.005;
valtozok.baselineVwindow_long=.1;
if nargin==1 | doitall==false
    doitall=false;
    selectedeventdatafield=get(handles.popupmenu3,'String');
    selectedeventdatafield=selectedeventdatafield{get(handles.popupmenu3,'Value')};
    % selectedeventdatavalues=get(handles.listbox1,'String');
    % selectedeventdatavalues=selectedeventdatavalues(get(handles.listbox1,'Value'));
    if ~ischar(handles.data.eventdata(1).(selectedeventdatafield))
        uniquevals=unique([handles.data.eventdata.(selectedeventdatafield)]);
    else
        uniquevals=unique({handles.data.eventdata.(selectedeventdatafield)});
    end
    selectedeventdatavalues=uniquevals(get(handles.listbox1,'Value'));
    neededevents=false(size(handles.data.eventdata));
    for selectedvali=1:length(selectedeventdatavalues)
        if isnumeric(handles.data.eventdata(1).(selectedeventdatafield))
            neededevents=neededevents | [handles.data.eventdata.(selectedeventdatafield)]==selectedeventdatavalues(selectedvali);
        else
            neededevents=neededevents | strcmp({handles.data.eventdata.(selectedeventdatafield)},selectedeventdatavalues{selectedvali});
        end
    end
    neededevents=find(neededevents);
else
    disp('all APs being checked.. please wait')
    neededevents=find(strcmpi({handles.data.eventdata.type},'ap'));
    apmaxtimes=[handles.data.eventdata(neededevents).maxtime];
    apprevisis=diff([handles.data.eventdata(neededevents(1)).maxtime,handles.data.eventdata(neededevents).maxtime]);
end

APwaves=struct;
eventdata=handles.data.eventdata;
stimdata=handles.data.stimdata;
bridgeddata=handles.data.bridgeddata;


% if ~isempty(eventdata) & ~isempty(fieldnames(eventdata))
apidxes=neededevents;
prevsweepnum=NaN;
% for eventi=1:length(eventdata)
%     eventdata(eventi).maxdv=0;
% end
for apii=1:length(apidxes)
    if doitall
        progressbar(apii/length(apidxes))
    end
    api=apidxes(apii);
    sweepnum=eventdata(api).sweepnum;
    if sweepnum~=prevsweepnum
        RS=stimdata(sweepnum).RS/10^6;
        y=bridgeddata(sweepnum).y;
        si=round(bridgeddata(sweepnum).si*10^6)/10^6;
        if valtozok.movingt>1
            if strcmp(valtozok.filter,'boxcar')
                valtozok.movingn=ceil(valtozok.movingt/si/1e6);
                yfiltered=moving(y,valtozok.movingn);
            elseif strcmp(valtozok.filter,'gauss')
                valtozok.movingn=ceil(valtozok.movingt/si*10/1e6);
                hgauss=fspecial('gaussian',[valtozok.movingn,1],valtozok.movingt/si/1e6);
                yfiltered=imfilter(y,hgauss','replicate')';
            end    
        else
            yfiltered=y';
        end
        dyfiltered=diff(yfiltered)/si;
        if valtozok.movingt>1
            if strcmp(valtozok.filter,'boxcar')
                ddyfiltered=diff(moving(dyfiltered,valtozok.movingn));
                dddyfiltered=diff(moving(ddyfiltered,valtozok.movingn));
                ddddyfiltered=diff(moving(dddyfiltered,valtozok.movingn));
            elseif strcmp(valtozok.filter,'gauss')
                ddyfiltered=diff(imfilter(dyfiltered,hgauss,'replicate'));
                dddyfiltered=diff(imfilter(ddyfiltered,hgauss,'replicate'));
                ddddyfiltered=diff(imfilter(dddyfiltered,hgauss,'replicate'));
            end
        else
            ddyfiltered=diff(dyfiltered);
            dddyfiltered=diff(ddyfiltered);
            ddddyfiltered=diff(dddyfiltered);
        end
        dyfiltered=nanmean([[NaN;dyfiltered],[dyfiltered;NaN]],2);
        ddyfiltered=[ddyfiltered(1);ddyfiltered;ddyfiltered(end)];
        dddyfiltered=nanmean([[NaN;NaN;NaN;dddyfiltered],[dddyfiltered;NaN;NaN;NaN]],2);
        ddddyfiltered=[ddddyfiltered(1);ddddyfiltered(1);ddddyfiltered;ddddyfiltered(end);ddddyfiltered(end)];
        prevsweepnum=sweepnum;
        stepback_forthresh_first=round(.0002/si);
        stepback=round(valtozok.timeback/si);
        stepforward=round(valtozok.timeforward/si);
        stepback_forthresh=round(valtozok.timebackforthreshold/si);
    end
    onseth=eventdata(api).onseth;
    %%
    maxh=eventdata(api).maxh;
    elore=min(length(yfiltered)-maxh,5);
    hatra=min(maxh,5);
    [~,maxhh]=max(yfiltered([-1*hatra:elore]+maxh));
    maxh=maxhh-1*hatra-1+maxh;
    %% 
    stepback_forthresh_orig=stepback_forthresh;
    [~,threshh]=max(dyfiltered(max(1,maxh-stepback_forthresh):maxh-round(3e-5/si)));
    threshh=threshh+max(1,maxh-stepback_forthresh);
    while ((dyfiltered(threshh)>30) | (max(dyfiltered(max(1,threshh-stepback_forthresh):threshh))>30 & max(yfiltered(max(1,threshh-stepback_forthresh):threshh))<=yfiltered(threshh)))  & threshh>3 %|    )
        threshh=threshh-1;
    end
    %% if there is a spike doublet
    stepback_forthresh=round(stepback_forthresh_orig/2);
    while ((dyfiltered(threshh)>30) | (max(dyfiltered(max(1,threshh-stepback_forthresh):threshh))>30 & max(yfiltered(max(1,threshh-stepback_forthresh):threshh))<=yfiltered(threshh)))  & threshh>3 %|    )
        threshh=threshh-1;
    end
    stepback_forthresh=round(stepback_forthresh_orig/4);
    while ((dyfiltered(threshh)>30) | (max(dyfiltered(max(1,threshh-stepback_forthresh):threshh))>30 & max(yfiltered(max(1,threshh-stepback_forthresh):threshh))<=yfiltered(threshh)))  & threshh>3 %|    )
        threshh=threshh-1;
    end
    stepback_forthresh=stepback_forthresh_orig;
    %% looking for biphasic AP rise - alternative universal solution
    
    v=yfiltered(threshh:maxh);
    dv=dyfiltered(threshh:maxh);
    ddv=ddyfiltered(threshh:maxh);
    dddv=dddyfiltered(threshh:maxh);
    ddddv=ddddyfiltered(threshh:maxh);
    maxdvh=length(v)-2;
    while maxdvh>3 & any(dv(maxdvh-3:maxdvh-1)>dv(maxdvh)) 
        maxdvh=maxdvh-1;
    end
    maxddvh=maxdvh;
    while  maxddvh>3 & any(ddv(maxddvh-3:maxddvh-1)>ddv(maxddvh))
        maxddvh=maxddvh-1;
    end
    maxdddvh=maxddvh;
    while maxdddvh>3 & any(dddv(maxdddvh-3:maxdddvh-1)>dddv(maxdddvh))
        maxdddvh=maxdddvh-1;
    end
    if maxdddvh==3
        maxdddvh=maxdvh;
        while maxdddvh>3 & any(dddv(maxdddvh-3:maxdddvh-1)<dddv(maxdddvh))
            maxdddvh=maxdddvh-1;
        end
        while maxdddvh>3 & any(dddv(maxdddvh-3:maxdddvh-1)>dddv(maxdddvh))
            maxdddvh=maxdddvh-1;
        end
    end
    maxddddvh=maxdddvh;
    while maxddddvh>3 & any(ddddv(maxddddvh-3:maxddddvh-1)>ddddv(maxddddvh))
        maxddddvh=maxddddvh-1;
    end
    if maxddddvh==3
        maxddddvh=maxdvh;
        while maxddddvh>3 & any(ddddv(maxddddvh-3:maxddddvh-1)<ddddv(maxddddvh))
        maxddddvh=maxddddvh-1;
        end
        while maxddddvh>3 & any(ddddv(maxddddvh-3:maxddddvh-1)>ddddv(maxddddvh))
        maxddddvh=maxddddvh-1;
        end
    end
    
    if maxddddvh>3
        biphasic=true;
        [~,firstmaxdvh]=max(dv(1:maxddddvh));
        maxdvhs=[firstmaxdvh,maxdvh]+threshh-1;
        if maxdddvh>3
            firstpeakh=maxdddvh+threshh-1;
        else
            firstpeakh=maxddddvh+threshh-1;
        end
        firstpeakval=yfiltered(firstpeakh);
        maxddvh=maxddvh+threshh-1;
        maxdddvh=maxdddvh+threshh-1;
        maxddddvh=maxddddvh+threshh-1;
    else
        biphasic=false;
        maxdvhs=maxdvh+threshh-1;
        firstpeakh=NaN;
        firstpeakval=eventdata(api).maxval;%yfiltered(threshh);
        maxddvh=threshh;
        maxdddvh=threshh;
        maxddddvh=threshh;
    end
    eventdata(api).phasedelay=diff(maxdvhs)*si;
    eventdata(api).phasevoltagediff=diff(yfiltered(maxdvhs));
    if isempty(eventdata(api).phasedelay) %| eventdata(api).phasedelay< si*5
        eventdata(api).phasedelay=0;
        eventdata(api).phasevoltagediff=0;
        
        biphasic=false;
        firstpeakh=NaN;
        firstpeakval=eventdata(api).maxval;%yfiltered(threshh);
        maxddvh=threshh;
        maxdddvh=threshh;
        maxddddvh=threshh;
    end
    
    
%     figure(3)
%     clf
%     subplot(5,1,1)
%     plot(v)
%     subplot(5,1,2)
%     plot(dv)
%     hold on
%     plot(maxdvh,dv(maxdvh),'ro')
%     subplot(5,1,3)
%     plot(ddv)
%     hold on
%     plot(maxddvh,ddv(maxddvh),'ro')
%     subplot(5,1,4)
%     plot(dddv)
%     hold on
%     plot(maxdddvh,dddv(maxdddvh),'ro')
%     subplot(5,1,5)
%     plot(ddddv)
%     hold on
%     plot(maxddddvh,ddddv(maxddddvh),'ro')
    
%     %% looking for biphasic AP rise
%     yfilterednow=yfiltered(threshh:maxh);
%     dyfilterednow=dyfiltered(threshh:maxh);
%     ddyfilterednow=ddyfiltered(threshh:maxh);
%     dddyfilterednow=dddyfiltered(threshh:maxh);
%     [pks,locs]=findpeaks(dyfilterednow);
%     needed=pks>30;
%     locs=locs(needed);
%     pks=pks(needed);
% %     locs=sort(locs,'ascend');
%     [pks,idx]=sort(pks,'descend');
%     locs=locs(idx);
%     if length(locs)==1
%         [~,minhely]=min(dyfilterednow(1:locs)-linspace(dyfilterednow(1),dyfilterednow(locs),locs)');
%         if minhely > 2
%             [pksnew,locsnew]=findpeaks(dyfilterednow(1:minhely)-linspace(dyfilterednow(1),dyfilterednow(minhely),minhely)');
%             locs=[locs;locsnew];
%             pks=dyfilterednow(locs);
%             needed=pks>30;
%             locs=locs(needed);
%             pks=pks(needed);
%             %         locs=sort(locs,'ascend');
%             [pks,idx]=sort(pks,'descend');
%             locs=locs(idx);
%         end
%     end
%     while any(abs(diff(locs))<=5)
%         if length(locs)>1
%             difik=abs(diff(locs));
%             neededlocs=[1;find(difik>5)+1];
%             locs=locs(neededlocs);
%             pks=pks(neededlocs);
%         end
%     end
% %     if abs(diff(locs))<5
% %         locs=locs(1);
% %         pks=pks(1);
% %     end
%     if length(locs)>1
%         biphasic=true;
%         locs=[locs(1),locs(2)];
%         maxdvhs=[min(locs),max(locs)]+threshh-1;
%         %%
%         dytoanal=dyfiltered(maxdvhs(1):maxdvhs(2));
%         dytoanal=dytoanal-linspace(dytoanal(1),dytoanal(end),length(dytoanal))';
%         [~,firstpeakh]=min(dytoanal);
%         firstpeakh=firstpeakh+maxdvhs(1)-1;
%         firstpeakval=yfiltered(firstpeakh);
%     else
%         biphasic=false;
%         maxdvhs=locs+threshh-1;
%         firstpeakh=NaN;
%         firstpeakval=eventdata(api).maxval;%yfiltered(threshh);
%     end
%     eventdata(api).phasedelay=diff(maxdvhs)*si;
%     eventdata(api).phasevoltagediff=diff(yfiltered(maxdvhs));
%     if isempty(eventdata(api).phasedelay) %| eventdata(api).phasedelay< si*5
%         eventdata(api).phasedelay=0;
%         eventdata(api).phasevoltagediff=0;
%         
%         biphasic=false;
%         firstpeakh=NaN;
%         firstpeakval=eventdata(api).maxval;%yfiltered(threshh);
%     end
    %%
    if doitall
        eventdata(api).prevISI=apprevisis(apii);
        eventdata(api).prevAPnum_500ms=sum(apmaxtimes>=eventdata(api).maxtime-.5 & apmaxtimes<eventdata(api).maxtime);
        eventdata(api).prevAPnum_100ms=sum(apmaxtimes>=eventdata(api).maxtime-.1 & apmaxtimes<eventdata(api).maxtime);
        eventdata(api).prevAPnum_50ms=sum(apmaxtimes>=eventdata(api).maxtime-.05 & apmaxtimes<eventdata(api).maxtime);
    end
    
    
%     
%     if isempty(eventdata(api).phasevoltagediff)
%         
%     end
    threshv=yfiltered(threshh);
    maxval=y(maxh);
    eventdata(api).threshv=threshv;
    eventdata(api).maxval=maxval;
    eventdata(api).APamplitude=eventdata(api).maxval-eventdata(api).threshv;
%     eventdata(api).firstpeakval=firstpeakval;
%     eventdata(api).firstpeakamplitude=eventdata(api).firstpeakval-eventdata(api).threshv;
%     eventdata(api).secondpeakamplitude=eventdata(api).APamplitude-eventdata(api).firstpeakamplitude;
    
    depolrateV=yfiltered(max(threshh-round(valtozok.depolratewindow/si),1):threshh-round(.0001/si));
    depolratet=(1:length(depolrateV))*si;
    p=polyfit(depolratet',depolrateV,1);
    eventdata(api).depolrate=p(1);
    eventdata(api).baselineV=mean(yfiltered(max(threshh-round(valtozok.baselineVwindow/si),1):threshh));
    eventdata(api).baselineV_long=mean(yfiltered(max(threshh-round(valtozok.baselineVwindow_long/si),1):threshh));
    
    dv=diff(y(threshh:maxh))/si;
    [maxdval,maxdvh]=max(dv);
    maxdvh=maxdvh+threshh-1;
    centerh=threshh;%maxh;%threshh;%maxdvh;
    
    t=[-stepback:stepforward]*si*1000;
    
    
    
    if length(y)-stepforward<centerh
        vegh=length(y);
        nanvegen=stepforward-(length(y)-centerh);
    else
        vegh=centerh+stepforward;
        nanvegen=0;
    end
    if stepback>=centerh
        kezdeth=1;
        nanelejen=stepback-centerh+1;
    else
        kezdeth=centerh-stepback;
        nanelejen=0;
    end
    
    v=[nan(nanelejen,1);y(kezdeth:vegh)';nan(nanvegen,1)]';
%     v=yfiltered(threshh:maxh);
    dv=[nan(nanelejen,1);dyfiltered(kezdeth:vegh);nan(nanvegen,1)]';%dyfiltered(threshh:maxh);
    ddv=[nan(nanelejen,1);ddyfiltered(kezdeth:vegh);nan(nanvegen,1)]';%ddyfiltered(threshh:maxh);
    dddv=[nan(nanelejen,1);dddyfiltered(kezdeth:vegh);nan(nanvegen,1)]';%dddyfiltered(threshh:maxh);
    ddddv=[nan(nanelejen,1);dddyfiltered(kezdeth:vegh);nan(nanvegen,1)]';%ddddyfiltered(threshh:maxh);
%     %             v=y([-stepback:stepforward]+centerh);
%     if valtozok.movingt>1
%         if strcmp(valtozok.filter,'boxcar')
%             dv=diff(moving(v,valtozok.movingn))'/si;
%             ddv=diff(moving(dv,valtozok.movingn))'/si;
%             dddv=diff(moving(ddv,valtozok.movingn))'/si;
%         elseif strcmp(valtozok.filter,'gauss')
%             dv=diff(imfilter(v,hgauss','replicate'))/si;
%             ddv=diff(imfilter(dv,hgauss','replicate'))/si;
%             dddv=diff(imfilter(ddv,hgauss','replicate'))/si;
%         end
%     else
%         dv=diff(v/si);
%         ddv=diff(dv/si);
%         dddv=diff(ddv/si);
%     end
%     vdv=mean([v(2:end);v(1:end-1)]);
%     vddv=mean([vdv(2:end);vdv(1:end-1)]);
%     vdddv=mean([vddv(2:end);vddv(1:end-1)]);
%     tdv=mean([t(2:end);t(1:end-1)]);
%     tddv=mean([tdv(2:end);tdv(1:end-1)]);
%     tdddv=mean([tddv(2:end);tddv(1:end-1)]);


%     vdv=v;
%     vddv=v;
%     vdddv=v;
%     tdv=t;
%     tddv=t;
%     tdddv=t;
    
    
    
    

    
    
    APwaves(apii).t=t';
    APwaves(apii).v=v'*1000;
    APwaves(apii).dv=dv';
    APwaves(apii).ddv=ddv';
    APwaves(apii).dddv=dddv';
    APwaves(apii).ddddv=ddddv';
%     APwaves(apii).tdv=tdv';
%     APwaves(apii).tddv=tddv';
%     APwaves(apii).tdddv=tdddv';
%     APwaves(apii).vdv=vdv'*1000;
%     APwaves(apii).vddv=vddv'*1000;
%     APwaves(apii).vdddv=vdddv'*1000;
    APwaves(apii).si=si;
    APwaves(apii).RS=RS;
    APwaves(apii).maxtime=eventdata(api).maxtime;
    APwaves(apii).stimulated=eventdata(api).stimulated;
    APwaves(apii).threshv=eventdata(api).threshv;
    APwaves(apii).depolrate=eventdata(api).depolrate;
    APwaves(apii).baselineV=eventdata(api).baselineV;
    APwaves(apii).maxdvhs=maxdvhs-centerh+stepback+1;
    APwaves(apii).spikeletpeakh=firstpeakh-centerh+stepback+1;
    APwaves(apii).biphasic=biphasic;
%     if isempty(maxdval)
%         pause
%     end
    APwaves(apii).maxdv=max(dyfiltered(maxdvhs));%maxdval;%max(dv);
    eventdata(api).maxdv=max(dyfiltered(maxdvhs));%maxdval;%max(dv);
    if length(maxdvhs)==1
        maxdvhs=[threshh,maxdvhs];
    end
    eventdata(api).threshh=threshh;
    eventdata(api).maxdv1_V=yfiltered(maxdvhs(1));
    eventdata(api).maxdv1_t=(maxdvhs(1)-threshh)*si;
    eventdata(api).maxdv2_V=yfiltered(maxdvhs(end));
    eventdata(api).maxdv2_t=(maxdvhs(end)-threshh)*si;
    eventdata(api).maxdv1_dV=dyfiltered(maxdvhs(1));
    eventdata(api).maxdv2_dV=dyfiltered(maxdvhs(end));
    eventdata(api).maxddv_V=yfiltered(maxddvh);
    eventdata(api).maxddv_ddV=ddyfiltered(maxddvh);
    eventdata(api).maxddv_t=(maxddvh-threshh)*si;
    eventdata(api).maxdddv_V=yfiltered(maxdddvh);
    eventdata(api).maxdddv_dddV=dddyfiltered(maxdddvh);
    eventdata(api).maxdddv_t=(maxdddvh-threshh)*si;
    eventdata(api).maxddddv_V=yfiltered(maxddddvh);
    eventdata(api).maxddddv_ddddV=ddddyfiltered(maxddddvh);
    eventdata(api).maxddddv_t=(maxddddvh-threshh)*si;
%     if biphasic==true
%         APwaves(apii).maxdv1=
%     else
%     end
%%

decay=yfiltered(maxh:min(maxh+100,length(yfiltered)));
decayt=[1:length(decay)]';

p=polyfit(decayt,decay,10);
% figure(33)
% clf
% plot(decayt,decay)
% hold on
% plot(decayt,polyval(p,decayt),'r-')
decay=decay-polyval(p,decayt);
% figure(34)
% plot(decay)
eventdata(api).oscillation=std(decay);
% std(decay)
end
handles.data.eventdata=eventdata;
handles.data.eventwaves=APwaves;
handles.data.selectedevents=neededevents;

function handles = plotthedata(handles)
if isfield(handles.data,'eventwaves') & ~isempty(fieldnames(handles.data.eventwaves))
    for ploti=1:4
%         Xvaltoplot=get(handles.popupmenu1,'String');
%         Xvaltoplot=Xvaltoplot{get(handles.popupmenu1,'Val')};
%         Yvaltoplot=get(handles.popupmenu2,'String');
%         Yvaltoplot=Yvaltoplot{get(handles.popupmenu2,'Val')};
        
        Xvaltoplot=handles.data.plotdata(ploti).X;
        Yvaltoplot=handles.data.plotdata(ploti).Y;
%         if any(strfind(Yvaltoplot,'d'))
%             Xvaltoplot=[Xvaltoplot,Yvaltoplot];
%         end
        axes(handles.(['axes',num2str(ploti)]));
        cla
        sis=unique([handles.data.eventwaves.si]);
        for i=1:length(sis)
            needed=find([handles.data.eventwaves.si]==sis(i));
            hold on
            X=[handles.data.eventwaves(needed).(Xvaltoplot)];
            Y=[handles.data.eventwaves(needed).(Yvaltoplot)];
            plot(X,Y);
            hold on
            
            
%             if any(strfind(Yvaltoplot,'dv'))
               for ii=1:length(needed) 
                   plot(X([handles.data.eventwaves(needed(ii)).maxdvhs],ii),Y([handles.data.eventwaves(needed(ii)).maxdvhs],ii),'ro')
                   if handles.data.eventwaves(needed(ii)).biphasic
                       plot(X([handles.data.eventwaves(needed(ii)).spikeletpeakh],ii),Y([handles.data.eventwaves(needed(ii)).spikeletpeakh],ii),'ko')
                   end
               end
%             end
            
            axis tight
        end
        if any(strfind(Yvaltoplot,'dddv'))
            ylabel('dddV/dt (mV/ms^3)')
        elseif any(strfind(Yvaltoplot,'ddv'))
            ylabel('ddV/dt (mV/ms^2)')
        elseif any(strfind(Yvaltoplot,'dv'))
            ylabel('dV/dt (mV/ms)')
        elseif any(strfind(Yvaltoplot,'v'))
            ylabel('Voltage (mV)')
        end
        if any(strfind(Xvaltoplot,'t'))
            xlabel('Time (ms)')
        elseif any(strfind(Xvaltoplot,'v'))
            xlabel('Voltage (mV)')
        end
        handles.(['axes',num2str(ploti)])=gca;
        %%
        % 
    end
%     disp('lol')
    whattoplot=get(handles.popupmenu8,'String');
    whattoplot=whattoplot{get(handles.popupmenu8,'Value')};
    if strcmp(whattoplot,'histogram')
        selectedeventdatafield=get(handles.popupmenu3,'String');
        selectedeventdatafield=selectedeventdatafield{get(handles.popupmenu3,'Value')};
        axes(handles.axes5);
        cla
        hist([handles.data.eventdata.(selectedeventdatafield)],100)
        xlabel(selectedeventdatafield)
        ylabel('count')
        handles.axes5=gca;
    elseif strcmp(whattoplot,'plot VS')
        selectedeventdatafield=get(handles.popupmenu3,'String');
        selectedeventdatafieldY=selectedeventdatafield{get(handles.popupmenu3,'Value')};
        selectedeventdatafield=get(handles.popupmenu7,'String');
        selectedeventdatafieldX=selectedeventdatafield{get(handles.popupmenu7,'Value')};
        X=[handles.data.eventdata.(selectedeventdatafieldX)];
        Y=[handles.data.eventdata.(selectedeventdatafieldY)];
        axes(handles.axes5);
        cla
        plot(X,Y,'ko')
        hold on
        if isfield(handles.data.eventdata,'axonalAP')
            aAPs=[handles.data.eventdata.axonalAP]==1;
            sAPs=[handles.data.eventdata.somaticAP]==1;
            if ~isempty(find(aAPs))
                plot(X(aAPs),Y(aAPs),'ro')
            end
            if ~isempty(find(sAPs))
                plot(X(sAPs),Y(sAPs),'bo')
            end
        end
        
        plot(X(handles.data.selectedevents),Y(handles.data.selectedevents),'kx','MarkerSize',16,'LineWidth',4)
        xlabel(selectedeventdatafieldX)
        ylabel(selectedeventdatafieldY)
        handles.axes5=gca;
    else
        disp(['dunno what to plot: ',whattoplot])
    end

% 


    
end
% UIWAIT makes aE_InspectEvents wait for user response (see UIRESUME)
% uiwait(handles.figure1);





% --- Outputs from this function are returned to the command line.
function varargout = aE_InspectEvents_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
axestoplot=get(handles.popupmenu4,'Value');
Xvaltoplot=get(handles.popupmenu1,'String');
Xvaltoplot=Xvaltoplot{get(handles.popupmenu1,'Val')};
Yvaltoplot=get(handles.popupmenu2,'String');
Yvaltoplot=Yvaltoplot{get(handles.popupmenu2,'Val')};
handles.data.plotdata(axestoplot).X=Xvaltoplot;
handles.data.plotdata(axestoplot).Y=Yvaltoplot;
handles = plotthedata(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be definhandles = plotthedata(handles);ed in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axestoplot=get(handles.popupmenu4,'Value');
Xvaltoplot=get(handles.popupmenu1,'String');
Xvaltoplot=Xvaltoplot{get(handles.popupmenu1,'Val')};
Yvaltoplot=get(handles.popupmenu2,'String');
Yvaltoplot=Yvaltoplot{get(handles.popupmenu2,'Val')};
handles.data.plotdata(axestoplot).X=Xvaltoplot;
handles.data.plotdata(axestoplot).Y=Yvaltoplot;
handles = plotthedata(handles);
guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
 handles = updateGUI(handles);
 guidata(hObject,handles);
% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = updateDATA(handles);
handles = plotthedata(handles);
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
axestoplot=get(handles.popupmenu4,'Value');
Xvalstoplot=get(handles.popupmenu1,'String');
set(handles.popupmenu1,'Val',find(strcmp(handles.data.plotdata(axestoplot).X,Xvalstoplot)));

Yvalstoplot=get(handles.popupmenu2,'String');
set(handles.popupmenu2,'Val',find(strcmp(handles.data.plotdata(axestoplot).Y,Yvalstoplot)));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection changehandles = updateGUI(handles) in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contenstringtoshow{i}=[handles.data.rules(i).selectedeventdatafield,' ',handles.data.rules(i).comparefunc,' ',handles.data.rules(i).value];ts of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%
selectedeventdatafield=get(handles.popupmenu5,'String');
selectedeventdatafield=selectedeventdatafield{get(handles.popupmenu5,'Value')};
modifyfunc=get(handles.popupmenu11,'String');
modifyfunc=modifyfunc{get(handles.popupmenu11,'Value')};
comparefunc=get(handles.popupmenu6,'String');
comparefunc=comparefunc{get(handles.popupmenu6,'Value')};

if isempty(fieldnames(handles.data.rules))
    next=1;
else
    next=length(handles.data.rules)+1;
end

handles.data.rules(next).selectedeventdatafield=selectedeventdatafield;
handles.data.rules(next).comparefunc=comparefunc;
handles.data.rules(next).modifyfunc=modifyfunc;

if ~ischar(handles.data.eventdata_orig(1).(selectedeventdatafield))    
    handles.data.rules(next).value=str2num(get(handles.edit1,'String'));
else
    handles.data.rules(next).value=get(handles.edit1,'String');
end
handles=applyrules(handles);
handles = updateGUI(handles);
guidata(hObject, handles);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selectedrule=get(handles.listbox2,'Value');
set(handles.listbox2,'Value',1);
if ~isempty(fieldnames(handles.data.rules))& length(handles.data.rules)>0
    handles.data.rules(selectedrule)=[];
end

handles=applyrules(handles);
handles = updateGUI(handles);
guidata(hObject, handles);

% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%



whattoplot=get(handles.popupmenu8,'String');
whattoplot=whattoplot{get(handles.popupmenu8,'Value')};
if strcmp(whattoplot,'histogram')
    disp('this function doesn''t work with histogram yet.. choose plot VS ')
elseif strcmp(whattoplot,'plot VS')
    axes(handles.axes5)
    [x,y]=ginput(1);
    selectedeventdatafield=get(handles.popupmenu3,'String');
    selectedeventdatafieldY=selectedeventdatafield{get(handles.popupmenu3,'Value')};
    selectedeventdatafield=get(handles.popupmenu7,'String');
    selectedeventdatafieldX=selectedeventdatafield{get(handles.popupmenu7,'Value')};
    X=[handles.data.eventdata.(selectedeventdatafieldX)];
    Y=[handles.data.eventdata.(selectedeventdatafieldY)];
    %%
    distmatrix=[x,X;y,Y]';
    xmu=nanmean(distmatrix);
    xsigma=nanstd(distmatrix);
    distmatrix=(distmatrix-repmat(xmu,length(distmatrix),1))./repmat(xsigma,length(distmatrix),1);
    %%
    distances=pdist(distmatrix);
    distances=distances(1:length(X));
    [~,idx]=nanmin(distances);
    idx=idx(1);
    
    
    if length(handles.data.eventdata)==0
        uniquevals='';
    else
        if ~ischar(handles.data.eventdata(1).(selectedeventdatafieldY))
            uniquevals=unique([handles.data.eventdata.(selectedeventdatafieldY)]);
        else
            uniquevals=unique({handles.data.eventdata.(selectedeventdatafieldY)});
        end
    end
    idxtoselect=find(uniquevals==Y(idx));
    set(handles.listbox1,'Value',idxtoselect);

    handles = updateDATA(handles);
    handles = plotthedata(handles);
    %%

end
guidata(hObject,handles);
%%
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%
eventdata_all=handles.data.eventdata_orig;
eventdataforclust=handles.data.eventdata;
eventdataforclust_normalized=[];
fieldnamek=fieldnames(eventdataforclust);
fieldnamestodel={'sweepnum','sweeptime','si','onseth','maxh','endh','h10','h90','type','maxtime','onsettime','maxdvalue','beselinesd_of_dy',...
    'injectedcurrentatpeak','injectedrelativecurrentatpeak','injectedcurrentatonset','injectedrelativecurrentatonset','stimulated','maxtimetosquarepulse',...
    'onsettimetosquarepulse','axonalAP','somaticAP','aAPfreq_in_window','sAPfreq_in_window','axonalAP_sporadic','axonalAP_persistent'};
fieldnamestodel=[fieldnamestodel,fieldnamek(~cellfun(@isempty,strfind(fieldnamek,'Princom')))'];
for i=1:length(fieldnamestodel)
    fieldnamek(strcmp(fieldnamek,fieldnamestodel{i}))=[];
end
reali=0;
for i=1:length(fieldnamek)
    if any(isnan([eventdataforclust.(fieldnamek{i})]))% | any(strcmp(fieldnamek{i},{'RS','si','fname','ID','class'}))
        eventdataforclust=rmfield(eventdataforclust,fieldnamek{i});
        eventdata_all=rmfield(eventdata_all,fieldnamek{i});
    else
        reali=reali+1;
%         DATAforclust_normalized.(fieldnamek{i})=zscore([DATAforclust.(fieldnamek{i})],1);
        [eventdataforclust_normalized.data(:,reali),eventdataforclust_normalized.mu(reali),eventdataforclust_normalized.sigma(reali)]=zscore([eventdataforclust.(fieldnamek{i})],1);
        eventdataforclust_normalized.alldata(:,reali)=([eventdata_all.(fieldnamek{i})]-eventdataforclust_normalized.mu(reali))/eventdataforclust_normalized.sigma(reali);
        eventdataforclust_normalized.header{reali}=fieldnamek{i};
    end
end
eventdataforclust_normalized.alldata(isinf(eventdataforclust_normalized.alldata))=0;
[indx,tf] = listdlg('ListString',eventdataforclust_normalized.header);
%%
if tf
    eventdataforclust_normalized.data=eventdataforclust_normalized.data(:,indx);
    eventdataforclust_normalized.alldata=eventdataforclust_normalized.alldata(:,indx);
    eventdataforclust_normalized.header=eventdataforclust_normalized.header(indx);
    eventdataforclust_normalized.mu=eventdataforclust_normalized.mu(indx);
    eventdataforclust_normalized.sigma=eventdataforclust_normalized.sigma(indx);
    
    [COEFF,SCORE,latent,tsquare] = princomp(eventdataforclust_normalized.data);
    %%
    figure(3)
    clf
    for i=1:min(7,size(SCORE,2))
        subplot(3,3,i)
        hist(SCORE(:,i),100)
    end
    subplot(3,3,8)
    plot(cumsum(latent)/sum(latent),'r-','LineWidth',3)
    subplot(3,3,9)
    plot(SCORE(:,1),SCORE(:,2),'ko')
    %%
%     nonnanidx=sum(isnan(eventdataforclust_normalized.alldata)')==0;
    princompstosave=10;
    scorecell=num2cell(SCORE);
    SCORE_all=[eventdataforclust_normalized.alldata*COEFF];
%     SCORE_all_nonnan=[eventdataforclust_normalized.alldata(nonnanidx,:)*COEFF];
    scorecell_all=num2cell(SCORE_all);
    for i=1:princompstosave
        [handles.data.eventdata.(['Princomp',num2str(i)])] = scorecell{:,i};
        [handles.data.eventdata_orig.(['Princomp',num2str(i)])] = scorecell_all{:,i};
    end
    
    % threshidx=find(strcmp('threshv',eventdataforclust_normalized.header));
    %
    % [~,threshVcoefforder]=sort(abs(COEFF(threshidx,:)),'descend');
    % for i=1:princompstosave
    %     [handles.data.eventdata.(['Princomp_threshV',num2str(i)])] = scorecell{:,threshVcoefforder(i)};
    % end
    
    
    % [handles.data.eventdata.Princomp1] = scorecell{:,1};
    % [handles.data.eventdata.Princomp2] = scorecell{:,2};
    % [handles.data.eventdata.Princomp3] = scorecell{:,3};
    %%
    figure(4)
    clf
    for i=1:min(7,size(SCORE,2))
        
        coeffnow=COEFF(:,i);
        sum(abs(coeffnow))
        subplot(3,3,i)
        bar(coeffnow)
        
        [coeffnow,idx]=sort(abs(COEFF(:,i)),'descend');
        coeffnow=coeffnow/sum(coeffnow);
        title([eventdataforclust_normalized.header(idx(1:3))])
    end
    
    
    %%
    set(handles.popupmenu7,'String',fieldnames(handles.data.eventdata));
     set(handles.popupmenu5,'String',fieldnames(handles.data.eventdata));
     set(handles.popupmenu3,'String',fieldnames(handles.data.eventdata));
    handles=updateeventdata(handles);
    guidata(hObject,handles)
end
% %%
% Y = pdist(eventdataforclust_normalized.data);
% % Y = pdist(SCORE(:,1:10));
% 
% Z = linkage(Y);
% [H,T,outperm]=dendrogram(Z,0,'Orientation','left');
% dendrogram(Z,0,'Orientation','left')%,'Label',{DATASUM(outperm).ID});


%%
% disp('pushbutton4 - no function now')

% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = plotthedata(handles);
guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7


% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 handles = updateDATA(handles,true);
 handles=updateeventdata(handles);
 guidata(hObject,handles);
 %%
%  disp('lol')


% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = plotthedata(handles);
guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8


% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu9.
function popupmenu9_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%
whattoplot=get(handles.popupmenu8,'String');
whattoplot=whattoplot{get(handles.popupmenu8,'Value')};
if strcmp(whattoplot,'plot VS')
    axes(handles.axes5)
    hAP=imfreehand;
    line=hAP.getPosition;
    %%
    selectedeventdatafield=get(handles.popupmenu3,'String');
    selectedeventdatafieldY=selectedeventdatafield{get(handles.popupmenu3,'Value')};
    selectedeventdatafield=get(handles.popupmenu7,'String');
    selectedeventdatafieldX=selectedeventdatafield{get(handles.popupmenu7,'Value')};
    X=[handles.data.eventdata.(selectedeventdatafieldX)];
    Y=[handles.data.eventdata.(selectedeventdatafieldY)];
    APk=inpolygon(X,Y,line(:,1),line(:,2));
    
    aptype=get(handles.popupmenu9,'String');
    aptype=aptype{get(handles.popupmenu9,'Value')};
    if any(strfind(aptype,'aAP'))
        [handles.data.eventdata(APk).axonalAP]=deal(true);
        [handles.data.eventdata(APk).somaticAP]=deal(false);
    elseif any(strfind(aptype,'sAP'))
        [handles.data.eventdata(APk).axonalAP]=deal(false);
        [handles.data.eventdata(APk).somaticAP]=deal(true);
    else
        [handles.data.eventdata(APk).axonalAP]=deal(false);
        [handles.data.eventdata(APk).somaticAP]=deal(false);
        %%
        APkf=find(APk);
        for i=1:length(APkf)
            handles.data.eventdata(APkf(i)).type='noise';
        end
        %%
%         for 
%         {handles.data.eventdata(APk).type}=deal('noise');
    end
    hAP.delete;
    handles=updateeventdata(handles);
end

handles = plotthedata(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popupmenu9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eventdata_orig=eventdata;
handles=updateeventdata(handles);
handles_orig=handles;
%%
h=findobj('Name','aE_InspectTraces');
handles = guidata(h);
selectedsamplenum=handles.data.selectedsamplenum;
handles.data.samples(selectedsamplenum).eventdata=handles_orig.data.eventdata_orig;
directory=handles.data.dirs.eventdir;
ID=handles.data.xlsdata(handles.data.samples(selectedsamplenum).loadedID-1).ID;
guidata(h,handles);
eventdata=handles_orig.data.eventdata_orig;
save([directory,'sorted/',ID],'eventdata');
handles=handles_orig;
eventdata=eventdata_orig;
guidata(hObject,handles);


% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10


% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu11.
function popupmenu11_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu11 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu11


% --- Executes during object creation, after setting all properties.
function popupmenu11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
