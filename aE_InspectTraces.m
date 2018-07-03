function varargout = aE_InspectTraces(varargin)
% With this script one can observe bridge balanced traces exported with the  script :analyseElphys_main.m
% USAGE: aE_InspectTraces(dirs,xlsdata) - after starting analyseElphys_main.m


%      AE_INSPECTTRACES, by itself, creates a new AE_INSPECTTRACES or raises the existing
%      singleton*.
%
%      H = AE_INSPECTTRACES returns the handle to a new AE_INSPECTTRACES or the handle to
%      the existing singleton*.
%
%      AE_INSPECTTRACES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AE_INSPECTTRACES.M with the given input arguments.
%
%      AE_INSPECTTRACES('Property','Value',...) creates a new AE_INSPECTTRACES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before aE_InspectTraces_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to aE_InspectTraces_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help aE_InspectTraces

% Last Modified by GUIDE v2.5 01-Jul-2018 17:04:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @aE_InspectTraces_OpeningFcn, ...
                   'gui_OutputFcn',  @aE_InspectTraces_OutputFcn, ...
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


% --- Executes just before aE_InspectTraces is made visible.
function aE_InspectTraces_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.0
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to aE_InspectTraces (see VARARGIN)

% Choose default command line output for aE_InspectTraces
handles.output = hObject;
handles.data.dirs=varargin{1};
handles.data.xlsdata=varargin{2};
set(handles.popupmenu6,'String',[{'none'};fieldnames(handles.data.xlsdata)]);


BrainStateProps.Statevalues={'Slow wave sleep', 'REM sleep', 'Quiet wakefulness','Active wakefulness','Delete'};
BrainStateProps.Statecolors={'blue','cyan','red','green'};
BrainStateProps.Statecolors_short={'b','c','r','g'};
handles.data.BrainStateProps=BrainStateProps;

% Update handles structure
guidata(hObject, handles);
resetdata(hObject,handles,'all')
% pushbutton4_Callback(hObject, eventdata, handles)
% UIWAIT makes aE_InspectTraces wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = aE_InspectTraces_OutputFcn(hObject, eventdata, handles) 
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
%%
samplenum=get(hObject,'Value');
%%
if handles.data.samples(samplenum).selectedID>1
    handles=loadthedata(handles);
    handles=updategui(handles);
    guidata(hObject, handles);
end
% handles=updatedatatoplot(handles);
% plotandupdate(handles)

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


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
selectedsamplenum=get(handles.popupmenu1,'Value');
handles.data.samples(selectedsamplenum).switch=get(hObject,'Value');
set(handles.listbox1,'Value',1)
handles=loadthedata(handles);
% handles=updatedatatoplot(handles);
guidata(hObject,handles);
% plotandupdate(handles);


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
selectedsamplenum=get(handles.popupmenu1,'Value');
IDs=handles.data.actualIDs;%get(hObject,'String');
selectedID=get(hObject,'Value');
handles.data.samples(selectedsamplenum).selectedID=find(strcmp(IDs{selectedID},handles.data.IDs));
if handles.data.samples(selectedsamplenum).selectedID>1
    load([handles.data.dirs.bridgeddir,handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID}],'lightdata');
    load([handles.data.dirs.bridgeddir,handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID}],'bridgeddata');
    handles.data.samples(selectedsamplenum).lightdata=lightdata;
    handles.data.samples(selectedsamplenum).bridgeddata=bridgeddata;
%     eventdataorig=eventdata;
    if ~isempty(dir([handles.data.dirs.eventdir,'sorted/',handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID},'.mat']))
        load([handles.data.dirs.eventdir,'sorted/',handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID}],'eventdata');
    elseif ~isempty(dir([handles.data.dirs.eventdir,handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID},'.mat']))
        load([handles.data.dirs.eventdir,handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID}],'eventdata');
        set(handles.checkbox3,'Value',0);
        set(handles.checkbox5,'Value',0);
    else
        eventdata=[];
    end
    handles.data.samples(selectedsamplenum).eventdata=eventdata;
%     eventdata=eventdataorig;
    handles.data.samples(selectedsamplenum).starttime=handles.data.samples(selectedsamplenum).lightdata(1).realtime;
    handles.data.samples(selectedsamplenum).endtime=handles.data.samples(selectedsamplenum).lightdata(end).realtime+handles.data.samples(selectedsamplenum).bridgeddata(end).si*length(handles.data.samples(selectedsamplenum).bridgeddata(end).y);
    
    handles.data.samples(selectedsamplenum).neededwaves=[];
    handles.data.samples(selectedsamplenum).maybeneededwaves=[];
    handles=loadthedata(handles);
    handles=updatedatatoplot(handles);
    guidata(hObject,handles);
    set(handles.listbox1,'Value',1)
else
    handles.data.samples(selectedsamplenum).lightdata=struct;
end
plotandupdate(handles);
handles=updategui(handles);
guidata(hObject,handles)
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
samplenow=get(handles.popupmenu1,'Value');
handles.data.samples(samplenow).selectedvariable=get(hObject,'Value');
set(handles.listbox1,'Value',1);
handles=loadthedata(handles);
handles=updatedatatoplot(handles);
guidata(hObject,handles)
% plotandupdate(handles)
handles=updategui(handles);
guidata(hObject,handles)


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

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selectedsamplenum=get(handles.popupmenu1,'Value');
selectedvariable=handles.data.samples(selectedsamplenum).selectedvariable;
selectedID=handles.data.samples(selectedsamplenum).selectedID;
selectedvariablename=handles.data.fieldnevek{selectedvariable};

if ischar(handles.data.samples(selectedsamplenum).lightdata(1).(selectedvariablename))
    reducednames=unique({handles.data.samples(selectedsamplenum).lightdata.(selectedvariablename)});
    ezszoveg=1;
else
    reducednames=unique([handles.data.samples(selectedsamplenum).lightdata.(selectedvariablename)]);
    ezszoveg=0;
end
valuesnow=get(handles.listbox1,'Value');
needed=[];
szovegmost=['-',selectedvariablename,':'];
for i=1:length(valuesnow)
    if ezszoveg==1
        needed=[needed,find(strcmp({handles.data.samples(selectedsamplenum).lightdata.(selectedvariablename)},reducednames{valuesnow(i)}))];
        szovegmost=[szovegmost,' ',reducednames{valuesnow(i)},','];
        if length(valuesnow)>5
            szovegmost=['-',selectedvariablename,': ',reducednames{valuesnow(1)},' ... ',reducednames{valuesnow(end)}];
        end
    else
        needed=[needed,find([handles.data.samples(selectedsamplenum).lightdata.(selectedvariablename)]==reducednames(valuesnow(i)))];
        szovegmost=[szovegmost,' ',num2str(reducednames(valuesnow(i))),','];
        if length(valuesnow)>5
            szovegmost=['-',selectedvariablename,': ',num2str(reducednames(valuesnow(1))),' ... ',num2str(reducednames(valuesnow(end)))];
        end
    end
end
handles.data.samples(selectedsamplenum).maybeneededwaves=sort(unique(needed));
handles.data.samples(selectedsamplenum).maybeneededwavespervariable{selectedvariable}=handles.data.samples(selectedsamplenum).neededwaves;
% timesofneededwaves=[handles.data.samples(selectedsamplenum).lightdata(handles.data.samples(selectedsamplenum).neededwaves).realtime];
% reallyneededwaves=timesofneededwaves>=str2num(get(handles.edit4,'String')) & timesofneededwaves<=str2num(get(handles.edit5,'String'));
% if reallyneededwaves(1)==0
%     reallyneededwaves(find(timesofneededwaves<=str2num(get(handles.edit4,'String')),1,'last'))=true;
% end
% handles.data.samples(selectedsamplenum).neededwaves=handles.data.samples(selectedsamplenum).neededwaves(reallyneededwaves);
handles.data.samples(selectedsamplenum).changes=szovegmost;
handles=loadthedata(handles);
handles=updatedatatoplot(handles);
guidata(hObject,handles)
% plotandupdate(handles)

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%
selectedsamplenum=get(handles.popupmenu1,'Value');
starttime=str2num(get(handles.edit4,'String'));
endtime=str2num(get(handles.edit5,'String'));
neededidx=find(handles.data.samples(selectedsamplenum).movementdata.time>=starttime & handles.data.samples(selectedsamplenum).movementdata.time<=endtime);
movie=handles.data.samples(selectedsamplenum).movementdata.video(:,:,neededidx);
% implay(movie,100);

%
vlc_wrapper cleanup

movieidx=handles.data.samples(selectedsamplenum).movementdata.videoidx(neededidx(1));
moviename=handles.data.samples(selectedsamplenum).movementdata.videofname{movieidx};
framenum=handles.data.samples(selectedsamplenum).movementdata.framenumber(neededidx(1));
endmovieidx=handles.data.samples(selectedsamplenum).movementdata.videoidx(neededidx(end));
%%
nextmoviestart=find(handles.data.samples(selectedsamplenum).movementdata.videoidx>movieidx,1,'first');
if isempty(nextmoviestart)
    maxframenum=handles.data.samples(selectedsamplenum).movementdata.framenumber(end);
else
    maxframenum=handles.data.samples(selectedsamplenum).movementdata.framenumber(nextmoviestart-1);
end
if movieidx==endmovieidx
    framenumend=handles.data.samples(selectedsamplenum).movementdata.framenumber(neededidx(end));
else
    framenumend=maxframenum;
end

%%
setupname=handles.data.xlsdata(handles.data.samples(selectedsamplenum).selectedID-1).setup;
locations=marcicucca_locations;
moviefiletoplay=[locations.tgtardir,'VIDEOdata/',setupname,'/',moviename];
a=dir(moviefiletoplay);
vh = vlc_wrapper('init');
vp = vlc_wrapper('open', vh, moviefiletoplay);
info=vlc_wrapper('info', vp);
vlcmatlabratio=info(2)/maxframenum;

%%
vlc_wrapper('frame', vp,round(framenum*vlcmatlabratio));
%%
vlc_wrapper('play', vp, 1);
while ~isempty(vlc_wrapper('frame', vp))& vlc_wrapper('frame', vp)>0 &  vlc_wrapper('frame', vp)<=framenumend*vlcmatlabratio 
    pause(1)
    disp(['VLC frame: ' num2str(vlc_wrapper('frame', vp)), ' real frame:',num2str(vlc_wrapper('frame', vp)/vlcmatlabratio)])
end
vlc_wrapper cleanup
%%
%%
% disp('a')



%old version of pushbutton2 when it was AND Select
% % % % % selectedsamplenum=get(handles.popupmenu1,'Value');
% % % % % selectedvariable=handles.data.samples(selectedsamplenum).selectedvariable;
% % % % % selectedID=handles.data.samples(selectedsamplenum).selectedID;
% % % % % selectedvariablename=handles.data.fieldnevek{selectedvariable};
% % % % % 
% % % % % if ischar(handles.data.samples(selectedsamplenum).lightdata(1).(selectedvariablename))
% % % % %     reducednames=unique({handles.data.samples(selectedsamplenum).lightdata.(selectedvariablename)});
% % % % %     ezszoveg=1;
% % % % % else
% % % % %     reducednames=unique([handles.data.samples(selectedsamplenum).lightdata.(selectedvariablename)]);
% % % % %     ezszoveg=0;
% % % % % end
% % % % % valuesnow=get(handles.listbox1,'Value');
% % % % % needed=[];
% % % % % szovegmost=['-',selectedvariablename,':'];
% % % % % for i=1:length(valuesnow)
% % % % %     if ezszoveg==1
% % % % %         needed=[needed,find(strcmp({handles.data.samples(selectedsamplenum).lightdata.(selectedvariablename)},reducednames{valuesnow(i)}))];
% % % % %         szovegmost=[szovegmost,' ',reducednames{valuesnow(i)},','];
% % % % %         if length(valuesnow)>5
% % % % %             szovegmost=['-',selectedvariablename,': ',reducednames{valuesnow(1)},' ... ',reducednames{valuesnow(end)}];
% % % % %         end
% % % % %     else
% % % % %         needed=[needed,find([handles.data.samples(selectedsamplenum).lightdata.(selectedvariablename)]==reducednames(valuesnow(i)))];
% % % % %         szovegmost=[szovegmost,' ',num2str(reducednames(valuesnow(i))),','];
% % % % %         if length(valuesnow)>5
% % % % %             szovegmost=['-',selectedvariablename,': ',num2str(reducednames(valuesnow(1))),' ... ',num2str(reducednames(valuesnow(end)))];
% % % % %         end
% % % % %     end
% % % % %     
% % % % % end
% % % % % needed=sort(unique(needed));
% % % % % common=ismember(handles.data.samples(selectedsamplenum).neededwaves,needed);
% % % % % handles.data.samples(selectedsamplenum).neededwaves=handles.data.samples(selectedsamplenum).neededwaves(common);
% % % % % handles.data.samples(selectedsamplenum).neededwavespervariable{selectedvariable}=handles.data.samples(selectedsamplenum).neededwaves;
% % % % % handles.data.samples(selectedsamplenum).changes=[handles.data.samples(selectedsamplenum).changes,' AND ',szovegmost];
% % % % % 
% % % % % % timesofneededwaves=[handles.data.samples(selectedsamplenum).lightdata(handles.data.samples(selectedsamplenum).neededwaves).realtime];
% % % % % % reallyneededwaves=timesofneededwaves>=str2num(get(handles.edit4,'String')) & timesofneededwaves<=str2num(get(handles.edit5,'String'));
% % % % % % handles.data.samples(selectedsamplenum).neededwaves=handles.data.samples(selectedsamplenum).neededwaves(reallyneededwaves);
% % % % % 
% % % % % handles=loadthedata(handles);
% % % % % handles=updatedatatoplot(handles);
% % % % % guidata(hObject,handles)
% % % % % %plotandupdate(handles)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resetdata(hObject,handles,'all')

function resetdata(hObject,handles,melyiket)
data=handles.data;
markers={'ko','rx','b^','gv','cs'};
if ischar(melyiket)
    data.samples=struct;
end
data.IDs=['no cell selected',{handles.data.xlsdata.ID}];
data.actualIDs=data.IDs;
temp=load([data.dirs.bridgeddir,data.IDs{2}],'lightdata');
data.fieldnevek=fieldnames(temp.lightdata);
for i=1:length(markers)
    if ischar(melyiket) | melyiket==i
        data.samples(i).marker=markers{i};
        data.samples(i).switch=0;
        data.samples(i).selectedvariable=1;
        data.samples(i).changes=[];
        data.samples(i).bridgeddata=struct;
        data.samples(i).stimdata=struct;
        data.samples(i).eventdata=[];
        data.samples(i).selectedID=1;
        data.samples(i).loadedID=1;
        data.samples(i).lightdata=[];
        data.samples(i).cutoffreq=0;
        data.samples(i).filterdeg=1;
        data.samples(i).plotdetails=struct;
%         for j=1:length(data.fieldnevek)
%             data.samples(i).neededwavespervariable{j}=1:length(data.WAVES);
%         end
        data.samples(i).neededwaves=[];%1:length(data.WAVES);
    end
end

handles.data=data;
if ischar(melyiket)
    axes(handles.axes1)
    cla reset
    handles.axes1=gca;
    set(handles.checkbox1,'Value',0);
    set(handles.popupmenu1,'String',{handles.data.samples.marker});
    set(handles.popupmenu1,'Value',1);
%     set(handles.popupmenu2,'String',handles.data.IDs);
    handles.data.actualIDs=handles.data.IDs;
    set(handles.popupmenu2,'Value',1);
    set(handles.text10,'String',[]);
    
    %%
    set(handles.popupmenu3,'String',handles.data.fieldnevek);
    set(handles.popupmenu4,'String',[{'stiumulus'},{handles.data.samples.marker},{'movement and pupil diameter'},{'PSD of field'}]);
    set(handles.popupmenu7,'String',[{'stiumulus'},{handles.data.samples.marker},{'movement and pupil diameter'},{'PSD of field'}]);
    set(handles.slider1,'Min',0);
    set(handles.slider1,'Max',1);
    set(handles.slider2,'Min',0);
    set(handles.slider2,'Max',1);
    set(handles.slider1,'Value',0);
    set(handles.slider2,'Value',1);
    set(handles.slider3,'Min',0);
    set(handles.slider3,'Max',1);
    set(handles.slider3,'Value',0.5);
    
    %%
    toget=[];
    xlsfieldek=(fieldnames(handles.data.xlsdata));
    for i=1:length(xlsfieldek)
        if ischar(handles.data.xlsdata(1).(xlsfieldek{i}))
            toget=[toget,i];
        end
    end
    set(handles.popupmenu5,'String',xlsfieldek(toget));
    set(handles.edit6,'String','');
end
handles=loadthedata(handles);
handles=updatedatatoplot(handles);
handles=updategui(handles);
%plotandupdate(handles)
guidata(hObject, handles);

function handles=loadthedata(handles)
for samplei=1:length(handles.data.samples)
    if ~(handles.data.samples(samplei).loadedID==handles.data.samples(samplei).selectedID) &handles.data.samples(samplei).switch==1 & handles.data.samples(samplei).selectedID>1
        load([handles.data.dirs.bridgeddir,handles.data.IDs{handles.data.samples(samplei).selectedID}],'stimdata');
        handles.data.samples(samplei).stimdata=stimdata;
        load([handles.data.dirs.bridgeddir,handles.data.IDs{handles.data.samples(samplei).selectedID}],'bridgeddata');
        handles.data.samples(samplei).bridgeddata=bridgeddata;
        handles.data.samples(samplei).loadedID=handles.data.samples(samplei).selectedID;
        if isfield(handles.data.dirs,'videodir')
            fname=handles.data.xlsdata(handles.data.samples(samplei).selectedID-1).HEKAfname;
            %%
%            medfiltlength=5;
           
            
            %movementdata
            a=dir([handles.data.dirs.videodir,'movement/',fname,'.mat']);
            handles.data.samples(samplei).movementnames={};
            if ~isempty(a)
                load([handles.data.dirs.videodir,'movement/',fname,'.mat']);
                %%
                movement=[];
                time=[];
                video=[];
                videofname={};
                framenumber=[];
                videoidx=[];
                for i=1:length(videodata)
                    video=cat(3,video,rawvideo(i).vid);
                    framenumber=[framenumber;[1:length(videodata(i).movement_selected)]'];
                    videoidx=[videoidx;ones(length(videodata(i).movement_selected),1)*i];
                    movement=[movement;videodata(i).movement_selected];
                    time=[time;videodata(i).time];
                    [time,ix]=sort(time);
                    movement=movement(ix);
                    video=video(:,:,ix);
                    videofname{i}=[videodata(i).timestamp,'.avi'];
                end
                si=mode(diff(time));
                movement=moving(movement,2);
                movement=medfilt1(movement,round(.5/si));
                a=dir([handles.data.dirs.videodir,'percentiles.mat']);
                if ~isempty(a)
                    load([handles.data.dirs.videodir,'percentiles.mat'])
                    minval=videopercentiles.movementpercentiles(5);
                    maxval=videopercentiles.movementpercentiles(95);
                else
                    movementsorted=sort(movement);
                    minval=movementsorted(round(length(movementsorted)*.05));
                    maxval=movementsorted(round(length(movementsorted)*.95));
                end
                movement=(movement-minval)/(maxval-minval);
                movement(movement>1)=1;
                movement(movement<0)=0;
                handles.data.samples(samplei).movementdata.movement=movement+1;
                handles.data.samples(samplei).movementdata.time=time;
                handles.data.samples(samplei).movementdata.video=video;
                handles.data.samples(samplei).movementdata.videofname=videofname;
                handles.data.samples(samplei).movementdata.videoidx=videoidx;
                handles.data.samples(samplei).movementdata.framenumber=framenumber;
                handles.data.samples(samplei).movementnames={'movement'};
            end
            %%
            a=dir([handles.data.dirs.videodir,'ROIs/',fname,'.mat']);
            if ~isempty(a)
                load([handles.data.dirs.videodir,'ROIs/',fname,'.mat']);
                for ROIi=1:length(ROIdata)
                    if isfield(ROIdata,'PCA') & isfield(ROIdata(ROIi).PCA,'bestPC')
                        ROIname=ROIdata(ROIi).ROIname;
                        movement=[];
                        time=[];
                        for i=1:length(videodata)
                            movement=[movement;ROIdata(ROIi).PCA(i).bestPC];
                            time=[time;videodata(i).time];
                            [time,ix]=sort(time);
                            movement=movement(ix);
                        end
                        si=mode(diff(time));
                        movement=moving(abs(movement),2);
                        movement=medfilt1(movement,round(.5/si));
                        movementsorted=sort(movement);
                        minval=movementsorted(round(length(movementsorted)*.05));
                        maxval=movementsorted(round(length(movementsorted)*.95));
                        movement=(movement-minval)/(maxval-minval);
                        movement(movement>1)=1;
                        movement(movement<0)=0;
                        handles.data.samples(samplei).movementdata.(ROIname)=movement+1;
                        handles.data.samples(samplei).movementnames=[handles.data.samples(samplei).movementnames,ROIname];
                        %%
                        if strcmp(ROIname,'Body') & isfield(ROIdata(ROIi).PCA,'breath')
                            ROIname='breath';
                            movement=[];
                            time=[];
                            for i=1:length(videodata)
                                movement=[movement;ROIdata(ROIi).PCA(i).breath];
                                time=[time;videodata(i).time];
                                [time,ix]=sort(time);
                                movement=movement(ix);
                            end
                            si=mode(diff(time));
                            movement=moving((movement),2);
%                             movement=medfilt1(movement,round(.5/si));
                            movementsorted=sort(movement);
                            minval=movementsorted(round(length(movementsorted)*.05));
                            maxval=movementsorted(round(length(movementsorted)*.95));
                            movement=(movement-minval)/(maxval-minval);
                            movement(movement>1)=1;
                            movement(movement<0)=0;
                            handles.data.samples(samplei).movementdata.(ROIname)=movement+1;
                            handles.data.samples(samplei).movementnames=[handles.data.samples(samplei).movementnames,ROIname];
                        end
                    end
                end
                
                % pupildata
            a=dir([handles.data.dirs.videodir,'eye/',fname,'.mat']);
            if ~isempty(a)
                load([handles.data.dirs.videodir,'eye/',fname]);
                handles.data.samples(samplei).pupildata=pupildata;
                diameter=handles.data.samples(samplei).pupildata.diameter;
                si=mode(diff(pupildata.time));
                if si==0
                    si=1/15;
                end
                diameter=moving(diameter,2);
                diameter=medfilt1(diameter,round(10/si));
                
                a=dir([handles.data.dirs.videodir,'percentiles.mat']);
                if ~isempty(a)
                    load([handles.data.dirs.videodir,'percentiles.mat'])
                    if exist('ROIdata','var') & any(strcmp({ROIdata.ROIname},'Eye'))
                        stats=regionprops(ROIdata(find(strcmp({ROIdata.ROIname},'Eye'))).mask,'MajorAxisLength');
                        minval=videopercentiles.pupilpercentiles(5)*stats.MajorAxisLength;
                        maxval=videopercentiles.pupilpercentiles(95)*stats.MajorAxisLength;
                    else
                        minval=videopercentiles.pupilpercentiles(5);
                        maxval=videopercentiles.pupilpercentiles(95);
                    end
                else
                    pupildatasorted=sort(diameter);
                    minval=pupildatasorted(round(length(pupildatasorted)*.05));
                    maxval=pupildatasorted(round(length(pupildatasorted)*.95));
                end
                diameter=(handles.data.samples(samplei).pupildata.diameter-minval)/(maxval-minval);
                diameter(diameter>1)=1;
                diameter(diameter<0)=0;
                
                handles.data.samples(samplei).pupildata.diameter=diameter;
            end

            end
            set(handles.popupmenu8,'Value',1);
            set(handles.popupmenu8,'String',handles.data.samples(samplei).movementnames);
        end
        if isfield(handles.data.dirs,'PSDdir')
            xlsnum=handles.data.samples(samplei).selectedID-1;
            fieldxlsnum=find(strcmp(handles.data.xlsdata(xlsnum).HEKAfname,{handles.data.xlsdata.HEKAfname}) & [handles.data.xlsdata.field]==1);
            if ~isempty(fieldxlsnum)
                load([handles.data.dirs.PSDdir,handles.data.xlsdata(fieldxlsnum).ID]);
                
                
%                 for sweepnum=1:length(PSDdata)
%                     y=PSDdata(sweepnum).y;
%                     if ~isempty(PSDdata(sweepnum).y)
%                         [b,a]=butter(1,.2/(1/si)/2,'high');
%                         ytofilt=[y(end:-1:1),y,y(end:-1:1)];
%                         ytofilt=filtfilt(b,a,ytofilt);
%                         y=ytofilt([1:length(y)]+length(y));
%                         y=(y-min(y));
%                         ysort=sort(y);
%                         y=y/ysort(round(length(y)*.99));
%                         %                         y=y/max(y);
%                         y(y>1)=1;
%                         PSDdata(sweepnum).y=y;
%                     end
%                 end
                handles.data.samples(samplei).PSDdata=PSDdata;
                set(handles.edit9,'String',max(PSDdata(1).powerMatrix(:)));
                set(handles.edit10,'String',8); %max(PSDdata(1).frequencyVector(:))
                set(handles.edit11,'String',min(PSDdata(1).frequencyVector(:)));
            else
                handles.data.samples(samplei).PSDdata=[];
            end

                
        end
        
        if isfield(handles.data.dirs,'brainstatedir')
            %%
            ID=handles.data.IDs{handles.data.samples(samplei).loadedID};
            a=dir([handles.data.dirs.brainstatedir,ID,'.mat']);
            if isempty(a)
                handles.data.samples(samplei).BrainStateData=struct;
            else
                load([handles.data.dirs.brainstatedir,ID],'BrainStateData');
                handles.data.samples(samplei).BrainStateData=BrainStateData;
            end
        end
    end
end

function plotandupdate(handles)

  starttime=str2num(get(handles.edit4,'String'));
        endtime=str2num(get(handles.edit5,'String'));
 handles=updatedatatoplot(handles);
selectedsamplenum=get(handles.popupmenu1,'Value');
set(handles.checkbox1,'Value',handles.data.samples(selectedsamplenum).switch);
selectedID=handles.data.samples(selectedsamplenum).selectedID;
selectedvariable=handles.data.samples(selectedsamplenum).selectedvariable;
% set(handles.popupmenu2,'Value',selectedID);
set(handles.popupmenu3,'Value',selectedvariable);
set(handles.edit1,'String',num2str(handles.data.samples(selectedsamplenum).cutoffreq));
set(handles.edit2,'String',num2str(handles.data.samples(selectedsamplenum).filterdeg));
selectedvariablename=handles.data.fieldnevek{selectedvariable};
if length(handles.data.samples(selectedsamplenum).lightdata)>1
    if ischar(handles.data.samples(selectedsamplenum).lightdata(1).(selectedvariablename))
        reducednames=unique({handles.data.samples(selectedsamplenum).lightdata.(selectedvariablename)});
        ezszoveg=1;
    else
        reducednames=unique([handles.data.samples(selectedsamplenum).lightdata.(selectedvariablename)]);
        ezszoveg=0;
    end
else
    reducednames=[];
    ezszoveg=0;
end
set(handles.listbox1,'String',reducednames);
set(handles.text1,'String',handles.data.samples(selectedsamplenum).changes)

axes(handles.axes1)
cla 
hold on
handles.axes1=gca;

axes(handles.axes2)
cla 
hold on
handles.axes2=gca;

axes(handles.axes3)
cla 
hold on
handles.axes3=gca;

for samplenum=1:length(handles.data.samples)
    y0=0;
    if handles.data.samples(samplenum).switch==1 & handles.data.samples(samplenum).selectedID>1 & ~isempty(handles.data.samples(samplenum).neededwaves)
        neededwaves=handles.data.samples(samplenum).neededwaves;
        marker=handles.data.samples(samplenum).marker;
        axes(handles.axes1)
      
        for sweepi=1:length(neededwaves)
            ettol=find(handles.data.samples(samplenum).datatoplot(sweepi).x>starttime,1,'first');
            eddig=find(handles.data.samples(samplenum).datatoplot(sweepi).x<endtime,1,'last');
            plot(handles.data.samples(samplenum).datatoplot(sweepi).x(ettol:eddig),handles.data.samples(samplenum).datatoplot(sweepi).yvoltage(ettol:eddig),marker(1),'LineWidth',2)
        end
        markevents=get(handles.checkbox3,'Value');
        if markevents==1 & ~isempty(handles.data.samples(samplenum).eventdata) %| isempty(fieldnames(handles.data.samples(samplenum).eventdata)))
            neededevents=[handles.data.samples(samplenum).eventdata.axonalAP]==1;%strcmp({handles.data.samples(samplenum).eventdata.type},'AP') & ~[handles.data.samples(samplenum).eventdata.stimulated];
            eventsnow=handles.data.samples(samplenum).eventdata(neededevents);
            neededevents=false(size(eventsnow));
              for sweepi=1:length(neededwaves)
                 sweepnum=neededwaves(sweepi);
                 neededevents([eventsnow.sweepnum]==sweepnum)=1;
              end
             plot([eventsnow(neededevents).maxtime],[eventsnow(neededevents).maxval],'ro','MarkerSize',8)
        end
        if get(handles.checkbox5,'Value')==1 & ~isempty(handles.data.samples(samplenum).eventdata) %| isempty(fieldnames(handles.data.samples(samplenum).eventdata)))
            neededevents=[handles.data.samples(samplenum).eventdata.somaticAP]==1;%strcmp({handles.data.samples(samplenum).eventdata.type},'AP') & ~[handles.data.samples(samplenum).eventdata.stimulated];
            eventsnow=handles.data.samples(samplenum).eventdata(neededevents);
            neededevents=false(size(eventsnow));
              for sweepi=1:length(neededwaves)
                 sweepnum=neededwaves(sweepi);
                 neededevents([eventsnow.sweepnum]==sweepnum)=1;
              end
             plot([eventsnow(neededevents).maxtime],[eventsnow(neededevents).maxval],'ko','MarkerSize',8)
        end
        if get(handles.checkbox6,'Value')==1 & ~isempty(handles.data.samples(samplenum).eventdata) %| isempty(fieldnames(handles.data.samples(samplenum).eventdata)))
            neededevents=strcmp({handles.data.samples(samplenum).eventdata.type},'ep');%strcmp({handles.data.samples(samplenum).eventdata.type},'AP') & ~[handles.data.samples(samplenum).eventdata.stimulated];
            eventsnow=handles.data.samples(samplenum).eventdata(neededevents);
            neededevents=false(size(eventsnow));
              for sweepi=1:length(neededwaves)
                 sweepnum=neededwaves(sweepi);
                 neededevents([eventsnow.sweepnum]==sweepnum)=1;
              end
             plot([eventsnow(neededevents).onsettime],[eventsnow(neededevents).baselineval],'rs','MarkerSize',8) 
             plot([eventsnow(neededevents).maxtime],[eventsnow(neededevents).maxval],'r^','MarkerSize',8)
        end
         if get(handles.checkbox7,'Value')==1 & ~isempty(handles.data.samples(samplenum).eventdata) %| isempty(fieldnames(handles.data.samples(samplenum).eventdata)))
            neededevents=strcmp({handles.data.samples(samplenum).eventdata.type},'ip');%strcmp({handles.data.samples(samplenum).eventdata.type},'AP') & ~[handles.data.samples(samplenum).eventdata.stimulated];
            eventsnow=handles.data.samples(samplenum).eventdata(neededevents);
            neededevents=false(size(eventsnow));
              for sweepi=1:length(neededwaves)
                 sweepnum=neededwaves(sweepi);
                 neededevents([eventsnow.sweepnum]==sweepnum)=1;
              end
             plot([eventsnow(neededevents).onsettime],[eventsnow(neededevents).baselineval],'bs','MarkerSize',8) 
             plot([eventsnow(neededevents).maxtime],[eventsnow(neededevents).maxval],'bv','MarkerSize',8)
        end
       markstates=get(handles.checkbox4,'Value');
       if markstates==1 & isfield(handles.data.dirs,'statedir')
           %%
           if ~isfield(handles.data.samples(samplenum),'statedata')
               load([handles.data.dirs.statedir,handles.data.IDs{handles.data.samples(samplenum).loadedID}],'statedata');
               handles.data.samples(samplenum).statedata=statedata;
           end
            for sweepi=1:length(neededwaves)
                 sweepnum=neededwaves(sweepi);
                 neededup=find([handles.data.samples(samplenum).statedata.UP.sweepnum]==sweepnum);
                 for upi=1:length(neededup)
                     onseth=handles.data.samples(samplenum).statedata.UP(neededup(upi)).onseth;
                     endh=handles.data.samples(samplenum).statedata.UP(neededup(upi)).endh;
                     voltagevalue=median(handles.data.samples(samplenum).datatoplot(sweepi).yvoltage(onseth:endh))-.001;
                     plot(handles.data.samples(samplenum).datatoplot(sweepi).x([onseth,endh]),[voltagevalue,voltagevalue],'r-','LineWidth',4)
                 end
                   neededdown=find([handles.data.samples(samplenum).statedata.DOWN.sweepnum]==sweepnum);
                 for downi=1:length(neededdown)
                     onseth=handles.data.samples(samplenum).statedata.DOWN(neededdown(downi)).onseth;
                     endh=handles.data.samples(samplenum).statedata.DOWN(neededdown(downi)).endh;
                     voltagevalue=median(handles.data.samples(samplenum).datatoplot(sweepi).yvoltage(onseth:endh))-.001;
                     plot(handles.data.samples(samplenum).datatoplot(sweepi).x([onseth,endh]),[voltagevalue,voltagevalue],'b-','LineWidth',4)
                 end
            end
       end
       markbrainstates=get(handles.checkbox8,'Value');
       if markbrainstates==1 & isfield(handles.data.dirs,'brainstatedir')
           %%
           BrainStateProps=handles.data.BrainStateProps;
           BrainStateData=handles.data.samples(samplenum).BrainStateData;
           if ~isempty(fieldnames(BrainStateData))
               ylimits=[min([handles.data.samples(samplenum).datatoplot.yvoltage]),max([handles.data.samples(samplenum).datatoplot.yvoltage])];
               for statei=1:length(BrainStateData)
                   if isempty(find(strcmp(BrainStateData(statei).name,BrainStateProps.Statevalues)))
                       color='k';
                   else
                       color=BrainStateProps.Statecolors_short{find(strcmp(BrainStateData(statei).name,BrainStateProps.Statevalues))};
                   end
                   rectangle('Position',[BrainStateData(statei).starttime,ylimits(1),BrainStateData(statei).endtime-BrainStateData(statei).starttime,diff(ylimits)],'EdgeColor',color,'LineWidth',2);
                   text(BrainStateData(statei).starttime,ylimits(1)+.95*diff(ylimits),BrainStateData(statei).name,'Color',color)
               end
           end
       end
       
%         plot([handles.data.samples(samplenum).datatoplot.x],[handles.data.samples(samplenum).datatoplot.yvoltage],marker(1),'LineWidth',2)
        y0=max([y0,max([handles.data.samples(samplenum).datatoplot.yvoltage])]);
        if length(neededwaves)>0
            prenum=handles.data.samples(samplenum).selectedID-1;
            if isfield(handles.data.xlsdata,'drugdata')
            for i=1:length(handles.data.xlsdata(prenum).drugdata)
                y0=y0+.005;
                startt=handles.data.xlsdata(prenum).drugdata(i).DrugWashinTime;
                origstartt=startt;
                endt=handles.data.xlsdata(prenum).drugdata(i).DrugWashoutTime;
                originalednt=endt;
                mintval=nanmin([handles.data.samples(samplenum).datatoplot.x]);
                maxtval=nanmax([handles.data.samples(samplenum).datatoplot.x]);
                if startt<mintval;
                    startt=mintval;
                end
                if isempty(endt)
                    endt=maxtval;
                end
                if endt<startt
                    endt=startt;
                end
                plot([startt, endt],[y0,y0],[marker(1),'-'],'LineWidth', 10)
                text(startt,y0,[handles.data.xlsdata(prenum).drugdata(i).DrugName,' - washin at ',num2str(round(origstartt))],'Color',[1 1 1])
            end
            end
        end
        handles.axes1=gca;
        axesname='axes2';
        popupname='popupmenu4';
        plotlowerpanels(marker,samplenum,starttime,endtime,neededwaves,axesname,popupname,handles);
        axesname='axes3';
        popupname='popupmenu7';
        plotlowerpanels(marker,samplenum,starttime,endtime,neededwaves,axesname,popupname,handles);
    end
end
axes(handles.axes3)
xlabel('Time (s)')
% axis tight
% xlim([starttime endtime])
handles.axes3=gca;

axes(handles.axes2)
% axis tight
xlabel([]);
set(gca,'xtick',[]);

handles.axes2=gca;

axes(handles.axes1)
axis tight
if starttime~=endtime
    xlim([starttime endtime])
end
xlabel([]);
set(gca,'xtick',[]);
ylabel('Voltage (V)')
handles.axes1=gca; 


linkaxes([handles.axes1,handles.axes2, handles.axes3],'x');
java.lang.System.gc()

function plotlowerpanels(marker,samplenum,starttime,endtime,neededwaves,axesname,popupname,handles)
axes(handles.(axesname))
if get(handles.(popupname),'Value')-1==0
    for sweepi=1:length(neededwaves)
        ettol=find(handles.data.samples(samplenum).datatoplot(sweepi).x>starttime,1,'first');
        eddig=find(handles.data.samples(samplenum).datatoplot(sweepi).x<endtime,1,'last');
        plot(handles.data.samples(samplenum).datatoplot(sweepi).x(ettol:eddig),handles.data.samples(samplenum).datatoplot(sweepi).ycurrent(ettol:eddig),marker(1),'LineWidth',2)
    end
%     xlabel('Time (s)')
    ylabel('Injected Current (pA)')
    axis tight
elseif samplenum==get(handles.(popupname),'Value')-1
    for sweepi=1:length(neededwaves)
        ettol=find(handles.data.samples(samplenum).datatoplot(sweepi).x>starttime,1,'first');
        eddig=find(handles.data.samples(samplenum).datatoplot(sweepi).x<endtime,1,'last');
        plot(handles.data.samples(samplenum).datatoplot(sweepi).x(ettol:eddig),handles.data.samples(samplenum).datatoplot(sweepi).yvoltage(ettol:eddig),marker(1),'LineWidth',2)
    end
    axis tight
    ylabel('Voltage (mV)')
elseif get(handles.(popupname),'Value')==length(handles.data.samples)+2
    %%
    hold on
    if isfield(handles.data.samples(samplenum),'pupildata')
        ettol=find(handles.data.samples(samplenum).pupildata.time>starttime,1,'first');
        eddig=find(handles.data.samples(samplenum).pupildata.time<endtime,1,'last');
        plot(handles.data.samples(samplenum).pupildata.time(ettol:eddig),handles.data.samples(samplenum).pupildata.diameter(ettol:eddig),'b-','LineWidth',2);
    end
    if isfield(handles.data.samples(samplenum),'movementdata')
        possiblefields=get(handles.popupmenu8,'String');
        movementfield=possiblefields{get(handles.popupmenu8,'Value')};
        ettol=find(handles.data.samples(samplenum).movementdata.time>starttime,1,'first');
        eddig=find(handles.data.samples(samplenum).movementdata.time<endtime,1,'last');
        plot(handles.data.samples(samplenum).movementdata.time(ettol:eddig),handles.data.samples(samplenum).movementdata.(movementfield)(ettol:eddig),'r-','LineWidth',2);
        
    end
    ylabel('\color{red}Movement \color{black}and \color{blue}pupil size \color{black}(AU)')
    ylim([0 2])
elseif get(handles.(popupname),'Value')==length(handles.data.samples)+3 & ~isempty(handles.data.samples(samplenum).PSDdatatoplot)
    
    imagesc(handles.data.samples(samplenum).PSDdatatoplot.time,handles.data.samples(samplenum).PSDdatatoplot.frequencyVector,handles.data.samples(samplenum).PSDdatatoplot.powerMatrix)
    
    
    caxis([0 str2num(get(handles.edit9,'String'))])
    set(gca,'YDir','normal');
    colormap linspecer
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    hold on
    if ~isempty(handles.data.samples(samplenum).PSDdatatoplot.frequencyVector)
        szorzo=max(handles.data.samples(samplenum).PSDdatatoplot.frequencyVector);
    else
        szorzo=1;
    end
    for sweepi=1:length(handles.data.samples(samplenum).PSDdatatoplot.trace)
        ettol=find(handles.data.samples(samplenum).PSDdatatoplot.trace(sweepi).x>starttime,1,'first');
        eddig=find(handles.data.samples(samplenum).PSDdatatoplot.trace(sweepi).x<endtime,1,'last');
        plot(handles.data.samples(samplenum).PSDdatatoplot.trace(sweepi).x(ettol:eddig),handles.data.samples(samplenum).PSDdatatoplot.trace(sweepi).y(ettol:eddig)*szorzo/3+2*szorzo/3,'Color',[0.5 .5 .5],'LineWidth',.5)
        
    end
%     disp('muahah')
    axis tight
end

markbrainstates=get(handles.checkbox8,'Value');
if markbrainstates==1 & isfield(handles.data.dirs,'brainstatedir')
    %%
    BrainStateProps=handles.data.BrainStateProps;
    BrainStateData=handles.data.samples(samplenum).BrainStateData;
    if ~isempty(fieldnames(BrainStateData))
%         ylimits=[min([handles.data.samples(samplenum).datatoplot.yvoltage]),max([handles.data.samples(samplenum).datatoplot.yvoltage])];
        ylimits=get(gca,'Ylim');
        for statei=1:length(BrainStateData)
            if isempty(find(strcmp(BrainStateData(statei).name,BrainStateProps.Statevalues)))
                color='k';
            else
                color=BrainStateProps.Statecolors_short{find(strcmp(BrainStateData(statei).name,BrainStateProps.Statevalues))};
            end
            rectangle('Position',[BrainStateData(statei).starttime,ylimits(1),BrainStateData(statei).endtime-BrainStateData(statei).starttime,diff(ylimits)],'EdgeColor',color,'LineWidth',2);
            text(BrainStateData(statei).starttime,ylimits(1)+.95*diff(ylimits),BrainStateData(statei).name,'Color',color)
        end
    end
end

%         plot([handles.data.samples(samplenum).datatoplot.x],[handles.data.samples(samplenum).datatoplot.ycurrent],marker(1),'LineWidth',2)
handles.(axesname)=gca;


function handles=updategui(handles)
selectedsamplenum=get(handles.popupmenu1,'Value');
set(handles.checkbox1,'Value',handles.data.samples(selectedsamplenum).switch);
% selectedID=handles.data.samples(selectedsamplenum).selectedID;
selectedvariable=handles.data.samples(selectedsamplenum).selectedvariable;
% set(handles.popupmenu2,'Value',selectedID);
set(handles.popupmenu3,'Value',selectedvariable);
set(handles.edit1,'String',num2str(handles.data.samples(selectedsamplenum).cutoffreq));
set(handles.edit2,'String',num2str(handles.data.samples(selectedsamplenum).filterdeg));
startval=get(handles.slider1,'Value');
endval=get(handles.slider2,'Value');
midval=mean([startval,endval]);
set(handles.slider3,'Value',midval);




if isfield(handles.data.samples,'starttime')
    starttime=min([handles.data.samples(selectedsamplenum).starttime]);
    endtime=max([handles.data.samples(selectedsamplenum).endtime]);
    timedifi=handles.data.samples(selectedsamplenum).endtime-handles.data.samples(selectedsamplenum).starttime;
    % starttime=min([handles.data.samples.starttime]);
    % endtime=max([handles.data.samples.endtime]);
    % timedifi=endtime-starttime;
    set(handles.edit4,'String',starttime+timedifi*startval);
    set(handles.edit5,'String',starttime+timedifi*endval);
    set(handles.edit7,'String',starttime+timedifi*midval);
    set(handles.edit8,'String',timedifi*(endval-startval));
end
selectedvariablename=handles.data.fieldnevek{selectedvariable};
if length(handles.data.samples(selectedsamplenum).lightdata)>1
    if ischar(handles.data.samples(selectedsamplenum).lightdata(1).(selectedvariablename))
        reducednames=unique({handles.data.samples(selectedsamplenum).lightdata.(selectedvariablename)});
        ezszoveg=1;
    else
        reducednames=unique([handles.data.samples(selectedsamplenum).lightdata.(selectedvariablename)]);
        ezszoveg=0;
    end
else
    reducednames=[];
    ezszoveg=0;
end
set(handles.listbox1,'String',reducednames);
set(handles.text1,'String',handles.data.samples(selectedsamplenum).changes)

    if ~isempty(handles.data.samples(selectedsamplenum).eventdata)
        stimapnum=length(find(strcmp({handles.data.samples(selectedsamplenum).eventdata.type},'AP') & [handles.data.samples(selectedsamplenum).eventdata.stimulated]));
        spontapnum=length(find(strcmp({handles.data.samples(selectedsamplenum).eventdata.type},'AP') & ~[handles.data.samples(selectedsamplenum).eventdata.stimulated]));
        if isfield(handles.data.samples(selectedsamplenum).eventdata,'axonalAP')
           stimaAPnum=length(find([handles.data.samples(selectedsamplenum).eventdata.axonalAP] & [handles.data.samples(selectedsamplenum).eventdata.stimulated]));
           stimsAPnum=length(find([handles.data.samples(selectedsamplenum).eventdata.somaticAP] & [handles.data.samples(selectedsamplenum).eventdata.stimulated]));
           spontaAPnum=length(find([handles.data.samples(selectedsamplenum).eventdata.axonalAP] & ~[handles.data.samples(selectedsamplenum).eventdata.stimulated]));
           spontsAPnum=length(find([handles.data.samples(selectedsamplenum).eventdata.somaticAP] & ~[handles.data.samples(selectedsamplenum).eventdata.stimulated]));
        end
    else
        stimapnum=0;
        spontapnum=0;
    end
    
    if handles.data.samples(selectedsamplenum).selectedID-1 > 0 & isfield(handles.data.xlsdata,'drugdata') & ~isempty(handles.data.xlsdata(handles.data.samples(selectedsamplenum).selectedID-1).drugdata)
        drugnames={handles.data.xlsdata(handles.data.samples(selectedsamplenum).selectedID-1).drugdata.DrugName};
        drugwashintimes=[handles.data.xlsdata(handles.data.samples(selectedsamplenum).selectedID-1).drugdata.DrugWashinTime];
        drugstring=' ';
        for i=1:length(drugnames)
            drugstring=[drugstring,drugnames{i},'-',num2str(drugwashintimes(i)),';  '];
        end
    else
        drugstring='   no drugs';
    end
    if isfield(handles.data.samples(selectedsamplenum).eventdata,'axonalAP')
        set(handles.text5,'String',[handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID},'  -  stimulated: ',num2str(stimaAPnum),' aAP and ',num2str(stimsAPnum),' sAP','  spontaneous: ',num2str(spontaAPnum),' aAP and ',num2str(spontsAPnum),' sAP','       ',drugstring]);
    else
        set(handles.text5,'String',[handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID},'  -  stimulated AP: ', num2str(stimapnum),'  spontaneous AP: ', num2str(spontapnum),'       ',drugstring]);
    end
    xlsdata=handles.data.xlsdata;
    fieldek=get(handles.popupmenu6,'String');
    neededfield=fieldek{get(handles.popupmenu6,'Value')};
    actualIDs=handles.data.actualIDs;
    actualxlsidxs=[];
    if get(handles.popupmenu6,'Value')>1
        for i=2:length(actualIDs)
            actualxlsidxs(i-1)=find(strcmp({xlsdata.ID},actualIDs{i}));
        end
        if isnumeric([xlsdata(1).(neededfield)])
            [~,neworder]=sort([xlsdata(actualxlsidxs).(neededfield)]);
            handles.data.actualIDs=[actualIDs{1},actualIDs(neworder+1)];
        else
            [~,neworder]=sort({xlsdata(actualxlsidxs).(neededfield)});
            handles.data.actualIDs=[actualIDs(1),actualIDs(neworder+1)];
        end
        
    end
    
    fieldstoadd={'anaesthesia','field'};
    IDstring={};
    actualIDs=handles.data.actualIDs;
    for i=1:length(actualIDs)
        xlsidx=find(strcmp({handles.data.xlsdata.ID},actualIDs{i}));
        IDstring{i}=actualIDs{i};
        for fieldi=1:length(fieldstoadd);
            if isfield(handles.data.xlsdata,fieldstoadd{fieldi})
                if i>1 & isnumeric(handles.data.xlsdata(xlsidx).(fieldstoadd{fieldi}))
                    if handles.data.xlsdata(xlsidx).(fieldstoadd{fieldi})==1
                        IDstring{i}=[IDstring{i} ,' - ', fieldstoadd{fieldi}];
                    elseif handles.data.xlsdata(xlsidx).(fieldstoadd{fieldi})==0
                    else
                        IDstring{i}=[IDstring{i} ,' - ', num2str(handles.data.xlsdata(xlsidx).(fieldstoadd{fieldi}))];
                    end
                else
                    IDstring{i}=[IDstring{i} ,' - ', handles.data.xlsdata(xlsidx).(fieldstoadd{fieldi})];
                end
            end
        end
    end
    set(handles.popupmenu2,'String',IDstring);
    




function handles=updatedatatoplot(handles)
% selectedsamplenum=get(handles.popupmenu1,'Value');
for selectedsamplenum=1:5
if handles.data.samples(selectedsamplenum).switch>0
    neededwaves=handles.data.samples(selectedsamplenum).maybeneededwaves;
    %
    if ~isempty(neededwaves)
        timesofneededwaves=[handles.data.samples(selectedsamplenum).lightdata(handles.data.samples(selectedsamplenum).maybeneededwaves).realtime];
        reallyneededwaves=timesofneededwaves>=str2num(get(handles.edit4,'String')) & timesofneededwaves<=str2num(get(handles.edit5,'String'));
        if reallyneededwaves(1)==0
            reallyneededwaves(find(timesofneededwaves<=str2num(get(handles.edit4,'String')),1,'last'))=true;
        end
        neededwaves=neededwaves(reallyneededwaves);
        handles.data.samples(selectedsamplenum).neededwaves=neededwaves;
    end
    %
    neededsamplenum=str2num(get(handles.edit3,'String'));
    cutoffreq=handles.data.samples(selectedsamplenum).cutoffreq;
    filterdeg=handles.data.samples(selectedsamplenum).filterdeg;
    shoulddownsample=get(handles.checkbox2,'Value');
    minfreq=str2num(get(handles.edit11,'string'));
    maxfreq=str2num(get(handles.edit10,'string'));
    
    plotdetails=handles.data.samples(selectedsamplenum).plotdetails;
    if ~isempty(neededwaves)
        if isempty(fieldnames(plotdetails)) | ~(length(neededwaves)==length(plotdetails.neededwaves)) | ~any(neededwaves==plotdetails.neededwaves) | ~(neededsamplenum==plotdetails.neededsamplenum) | ~all(cutoffreq==plotdetails.cutoffreq) | ~(filterdeg==plotdetails.filterdeg) | ~(shoulddownsample==plotdetails.shoulddownsample) | ~(minfreq==plotdetails.minfreq) | ~(maxfreq==plotdetails.maxfreq)
            datatoplot=struct;
            PSDdatatoplot=struct;
            PSDdatatoplot.trace=struct;
            plotdetails.neededwaves=neededwaves;
            plotdetails.neededsamplenum=neededsamplenum;
            plotdetails.cutoffreq=cutoffreq;
            plotdetails.filterdeg=filterdeg;
            plotdetails.shoulddownsample=shoulddownsample;
            plotdetails.minfreq=minfreq;
            plotdetails.maxfreq=maxfreq;
            bigmatrix=[];
            for sweepi=1:length(neededwaves)
                sweepnum=neededwaves(sweepi);
                datatoplot(sweepi).ycurrent=handles.data.samples(selectedsamplenum).stimdata(sweepnum).y;
                hossz=length(datatoplot(sweepi).ycurrent);
                start=handles.data.samples(selectedsamplenum).bridgeddata(sweepnum).realtime;
                si=handles.data.samples(selectedsamplenum).bridgeddata(sweepnum).si;
                datatoplot(sweepi).x=[0:si:si*(hossz-1)]+start;
                if any(cutoffreq>0) & filterdeg>0
                    if length(cutoffreq)==1
                        [b,a]=butter(handles.data.samples(selectedsamplenum).filterdeg,handles.data.samples(selectedsamplenum).cutoffreq/(1/si)/2,'low');
                    else
                        [b,a]=butter(handles.data.samples(selectedsamplenum).filterdeg,handles.data.samples(selectedsamplenum).cutoffreq/(1/si)/2,'bandpass');
                    end
                    y=filtfilt(b,a,handles.data.samples(selectedsamplenum).bridgeddata(sweepnum).y);
                else
                    y=handles.data.samples(selectedsamplenum).bridgeddata(sweepnum).y;
                end
                datatoplot(sweepi).yvoltage=y;
            end
            if ~isempty(fieldnames(datatoplot))
                samplenumnow=length([datatoplot.yvoltage]);
                if neededsamplenum<samplenumnow & neededsamplenum>0 & shoulddownsample==1
                    %                 disp('downsampling started')
                    ratio=round(samplenumnow/neededsamplenum);
                    for sweepi=1:length(neededwaves)
                        datatoplot(sweepi).yvoltage=downsample(datatoplot(sweepi).yvoltage,ratio);
                        datatoplot(sweepi).ycurrent=downsample(datatoplot(sweepi).ycurrent,ratio);
                        datatoplot(sweepi).x=downsample(datatoplot(sweepi).x,ratio);
                    end
                    %                 disp('downsampling finished')
                end
            end
            handles.data.samples(selectedsamplenum).datatoplot=datatoplot;
            handles.data.samples(selectedsamplenum).plotdetails=plotdetails;
            %%
            if ~isempty(handles.data.samples(selectedsamplenum).PSDdata) & (get(handles.popupmenu4,'Value')==length(handles.data.samples)+3 | get(handles.popupmenu7,'Value')==length(handles.data.samples)+3)
                timesofPSDwaves=[handles.data.samples(selectedsamplenum).PSDdata.realtime];
                reallyneededwaves=timesofPSDwaves>=str2num(get(handles.edit4,'String')) & timesofPSDwaves<=str2num(get(handles.edit5,'String'));
                if reallyneededwaves(1)==0
                    reallyneededwaves(find(timesofPSDwaves<=str2num(get(handles.edit4,'String')),1,'last'))=true;
                end
                neededwaves=find(reallyneededwaves);
                minfreq=str2num(get(handles.edit11,'string'));
                maxfreq=str2num(get(handles.edit10,'string'));
                
                
                for sweepi=1:length(neededwaves)
                    sweepnum=neededwaves(sweepi);
                    if isempty(bigmatrix)
                        frequencyVector=handles.data.samples(selectedsamplenum).PSDdata(sweepnum).frequencyVector;
                        neededfrequencies=frequencyVector>=minfreq & frequencyVector<=maxfreq;
                        frequencyVector=frequencyVector(neededfrequencies);
                        starttime=handles.data.samples(selectedsamplenum).PSDdata(sweepnum).realtime;
                        si=handles.data.samples(selectedsamplenum).PSDdata(sweepnum).si_powerMatrix;
                        bigmatrix=handles.data.samples(selectedsamplenum).PSDdata(sweepnum).powerMatrix(neededfrequencies,:);
                        
                    else
                        timeskipped=(handles.data.samples(selectedsamplenum).PSDdata(sweepnum).realtime-(starttime+size(bigmatrix,2)*si));
                        colstoadd=round(timeskipped/si);
                        if ~isempty(handles.data.samples(selectedsamplenum).PSDdata(sweepnum).powerMatrix)
                            bigmatrix=[bigmatrix,zeros(size(bigmatrix,1),colstoadd),handles.data.samples(selectedsamplenum).PSDdata(sweepnum).powerMatrix(neededfrequencies,:)];
                        end
                        %             starttime+size(bigmatrix,2)*si
                        %             return
                    end
                    
                    PSDdatatoplot.trace(sweepi).x=[0:si:si*(length(handles.data.samples(selectedsamplenum).PSDdata(sweepnum).y)-1)]+handles.data.samples(selectedsamplenum).PSDdata(sweepnum).realtime;
                    
                    if ~isempty(handles.data.samples(selectedsamplenum).PSDdata(sweepnum).y)
                        y=handles.data.samples(selectedsamplenum).PSDdata(sweepnum).y;
                        dofilt=0;
                        if dofilt==1
                            [b,a]=butter(1,.1/(1/si)/2,'high');
                            ytofilt=[y(end:-1:1),y,y(end:-1:1)];
                            ytofilt=filtfilt(b,a,ytofilt);
                            y=ytofilt([1:length(y)]+length(y));
                        end
                        y=(y-min(y));
                        ysort=sort(y);
                        y=y/ysort(round(length(y)*.99));
                        %                         y=y/max(y);
                        y(y>1)=1;
                        PSDdatatoplot.trace(sweepi).y=y;
                        PSDdatatoplot.trace(sweepi).realtime=handles.data.samples(selectedsamplenum).PSDdata(sweepnum).realtime;
                        PSDdatatoplot.trace(sweepi).si=handles.data.samples(selectedsamplenum).PSDdata(sweepnum).si_powerMatrix;
                    end
                end
                
                
                
                %%
                if ~isempty(fieldnames(PSDdatatoplot.trace))
                    samplenumnow=length([PSDdatatoplot.trace.y]);
                    if neededsamplenum<samplenumnow & neededsamplenum>0 & shoulddownsample==1
                        %                 disp('downsampling started')
                        ratio=round(samplenumnow/neededsamplenum);
                        for sweepi=1:length(neededwaves)
                            PSDdatatoplot.trace(sweepi).x=downsample(PSDdatatoplot.trace(sweepi).x,ratio);
                            PSDdatatoplot.trace(sweepi).y=downsample(PSDdatatoplot.trace(sweepi).y,ratio);
                        end
                        %                 disp('downsampling finished')
                    end
                end
               %%     
               
               
               
               
               bigtime=(1:size(bigmatrix,2))*si+starttime;
                if ~isempty(bigtime)
                    samplenumnow=length(bigtime);
                    if neededsamplenum<samplenumnow & neededsamplenum>0 & shoulddownsample==1
                        %                 disp('downsampling started')
                        ratio=round(samplenumnow/neededsamplenum);
                        bigmatrix=downsample(bigmatrix',ratio)';
                        bigtime=downsample(bigtime,ratio);
                      %                 disp('downsampling finished')
                    end
                end
                
               %% 
                PSDdatatoplot.powerMatrix=bigmatrix;
                PSDdatatoplot.frequencyVector=frequencyVector;
                PSDdatatoplot.si_powerMatrix=si;
                PSDdatatoplot.time=bigtime;
                if ~isempty(PSDdatatoplot.powerMatrix)
                    values=sort(PSDdatatoplot.powerMatrix(:));
                    PSDdatatoplot.intensitypercentiles=values(ceil([.01:.01:1]*length(values)));
                else
                    PSDdatatoplot.intensitypercentiles=ones(size([.01:.01:1]));
                end
                handles.data.samples(selectedsamplenum).PSDdatatoplot=PSDdatatoplot;
%                 set(handles.edit9,'String',handles.data.samples(selectedsamplenum).PSDdatatoplot.intensitypercentiles(99));
            else
                handles.data.samples(selectedsamplenum).PSDdatatoplot=[];
            end
        end
    end
    
end
end




function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
selectedsamplenum=get(handles.popupmenu1,'Value');
freqstring=get(hObject,'String');
if any(strfind(freqstring,','))
    idx=strfind(freqstring,',');
    freq(1)=(str2num(freqstring(1:idx-1)));
    freq(2)=(str2num(freqstring(idx+1:end)));
    handles.data.samples(selectedsamplenum).cutoffreq=[min(freq),max(freq)];
else
    handles.data.samples(selectedsamplenum).cutoffreq=(str2num(freqstring));
end
handles=loadthedata(handles);
handles=updatedatatoplot(handles);
guidata(hObject,handles)
%plotandupdate(handles)

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



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
selectedsamplenum=get(handles.popupmenu1,'Value');
handles.data.samples(selectedsamplenum).filterdeg=abs(round(str2num(get(hObject,'String'))));
handles=loadthedata(handles);
handles=updatedatatoplot(handles);
guidata(hObject,handles)
%plotandupdate(handles)

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
handles=loadthedata(handles);
handles=updatedatatoplot(handles);
%plotandupdate(handles)

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


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=updatedatatoplot(handles);
plotandupdate(handles)
guidata(hObject, handles)


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if get(hObject,'Value')>get(handles.slider2,'Value')
    set(hObject,'Value',get(handles.slider2,'Value'))
end
handles=updategui(handles);
handles=updatedatatoplot(handles);
plotandupdate(handles)
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if get(hObject,'Value')<get(handles.slider1,'Value')
    set(hObject,'Value',get(handles.slider1,'Value'))
end
handles=updategui(handles);
handles=updatedatatoplot(handles);
plotandupdate(handles)
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
minval=get(handles.slider1,'Value');
maxval=get(handles.slider2,'Value');
oldval=mean([minval,maxval]);
newval=get(handles.slider3,'Value' );
difi=newval-oldval;
minval=max([0,minval+difi]);
maxval=min([1,maxval+difi]);
set(handles.slider1,'Value',minval)
set(handles.slider2,'Value',maxval)
handles=updategui(handles);
handles=updatedatatoplot(handles);
plotandupdate(handles)
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
startval=get(handles.slider1,'Value');
endval=get(handles.slider2,'Value');
selectedsamplenum=get(handles.popupmenu1,'Value');
timedifi=handles.data.samples(selectedsamplenum).endtime-handles.data.samples(selectedsamplenum).starttime;
% timedifi=max([handles.data.samples.endtime])-min([handles.data.samples.starttime]);
timeentered=str2num(get(hObject,'String'));
startvalneeded=(timeentered-handles.data.samples(selectedsamplenum).starttime)/timedifi;
if startvalneeded<get(handles.slider1,'Min')
    startvalset=get(handles.slider1,'Min');
elseif startvalneeded>endval
    startvalset=endval;
elseif startvalneeded>get(handles.slider1,'Max')
    startvalset=get(handles.slider1,'Max');
else
    startvalset=startvalneeded;
end
set(handles.slider1,'Value',startvalset);
handles=updategui(handles);
guidata(hObject, handles)
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



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
startval=get(handles.slider1,'Value');
endval=get(handles.slider2,'Value');
selectedsamplenum=get(handles.popupmenu1,'Value');
timedifi=handles.data.samples(selectedsamplenum).endtime-handles.data.samples(selectedsamplenum).starttime;
% timedifi=max([handles.data.samples.endtime])-min([handles.data.samples.starttime]);
timeentered=str2num(get(hObject,'String'));
endvalneeded=(timeentered-handles.data.samples(selectedsamplenum).starttime)/timedifi;
if endvalneeded>get(handles.slider2,'Max')
    endvalset=get(handles.slider2,'Max');
elseif endvalneeded<startval
    endvalset=startval;
elseif endvalneeded<get(handles.slider2,'Min')
    endvalset=get(handles.slider2,'Min');
else
    endvalset=endvalneeded;
end
set(handles.slider2,'Value',endvalset);
handles=updategui(handles);
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
startval=get(handles.slider1,'Value');
endval=get(handles.slider2,'Value');
selectedsamplenum=get(handles.popupmenu1,'Value');
timedifi=handles.data.samples(selectedsamplenum).endtime-handles.data.samples(selectedsamplenum).starttime;
% timedifi=max([handles.data.samples.endtime])-min([handles.data.samples.starttime]);
timeentered=str2num(get(hObject,'String'));
valneeded=(timeentered-min([handles.data.samples.starttime]))/timedifi;
if valneeded>get(handles.slider3,'Max')
    valset=get(handles.slider3,'Max');
elseif valneeded<startval
    valset=startval;
else
    valset=valneeded;
end
set(handles.slider3,'Value',valset);
slider3_Callback(hObject,eventdata,handles);
handles=updategui(handles);
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)

% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')==true
    handles=checkforaAPs(handles);
    guidata(hObject, handles)
end
% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
handles=updategui(handles);
guidata(hObject, handles)

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


% --- Executes on selection change in popupmenu5.
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



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
xlsdata=handles.data.xlsdata;
searchstr=get(hObject,'String');
searchfieldnum=get(handles.popupmenu5,'Value');
searchfield=get(handles.popupmenu5,'String');
searchfield=searchfield{searchfieldnum};

searchtext=get(handles.text10,'String');
searchtext=[searchtext;{[searchfield,':.:',searchstr]}];

if ~isempty(searchstr')
    ezaz=true(size(xlsdata));
    for si=1:length(searchtext)
        searchtextline=searchtext{si};
        idx=strfind(searchtextline,':.:');
        searchfield=searchtextline(1:idx-1);
        searchstr=searchtextline(idx+3:end);
        for i=1:length(xlsdata)
            if ~any(regexp(xlsdata(i).(searchfield),searchstr))
                ezaz(i)=false;
            end
        end
        ezaz=find(ezaz);
        
        handles.data.actualIDs=[{'no cell selected'},{xlsdata(ezaz).ID}];
%         else
%             set(handles.popupmenu2,'String',[{'no cell selected'},{xlsdata.ID}]);
    end
    set(handles.text10,'String',searchtext);
end
set(handles.popupmenu2,'Value',1);
handles=updategui(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xlsdata=handles.data.xlsdata;
%  set(handles.popupmenu2,'String',[{'no cell selected'},{xlsdata.ID}]);
 set(handles.text10,'String','');
 handles.data.actualIDs=handles.data.IDs;
 handles=updategui(handles);
 guidata(hObject,handles);


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dirs=handles.data.dirs;
xlsdata=handles.data.xlsdata;
selectedsamplenum=get(handles.popupmenu1,'Value');
ID=handles.data.IDs{handles.data.samples(selectedsamplenum).loadedID};
icxlsnum=find(strcmp(ID,{xlsdata.ID}));
timeborders=[str2num(get(handles.edit4,'String')),str2num(get(handles.edit5,'String'))];
additionaldata=struct;
additionaldata.eventdata=handles.data.samples(selectedsamplenum).eventdata;
if isfield(handles.data.samples,'BrainStateData')
    additionaldata.BrainStateData=handles.data.samples(selectedsamplenum).BrainStateData;
end

dataout=aE_ic_field_analysis(dirs,xlsdata,icxlsnum,timeborders,'trough',additionaldata); %peak
valtozok_fieldplot.timebefore=.5;
valtozok_fieldplot.timestep=.05;
valtozok_fieldplot.timeafter=.5;
% valtozok_fieldplot.ICylim
% valtozok_fieldplot.Fieldylim
valtozok_fieldplot.highlightaAP=1;
valtozok_fieldplot.highlightaAP_timeback=.002;
valtozok_fieldplot.highlightaAP_timeforward=.02;

if isfield(additionaldata,'BrainStateData')
    statestodo=[unique({additionaldata.BrainStateData.name}),'All'];
else
    statestodo={'All'};
end
for statei=1:length(statestodo)
    statename=statestodo{statei};
    %%
    FieldData=dataout.FieldData;
    if ~strcmp(statename,'All')
        FieldData=FieldData(strcmp({FieldData.brainstatename},statename));
    end
    valtozok_fieldplot.figurenum=10+statei;
    persistent_generatefigures_2017_fieldanal(FieldData,valtozok_fieldplot)
    subplot(4,1,1);
    title(statename)
end





% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%

handles=updategui(handles);;
guidata(hObject,handles);
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


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
if get(hObject,'Value')==true
    handles=checkforaAPs(handles);
    guidata(hObject, handles)
end



% 
% if get(hObject,'Value')==true
%     samplenum=get(handles.popupmenu1,'Value');
%     button=2;
%     if isfield(handles.data.samples(samplenum).eventdata,'axonalAP')
%         button = questdlg('Do you want to rerun the aAP detection','aAP detection','Yes','No','Yes');
%     end
%     
%     if ~isfield(handles.data.samples(samplenum).eventdata,'axonalAP') | button==1
%         
%         timeborders=[str2num(get(handles.edit4,'String')),str2num(get(handles.edit5,'String'))];
%         %     for samplenum=1:length(handles.data.samples)
%         if handles.data.samples(samplenum).loadedID-1>0 & samplenum==get(handles.popupmenu1,'Value');
%             valtozok.timeborders=timeborders;
%             [dataout]=aE_sort_sAP_aAP(handles.data.dirs,handles.data.xlsdata,handles.data.samples(samplenum).loadedID-1,valtozok);
%             handles.data.samples(samplenum).eventdata=dataout.eventdata;
%         end
%         %     end
%         guidata(hObject, handles)
%     end
%     
% end

function handles=checkforaAPs(handles)
samplenum=get(handles.popupmenu1,'Value');
    button='No';
    if isfield(handles.data.samples(samplenum).eventdata,'axonalAP')
        button = questdlg('Do you want to rerun the aAP detection','aAP detection','Yes','No','Yes');
    end
    
    if ~isfield(handles.data.samples(samplenum).eventdata,'axonalAP') | strcmp(button,'Yes')
        
        timeborders=[str2num(get(handles.edit4,'String')),str2num(get(handles.edit5,'String'))];
        %     for samplenum=1:length(handles.data.samples)
        if handles.data.samples(samplenum).loadedID-1>0 & samplenum==get(handles.popupmenu1,'Value');
            valtozok.timeborders=timeborders;
            [dataout]=aE_sort_sAP_aAP(handles.data.dirs,handles.data.xlsdata,handles.data.samples(samplenum).loadedID-1,valtozok);
            %%
%             eventdata_GUI=eventdata;
            eventdata=dataout.eventdata;
            APwaves=dataout.APwaves;
            save([handles.data.dirs.eventdir,'sorted/',handles.data.xlsdata(handles.data.samples(samplenum).loadedID-1).ID],'eventdata','APwaves');
%             eventdata=eventdata_GUI;
            handles.data.samples(samplenum).eventdata=dataout.eventdata;
        end
        %     end
        
    end



% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(handles.d)
%%
selectedsamplenum=get(handles.popupmenu1,'Value');
set(handles.edit9,'String',handles.data.samples(selectedsamplenum).PSDdatatoplot.intensitypercentiles(99));
%%
plotandupdate(handles)
% assignin('base', 'data_InspectTraces', handles.data);


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotandupdate(handles)
% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes durinPush Buttong object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selectedsamplenum=get(handles.popupmenu1,'Value');

Statevalues=handles.data.BrainStateProps.Statevalues;

[selection,ok]=listdlg('ListString',Statevalues);
if ok==1
    rect=getrect(handles.axes1);
    %%
    starttime=rect(1);
    endtime=rect(1)+rect(3);
    
    BrainStateData=handles.data.samples(selectedsamplenum).BrainStateData;
    if isempty(fieldnames(BrainStateData))
        NEXT=1;
    else
        NEXT=length(BrainStateData)+1;
    end
    BrainStateData(NEXT).name=Statevalues{selection};
    BrainStateData(NEXT).starttime=starttime;
    BrainStateData(NEXT).endtime=endtime;
    % Ide olyan kell még, hogy ha ugyanabban az időben van más is,
    % akkor abból ki kell vágni ezt a részt.
    %%
    statestodel=zeros(size(BrainStateData));
    if NEXT>1
        for statei=1:NEXT-1
            if (BrainStateData(NEXT).starttime>BrainStateData(statei).starttime && BrainStateData(NEXT).starttime<BrainStateData(statei).endtime) || (BrainStateData(NEXT).endtime>BrainStateData(statei).starttime && BrainStateData(NEXT).endtime<BrainStateData(statei).endtime) || (BrainStateData(NEXT).endtime>BrainStateData(statei).endtime && BrainStateData(NEXT).starttime<BrainStateData(statei).starttime)
                %there is an intersect
                if strcmp(BrainStateData(NEXT).name,BrainStateData(statei).name);
                    % same state - merging them
                    BrainStateData(NEXT).starttime=min(BrainStateData(NEXT).starttime,BrainStateData(statei).starttime) ;
                    BrainStateData(NEXT).endtime=max(BrainStateData(NEXT).endtime,BrainStateData(statei).endtime) ;
                    statestodel(statei)=1;
                else
                    % different state
                    if BrainStateData(NEXT).starttime>BrainStateData(statei).starttime && BrainStateData(NEXT).starttime<BrainStateData(statei).endtime && BrainStateData(NEXT).starttime<BrainStateData(statei).endtime && BrainStateData(NEXT).endtime<BrainStateData(statei).endtime
                        % new state is completely within an old state
                        BrainStateData(NEXT+1)=BrainStateData(statei);
                        BrainStateData(NEXT+1).starttime=BrainStateData(NEXT).endtime;
                        BrainStateData(statei).endtime=BrainStateData(NEXT).starttime;
                    elseif BrainStateData(NEXT).starttime>BrainStateData(statei).starttime && BrainStateData(NEXT).starttime<BrainStateData(statei).endtime
                        % new state starts later
                        BrainStateData(statei).endtime=BrainStateData(NEXT).starttime;
                    elseif BrainStateData(NEXT).endtime>BrainStateData(statei).starttime && BrainStateData(NEXT).endtime<BrainStateData(statei).endtime
                        % new state starts earlier
                        BrainStateData(statei).starttime=BrainStateData(NEXT).endtime;
                    else
                        statestodel(statei)=1;
                        % old state is within the new state
                    end
                    if BrainStateData(statei).endtime <= BrainStateData(statei).starttime
                        statestodel(statei)=1;
                    end
                end
            end
        end
        if strcmp(Statevalues{selection},'Delete')
            statestodel(NEXT)=1;
        end
        BrainStateData(find(statestodel))=[];
    end
    %%
    
    
    %%
    
    handles.data.samples(selectedsamplenum).BrainStateData=BrainStateData;
end
plotandupdate(handles);

save([handles.data.dirs.brainstatedir,handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID}],'BrainStateData');
guidata(hObject,handles);


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selectedsamplenum=get(handles.popupmenu1,'Value');
axes(handles.axes1)
rect = getrect;
handles.axes1=gca;
%%
starttime=min([handles.data.samples(selectedsamplenum).starttime]);
endtime=max([handles.data.samples(selectedsamplenum).endtime]);
timedifi=handles.data.samples(selectedsamplenum).endtime-handles.data.samples(selectedsamplenum).starttime;
startval=(rect(1)-starttime)/timedifi;
endval=(rect(1)+rect(3)-starttime)/timedifi;
midval=mean([startval,endval]);

set(handles.slider1,'Value',startval);
set(handles.slider2,'Value',endval);
set(handles.slider3,'Value',midval);
%%

    
% sliderek beállítása
handles=updategui(handles);
handles=updatedatatoplot(handles);
plotandupdate(handles)
guidata(hObject, handles)


% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
valtozok=struct;
valtozok.movingvindowsize=5;
valtozok.movingvindowstep=.5; %seconds for downsampling and median filtering
valtozok.timeborders=[str2num(get(handles.edit4,'String')),str2num(get(handles.edit5,'String'))];
valtozok.frequencyrange=[1 4];
valtozok.PSDonfield=true;
valtozok.minsweeplength=valtozok.movingvindowsize/2;

selectedsamplenum=get(handles.popupmenu1,'Value');

aE_V0_vs_PSD(handles.data.dirs,handles.data.xlsdata,handles.data.samples(selectedsamplenum).selectedID-1,valtozok)