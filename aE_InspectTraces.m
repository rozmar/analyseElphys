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

% Last Modified by GUIDE v2.5 21-Jan-2019 14:24:10

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
set(handles.popupmenu9,'String',{'All','sporadic','persistent'})
set(handles.popupmenu10,'String',{'upper','lower','both'})

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

handles=updategui(handles);
guidata(hObject,handles)
plotandupdate(handles);
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
    set(handles.popupmenu4,'String',[{'stiumulus'},{handles.data.samples.marker},{'movement and pupil diameter'},{'PSD of field'},{'Behaviour'},{'Breathing'},{'PSD of Breathing'},{'Breathing field coherence'},{'PSD of IC'}]);
    set(handles.popupmenu7,'String',[{'stiumulus'},{handles.data.samples.marker},{'movement and pupil diameter'},{'PSD of field'},{'Behaviour'},{'Breathing'},{'PSD of Breathing'},{'Breathing field coherence'},{'PSD of IC'}]);
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
        disp('loading elphys')
        load([handles.data.dirs.bridgeddir,handles.data.IDs{handles.data.samples(samplei).selectedID}],'stimdata');
        handles.data.samples(samplei).stimdata=stimdata;
        load([handles.data.dirs.bridgeddir,handles.data.IDs{handles.data.samples(samplei).selectedID}],'bridgeddata');
        handles.data.samples(samplei).bridgeddata=bridgeddata;
        handles.data.samples(samplei).loadedID=handles.data.samples(samplei).selectedID;
        xlsnum=handles.data.samples(samplei).selectedID-1;
        if isfield(handles.data.xlsdata,'Video_time_offset')
            videotimeoffset=handles.data.xlsdata(xlsnum).Video_time_offset;
            %%
            if ischar(videotimeoffset)
                videotimeoffset=0;
                calculatetimeoffset=true;
            else
                calculatetimeoffset=false;
            end
        else
            calculatetimeoffset=false;
        end
        %%
        if isfield(handles.data.dirs,'videodir')
            fname=handles.data.xlsdata(handles.data.samples(samplei).selectedID-1).HEKAfname;

%            medfiltlength=5;
           
            
            %movementdata
            a=dir([handles.data.dirs.videodir,'movement/',fname,'.mat']);
            handles.data.samples(samplei).movementnames={};
            %%
            if ~isempty(a)
                disp('loading video')
                load([handles.data.dirs.videodir,'movement/',fname,'.mat']);
                if isfield(videodata,'time_behaviour')
                    exportbehaviour=true;
                else
                    exportbehaviour=false;
                end
                %%
                movement=[];
                time=[];
                videoo=[];
                videofname={};
                framenumber=[];
                videoidx=[];
                
                time_behaviour=[];
                tonetype=[];
                sound=[];
                lick=[];
                water=[];
                trigger=[];
                for i=1:length(videodata)
                    if isempty(videoo)
                        videoo=rawvideo(i).vid;
                    else
                        videoo=cat(3,videoo,rawvideo(i).vid);
                    end
                    framenumber=[framenumber;[1:length(videodata(i).movement_selected)]'];
                    videoidx=[videoidx;ones(length(videodata(i).movement_selected),1)*i];
                    movement=[movement;videodata(i).movement_selected];
                    time=[time;videodata(i).time];
                    [time,ix]=sort(time);
                    movement=movement(ix);
                    videoo=videoo(:,:,ix);
                    videofname{i}=[videodata(i).timestamp,'.avi'];
                    if exportbehaviour==true
                        time_behaviour=[time_behaviour;videodata(i).time_behaviour];
                        tonetype=[tonetype;videodata(i).tonetype];
                        sound=[sound;videodata(i).sound];
                        lick=[lick;videodata(i).lick];
                        water=[water;videodata(i).water];
                        
                        [time_behaviour,ix]=sort(time_behaviour);
                        tonetype=tonetype(ix);
                        sound=sound(ix);
                        lick=lick(ix);
                        water=water(ix);
                        if isfield(videodata,'trigger')
                            trigger=[trigger;videodata(i).trigger];
                            trigger=trigger(ix);
                        end
                        %% for some reason, the sound stimuli was not recorded for the second sound..
                        tonetypechangeidx=find(diff(tonetype)<0);
                        sound(tonetypechangeidx+1)=1;
                    end
                    
                end
                
                if calculatetimeoffset == true
                    %%
                    trigtimes=time_behaviour(find(diff(trigger)>0));
                    sweeptimes=[handles.data.samples(samplei).bridgeddata.realtime];
                    alldistances=[];
                    for trigi = 1:length(trigtimes)
                        distances=sweeptimes-trigtimes(trigi);
                        alldistances=[alldistances,distances];
                    end
                    alldistances=sort(alldistances);
                    [values,bins]=hist(alldistances,[-600:.1:600]);
                    [~,maxidx]=max(values(2:length(values)-1));
                    videotimeoffset=bins(maxidx+1);
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
                handles.data.samples(samplei).movementdata.time=time+videotimeoffset;
                
                if exportbehaviour==true
                    handles.data.samples(samplei).movementdata.time_behaviour=time_behaviour+videotimeoffset;
                    handles.data.samples(samplei).movementdata.tonetype=tonetype;
                    handles.data.samples(samplei).movementdata.sound=sound;
                    handles.data.samples(samplei).movementdata.lick=lick;
                    handles.data.samples(samplei).movementdata.water=water;
                    handles.data.samples(samplei).movementdata.trigger=trigger;
                end
                handles.data.samples(samplei).movementdata.video=videoo;
                handles.data.samples(samplei).movementdata.videofname=videofname;
                handles.data.samples(samplei).movementdata.videoidx=videoidx;
                handles.data.samples(samplei).movementdata.framenumber=framenumber;
                handles.data.samples(samplei).movementnames={'movement'};
                
                
                    
            end
            %%
            a=dir([handles.data.dirs.videodir,'ROIs/',fname,'.mat']);
            if ~isempty(a)
                disp('loading video analysis')
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
                disp('loading pupil diameters')
                load([handles.data.dirs.videodir,'eye/',fname]);
                pupildata.time=pupildata.time+videotimeoffset;
                handles.data.samples(samplei).pupildata=pupildata;
                diameter=handles.data.samples(samplei).pupildata.diameter;
                si=mode(diff(pupildata.time));
                if si==0
                    si=1/15;
                end
                diameter=moving(diameter,2);
                diameter=medfilt1(diameter,round(10/si));
                
                %                 a=dir([handles.data.dirs.videodir,'percentiles.mat']);
                %                 if ~isempty(a)
                %                     load([handles.data.dirs.videodir,'percentiles.mat'])
                %                     if exist('ROIdata','var') & any(strcmp({ROIdata.ROIname},'Eye'))
                %                         stats=regionprops(ROIdata(find(strcmp({ROIdata.ROIname},'Eye'))).mask,'MajorAxisLength');
                %                         minval=videopercentiles.pupilpercentiles(5)*stats.MajorAxisLength;
                %                         maxval=videopercentiles.pupilpercentiles(95)*stats.MajorAxisLength;
                %                     else
                %                         minval=videopercentiles.pupilpercentiles(5);
                %                         maxval=videopercentiles.pupilpercentiles(95);
                %                     end
                %                 else
                %                     pupildatasorted=sort(diameter);
                %                     minval=pupildatasorted(round(length(pupildatasorted)*.05));
                %                     maxval=pupildatasorted(round(length(pupildatasorted)*.95));
                %                 end
                %                 diameter=(handles.data.samples(samplei).pupildata.diameter-minval)/(maxval-minval);
                
                
                %% pontosabb módszer
                moviename=handles.data.samples(samplei).movementdata.videofname{1};
                setupname=handles.data.xlsdata(handles.data.samples(samplei).selectedID-1).setup;
                
                locations=marcicucca_locations;
                moviefiletoplay=[locations.tgtardir,'VIDEOdata/',setupname,'/',moviename];
                movieobj=VideoReader([moviefiletoplay]);
                szorzo=movieobj.Height/size(videodata(1).originalpic_all,1);
                if exist('ROIdata','var') & any(strcmp({ROIdata.ROIname},'Eye'))
                    stats=regionprops(ROIdata(find(strcmp({ROIdata.ROIname},'Eye'))).mask,'MajorAxisLength','MinorAxisLength');
                else
                    disp('No eye ROI!!!')
                    stats.MajorAxisLength=1;
                    stats.MinorAxisLength=1;
                end
                %%

                diameter=handles.data.samples(samplei).pupildata.diameter/szorzo/stats.MinorAxisLength*2;
                
                pupildiametersnow=sort(diameter);
                percentile95=pupildiametersnow(round(.95*length(pupildiametersnow)));
                if percentile95>1
                    diameter=diameter/percentile95;
                end
                
                
                
                diameter(diameter>1)=1;
                diameter(diameter<0)=0;
                
                handles.data.samples(samplei).pupildata.diameter=diameter;
            end

            end
            set(handles.popupmenu8,'Value',1);
            set(handles.popupmenu8,'String',handles.data.samples(samplei).movementnames);
        end
        
        
        if isfield(handles.data.dirs,'PSDdir')
            disp('loading PSD of IC')
            xlsnum=handles.data.samples(samplei).selectedID-1;
%             fieldxlsnum=find(strcmp(handles.data.xlsdata(xlsnum).HEKAfname,{handles.data.xlsdata.HEKAfname}) & [handles.data.xlsdata.field]==1);
            if ~isempty(xlsnum)
                if isfield(handles.data.dirs,'PSDdir_high') | isfield(handles.data.dirs,'PSDdir_log')
                    
                    answerstring=[',''Low PSD (0.5-4 Hz)'''];
                    
                    if isfield(handles.data.dirs,'PSDdir_high')
                        directory=handles.data.dirs.PSDdir_high;
                        a=dir([directory,handles.data.xlsdata(xlsnum).ID,'.mat']);
                        if ~isempty(a)
                            answerstring=[answerstring,',''High PSD (1-220 Hz)'''];
                        end
                    end
                    if isfield(handles.data.dirs,'PSDdir_log')
                        directory=handles.data.dirs.PSDdir_log;
                        a=dir([directory,handles.data.xlsdata(xlsnum).ID,'.mat']);
                        if ~isempty(a)
                            answerstring=[answerstring,',''LOG PSD (0.5-300 Hz)'''];
                        end
                    end
                    eval(['answer = questdlg(''Low or high frequency PSD should be loaded (IC)'',''Low or high PSD''',answerstring,',''Low PSD (0.5-4 Hz)'');']);
                    switch answer
                        case 'Low PSD (0.5-4 Hz)'
                            directory=handles.data.dirs.PSDdir;
                        case 'High PSD (1-220 Hz)'
                            directory=handles.data.dirs.PSDdir_high;
                        case 'LOG PSD (0.5-300 Hz)'
                            directory=handles.data.dirs.PSDdir_log;    
                    end
                else
                    directory=handles.data.dirs.PSDdir;
                end
                load([directory,handles.data.xlsdata(xlsnum).ID]);
                %%
                if isfield(PSDdata,'compress_offset')
                    for sweepi=1:length(PSDdata)
                        PSDdata(sweepi).powerMatrix=double(PSDdata(sweepi).powerMatrix)*PSDdata(sweepi).compress_multiplier+PSDdata(sweepi).compress_offset;
                    end
                end
                %%
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
                handles.data.samples(samplei).PSDdata_ic=PSDdata;
                set(handles.edit9,'String',max(PSDdata(1).powerMatrix(:)));
                set(handles.edit10,'String',8); %max(PSDdata(1).frequencyVector(:))
                set(handles.edit11,'String',min(PSDdata(1).frequencyVector(:)));
            else
                handles.data.samples(samplei).PSDdata_ic=[];
            end                
        end
        
        
        if isfield(handles.data.dirs,'PSDdir')
            disp('loading PSD of field')
            xlsnum=handles.data.samples(samplei).selectedID-1;
            fieldxlsnum=find(strcmp(handles.data.xlsdata(xlsnum).HEKAfname,{handles.data.xlsdata.HEKAfname}) & [handles.data.xlsdata.field]==1);
            if length(fieldxlsnum)>1
                disp('error, more than 1 field files found.. user should select..')
                ok=0;
                Selection=1;
                while ok==0 | length(Selection)~=1
                    [Selection,ok] = listdlg('PromptString','Select the corresponding field file:','ListString',{handles.data.xlsdata(fieldxlsnum).ID});
                end
                fieldxlsnum=fieldxlsnum(Selection);
            end
            if ~isempty(fieldxlsnum)
                if isfield(handles.data.dirs,'PSDdir_high') | isfield(handles.data.dirs,'PSDdir_log')
                    
                    answerstring=[',''Low PSD (0.5-4 Hz)'''];
                    if isfield(handles.data.dirs,'PSDdir_high')
                        directory=handles.data.dirs.PSDdir_high;
                        a=dir([directory,handles.data.xlsdata(fieldxlsnum).ID,'.mat']);
                        if ~isempty(a)
                            answerstring=[answerstring,',''High PSD (1-220 Hz)'''];
                        end
                    end
                    if isfield(handles.data.dirs,'PSDdir_log')
                        directory=handles.data.dirs.PSDdir_log;
                        a=dir([directory,handles.data.xlsdata(fieldxlsnum).ID,'.mat']);
                        if ~isempty(a)
                            answerstring=[answerstring,',''LOG PSD (0.5-300 Hz)'''];
                        end
                    end
                    eval(['answer = questdlg(''Low or high frequency PSD should be loaded (field)'',''Low or high PSD''',answerstring,',''Low PSD (0.5-4 Hz)'');']);
                    switch answer
                        case 'Low PSD (0.5-4 Hz)'
                            directory=handles.data.dirs.PSDdir;
                        case 'High PSD (1-220 Hz)'
                            directory=handles.data.dirs.PSDdir_high;
                        case 'LOG PSD (0.5-300 Hz)'
                            directory=handles.data.dirs.PSDdir_log;
                    end
                    %                     answer = questdlg('Low or high frequency PSD should be loaded','Low or high PSD','Low PSD (0.5-4 Hz)','High PSD (1-220 Hz)','Low PSD (0.5-4 Hz)');
                    %                     switch answer
                    %                         case 'Low PSD (0.5-4 Hz)'
                    %                             directory=handles.data.dirs.PSDdir;
                    %                         case 'High PSD (1-220 Hz)'
                    %                             directory=handles.data.dirs.PSDdir_high;
                    
                else
                    directory=handles.data.dirs.PSDdir;
                end
                load([directory,handles.data.xlsdata(fieldxlsnum).ID]);
                %%
                if isfield(PSDdata,'compress_offset')
                    for sweepi=1:length(PSDdata)
                        PSDdata(sweepi).powerMatrix=double(PSDdata(sweepi).powerMatrix)*PSDdata(sweepi).compress_multiplier+PSDdata(sweepi).compress_offset;
                    end
                end
                %%
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
                handles.data.samples(samplei).PSDdata_field=PSDdata;
                set(handles.edit9,'String',max(PSDdata(1).powerMatrix(:)));
                set(handles.edit10,'String',8); %max(PSDdata(1).frequencyVector(:))
                set(handles.edit11,'String',min(PSDdata(1).frequencyVector(:)));
            else
                handles.data.samples(samplei).PSDdata_field=[];
            end                
        end
        
        if isfield(handles.data.dirs,'brainstatedir')
            disp('loading brainstates')
            %%
            ID=handles.data.IDs{handles.data.samples(samplei).loadedID};
            a=dir([handles.data.dirs.brainstatedir,ID,'.mat']);
            if isempty(a)
                handles.data.samples(samplei).BrainStateData=struct;
            else
                temp=load([handles.data.dirs.brainstatedir,ID],'BrainStateData');
                handles.data.samples(samplei).BrainStateData=temp.BrainStateData;
                if isfield(temp,'BrainStateRules')
                    handles.data.samples(samplei).BrainStateRules=temp.BrainStateRules;
                else
                    handles.data.samples(samplei).BrainStateRules=struct;
                end
                
            end
        end
        %%
%         if isfield(handles.data.dirs,'breathingdir')
%             disp('loading breathing')
%             xlsnum=handles.data.samples(samplei).selectedID-1;
%             fieldxlsnum=find(strcmp(handles.data.xlsdata(xlsnum).HEKAfname,{handles.data.xlsdata.HEKAfname}) & [handles.data.xlsdata.field]==1);
%             if ~isempty(fieldxlsnum)
%                 a=dir([handles.data.dirs.breathingdir,handles.data.xlsdata(fieldxlsnum).ID,'.mat']);
%                 if ~isempty(a)
%                     load([handles.data.dirs.breathingdir,handles.data.xlsdata(fieldxlsnum).ID,'.mat']);
%                     %% filter and downsample breathing
%                     breathing=struct;
%                     breathing.y=[];
%                     breathing.time=[];
%                     for sweepi=1:length(rawdata)
%                         y=rawdata(sweepi).y;
%                         si=rawdata(sweepi).si;
%                         d{1} = designfilt('lowpassiir','PassbandFrequency',4,'StopbandFrequency',8,'SampleRate',1/si,'DesignMethod','butter');
%                         d{2} = designfilt('highpassiir','PassbandFrequency',1,'StopbandFrequency',.1,'SampleRate',1/si,'DesignMethod','butter');
% %                         [b,a]=butter(1,[.5,8]/(1/si)/2,'bandpass'); % data is filtered between 0.5 and 8 Hz
% %                         y=filtfilt(b,a,y);
%                         y=filtfilt(d{1},y);
%                         y=filtfilt(d{2},y);
%                         time=[si:si:si*length(y)]-si+rawdata(sweepi).realtime;
%                         downsamplenum=round(.001/si); %1kHz sampling rate is enough
%                         y=downsample(y,downsamplenum);
%                         time=downsample(time,downsamplenum);
%                         breathing.y=[breathing.y,y];
%                         breathing.time=[breathing.time,time];
% 
%                     end
%                     handles.data.samples(samplei).BreathingData=breathing;
%                 end
%                 
%             end
%         end
        
        if isfield(handles.data.dirs,'breathingdir_psd')
            disp('loading breathing PSD')
            xlsnum=handles.data.samples(samplei).selectedID-1;
            fieldxlsnum=find(strcmp(handles.data.xlsdata(xlsnum).HEKAfname,{handles.data.xlsdata.HEKAfname}) & [handles.data.xlsdata.field]==1);
            if ~isempty(fieldxlsnum)
                a=dir([handles.data.dirs.breathingdir,handles.data.xlsdata(fieldxlsnum).ID,'.mat']);
                if ~isempty(a)
                    load([handles.data.dirs.breathingdir_psd,handles.data.xlsdata(fieldxlsnum).ID,'.mat']);
                    %%
                    if isfield(PSDdata,'compress_offset')
                        for sweepi=1:length(PSDdata)
                            PSDdata(sweepi).powerMatrix=double(PSDdata(sweepi).powerMatrix)*PSDdata(sweepi).compress_multiplier+PSDdata(sweepi).compress_offset;
                        end
                    end
                    handles.data.samples(samplei).BreathingData_PSD=PSDdata;
                end
                
            end
        end
%         if isfield(handles.data.dirs,'cross_spectrum_breathing')
%             disp('loading cross spectrum breathing')
%             xlsnum=handles.data.samples(samplei).selectedID-1;
%             fieldxlsnum=find(strcmp(handles.data.xlsdata(xlsnum).HEKAfname,{handles.data.xlsdata.HEKAfname}) & [handles.data.xlsdata.field]==1);
%             if ~isempty(fieldxlsnum)
%                 a=dir([handles.data.dirs.breathingdir,handles.data.xlsdata(fieldxlsnum).ID,'.mat']);
%                 if ~isempty(a)
%                     load([handles.data.dirs.cross_spectrum_breathing,handles.data.xlsdata(fieldxlsnum).ID,'.mat']);
%                     %%
%                     if isfield(CROSSdata,'compress_offset')
%                         for sweepi=1:length(CROSSdata)
%                             CROSSdata(sweepi).coherence=double(CROSSdata(sweepi).coherence)*CROSSdata(sweepi).compress_multiplier+CROSSdata(sweepi).compress_offset;
%                         end
%                     end
%                     handles.data.samples(samplei).Breathing_field_CROSS=CROSSdata;
%                 end
%                 
%             end
%         end
    end
    %%
    if exist('xlsnum','var') & isfield(handles.data.xlsdata,'Cranio_Lat')
    cellcoord_lat=[handles.data.xlsdata(xlsnum).Cranio_Lat]+(([handles.data.xlsdata(xlsnum).locationX]-[handles.data.xlsdata(xlsnum).Cranio_center_X])/1000).*[handles.data.xlsdata(xlsnum).Lateral_dir_X];
cellcoord_AP=[handles.data.xlsdata(xlsnum).Cranio_AP]+(([handles.data.xlsdata(xlsnum).locationy]-[handles.data.xlsdata(xlsnum).Cranio_center_Y])/1000).*[handles.data.xlsdata(xlsnum).Rostral_dir_Y];
cellcoord_z=handles.data.xlsdata(xlsnum).locationz;
    fieldcoord_lat=[handles.data.xlsdata(fieldxlsnum).Cranio_Lat]+(([handles.data.xlsdata(fieldxlsnum).locationX]-[handles.data.xlsdata(fieldxlsnum).Cranio_center_X])/1000).*[handles.data.xlsdata(fieldxlsnum).Lateral_dir_X];
fieldcoord_AP=[handles.data.xlsdata(fieldxlsnum).Cranio_AP]+(([handles.data.xlsdata(fieldxlsnum).locationy]-[handles.data.xlsdata(fieldxlsnum).Cranio_center_Y])/1000).*[handles.data.xlsdata(fieldxlsnum).Rostral_dir_Y];
fieldcoord_z=handles.data.xlsdata(fieldxlsnum).locationz;
craniocoord=[handles.data.xlsdata(xlsnum).Cranio_AP,handles.data.xlsdata(xlsnum).Cranio_Lat];
cellcoord=[cellcoord_AP,cellcoord_lat,cellcoord_z];
fieldcoord=[fieldcoord_AP,fieldcoord_lat,fieldcoord_z];
    disp(['cranio coord: AP - ',num2str(craniocoord(1)),'   Lateral - ',num2str(craniocoord(2))])
    disp(['cell coord: AP - ',num2str(cellcoord(1)),'   Lateral - ',num2str(cellcoord(2)),'   depth - ',num2str(cellcoord(3))])
    disp(['field coord: AP - ',num2str(fieldcoord(1)),'   Lateral - ',num2str(fieldcoord(2)),'   depth - ',num2str(fieldcoord(3))])
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


axes(handles.axes2)
cla 
hold on
handles.axes2=gca;

axes(handles.axes3)
cla 
hold on
handles.axes3=gca;

axes(handles.axes1)
cla 
hold on
handles.axes1=gca;


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
            neededevents=[handles.data.samples(samplenum).eventdata.axonalAP]==1;
            eventsnow=handles.data.samples(samplenum).eventdata(neededevents);
            neededevents=false(size(eventsnow));
              for sweepi=1:length(neededwaves)
                 sweepnum=neededwaves(sweepi);
                 neededevents([eventsnow.sweepnum]==sweepnum)=1;
              end
             plot([eventsnow(neededevents).maxtime],[eventsnow(neededevents).maxval],'ro','MarkerSize',8)
             
             handles.data.samples(samplenum).eventdata=persistent_sort_sporadic_persistent_aAPs(handles.data.samples(samplenum).eventdata);
             
             neededevents=[handles.data.samples(samplenum).eventdata.axonalAP_persistent]==1;
            eventsnow=handles.data.samples(samplenum).eventdata(neededevents);
            neededevents=false(size(eventsnow));
              for sweepi=1:length(neededwaves)
                 sweepnum=neededwaves(sweepi);
                 neededevents([eventsnow.sweepnum]==sweepnum)=1;
              end
             plot([eventsnow(neededevents).maxtime],[eventsnow(neededevents).maxval],'rx','MarkerSize',8)
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
hObject=findall(gcf,'Name','aE_InspectTraces');
guidata(hObject,handles);

function plotlowerpanels(marker,samplenum,starttime,endtime,neededwaves,axesname,popupname,handles)
axes(handles.(axesname))
cla
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
    if isfield(handles.data.samples(samplenum),'movementdata') & ~isempty(get(handles.popupmenu8,'String'))
        possiblefields=get(handles.popupmenu8,'String');
        movementfield=possiblefields{get(handles.popupmenu8,'Value')};
        ettol=find(handles.data.samples(samplenum).movementdata.time>starttime,1,'first');
        eddig=find(handles.data.samples(samplenum).movementdata.time<endtime,1,'last');
        plot(handles.data.samples(samplenum).movementdata.time(ettol:eddig),handles.data.samples(samplenum).movementdata.(movementfield)(ettol:eddig),'r-','LineWidth',2);
        
    end
    ylabel('\color{red}Movement \color{black}and \color{blue}pupil size \color{black}(AU)')
    ylim([0 2])
elseif get(handles.(popupname),'Value')==length(handles.data.samples)+3 & ~isempty(handles.data.samples(samplenum).PSDdata_fieldtoplot)
    
    imagesc(handles.data.samples(samplenum).PSDdata_fieldtoplot.time,handles.data.samples(samplenum).PSDdata_fieldtoplot.frequencyVector,handles.data.samples(samplenum).PSDdata_fieldtoplot.powerMatrix)
%     contourf(handles.data.samples(samplenum).PSDdata_fieldtoplot.time,handles.data.samples(samplenum).PSDdata_fieldtoplot.frequencyVector,handles.data.samples(samplenum).PSDdata_fieldtoplot.powerMatrix)
    
    if strcmp(axesname,'axes2')
        cmaxval=handles.data.samples(samplenum).plotdetails.PSD.cmax(1);
    elseif strcmp(axesname,'axes3')
        cmaxval=handles.data.samples(samplenum).plotdetails.PSD.cmax(2);
    end
%     caxis([0 str2num(get(handles.edit9,'String'))])
    caxis([0 cmaxval])
    
    set(gca,'YDir','normal');
    colormap linspecer
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    %%
    set(gca, 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
    ytickidxs=[];
    yticknow=get(gca,'YTick');
    ylimnow=get(gca,'Ylim');
   
    linearyticks = linspace(ylimnow(1),ylimnow(end),length(handles.data.samples(samplenum).PSDdata_fieldtoplot.frequencyVector));
    for ticki=1:length(yticknow)
        [~,ytickidxs(ticki)]=min(abs(linearyticks-yticknow(ticki)));
    end
    yticklabelnow=handles.data.samples(samplenum).PSDdata_fieldtoplot.frequencyVector(ytickidxs);
    set(gca,'Yticklabel',round(yticklabelnow*100)/100)
    %%
%     set(gca,'YTick',round(handles.data.samples(samplenum).PSDdata_fieldtoplot.frequencyVector))
    hold on
    if ~isempty(handles.data.samples(samplenum).PSDdata_fieldtoplot.frequencyVector)
        szorzo=max(handles.data.samples(samplenum).PSDdata_fieldtoplot.frequencyVector);
    else
        szorzo=1;
    end
    for sweepi=1:length(handles.data.samples(samplenum).PSDdata_fieldtoplot.trace)
        ettol=find(handles.data.samples(samplenum).PSDdata_fieldtoplot.trace(sweepi).x>starttime,1,'first');
        eddig=find(handles.data.samples(samplenum).PSDdata_fieldtoplot.trace(sweepi).x<endtime,1,'last');
        plot(handles.data.samples(samplenum).PSDdata_fieldtoplot.trace(sweepi).x(ettol:eddig),handles.data.samples(samplenum).PSDdata_fieldtoplot.trace(sweepi).y(ettol:eddig)*szorzo/3+2*szorzo/3,'Color',[0.5 .5 .5],'LineWidth',.5)
        
    end
    %     disp('muahah')
    axis tight
elseif get(handles.(popupname),'Value')==length(handles.data.samples)+4 & isfield(handles.data.samples(samplenum).movementdata,'time_behaviour')
    %%
    hold on
    ettol=find(handles.data.samples(samplenum).movementdata.time_behaviour>starttime,1,'first');
    eddig=find(handles.data.samples(samplenum).movementdata.time_behaviour<endtime,1,'last');
    
    plot(handles.data.samples(samplenum).movementdata.time_behaviour(ettol:eddig),(handles.data.samples(samplenum).movementdata.tonetype(ettol:eddig)+1).*handles.data.samples(samplenum).movementdata.sound(ettol:eddig),'k-','LineWidth',1)
%     plot(handles.data.samples(samplenum).movementdata.time_behaviour(ettol:eddig-1),diff(handles.data.samples(samplenum).movementdata.tonetype(ettol:eddig)),'c-')
    plot(handles.data.samples(samplenum).movementdata.time_behaviour(ettol:eddig),handles.data.samples(samplenum).movementdata.lick(ettol:eddig),'r-','LineWidth',1)
    plot(handles.data.samples(samplenum).movementdata.time_behaviour(ettol:eddig),handles.data.samples(samplenum).movementdata.water(ettol:eddig),'b-','LineWidth',1)
%     plot(handles.data.samples(samplenum).movementdata.time(ettol:eddig),handles.data.samples(samplenum).movementdata.(movementfield)(ettol:eddig),'r-','LineWidth',2);
elseif get(handles.(popupname),'Value')==length(handles.data.samples)+5 & isfield(handles.data.samples(samplenum),'BreathingData')
    %%
    ettol=find(handles.data.samples(samplenum).BreathingData.time>=starttime,1,'first');
    if isempty(ettol)
        ettol=1;
    end
    eddig=find(handles.data.samples(samplenum).BreathingData.time<=endtime,1,'last');
    if isempty(eddig)
        eddig=length(handles.data.samples(samplenum).BreathingData.time);
    end
%     x=handles.data.samples(samplenum).BreathingData.time(ettol:eddig-1);
%     y=diff(handles.data.samples(samplenum).BreathingData.y(ettol:eddig));
    x=handles.data.samples(samplenum).BreathingData.time(ettol:eddig);
    y=handles.data.samples(samplenum).BreathingData.y(ettol:eddig);
    plot(x,y,'k-','LineWidth',1)
    ylim([min(y), max(y)])
elseif get(handles.(popupname),'Value')==length(handles.data.samples)+6 & isfield(handles.data.samples(samplenum),'BreathingPSDdatatoplot')
    imagesc(handles.data.samples(samplenum).BreathingPSDdatatoplot.time,handles.data.samples(samplenum).BreathingPSDdatatoplot.frequencyVector,handles.data.samples(samplenum).BreathingPSDdatatoplot.powerMatrix)
    
    
%     caxis([0 str2num(get(handles.edit9,'String'))])
    
    if strcmp(axesname,'axes2')
        cmaxval=handles.data.samples(samplenum).plotdetails.PSD.cmax(1);
    elseif strcmp(axesname,'axes3')
        cmaxval=handles.data.samples(samplenum).plotdetails.PSD.cmax(2);
    end
%     caxis([0 str2num(get(handles.edit9,'String'))])
    caxis([0 cmaxval])
    
    set(gca,'YDir','normal');
    colormap linspecer
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    hold on
    if ~isempty(handles.data.samples(samplenum).BreathingPSDdatatoplot.frequencyVector)
        szorzo=max(handles.data.samples(samplenum).BreathingPSDdatatoplot.frequencyVector);
    else
        szorzo=1;
    end
    for sweepi=1:length(handles.data.samples(samplenum).BreathingPSDdatatoplot.trace)
        ettol=find(handles.data.samples(samplenum).BreathingPSDdatatoplot.trace(sweepi).x>starttime,1,'first');
        eddig=find(handles.data.samples(samplenum).BreathingPSDdatatoplot.trace(sweepi).x<endtime,1,'last');
        plot(handles.data.samples(samplenum).BreathingPSDdatatoplot.trace(sweepi).x(ettol:eddig),handles.data.samples(samplenum).BreathingPSDdatatoplot.trace(sweepi).y(ettol:eddig)*szorzo/3+2*szorzo/3,'Color',[0.5 .5 .5],'LineWidth',.5)
        
    end
    %     disp('muahah')
    axis tight
elseif get(handles.(popupname),'Value')==length(handles.data.samples)+7 & isfield(handles.data.samples(samplenum),'Breathing_Field_CROSSdatatoplot')
    imagesc(handles.data.samples(samplenum).Breathing_Field_CROSSdatatoplot.time,handles.data.samples(samplenum).Breathing_Field_CROSSdatatoplot.frequencyVector,handles.data.samples(samplenum).Breathing_Field_CROSSdatatoplot.coherence)
    
    if strcmp(axesname,'axes2')
        cmaxval=handles.data.samples(samplenum).plotdetails.PSD.cmax(1);
    elseif strcmp(axesname,'axes3')
        cmaxval=handles.data.samples(samplenum).plotdetails.PSD.cmax(2);
    end
%     caxis([0 str2num(get(handles.edit9,'String'))])
    caxis([0 cmaxval])
    
%     caxis([0 str2num(get(handles.edit9,'String'))])
    set(gca,'YDir','normal');
    colormap linspecer
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    hold on
%     if ~isempty(handles.data.samples(samplenum).Breathing_Field_CROSSdatatoplot.frequencyVector)
%         szorzo=max(handles.data.samples(samplenum).Breathing_Field_CROSSdatatoplot.frequencyVector);
%     else
%         szorzo=1;
%     end
%     for sweepi=1:length(handles.data.samples(samplenum).BreathingPSDdatatoplot.trace)
%         ettol=find(handles.data.samples(samplenum).BreathingPSDdatatoplot.trace(sweepi).x>starttime,1,'first');
%         eddig=find(handles.data.samples(samplenum).BreathingPSDdatatoplot.trace(sweepi).x<endtime,1,'last');
%         plot(handles.data.samples(samplenum).BreathingPSDdatatoplot.trace(sweepi).x(ettol:eddig),handles.data.samples(samplenum).BreathingPSDdatatoplot.trace(sweepi).y(ettol:eddig)*szorzo/3+2*szorzo/3,'Color',[0.5 .5 .5],'LineWidth',.5)
%         
%     end
    %     disp('muahah')
    axis tight
  
elseif get(handles.(popupname),'Value')==length(handles.data.samples)+8 & ~isempty(handles.data.samples(samplenum).PSDdata_ictoplot)
    
    imagesc(handles.data.samples(samplenum).PSDdata_ictoplot.time,handles.data.samples(samplenum).PSDdata_ictoplot.frequencyVector,handles.data.samples(samplenum).PSDdata_ictoplot.powerMatrix)
    
    if strcmp(axesname,'axes2')
        cmaxval=handles.data.samples(samplenum).plotdetails.PSD.cmax(1);
    elseif strcmp(axesname,'axes3')
        cmaxval=handles.data.samples(samplenum).plotdetails.PSD.cmax(2);
    end
%     caxis([0 str2num(get(handles.edit9,'String'))])
    caxis([0 cmaxval])
%     caxis([0 str2num(get(handles.edit9,'String'))])
    set(gca,'YDir','normal');
    colormap linspecer
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    hold on
    if ~isempty(handles.data.samples(samplenum).PSDdata_ictoplot.frequencyVector)
        szorzo=max(handles.data.samples(samplenum).PSDdata_ictoplot.frequencyVector);
    else
        szorzo=1;
    end
    for sweepi=1:length(handles.data.samples(samplenum).PSDdata_ictoplot.trace)
        ettol=find(handles.data.samples(samplenum).PSDdata_ictoplot.trace(sweepi).x>starttime,1,'first');
        eddig=find(handles.data.samples(samplenum).PSDdata_ictoplot.trace(sweepi).x<endtime,1,'last');
        plot(handles.data.samples(samplenum).PSDdata_ictoplot.trace(sweepi).x(ettol:eddig),handles.data.samples(samplenum).PSDdata_ictoplot.trace(sweepi).y(ettol:eddig)*szorzo/3+2*szorzo/3,'Color',[0.5 .5 .5],'LineWidth',.5)
        
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
    if isfield(handles.data.xlsdata,'axonalAPnum')
        fieldstoadd=[fieldstoadd,{'axonalAPnum'}];
    end
    IDstring={};
    actualIDs=handles.data.actualIDs;
    for i=1:length(actualIDs)
        xlsidx=find(strcmp({handles.data.xlsdata.ID},actualIDs{i}));
        IDstring{i}=actualIDs{i};
        for fieldi=1:length(fieldstoadd);
            if i>1 & isfield(handles.data.xlsdata,fieldstoadd{fieldi})
                if strcmp(fieldstoadd{fieldi},'field') & isnumeric(handles.data.xlsdata(xlsidx).(fieldstoadd{fieldi}))
                    if handles.data.xlsdata(xlsidx).(fieldstoadd{fieldi})==1
                        IDstring{i}=[IDstring{i} ,' - ', fieldstoadd{fieldi}];
%                     elseif handles.data.xlsdata(xlsidx).(fieldstoadd{fieldi})==0
%                     else
                        
                    end
                elseif isnumeric(handles.data.xlsdata(xlsidx).(fieldstoadd{fieldi}))
                    IDstring{i}=[IDstring{i} ,' - ', num2str(handles.data.xlsdata(xlsidx).(fieldstoadd{fieldi}))];
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
        PSD_cmax=str2num(get(handles.edit9,'string'));
        panelnames=get(handles.popupmenu10,'string');
        selectedpanel=panelnames{get(handles.popupmenu10,'Value')};
        plotdetails=handles.data.samples(selectedsamplenum).plotdetails;
        if ~isempty(fieldnames(plotdetails))
        if strcmp(selectedpanel,'upper')
            mehetPSD=~(plotdetails.PSD.cmax(1)==PSD_cmax);
        elseif strcmp(selectedpanel,'lower')
            mehetPSD=~(plotdetails.PSD.cmax(2)==PSD_cmax);
        elseif strcmp(selectedpanel,'both')
            mehetPSD=~(plotdetails.PSD.cmax(2)==PSD_cmax & plotdetails.PSD.cmax(1)==PSD_cmax);    
        end
        end
        if ~isempty(neededwaves)
            if isempty(fieldnames(plotdetails)) | ~(length(neededwaves)==length(plotdetails.neededwaves)) | ~any(neededwaves==plotdetails.neededwaves) | ~(neededsamplenum==plotdetails.neededsamplenum) | ~all(cutoffreq==plotdetails.cutoffreq) | ~(filterdeg==plotdetails.filterdeg) | ~(shoulddownsample==plotdetails.shoulddownsample) | ~(minfreq==plotdetails.minfreq) | ~(maxfreq==plotdetails.maxfreq)  | mehetPSD
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
                if ~isfield(plotdetails,'PSD') | strcmp(selectedpanel,'both')
                    plotdetails.PSD.cmax(1)=PSD_cmax;
                    plotdetails.PSD.cmax(2)=PSD_cmax;
                elseif strcmp(selectedpanel,'upper')
                    plotdetails.PSD.cmax(1)=PSD_cmax;
                elseif strcmp(selectedpanel,'lower')
                    plotdetails.PSD.cmax(2)=PSD_cmax;
                end
                
                
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
                    if isfield(handles.data.samples(selectedsamplenum).bridgeddata,'offset') & get(handles.checkbox9,'Value')==true
                        datatoplot(sweepi).yvoltage=y-handles.data.samples(selectedsamplenum).bridgeddata(sweepnum).offset;
                    else
                        datatoplot(sweepi).yvoltage=y;
                    end
                    
                    
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
                
                 if ~isempty(handles.data.samples(selectedsamplenum).PSDdata_ic) & (get(handles.popupmenu4,'Value')==length(handles.data.samples)+8 | get(handles.popupmenu7,'Value')==length(handles.data.samples)+8)
                    PSDdatatoplot=preparePSDdataforplotting(handles.data.samples(selectedsamplenum).PSDdata_ic,handles,neededsamplenum,shoulddownsample);
                    handles.data.samples(selectedsamplenum).PSDdata_ictoplot=PSDdatatoplot;
                    %                 set(handles.edit9,'String',handles.data.samples(selectedsamplenum).PSDdata_ictoplot.intensitypercentiles(99));
                else
                    handles.data.samples(selectedsamplenum).PSDdata_ictoplot=[];
                end
                
                if ~isempty(handles.data.samples(selectedsamplenum).PSDdata_field) & (get(handles.popupmenu4,'Value')==length(handles.data.samples)+3 | get(handles.popupmenu7,'Value')==length(handles.data.samples)+3)
                    PSDdatatoplot=preparePSDdataforplotting(handles.data.samples(selectedsamplenum).PSDdata_field,handles,neededsamplenum,shoulddownsample);
                    handles.data.samples(selectedsamplenum).PSDdata_fieldtoplot=PSDdatatoplot;
                    %                 set(handles.edit9,'String',handles.data.samples(selectedsamplenum).PSDdata_fieldtoplot.intensitypercentiles(99));
                else
                    handles.data.samples(selectedsamplenum).PSDdata_fieldtoplot=[];
                end
                
                
                if isfield(handles.data.samples(selectedsamplenum),'BreathingData_PSD') & ~isempty(handles.data.samples(selectedsamplenum).BreathingData_PSD) & (get(handles.popupmenu4,'Value')==length(handles.data.samples)+6 | get(handles.popupmenu7,'Value')==length(handles.data.samples)+6)
                    PSDdatatoplot=preparePSDdataforplotting(handles.data.samples(selectedsamplenum).BreathingData_PSD,handles,neededsamplenum,shoulddownsample);
                    handles.data.samples(selectedsamplenum).BreathingPSDdatatoplot=PSDdatatoplot;
                    %                 set(handles.edit9,'String',handles.data.samples(selectedsamplenum).PSDdata_fieldtoplot.intensitypercentiles(99));
                else
                    handles.data.samples(selectedsamplenum).BreathingPSDdatatoplot=[];
                end
                
                if isfield(handles.data.samples(selectedsamplenum),'Breathing_field_CROSS') & ~isempty(handles.data.samples(selectedsamplenum).Breathing_field_CROSS) & (get(handles.popupmenu4,'Value')==length(handles.data.samples)+7 | get(handles.popupmenu7,'Value')==length(handles.data.samples)+7)
                    CROSSdatatoplot=prepareCROSSdataforplotting(handles.data.samples(selectedsamplenum).Breathing_field_CROSS,handles,neededsamplenum,shoulddownsample);
                    handles.data.samples(selectedsamplenum).Breathing_Field_CROSSdatatoplot=CROSSdatatoplot;
                    %                 set(handles.edit9,'String',handles.data.samples(selectedsamplenum).PSDdata_fieldtoplot.intensitypercentiles(99));
                else
                    handles.data.samples(selectedsamplenum).Breathing_Field_CROSSdatatoplot=[];
                end
                
            end
        end
        
    end
end

function PSDdatatoplot=preparePSDdataforplotting(PSDdata,handles,neededsamplenum,shoulddownsample)
hosszofPSDwaves=zeros(1,length(PSDdata));
for sweepi=1:length(PSDdata)
    if ~isempty(PSDdata(sweepi).y)
        hosszofPSDwaves(1,sweepi)=length(PSDdata(sweepi).y)*PSDdata(sweepi).si_powerMatrix;
    else
        hosszofPSDwaves(1,sweepi)=0;
    end
end
bigmatrix=[];
timesofPSDwaves=[PSDdata.realtime];
reallyneededwaves=timesofPSDwaves>=str2num(get(handles.edit4,'String')) & timesofPSDwaves<=str2num(get(handles.edit5,'String')) & hosszofPSDwaves > 3;
if reallyneededwaves(1)==0
    reallyneededwaves(find(timesofPSDwaves<=str2num(get(handles.edit4,'String')),1,'last'))=true;
end
neededwaves=find(reallyneededwaves);
minfreq=str2num(get(handles.edit11,'string'));
maxfreq=str2num(get(handles.edit10,'string'));

for sweepi=1:length(neededwaves)
    sweepnum=neededwaves(sweepi);
    if isempty(bigmatrix)
        frequencyVector=PSDdata(sweepnum).frequencyVector;
        neededfrequencies=frequencyVector>=minfreq & frequencyVector<=maxfreq;
        frequencyVector=frequencyVector(neededfrequencies);
        starttime=PSDdata(sweepnum).realtime;
        si=PSDdata(sweepnum).si_powerMatrix;
        bigmatrix=PSDdata(sweepnum).powerMatrix(neededfrequencies,:);
        
    else
        timeskipped=(PSDdata(sweepnum).realtime-(starttime+size(bigmatrix,2)*si));
        colstoadd=round(timeskipped/si);
        if ~isempty(PSDdata(sweepnum).powerMatrix)
            bigmatrix=[bigmatrix,zeros(size(bigmatrix,1),colstoadd),PSDdata(sweepnum).powerMatrix(neededfrequencies,:)];
        end
        %             starttime+size(bigmatrix,2)*si
        %             return
    end
    
    PSDdatatoplot.trace(sweepi).x=[0:si:si*(length(PSDdata(sweepnum).y)-1)]+PSDdata(sweepnum).realtime;
    
    if ~isempty(PSDdata(sweepnum).y)
        y=PSDdata(sweepnum).y;
%         dofilt=0;
%         if dofilt==1
%             [b,a]=butter(1,.1/(1/si)/2,'high');
%             ytofilt=[y(end:-1:1),y,y(end:-1:1)];
%             ytofilt=filtfilt(b,a,ytofilt);
%             y=ytofilt([1:length(y)]+length(y));
%         end
        y=(y-min(y));
        ysort=sort(y);
        y=y/ysort(round(length(y)*.99));
        %                         y=y/max(y);
        y(y>1)=1;
        PSDdatatoplot.trace(sweepi).y=y;
        PSDdatatoplot.trace(sweepi).realtime=PSDdata(sweepnum).realtime;
        PSDdatatoplot.trace(sweepi).si=PSDdata(sweepnum).si_powerMatrix;
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

%% normalize 
% bigmatrix=bigmatrix.*repmat(frequencyVector',1,size(bigmatrix,2));

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


function CROSSdatatoplot=prepareCROSSdataforplotting(CROSSdata,handles,neededsamplenum,shoulddownsample)
hosszofPSDwaves=zeros(1,length(CROSSdata));
for sweepi=1:length(CROSSdata)
    if ~isempty(CROSSdata(sweepi).y1)
        hosszofPSDwaves(1,sweepi)=length(CROSSdata(sweepi).y1)*CROSSdata(sweepi).si_cross;
    else
        hosszofPSDwaves(1,sweepi)=0;
    end
end
bigmatrix=[];
timesofPSDwaves=[CROSSdata.realtime];
reallyneededwaves=timesofPSDwaves>=str2num(get(handles.edit4,'String')) & timesofPSDwaves<=str2num(get(handles.edit5,'String')) & hosszofPSDwaves > 3;
if reallyneededwaves(1)==0
    reallyneededwaves(find(timesofPSDwaves<=str2num(get(handles.edit4,'String')),1,'last'))=true;
end
neededwaves=find(reallyneededwaves);
minfreq=str2num(get(handles.edit11,'string'));
maxfreq=str2num(get(handles.edit10,'string'));

for sweepi=1:length(neededwaves)
    sweepnum=neededwaves(sweepi);
    if isempty(bigmatrix)
        frequencyVector=CROSSdata(sweepnum).frequencyVector;
        neededfrequencies=frequencyVector>=minfreq & frequencyVector<=maxfreq;
        frequencyVector=frequencyVector(neededfrequencies);
        starttime=CROSSdata(sweepnum).realtime;
        si=CROSSdata(sweepnum).si_cross;
        bigmatrix=CROSSdata(sweepnum).coherence(neededfrequencies,:);
        
    else
        timeskipped=(CROSSdata(sweepnum).realtime-(starttime+size(bigmatrix,2)*si));
        colstoadd=round(timeskipped/si);
        if ~isempty(CROSSdata(sweepnum).coherence)
            bigmatrix=[bigmatrix,zeros(size(bigmatrix,1),colstoadd),CROSSdata(sweepnum).coherence(neededfrequencies,:)];
        end
        %             starttime+size(bigmatrix,2)*si
        %             return
    end
    
    CROSSdatatoplot.trace(sweepi).x=[0:si:si*(length(CROSSdata(sweepnum).y1)-1)]+CROSSdata(sweepnum).realtime;
    
    if ~isempty(CROSSdata(sweepnum).y1)
        y=CROSSdata(sweepnum).y1;
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
        CROSSdatatoplot.trace(sweepi).y1=y;
        CROSSdatatoplot.trace(sweepi).realtime=CROSSdata(sweepnum).realtime;
        CROSSdatatoplot.trace(sweepi).si=CROSSdata(sweepnum).si_cross;
    end
end



%%
if ~isempty(fieldnames(CROSSdatatoplot.trace))
    samplenumnow=length([CROSSdatatoplot.trace.y1]);
    if neededsamplenum<samplenumnow & neededsamplenum>0 & shoulddownsample==1
        %                 disp('downsampling started')
        ratio=round(samplenumnow/neededsamplenum);
        for sweepi=1:length(neededwaves)
            CROSSdatatoplot.trace(sweepi).x=downsample(CROSSdatatoplot.trace(sweepi).x,ratio);
            CROSSdatatoplot.trace(sweepi).y1=downsample(CROSSdatatoplot.trace(sweepi).y1,ratio);
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

%% normalize 
% bigmatrix=bigmatrix.*repmat(frequencyVector',1,size(bigmatrix,2));

%%
CROSSdatatoplot.coherence=bigmatrix;
CROSSdatatoplot.frequencyVector=frequencyVector;
CROSSdatatoplot.si_cross=si;
CROSSdatatoplot.time=bigtime;
if ~isempty(CROSSdatatoplot.coherence)
    values=sort(CROSSdatatoplot.coherence(:));
    CROSSdatatoplot.intensitypercentiles=values(ceil([.01:.01:1]*length(values)));
else
    CROSSdatatoplot.intensitypercentiles=ones(size([.01:.01:1]));
end

% disp('lol')


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
guidata(hObject, handles)
plotandupdate(handles)


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
guidata(hObject, handles)
plotandupdate(handles)


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
guidata(hObject, handles)
plotandupdate(handles)


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
guidata(hObject, handles)
plotandupdate(handles)


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
        
%         else
%             set(handles.popupmenu2,'String',[{'no cell selected'},{xlsdata.ID}]);
    end
    ezaz=find(ezaz);
        
    handles.data.actualIDs=[{'no cell selected'},{xlsdata(ezaz).ID}];
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

data=struct;
data.eventdata=handles.data.samples(selectedsamplenum).eventdata;
data.bridgeddata=handles.data.samples(selectedsamplenum).bridgeddata;
if isfield(handles.data.samples,'BrainStateData')
    data.BrainStateData=handles.data.samples(selectedsamplenum).BrainStateData;
end
data.starttime=handles.data.samples(selectedsamplenum).starttime;
data.endtime=handles.data.samples(selectedsamplenum).endtime;
data.dirs=dirs;
data.xlsdata=xlsdata;
data.icxlsnum=icxlsnum;
data.fieldxlsnum=find(strcmp(handles.data.xlsdata(icxlsnum).HEKAfname,{handles.data.xlsdata.HEKAfname}) & [handles.data.xlsdata.field]==1);
aE_ic_field_analysis_GUI(data)



% additionaldata=struct;
% additionaldata.eventdata=handles.data.samples(selectedsamplenum).eventdata;
% if isfield(handles.data.samples,'BrainStateData')
%     additionaldata.BrainStateData=handles.data.samples(selectedsamplenum).BrainStateData;
% end
% 
% dataout=aE_ic_field_analysis(dirs,xlsdata,icxlsnum,timeborders,'trough',additionaldata); %peak
% valtozok_fieldplot.timebefore=.5;%1;%.5;
% valtozok_fieldplot.timestep=.025;
% valtozok_fieldplot.timeafter=.5;%;1;%.5;
% % valtozok_fieldplot.ICylim
% % valtozok_fieldplot.Fieldylim
% valtozok_fieldplot.highlightaAP=1;
% valtozok_fieldplot.highlightaAP_timeback=.002;
% % valtozok_fieldplot.highlightaAP_timeforward=.02;
% 
% if isfield(additionaldata,'BrainStateData')
%     statestodo=[unique({additionaldata.BrainStateData.name}),'All'];
% else
%     statestodo={'All'};
% end
% for statei=1:length(statestodo)
%     statename=statestodo{statei};
%     %%
%     FieldData=dataout.FieldData;
%     if ~strcmp(statename,'All')
%         FieldData=FieldData(strcmp({FieldData.brainstatename},statename));
%     end
%     valtozok_fieldplot.figurenum=10+statei;
%     persistent_generatefigures_2017_fieldanal(FieldData,valtozok_fieldplot)
%     subplot(4,1,1);
%     title(statename)
% end





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
            valtozok.subtractoffset=get(handles.checkbox9,'Value');
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

panelnames=get(handles.popupmenu10,'string');
selectedpanel=panelnames{get(handles.popupmenu10,'Value')};

if strcmp(selectedpanel,'upper')
    popupname='popupmenu4';
elseif strcmp(selectedpanel,'lower')
    popupname='popupmenu7';
end

if get(handles.(popupname),'Value')==length(handles.data.samples)+3
    set(handles.edit9,'String',handles.data.samples(selectedsamplenum).PSDdata_fieldtoplot.intensitypercentiles(99));
elseif  get(handles.(popupname),'Value')==length(handles.data.samples)+6
    set(handles.edit9,'String',handles.data.samples(selectedsamplenum).BreathingPSDdatatoplot.intensitypercentiles(99));
    elseif  get(handles.(popupname),'Value')==length(handles.data.samples)+8
        set(handles.edit9,'String',handles.data.samples(selectedsamplenum).PSDdata_ictoplot.intensitypercentiles(99));
end
    
        
        
        

% plotdetails=handles.data.samples(selectedsamplenum).plotdetails;


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

BrainStateRules=handles.data.samples(selectedsamplenum).BrainStateRules;
save([handles.data.dirs.brainstatedir,handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID}],'BrainStateData','BrainStateRules');
guidata(hObject,handles);
plotandupdate(handles);

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
guidata(hObject, handles)
plotandupdate(handles)



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


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%
selectedsamplenum=get(handles.popupmenu1,'Value');
handles.data.selectedsamplenum=selectedsamplenum;
aE_BrainState_annotator(handles.data.samples(selectedsamplenum),handles.data.BrainStateProps);
% pupildata=handles.data.samples(selectedsamplenum).pupildata;
% movementdata=handles.data.samples(selectedsamplenum).movementdata;
guidata(hObject,handles)


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selectedsamplenum=get(handles.popupmenu1,'Value');
bridgeddata=handles.data.samples(selectedsamplenum).bridgeddata;
eventdata=handles.data.samples(selectedsamplenum).eventdata;
stimdata=handles.data.samples(selectedsamplenum).stimdata;

[MembPot, countPersistentStim]=persistentDiffMPFunction(eventdata,stimdata,bridgeddata);
sucStim=size(MembPot,2);
if MembPot(1).Before==0
    sucStim=0;
end
if MembPot(1).Before~=0 && MembPot(2).Before==0
   sucStim=1; 
end    

figure('Name', 'Membrane Potential Before And After Of Peristent Stimulation')
scatter(1:size(MembPot,2),[MembPot.Before])
hold on
scatter(1:size(MembPot,2),[MembPot.After], 'filled')
str=sprintf('Count Of Persistent Stim= %d \n Successful stim= %d',countPersistentStim, sucStim);
title(str)
xlabel('Number Of Events, Where were Evoked Spikes')
ylabel('Membrane Potential mV*10^(-2)')
legend({'Before','After'},'Location','northwest')
%legend(strcat('Count Of Persistent Stim=', num2str(countPersistentStim)))
hold off


% disp('lol')


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9
%%
% disp('lol')
%%
if get(hObject,'Value')==true
    hossz=30;%seconds for lowpass filter and averaging
    selectedsamplenum=get(handles.popupmenu1,'Value');
    xlsnum=handles.data.samples(selectedsamplenum).selectedID-1;
    fieldxlsnum=find(strcmp(handles.data.xlsdata(xlsnum).HEKAfname,{handles.data.xlsdata.HEKAfname}) & [handles.data.xlsdata.field]==1);
    if ~isempty(fieldxlsnum)
        a=dir([handles.data.dirs.offsetdir,handles.data.xlsdata(xlsnum).ID,'.mat']);
        bridgeddata=handles.data.samples(selectedsamplenum).bridgeddata;
        if isempty(a)
            load([handles.data.dirs.rawexporteddir,handles.data.xlsdata(fieldxlsnum).ID],'rawdata');
            fielddata=rawdata;
            %%
            
            offsetdata=struct;
            startvalue=[];
            for sweep=1:length(bridgeddata)
                fieldsweep=find([fielddata.realtime]==bridgeddata(sweep).realtime);
                offsetdata(sweep).realtime=bridgeddata(sweep).realtime;
                si=bridgeddata(sweep).si;
                y=fielddata(fieldsweep).y;
                step=round(hossz/si);
                offsetdata(sweep).hossz=length(y)*si;
                offsetdata(sweep).startvalue=median(y(1:min(length(y),step)));
                offsetdata(sweep).endvalue=median(y(max(1,length(y)-step):end));
                if offsetdata(sweep).hossz>hossz
                    y_filt=[y(step:-1:1),y,y(end:-1:end-step+1)];
                    ninq=(1/si)/2;
                    [b,a]=butter(1,(1/hossz)/ninq,'low');
                    y_filt=filtfilt(b,a,y_filt);
                    y_filt=y_filt(step+1:end-step);
                    offsetdata(sweep).y=y_filt;
                else
                    offsetdata(sweep).y=ones(size(y))*mean([offsetdata(sweep).startvalue,offsetdata(sweep).endvalue]);
                end
                if isempty(startvalue)
                    startvalue=offsetdata(sweep).startvalue; % define the starting point
                end
                offsetdata(sweep).y=offsetdata(sweep).y-startvalue;
            end
            save([handles.data.dirs.offsetdir,handles.data.xlsdata(xlsnum).ID],'offsetdata');
        else
            load([handles.data.dirs.offsetdir,handles.data.xlsdata(xlsnum).ID],'offsetdata');
        end
        for sweep=1:length(bridgeddata)
            bridgeddata(sweep).offset=offsetdata(sweep).y;
        end
        handles.data.samples(selectedsamplenum).bridgeddata=bridgeddata;
        %%
        %     events=handles.data.samples(selectedsamplenum).eventdata;
        
    else
        disp('there is no field so can''t cancel out offset..')
        set(hObject,'Value',false);
    end
    guidata(hObject,handles)
end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selectedsamplenum=get(handles.popupmenu1,'Value');
starttime=str2num(get(handles.edit4,'String'));
endtime=str2num(get(handles.edit5,'String'));
bridgeddata=handles.data.samples(selectedsamplenum).bridgeddata;
xlsnum=handles.data.samples(selectedsamplenum).selectedID-1;
fieldxlsnum=find(strcmp(handles.data.xlsdata(xlsnum).HEKAfname,{handles.data.xlsdata.HEKAfname}) & [handles.data.xlsdata.field]==1);

load([handles.data.dirs.rawexporteddir,handles.data.xlsdata(fieldxlsnum).ID],'rawdata');
fielddata=rawdata;
% load([handles.data.dirs.breathingdir,handles.data.xlsdata(fieldxlsnum).ID],'rawdata');
% bridgeddata=rawdata;
apdata=handles.data.samples(selectedsamplenum).eventdata;
apdata=apdata(strcmp({apdata.type},'AP'));
%%
neededicsweeps=find([bridgeddata.realtime]>starttime&[bridgeddata.realtime]<endtime);
for sweepi=1:length(neededicsweeps)
    icsweepnum=neededicsweeps(sweepi);
    fieldsweepnum=find([fielddata.realtime]==bridgeddata(icsweepnum).realtime);
    if ~isempty(fieldsweepnum)
        
        si=bridgeddata(icsweepnum).si;
        ic=bridgeddata(icsweepnum).y';

%%
        field=-1*fielddata(fieldsweepnum).y';
        time=[1:length(ic)]*si;
        hossz=length(ic)*si;
        
        if hossz>3
                    %removing APs
%         apdatanow=apdata([apdata.sweepnum]==icsweepnum);
%         for api=length(apdatanow):-1:1
%             ic(apdatanow(api).onseth:round(apdatanow(api).halfwidth/si/1000*2)+apdatanow(api).onseth)=apdatanow(api).baselineval;
%         end

            %%
%             d{1} = designfilt(['lowpassiir'],'PassbandFrequency',10,'StopbandFrequency',20,'SampleRate',1/si,'DesignMethod','butter');
        d{2} = designfilt(['highpassiir'],'PassbandFrequency',1   ,'StopbandFrequency',.1,'SampleRate',1/si,'DesignMethod','butter');
%         ic=filtfilt(d{1},ic);
        ic=filtfilt(d{2},ic);
%         field=filtfilt(d{1},field);
        field=filtfilt(d{2},field);

        %%
        downratio=round(.001/si);
        ic=downsample(ic,downratio);
        field=downsample(field,downratio);
        time=downsample(time,downratio);
        si=si*downratio;
            movingwin=[5, .5];
            bandwidth=1;%Hz
            params=struct;
%             params.tapers=[.5, movingwin(1) ,1];
            params.tapers=[movingwin(1)*bandwidth movingwin(1)*bandwidth*2-1]
            params.Fs=1/si;%sampling rate
            params.fpass=[0.5 10];%frequency range
            params.pad=2;%frequency resolution
            params.err=[2 .05];
            [S1d,f1d] = mtspectrumc( [ic,field], params );
            [C1d,phi1d,S12,S1,S2,f1dc]=coherencyc(field,ic,params);

            
%             [Sfield,t,f] = mtspecgramc(field, movingwin, params );
%             [Sic,t,f] = mtspecgramc(ic, movingwin, params );
            
            [C,phi,Sicfield,Sfield,Sic,t,f,confC,phistd,Cerr]=cohgramc(field,ic,movingwin,params);
            %%
            figure(5)
            clf
            subplot(5,2,1)
            plot(time,ic)
            subplot(5,2,3)
            plot(time,field)
            subplot(5,2,2)
            imagesc(t,f,Sic')
            set(gca,'YDir','normal')
            colormap linspecer
            subplot(5,2,4)
            imagesc(t,f,Sfield')
            set(gca,'YDir','normal')
            colormap linspecer
            subplot(5,2,6)
            imagesc(t,f,C')
            set(gca,'YDir','normal')
            colormap linspecer
            title('coherence')
            subplot(5,2,5)
            imagesc(t,f,phi')
            set(gca,'YDir','normal')
            colormap linspecer
            title('phase')
            subplot(5,2,7)
            imagesc(t,f,phistd')
            set(gca,'YDir','normal')
            colormap linspecer
            title('phase std')
            subplot(5,2,8)
            imagesc(t,f,abs(Sicfield)')
            set(gca,'YDir','normal')
            colormap linspecer
            title('cross spectrum')
            subplot(5,2,9)
            plot(f1d,S1d./repmat(max(S1d),size(S1d,1),1))
            legend('ic','field')
             subplot(5,2,10)
            plot(f1dc,C1d)
%             pause
        end
    end
end
disp('lol')
%%

reclength=handles.data.samples(selectedsamplenum).bridgeddata(end).realtime+handles.data.samples(selectedsamplenum).bridgeddata(end).si*length(handles.data.samples(selectedsamplenum).bridgeddata(end).y)-handles.data.samples(selectedsamplenum).bridgeddata(1).realtime;
stepsize=1;
window=30;
time=0:stepsize:reclength;
aAPnum=zeros(size(time));
aAPfreq=zeros(size(time));
aAPidxes=find([handles.data.samples(selectedsamplenum).eventdata.axonalAP]);%&~ [handles.data.samples(selectedsamplenum).eventdata.stimulated]
aAPtimes=[handles.data.samples(selectedsamplenum).eventdata(aAPidxes).maxtime]-handles.data.samples(selectedsamplenum).bridgeddata(1).realtime;
sAPidxes=find([handles.data.samples(selectedsamplenum).eventdata.somaticAP] );%&~ [handles.data.samples(selectedsamplenum).eventdata.stimulated]
sAPtimes=[handles.data.samples(selectedsamplenum).eventdata(sAPidxes).maxtime]-handles.data.samples(selectedsamplenum).bridgeddata(1).realtime;
for i=1:length(time)
    timestart=max(0,time(i)-window/2);
    timeend=min(time(end),time(i)+window/2);
    timeinterval=timeend-timestart;
    aAPnum(i)=sum(aAPtimes>time(i)-window/2 & aAPtimes<=time(i)+window/2);
    aAPfreq(i)=aAPnum(i)/timeinterval;
    sAPnum(i)=sum(sAPtimes>time(i)-window/2 & sAPtimes<=time(i)+window/2);
    sAPfreq(i)=sAPnum(i)/timeinterval;
end
%%
figure(33)
clf

plot(time,aAPfreq,'r-')
hold on
plot(time,sAPfreq,'k-')
xlabel('time since recording (sec)')
ylabel('AP freqency')
legend('axonal','somatic')


% --- Executes on selection change in popupmenu9.
function popupmenu9_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu9 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu9


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


% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10
selectedsamplenum=get(handles.popupmenu1,'Value');
set(handles.edit9,'String', handles.data.samples(selectedsamplenum).plotdetails.PSD.cmax(get(hObject,'Value')));


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


% --- Executes on button press in AxSomLDA.
function AxSomLDA_Callback(hObject, eventdata, handles)
% hObject    handle to AxSomLDA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

   %Initialize set of axonal and somatic spikes
selectedsamplenum=get(handles.popupmenu1,'Value');
eventdata=handles.data.samples(selectedsamplenum).eventdata;

%call ax_somClassifier function, return: sorted spikes, classes: axonal,somatic
[aAPsMatrix,sAPsMatrix]=ax_somClassifier(eventdata);

%Lda Viewer
figure;
  %Plot classified points
 X=[aAPsMatrix;sAPsMatrix];
 
plot(X(1:length(X)/2, 1), X(1:length(X)/2, 2),'rx', X((length(X)/2)+1:end, 1), X((length(X)/2)+1:end, 2),'bo');
title('Axonal-Somatic LDA')
hold on;

minVal = min(min(X)) - 1;
maxVal = max(max(X)) + 1;
axis([minVal maxVal minVal maxVal]);
axis square;

fprintf('Program paused. Press enter to continue.\n');
pause;

 %=============== Part 2: Performing Linear Discriminant Analysis ===============
fprintf('\nRunning LDA on example dataset.\n\n');
X1=X(1:length(X)/2,:);
X2=X((length(X)/2)+1:end,:);
m1 = mean(X1);
m2 = mean(X2);
%fprintf('Centroid of class 1 is at [%f %f] \n', m1(1), m1(2), m1(3), m1(4));
%fprintf('Centroid of class 2 is at [%f %f] \n', m2(1), m2(2), m1(3), m1(4));

Sb = (m1-m2)'*(m1-m2); % this is the rank-1 between-class covariance matrix

% these are the class-specific centralized data matrices
Xc1 = bsxfun(@minus, X1, m1);
Xc2 = bsxfun(@minus, X2, m2);

S1 = Xc1'*Xc1; S2=Xc2'*Xc2; % these are the class-specific within-class covariance matrices
Sw = S1+S2; % this is the aggregated within-class covariance matrix

% Performing LDA
[U, S] = eig(inv(Sw)*Sb);

[v, k] = max(diag(S));
w = U(:, k);
% an equivalent form would be w = inv(Sw)*(m1-m2)'; for maximizing the fraction and w = (m1-m2)'; for maximizing the nominator

drawLine = @(p1, p2, varargin) plot([p1(1) p2(1)], [p1(2) p2(2)], varargin{:});
% Draw the direction of the optimal projection w.
multipier = max(maxVal/abs(w(1)), maxVal/abs(w(1))); % this is just a stretching factor for w to be drawn
drawLine(-multipier * w, multipier * w, '-k', 'LineWidth', 2);

fprintf('w = [%f %f] \n', w(1), w(2));

fprintf('Program paused. Press enter to continue.\n');
pause;

% =================== Part 3: Dimension Reduction ===================
%  The code here plots the data points into the reduced-dimensional subspace.
%
fprintf('\nDimension reduction on example dataset.\n\n');
fprintf('The original first example: [%f %f]\n', X(1, 1), X(1, 2));
Z = X * w;
fprintf('Projection of the first example: %f\n', Z(1));

X_rec = Z * w';
fprintf('Approximation of the first example: [%f %f]\n', X_rec(1, 1), X_rec(1, 2));
rm1 = mean(X_rec(1:length(X)/2,:));
rm2 = mean(X_rec((length(X)/2)+1:end,:));
rm3 = m1*w*w';
rm4 = m2*w*w';
fprintf('\nThe reconstructed centroid of class 1 at = [%f %f] \n', rm3(1), rm3(2));
fprintf('The reconstructed centroid of class 2 at = [%f %f] \n\n', rm4(1), rm4(2));

%  Draw lines connecting the projected points to the original points
plot(X_rec(1:length(X)/2, 1), X_rec(1:length(X)/2, 2),'rx', X_rec((length(X)/2)+1:end, 1), X_rec((length(X)/2)+1:end, 2),'bo');


for i = 1:size(X, 1)
  drawLine(X(i,:), X_rec(i,:), '--k', 'LineWidth', 1);
end
fprintf('The optimal value of w^T*S_B*w is %f\n', (w'*Sb*w)); % this equals the squared l2-distance of the projected class means
fprintf('The optimal value of w^T*S_W*w is %f\n', (w'*Sw*w));
fprintf('The optimal value of w^T*S_B*w/w^T*S_W*w is %f\n', (w'*Sb*w)/(w'*Sw*w));
%distance,similarity to spike vector
axsomdistancesimilarity(aAPsMatrix,sAPsMatrix);
