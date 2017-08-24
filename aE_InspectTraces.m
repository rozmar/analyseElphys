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

% Last Modified by GUIDE v2.5 02-Jun-2017 16:05:49

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
handles=loadthedata(handles);
updategui(handles)
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
IDs=get(hObject,'String');
selectedID=get(hObject,'Value');
handles.data.samples(selectedsamplenum).selectedID=find(strcmp(IDs{selectedID},handles.data.IDs));
if handles.data.samples(selectedsamplenum).selectedID>1
    load([handles.data.dirs.bridgeddir,handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID}],'lightdata');
    load([handles.data.dirs.bridgeddir,handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID}],'bridgeddata');
    handles.data.samples(selectedsamplenum).lightdata=lightdata;
    handles.data.samples(selectedsamplenum).bridgeddata=bridgeddata;
%     eventdataorig=eventdata;
    if ~isempty(dir([handles.data.dirs.eventdir,handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID},'.mat']))
        load([handles.data.dirs.eventdir,handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID}],'eventdata');
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
    if ~isempty(eventdata)
        stimapnum=length(find(strcmp({handles.data.samples(selectedsamplenum).eventdata.type},'AP') & [handles.data.samples(selectedsamplenum).eventdata.stimulated]));
        persapnum=length(find(strcmp({handles.data.samples(selectedsamplenum).eventdata.type},'AP') & ~[handles.data.samples(selectedsamplenum).eventdata.stimulated]));
    else
        stimapnum=0;
        persapnum=0;
    end
    
    if isfield(handles.data.xlsdata,'drugdata') & ~isempty(handles.data.xlsdata(handles.data.samples(selectedsamplenum).selectedID-1).drugdata)
        drugnames={handles.data.xlsdata(handles.data.samples(selectedsamplenum).selectedID-1).drugdata.DrugName};
        drugwashintimes=[handles.data.xlsdata(handles.data.samples(selectedsamplenum).selectedID-1).drugdata.DrugWashinTime];
        drugstring=' ';
        for i=1:length(drugnames)
            drugstring=[drugstring,drugnames{i},'-',num2str(drugwashintimes(i)),';  '];
        end
    else
        drugstring='   no drugs';
    end
    set(handles.text5,'String',['stimulated AP: ', num2str(stimapnum),'  persistent AP: ', num2str(persapnum),'       ',drugstring]);
else
    handles.data.samples(selectedsamplenum).lightdata=struct;
end
plotandupdate(handles);
updategui(handles)

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
updategui(handles)


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
needed=sort(unique(needed));
common=ismember(handles.data.samples(selectedsamplenum).neededwaves,needed);
handles.data.samples(selectedsamplenum).neededwaves=handles.data.samples(selectedsamplenum).neededwaves(common);
handles.data.samples(selectedsamplenum).neededwavespervariable{selectedvariable}=handles.data.samples(selectedsamplenum).neededwaves;
handles.data.samples(selectedsamplenum).changes=[handles.data.samples(selectedsamplenum).changes,' AND ',szovegmost];

% timesofneededwaves=[handles.data.samples(selectedsamplenum).lightdata(handles.data.samples(selectedsamplenum).neededwaves).realtime];
% reallyneededwaves=timesofneededwaves>=str2num(get(handles.edit4,'String')) & timesofneededwaves<=str2num(get(handles.edit5,'String'));
% handles.data.samples(selectedsamplenum).neededwaves=handles.data.samples(selectedsamplenum).neededwaves(reallyneededwaves);

handles=loadthedata(handles);
handles=updatedatatoplot(handles);
guidata(hObject,handles)
%plotandupdate(handles)


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
        data.samples(i).eventdata=struct;
        data.samples(i).selectedID=1;
        data.samples(i).loadedID=1;
        data.samples(i).lightdata=[];
        data.samples(i).cutoffreq=0;
        data.samples(i).filterdeg=3;
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
    set(handles.popupmenu2,'String',handles.data.IDs);
    set(handles.popupmenu2,'Value',1);
    set(handles.text10,'String',[]);
    
    %%
    set(handles.popupmenu3,'String',handles.data.fieldnevek);
    set(handles.popupmenu4,'String',[{'stiumulus'},{handles.data.samples.marker}]);
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
        axes(handles.axes2)
        if get(handles.popupmenu4,'Value')-1==0
            for sweepi=1:length(neededwaves)
                ettol=find(handles.data.samples(samplenum).datatoplot(sweepi).x>starttime,1,'first');
            eddig=find(handles.data.samples(samplenum).datatoplot(sweepi).x<endtime,1,'last');
                plot(handles.data.samples(samplenum).datatoplot(sweepi).x(ettol:eddig),handles.data.samples(samplenum).datatoplot(sweepi).ycurrent(ettol:eddig),marker(1),'LineWidth',2)
            end
        elseif samplenum==get(handles.popupmenu4,'Value')-1
            for sweepi=1:length(neededwaves)
                ettol=find(handles.data.samples(samplenum).datatoplot(sweepi).x>starttime,1,'first');
                 eddig=find(handles.data.samples(samplenum).datatoplot(sweepi).x<endtime,1,'last');
                plot(handles.data.samples(samplenum).datatoplot(sweepi).x(ettol:eddig),handles.data.samples(samplenum).datatoplot(sweepi).yvoltage(ettol:eddig),marker(1),'LineWidth',2)
            end
        end
        %         plot([handles.data.samples(samplenum).datatoplot.x],[handles.data.samples(samplenum).datatoplot.ycurrent],marker(1),'LineWidth',2)
        handles.axes2=gca;
    end
    
end
axes(handles.axes1)
axis tight
xlabel([]);
set(gca,'xtick',[]);
ylabel('Voltage (V)')
handles.axes1=gca;
axes(handles.axes2)
xlabel('Time (s)')
ylabel('Injected Current (pA)')
axis tight
handles.axes2=gca;
linkaxes([handles.axes1,handles.axes2],'x');
java.lang.System.gc()

function updategui(handles)
selectedsamplenum=get(handles.popupmenu1,'Value');
set(handles.checkbox1,'Value',handles.data.samples(selectedsamplenum).switch);
selectedID=handles.data.samples(selectedsamplenum).selectedID;
selectedvariable=handles.data.samples(selectedsamplenum).selectedvariable;
% set(handles.popupmenu2,'Value',selectedID);
set(handles.popupmenu3,'Value',selectedvariable);
set(handles.edit1,'String',num2str(handles.data.samples(selectedsamplenum).cutoffreq));
set(handles.edit2,'String',num2str(handles.data.samples(selectedsamplenum).filterdeg));
startval=get(handles.slider1,'Value');
endval=get(handles.slider2,'Value');
midval=mean([startval,endval]);
set(handles.slider3,'Value',midval);

% timedifi=handles.data.samples(selectedsamplenum).endtime-handles.data.samples(selectedsamplenum).starttime;
starttime=min([handles.data.samples.starttime]);
endtime=max([handles.data.samples.endtime]);
timedifi=endtime-starttime;
set(handles.edit4,'String',starttime+timedifi*startval);
set(handles.edit5,'String',starttime+timedifi*endval);
set(handles.edit7,'String',starttime+timedifi*midval);
set(handles.edit8,'String',timedifi*(endval-startval));
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
    plotdetails=handles.data.samples(selectedsamplenum).plotdetails;
    if ~isempty(neededwaves)
    if isempty(fieldnames(plotdetails)) | ~(length(neededwaves)==length(plotdetails.neededwaves)) |~(neededwaves==plotdetails.neededwaves) | ~(neededsamplenum==plotdetails.neededsamplenum) | ~all(cutoffreq==plotdetails.cutoffreq) | ~(filterdeg==plotdetails.filterdeg) | ~(shoulddownsample==plotdetails.shoulddownsample)
        datatoplot=struct;
        plotdetails.neededwaves=neededwaves;
        plotdetails.neededsamplenum=neededsamplenum;
        plotdetails.cutoffreq=cutoffreq;
        plotdetails.filterdeg=filterdeg;
        plotdetails.shoulddownsample=shoulddownsample;
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
updategui(handles)
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
updategui(handles)
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
updategui(handles)
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
% timedifi=handles.data.samples(selectedsamplenum).endtime-handles.data.samples(selectedsamplenum).starttime;
timedifi=max([handles.data.samples.endtime])-min([handles.data.samples.starttime]);
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
updategui(handles)
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
% timedifi=handles.data.samples(selectedsamplenum).endtime-handles.data.samples(selectedsamplenum).starttime;
timedifi=max([handles.data.samples.endtime])-min([handles.data.samples.starttime]);
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
updategui(handles)
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
% timedifi=handles.data.samples(selectedsamplenum).endtime-handles.data.samples(selectedsamplenum).starttime;
timedifi=max([handles.data.samples.endtime])-min([handles.data.samples.starttime]);
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
updategui(handles)
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
figure(1)
close(1)
figure(2)
close(2)
figure(3)
close(3)
figure(4)
close(4)
figure(5)
close(5)
figure(6)
close(6)
figure(33)
close(33)
timeback=.002;
timeforward=.002;
depolratewindow=.0015;
baselineVwindow=.005;
APwaves=struct;
if get(hObject,'Value')==true
    timeborders=[str2num(get(handles.edit4,'String')),str2num(get(handles.edit5,'String'))];
    for samplenum=1:length(handles.data.samples)
        eventdata=handles.data.samples(samplenum).eventdata;
        bridgeddata=handles.data.samples(samplenum).bridgeddata;
        stimdata=handles.data.samples(samplenum).stimdata;
        if ~isempty(eventdata) & ~isempty(fieldnames(eventdata))
            apidxes=find(strcmp({eventdata.type},'AP'));
            prevsweepnum=NaN;
            for eventi=1:length(eventdata)
                eventdata(eventi).maxdv=0;
            end
            for apii=1:length(apidxes)
                api=apidxes(apii);
                sweepnum=eventdata(api).sweepnum;
                if sweepnum~=prevsweepnum
                    RS=stimdata(sweepnum).RS/10^6;
                    y=bridgeddata(sweepnum).y;
                    si=round(bridgeddata(sweepnum).si*10^6)/10^6;
                    yfiltered=moving(y,2);
                    dyfiltered=diff(yfiltered)/si;
%                     dyfiltered_longfilt=moving(y,round(.0005/si));
                    prevsweepnum=sweepnum;
                    stepbackforthresh=round(.0002/si);
                    stepback=round(timeback/si);
                stepforward=round(timeforward/si);
                end
                onseth=eventdata(api).onseth;
                maxh=eventdata(api).maxh;
                threshh=maxh-stepbackforthresh;
                while dyfiltered(threshh)<max(dyfiltered(max(1,threshh-5):threshh)) & threshh>stepback+3
                    threshh=threshh-1;
                end
                while ~(dyfiltered(threshh)<20 & mean(dyfiltered(max(1,threshh-round(.0001/si)):threshh))<20) & threshh>stepback+3
                    threshh=threshh-2;
                end
                threshv=yfiltered(threshh);
                eventdata(api).threshv=threshv;
                eventdata(api).APamplitude=eventdata(api).maxval-eventdata(api).threshv;
                depolrateV=yfiltered(max(threshh-round(depolratewindow/si),1):threshh);
                depolratet=(1:length(depolrateV))*si;
                p=polyfit(depolratet',depolrateV,1);
                eventdata(api).depolrate=p(1);
                eventdata(api).baselineV=mean(yfiltered(max(threshh-round(baselineVwindow/si),1):threshh));
                
                dv=diff(y(threshh:maxh))/si;
                [maxdval,maxdvh]=max(dv);
                maxdvh=maxdvh+threshh-1;
                centerh=threshh;%maxdvh;
                
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
                movingn=3;
                %             v=y([-stepback:stepforward]+centerh);
                dv=diff(moving(v,movingn))'/si;
                ddv=diff(moving(dv,movingn))'/si;
                dddv=diff(moving(ddv,movingn))'/si;
                vdv=mean([v(2:end);v(1:end-1)]);
                vddv=mean([vdv(2:end);vdv(1:end-1)]);
                vdddv=mean([vddv(2:end);vddv(1:end-1)]);
                tdv=mean([t(2:end);t(1:end-1)]);
                tddv=mean([tdv(2:end);tdv(1:end-1)]);
                tdddv=mean([tddv(2:end);tddv(1:end-1)]);
                APwaves(apii).t=t';
                APwaves(apii).v=v'*1000;
                APwaves(apii).dv=dv';
                APwaves(apii).ddv=ddv';
                APwaves(apii).dddv=dddv';
                APwaves(apii).tdv=tdv';
                APwaves(apii).tddv=tddv';
                APwaves(apii).tdddv=tdddv';
                APwaves(apii).vdv=vdv'*1000;
                APwaves(apii).vddv=vddv'*1000;
                APwaves(apii).vdddv=vdddv'*1000;
                APwaves(apii).si=si;
                APwaves(apii).RS=RS;
                APwaves(apii).maxtime=eventdata(api).maxtime;
                APwaves(apii).stimulated=eventdata(api).stimulated;
                APwaves(apii).maxdv=max(dv);
                eventdata(api).maxdv=max(dv);
                % elso derivalt helyi minimumanak keresese a masodik derivalton
                difi=maxdvh-threshh;
                ddvsegment=ddv(stepback-difi:stepback);
                vege=length(ddvsegment);
                while vege>1 & (ddvsegment(vege)<ddvsegment(vege-1) | ddvsegment(vege)<0)
                    vege=vege-1;
                end
                eleje=1;
                while vege>eleje & (ddvsegment(eleje)<ddvsegment(eleje+1) | ddvsegment(eleje)<0)
                    eleje=eleje+1;
                end
%                 if any(moving(ddvsegment(eleje:vege),2,'max')<0)
%                     %                 figure(1)
%                     %                 clf
%                     %                 plot(ddvsegment,'k-')
%                     %                 hold on
%                     %                 plot(ddvsegment(1:vege),'r-')
%                     %                 pause
%                     eventdata(api).axonalAP=true;
%                     eventdata(api).somaticAP=false;
%                     APwaves(apii).axonalAP=true;
%                     APwaves(apii).somaticAP=false;
%                 else
%                     eventdata(api).axonalAP=false;
%                     eventdata(api).somaticAP=true;
%                     APwaves(apii).axonalAP=false;
%                     APwaves(apii).somaticAP=true;
%                 end
            end

           figure(6)
           clf
            needed=find([eventdata.maxtime]>=timeborders(1) &[eventdata.maxtime]<=timeborders(2) &strcmp({eventdata.type},'AP'));
            plot([eventdata(needed).APamplitude],[eventdata(needed).maxdv],'ko')
            ylabel('maxdV/dt (mv/ms)')
            xlabel('amplitude');
            [~,mindvperdt]=ginput(2);
            figure(3)
            clf
            hold on
            % hist([apdata.threshv],100)
%             needed=find([eventdata.maxtime]>=timeborders(1) &[eventdata.maxtime]<=timeborders(2) & [eventdata.axonalAP] );
%             plot([eventdata(needed).APamplitude],[eventdata(needed).threshv],'ro')
            needed=find([eventdata.maxtime]>=timeborders(1) &[eventdata.maxtime]<=timeborders(2) & [eventdata.maxdv]>min(mindvperdt) &[eventdata.maxdv]<max(mindvperdt));
            plot([eventdata(needed).APamplitude],[eventdata(needed).threshv]*1000,'ko')
            needed=find([eventdata.maxtime]>=timeborders(1) &[eventdata.maxtime]<=timeborders(2) & [eventdata.stimulated]==1 & [eventdata.maxdv]>min(mindvperdt) &[eventdata.maxdv]<max(mindvperdt));
            plot([eventdata(needed).APamplitude],[eventdata(needed).threshv]*1000,'bo')
            ylabel('AP threshold (V)')
            xlabel('AP amplitude (mV)')
            
            figure(4)
            clf
            hold on
            % hist([apdata.threshv],100)
%             needed=find([eventdata.maxtime]>=timeborders(1) &[eventdata.maxtime]<=timeborders(2) & [eventdata.axonalAP] );
%             plot([eventdata(needed).baselineV],[eventdata(needed).threshv],'ro')
            needed=find([eventdata.maxtime]>=timeborders(1) &[eventdata.maxtime]<=timeborders(2)& [eventdata.maxdv]>min(mindvperdt) &[eventdata.maxdv]<max(mindvperdt));
            plot([eventdata(needed).baselineV]*1000,[eventdata(needed).threshv]*1000,'ko')
            needed=find([eventdata.maxtime]>=timeborders(1) &[eventdata.maxtime]<=timeborders(2) & [eventdata.stimulated]==1& [eventdata.maxdv]>min(mindvperdt) &[eventdata.maxdv]<max(mindvperdt));
            plot([eventdata(needed).baselineV]*1000,[eventdata(needed).threshv]*1000,'bo')
            ylabel('AP threshold (mV)')
            xlabel('baseline V0 (mV)')

            figure(5)
            clf
            hold on
            % hist([apdata.threshv],100)
%             needed=find([eventdata.maxtime]>=timeborders(1) &[eventdata.maxtime]<=timeborders(2) & [eventdata.axonalAP] );
%             plot([eventdata(needed).depolrate],[eventdata(needed).threshv],'ro')
            needed=find([eventdata.maxtime]>=timeborders(1) &[eventdata.maxtime]<=timeborders(2) & [eventdata.maxdv]>min(mindvperdt) &[eventdata.maxdv]<max(mindvperdt));
            plot([eventdata(needed).depolrate],[eventdata(needed).threshv]*1000,'ko','LineWidth',3,'MarkerFaceColor',[0 0 0])
            
            needed=find([eventdata.maxtime]>=timeborders(1) &[eventdata.maxtime]<=timeborders(2) & [eventdata.stimulated]==1& [eventdata.maxdv]>min(mindvperdt) &[eventdata.maxdv]<max(mindvperdt));
            plot([eventdata(needed).depolrate],[eventdata(needed).threshv]*1000,'bo','LineWidth',3,'MarkerFaceColor',[0 0 1])
            ylabel('AP threshold (mV)')
            xlabel('rate of depolarization (mV/ms)')
            pause
            [x,y]=ginput(1);
            needed=find([eventdata.maxtime]>=timeborders(1) &[eventdata.maxtime]<=timeborders(2) &[eventdata.maxdv]>min(mindvperdt) &[eventdata.maxdv]<max(mindvperdt));
            plot([x,x],[min([eventdata(needed).threshv]),max([eventdata(needed).threshv])]*1000,'k-','LineWidth',2)
            plot([min([eventdata(needed).depolrate]),max([eventdata(needed).depolrate])],[y,y],'k-','LineWidth',2)
            for api=1:length(eventdata)
                if any(api==apidxes);
                    apii=find(api==apidxes);
                    if [eventdata(api).maxdv]<min(mindvperdt) | [eventdata(api).maxdv]>max(mindvperdt)
                        eventdata(api).axonalAP=false;
                        eventdata(api).somaticAP=false;
                         APwaves(apii).axonalAP=false;
                        APwaves(apii).somaticAP=false;
                        eventdata(api).type='ep';
                    elseif eventdata(api).depolrate<=x & eventdata(api).threshv<=y/1000
                        eventdata(api).axonalAP=true;
                        eventdata(api).somaticAP=false;
                        APwaves(apii).axonalAP=true;
                        APwaves(apii).somaticAP=false;
                    else
                        eventdata(api).axonalAP=false;
                        eventdata(api).somaticAP=true;
                         APwaves(apii).axonalAP=false;
                        APwaves(apii).somaticAP=true;
                    end
                else
                    eventdata(api).axonalAP=false;
                    eventdata(api).somaticAP=false;
                end
            end
            
            for ii=1:2
                figure(ii)
                clf
                sis=unique([APwaves.si]);
                for i=1:length(sis)
                    if ii==1
                        needed=find([APwaves.si]==sis(i) & [APwaves.axonalAP]==true & [APwaves.maxtime]>=timeborders(1) &[APwaves.maxtime]<=timeborders(2) &[APwaves.stimulated]==false);
                    else
                        needed=find([APwaves.si]==sis(i) & [APwaves.somaticAP]==true & [APwaves.maxtime]>=timeborders(1) &[APwaves.maxtime]<=timeborders(2) &[APwaves.stimulated]==false);
                    end
                    subplot(2,2,1)
                    hold on
                    plot([APwaves(needed).t],[APwaves(needed).v]);
                    axis tight
                    ylabel('Voltage (mV)')
                    xlabel('Time (ms)')
                    subplot(2,2,2)
                    hold on
                    plot([APwaves(needed).tdv],[APwaves(needed).dv]);
                    axis tight
                    ylabel('dV/dt (mV/ms)')
                    xlabel('Time (ms)')
                    subplot(2,2,3)
                    hist([APwaves(needed).RS]);
                    xlabel('RS (MOhm)')
                    ylabel('count')
                    subplot(2,2,4)
                    hold on
                    plot([APwaves(needed).vdv],[APwaves(needed).dv]);
                    axis tight
                     ylabel('dV/dt (mV/ms)')
                    xlabel('Voltage (mV)')
                end
            end
            
            
            
            handles.data.samples(samplenum).eventdata=eventdata;
            figure(33)
            clf
            hist([APwaves.RS]);
            xlabel('RS (MOhm)')
            axonalapnum=length(find([APwaves.axonalAP]==true&[APwaves.stimulated]==false& [APwaves.maxtime]>=timeborders(1) &[APwaves.maxtime]<=timeborders(2)  ));
            somaticaptnum=length(find([APwaves.somaticAP]==true & [APwaves.stimulated]==false& [APwaves.maxtime]>=timeborders(1) &[APwaves.maxtime]<=timeborders(2) ));
            stimulatedapnum=length(find([APwaves.somaticAP]==true&[APwaves.stimulated]==true& [APwaves.maxtime]>=timeborders(1) &[APwaves.maxtime]<=timeborders(2) ));
            disp([num2str(axonalapnum),' axonal APs and ',num2str(somaticaptnum),' somatic APs found - ' ,num2str(stimulatedapnum),'AP by somatic current injection during the inspected ',num2str(diff(timeborders)),' seconds'])
            
            ID=handles.data.IDs{handles.data.samples(samplenum).loadedID};
            
            figure(5)
            saveas(gcf,[handles.data.dirs.figuresdir,ID,'_threshold_rateofdepol'],'pdf')
            saveas(gcf,[handles.data.dirs.figuresdir,ID,'_threshold_rateofdepol'],'jpg')
            figure(1)
            saveas(gcf,[handles.data.dirs.figuresdir,ID,'_axonal_spikes'],'pdf')
            saveas(gcf,[handles.data.dirs.figuresdir,ID,'_axonal_spikes'],'jpg')
            figure(2)
            saveas(gcf,[handles.data.dirs.figuresdir,ID,'_somatic_spikes'],'pdf')
            saveas(gcf,[handles.data.dirs.figuresdir,ID,'_somatic_spikes'],'jpg')
        end
    end
    
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
updategui(handles)
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
        set(handles.popupmenu2,'String',[{'no cell selected'},{xlsdata(ezaz).ID}]);
%         else
%             set(handles.popupmenu2,'String',[{'no cell selected'},{xlsdata.ID}]);
    end
    set(handles.text10,'String',searchtext);
end

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
 set(handles.popupmenu2,'String',[{'no cell selected'},{xlsdata.ID}]);
 set(handles.text10,'String','');


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
aE_ic_field_analysis(dirs,xlsdata,icxlsnum,timeborders,'trough',additionaldata)