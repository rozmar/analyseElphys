function varargout = aE_InspectTraces(varargin)
% AE_INSPECTTRACES MATLAB code for aE_InspectTraces.fig

%USAGE: aE_InspectTraces(dirs,xlsdata) - after starting
%analyseElphys_main.m


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

% Last Modified by GUIDE v2.5 03-Jan-2017 16:18:40

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
handles=updatedatatoplot(handles);
guidata(hObject,handles);
plotandupdate(handles);


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
selectedsamplenum=get(handles.popupmenu1,'Value');
handles.data.samples(selectedsamplenum).selectedID=get(hObject,'Value');
if handles.data.samples(selectedsamplenum).selectedID>1
    load([handles.data.dirs.bridgeddir,handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID}],'lightdata');
    handles.data.samples(selectedsamplenum).lightdata=lightdata;
    eventdataorig=eventdata;
    load([handles.data.dirs.eventdir,handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID}],'eventdata');
    handles.data.samples(selectedsamplenum).eventdata=eventdata;
    eventdata=eventdataorig;
    handles.data.samples(selectedsamplenum).starttime=handles.data.samples(selectedsamplenum).lightdata(1).realtime;
    handles.data.samples(selectedsamplenum).endtime=handles.data.samples(selectedsamplenum).lightdata(end).realtime;
else
    handles.data.samples(selectedsamplenum).lightdata=struct;
end
handles.data.samples(selectedsamplenum).neededwaves=[];
handles=loadthedata(handles);
handles=updatedatatoplot(handles);
guidata(hObject,handles);
set(handles.listbox1,'Value',1)
stimapnum=length(find(strcmp({handles.data.samples(selectedsamplenum).eventdata.type},'AP') & [handles.data.samples(selectedsamplenum).eventdata.stimulated]));
persapnum=length(find(strcmp({handles.data.samples(selectedsamplenum).eventdata.type},'AP') & ~[handles.data.samples(selectedsamplenum).eventdata.stimulated]));

if ~isempty(handles.data.xlsdata(handles.data.samples(selectedsamplenum).selectedID-1).drugdata)
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
handles.data.samples(selectedsamplenum).neededwaves=sort(unique(needed));
handles.data.samples(selectedsamplenum).neededwavespervariable{selectedvariable}=handles.data.samples(selectedsamplenum).neededwaves;
timesofneededwaves=[handles.data.samples(selectedsamplenum).lightdata(handles.data.samples(selectedsamplenum).neededwaves).realtime];
reallyneededwaves=timesofneededwaves>=str2num(get(handles.edit4,'String')) & timesofneededwaves<=str2num(get(handles.edit5,'String'));
handles.data.samples(selectedsamplenum).neededwaves=handles.data.samples(selectedsamplenum).neededwaves(reallyneededwaves);
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

timesofneededwaves=[handles.data.samples(selectedsamplenum).lightdata(handles.data.samples(selectedsamplenum).neededwaves).realtime];
reallyneededwaves=timesofneededwaves>=str2num(get(handles.edit4,'String')) & timesofneededwaves<=str2num(get(handles.edit5,'String'));
handles.data.samples(selectedsamplenum).neededwaves=handles.data.samples(selectedsamplenum).neededwaves(reallyneededwaves);

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
    %%
    set(handles.popupmenu3,'String',handles.data.fieldnevek);
    set(handles.slider1,'Min',0);
    set(handles.slider1,'Max',1);
    set(handles.slider2,'Min',0);
    set(handles.slider2,'Max',1);
    set(handles.slider1,'Value',0);
    set(handles.slider2,'Value',1);
    %%
    % set(handles.popupmenu3,'Value',1);
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
 handles=updatedatatoplot(handles);
selectedsamplenum=get(handles.popupmenu1,'Value');
set(handles.checkbox1,'Value',handles.data.samples(selectedsamplenum).switch);
selectedID=handles.data.samples(selectedsamplenum).selectedID;
selectedvariable=handles.data.samples(selectedsamplenum).selectedvariable;
set(handles.popupmenu2,'Value',selectedID);
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
cla reset
hold on
handles.axes1=gca;

axes(handles.axes2)
cla reset
hold on
handles.axes2=gca;
for samplenum=1:length(handles.data.samples)
    y0=0;
    if handles.data.samples(samplenum).switch==1 & handles.data.samples(samplenum).selectedID>1 & ~isempty(handles.data.samples(samplenum).neededwaves)
        neededwaves=handles.data.samples(samplenum).neededwaves;
        marker=handles.data.samples(samplenum).marker;
        axes(handles.axes1)
        
        
        for sweepi=1:length(neededwaves)
            plot(handles.data.samples(samplenum).datatoplot(sweepi).x,handles.data.samples(samplenum).datatoplot(sweepi).yvoltage,marker(1),'LineWidth',2)
        end
        markevents=get(handles.checkbox3,'Value');
        if markevents==1
            neededevents=strcmp({handles.data.samples(samplenum).eventdata.type},'AP') & ~[handles.data.samples(samplenum).eventdata.stimulated];
            eventsnow=handles.data.samples(samplenum).eventdata(neededevents);
            neededevents=false(size(eventsnow));
              for sweepi=1:length(neededwaves)
                 sweepnum=neededwaves(sweepi);
                 neededevents([eventsnow.sweepnum]==sweepnum)=1;
              end
             plot([eventsnow(neededevents).maxtime],[eventsnow(neededevents).maxval],'ro','MarkerSize',8)
       end
%         plot([handles.data.samples(samplenum).datatoplot.x],[handles.data.samples(samplenum).datatoplot.yvoltage],marker(1),'LineWidth',2)
        y0=max([y0,max([handles.data.samples(samplenum).datatoplot.yvoltage])]);
        if length(neededwaves)>0
            prenum=handles.data.samples(samplenum).selectedID-1;
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
        handles.axes1=gca;
        axes(handles.axes2)
        for sweepi=1:length(neededwaves)
            plot(handles.data.samples(samplenum).datatoplot(sweepi).x,handles.data.samples(samplenum).datatoplot(sweepi).ycurrent,marker(1),'LineWidth',2)
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
set(handles.popupmenu2,'Value',selectedID);
set(handles.popupmenu3,'Value',selectedvariable);
set(handles.edit1,'String',num2str(handles.data.samples(selectedsamplenum).cutoffreq));
set(handles.edit2,'String',num2str(handles.data.samples(selectedsamplenum).filterdeg));
startval=get(handles.slider1,'Value');
endval=get(handles.slider2,'Value');
timedifi=handles.data.samples(selectedsamplenum).endtime-handles.data.samples(selectedsamplenum).starttime;
set(handles.edit4,'String',handles.data.samples(selectedsamplenum).starttime+timedifi*startval);
set(handles.edit5,'String',handles.data.samples(selectedsamplenum).starttime+timedifi*endval);
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
selectedsamplenum=get(handles.popupmenu1,'Value');
if handles.data.samples(selectedsamplenum).switch>0
    neededwaves=handles.data.samples(selectedsamplenum).neededwaves;
    neededsamplenum=str2num(get(handles.edit3,'String'));
    cutoffreq=handles.data.samples(selectedsamplenum).cutoffreq;
    filterdeg=handles.data.samples(selectedsamplenum).filterdeg;
    shoulddownsample=get(handles.checkbox2,'Value');
    plotdetails=handles.data.samples(selectedsamplenum).plotdetails;
    
    if isempty(fieldnames(plotdetails)) | ~(length(neededwaves)==length(plotdetails.neededwaves)) |~(neededwaves==plotdetails.neededwaves) | ~(neededsamplenum==plotdetails.neededsamplenum) | ~(cutoffreq==plotdetails.cutoffreq) | ~(filterdeg==plotdetails.filterdeg) | ~(shoulddownsample==plotdetails.shoulddownsample)
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
            if cutoffreq>0 & filterdeg>0
                [b,a]=butter(handles.data.samples(selectedsamplenum).filterdeg,handles.data.samples(selectedsamplenum).cutoffreq/(1/si)/2,'low');
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




function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
selectedsamplenum=get(handles.popupmenu1,'Value');
handles.data.samples(selectedsamplenum).cutoffreq=round(str2num(get(hObject,'String')));
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



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
startval=get(handles.slider1,'Value');
endval=get(handles.slider2,'Value');
selectedsamplenum=get(handles.popupmenu1,'Value');
timedifi=handles.data.samples(selectedsamplenum).endtime-handles.data.samples(selectedsamplenum).starttime;
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
timedifi=handles.data.samples(selectedsamplenum).endtime-handles.data.samples(selectedsamplenum).starttime;
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

% Hint: get(hObject,'Value') returns toggle state of checkbox3
