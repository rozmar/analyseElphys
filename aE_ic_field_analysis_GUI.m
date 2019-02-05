function varargout = aE_ic_field_analysis_GUI(varargin)
% AE_IC_FIELD_ANALYSIS_GUI MATLAB code for aE_ic_field_analysis_GUI.fig
%      AE_IC_FIELD_ANALYSIS_GUI, by itself, creates a new AE_IC_FIELD_ANALYSIS_GUI or raises the existing
%      singleton*.
%
%      H = AE_IC_FIELD_ANALYSIS_GUI returns the handle to a new AE_IC_FIELD_ANALYSIS_GUI or the handle to
%      the existing singleton*.
%
%      AE_IC_FIELD_ANALYSIS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AE_IC_FIELD_ANALYSIS_GUI.M with the given input arguments.
%
%      AE_IC_FIELD_ANALYSIS_GUI('Property','Value',...) creates a new AE_IC_FIELD_ANALYSIS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before aE_ic_field_analysis_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to aE_ic_field_analysis_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help aE_ic_field_analysis_GUI

% Last Modified by GUIDE v2.5 04-Jan-2019 17:02:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @aE_ic_field_analysis_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @aE_ic_field_analysis_GUI_OutputFcn, ...
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


% --- Executes just before aE_ic_field_analysis_GUI is made visible.
function aE_ic_field_analysis_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to aE_ic_field_analysis_GUI (see VARARGIN)

% Choose default command line output for aE_ic_field_analysis_GUI
handles.data=varargin{1};
%%
set(handles.popupmenu3,'String',{'trough','peak'});
sourcenames={'intracell'};
a=dir([handles.data.dirs.rawexporteddir,handles.data.xlsdata(handles.data.fieldxlsnum).ID,'.mat']);
if ~isempty(a) 
    sourcenames=[sourcenames,{'field'}];
end
if isfield(handles.data.dirs,'breathingdir')
    a=dir([handles.data.dirs.breathingdir,handles.data.xlsdata(handles.data.fieldxlsnum).ID,'.mat']);
    if ~isempty(a)
        sourcenames=[sourcenames,{'breathing'}];
    end
end
%%
set(handles.popupmenu2,'String',sourcenames);
set(handles.popupmenu4,'String',{'all','sporadic','persistent'});
handles=loadthedata(handles);
handles.output = hObject;
% Update handles structurehandles.data.fielddata=load([handles.data.dirs.rawexporteddir,handles.data.xlsdata(handles.data.fieldxlsnum).ID],'rawdata');
guidata(hObject, handles);

% UIWAIT makes aE_ic_field_analysis_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = aE_ic_field_analysis_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
minval=get(handles.slider1,'Value');
maxval=get(handles.slider3,'Value');
if maxval <= minval
    minval = maxval - 1;
    set(handles.slider1,'Value',minval);
end
medval=mean([minval,maxval]);
set(handles.slider2,'Value',medval);
plotthedata(handles)



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
minval=get(handles.slider1,'Value');
maxval=get(handles.slider3,'Value');
medval=get(handles.slider2,'Value');
radius=diff([minval,maxval]/2);
radius=min([medval-get(handles.slider2,'Min'),radius]);
radius=min([get(handles.slider2,'Max')-medval,radius]);
set(handles.slider1,'Value',medval-radius);
set(handles.slider3,'Value',medval+radius);
plotthedata(handles)
    

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.l rejection via a vomeronasal r
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
minval=get(handles.slider1,'Value');
maxval=get(handles.slider3,'Value');
if maxval <= minval
    maxval = minval + 1;
    set(handles.slider3,'Value',maxval);
end
medval=mean([minval,maxval]);
set(handles.slider2,'Value',medval);
plotthedata(handles)

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
string=get(hObject,'String');
if isnan(str2double(string))
    set(hObject,'String','0')
else
    set(hObject,'String',num2str(str2double(string)))
end

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
axes(handles.axes1)
rect = getrect;
handles.axes1=gca;
startval=max([rect(1),get(handles.slider2,'Min')]);
endval=min([rect(1)+rect(3),get(handles.slider2,'Max')]);
midval=mean([startval,endval]);
set(handles.slider1,'Value',startval);
set(handles.slider2,'Value',midval);
set(handles.slider3,'Value',endval);
guidata(hObject, handles)
plotthedata(handles)
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotthedata(handles)

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
checkfiltershouldbedesigned(handles)

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
checkfiltershouldbedesigned(handles)

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
checkfiltershouldbedesigned(handles)

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


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
checkfiltershouldbedesigned(handles)

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



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
checkfiltershouldbedesigned(handles)

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

function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


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

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filters=get(handles.popupmenu1,'String');
filter=filters{get(handles.popupmenu1,'Value')};
if ~strcmp(filter,'Raw')
    set(handles.pushbutton3,'String','Working...')
    pause(.1)
    [filter,filterdata]=testfilter(handles);
    handles.data.filter=filter;
    handles.data.filterdata=filterdata;
    set(handles.pushbutton3,'Backgroundcolor',[.94,.94,.94]);
    if ~isfield(handles.data,'filter_applied') | (isfield(handles.data,'filterdata_applied') & ~isequal(handles.data.filterdata,handles.data.filterdata_applied))
        set(handles.pushbutton4,'Backgroundcolor','r');
    else
        set(handles.pushbutton4,'Backgroundcolor',[.94,.94,.94]);
    end
    set(handles.pushbutton3,'String','Design filter')
    guidata(hObject,handles);
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tapertime=str2double(get(handles.edit6,'String'));

if isfield(handles.data,'filter') &  (~isfield(handles.data,'filterdata_applied') | (isfield(handles.data,'filterdata_applied') & ~isequal(handles.data.filterdata,handles.data.filterdata_applied)))
    set(handles.pushbutton4,'String','Working...')
    pause(.1)
%     pause
    for sweepnum=1:length(handles.data.fielddata)
        si=handles.data.fielddata(sweepnum).si;
        filteridx=find([handles.data.filter.si]==si);
        y=handles.data.fielddata(sweepnum).y;
        %         [b,a]=tf(handles.data.filter(filteridx).d);
        %         y=filtfilt(b,a,y);
        %             y=filtfilt(handles.data.filter(filteridx).d,y);
        if tapertime > 0
            taperlength = min(round(tapertime/si),length(y));
            y=[y(taperlength:-1:1),y,y(end:-1:length(y)-taperlength+1)];
        end
        for filti=1:length(handles.data.filter(filteridx).d)
            y=handles.data.filter(filteridx).d{filti}.filtfilt(y);
        end
        if tapertime > 0
            y=[y(taperlength+1:end-taperlength)];
        end
        handles.data.fielddata(sweepnum).y_filt=y;
    end
    handles.data.filter_applied=handles.data.filter;
    handles.data.filterdata_applied=handles.data.filterdata;
    set(handles.pushbutton4,'Backgroundcolor',[.94,.94,.94]);
    set(handles.pushbutton4,'String','Apply filter')
end
guidata(hObject,handles);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotthedata(handles,true)

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton6,'String','Working...')
pause(.1)

if ~isfield(handles.data.eventdata,'axonalAP_sporadic')
    handles.data.eventdata=persistent_sort_sporadic_persistent_aAPs(handles.data.eventdata);
end
timebefore=.5;%1;%.5;
timeafter=.5;%1;%.5;
localextremumwin=.1;
types=get(handles.popupmenu3,'String');
selected=get(handles.popupmenu3,'Value');
type=types{selected};
apdata=handles.data.eventdata(strcmp({handles.data.eventdata.type},'AP') & [handles.data.eventdata.stimulated]==0);%& [handles.data.eventdata.baselineval]<-.05
aapdata=handles.data.eventdata(strcmp({handles.data.eventdata.type},'AP') & [handles.data.eventdata.axonalAP] & [handles.data.eventdata.stimulated]==0);%& [handles.data.eventdata.baselineval]<-.05
sapdata=handles.data.eventdata(strcmp({handles.data.eventdata.type},'AP') & [handles.data.eventdata.somaticAP] & [handles.data.eventdata.stimulated]==0);%& [handles.data.eventdata.baselineval]<-.05
epdata=handles.data.eventdata(strcmp({handles.data.eventdata.type},'ep'));
ipdata=handles.data.eventdata(strcmp({handles.data.eventdata.type},'ip'));
FieldData=struct;
NEXT=0;
for fieldsweepnum= 1:length(handles.data.fielddata)
    %     progressbar(fieldsweepnum/ length(field.bridgeddata));
    hossz=length(handles.data.fielddata(fieldsweepnum).y);
    si=handles.data.fielddata(fieldsweepnum).si;
    sweepstart=handles.data.fielddata(fieldsweepnum).realtime;
    sweepend=handles.data.fielddata(fieldsweepnum).realtime+hossz*si;
    %     if isempty(timeborders) | (sweepstart>=timeborders(1) & sweepstart<=timeborders(2)) | (sweepend>=timeborders(1) & sweepend<=timeborders(2)) | (sweepstart<=timeborders(1) & sweepend>=timeborders(2))
    
    
    yfield=handles.data.fielddata(fieldsweepnum).y_filt;
    
%     yfieldorig=yfield;
%     yfieldtoshow=yfield;
    stepback=round(timebefore/si);
    stepforward=round(timeafter/si);
    sweepnum=find([handles.data.bridgeddata.realtime]==handles.data.fielddata(fieldsweepnum).realtime);
    lestep=round(localextremumwin/si);
    minh2=stepback;
    maxh=stepback;
    maxh2=stepback;
    if ~isempty(sweepnum)%% analysis relative to the Field
        if strcmp(type,'peak') % the waveform is upside-down if the recording is intracellular
            yfield=yfield*-1;
        end
        icy=handles.data.bridgeddata(sweepnum).y;
        sweeptime=[0:length(icy)-1]*si+handles.data.bridgeddata(sweepnum).realtime;
        while maxh2<length(yfield)-stepforward-stepback
            maxh=maxh2;
            while maxh<length(yfield)-lestep & max(yfield(maxh:maxh+lestep))>max(yfield(maxh-lestep:maxh))
                maxh=maxh+lestep;
            end
            if maxh+2*lestep<length(yfield)
                minh2=maxh+lestep;
            end
            while minh2<length(yfield)-2*lestep & minh2<length(yfield)-2*stepforward& min(yfield(minh2:minh2+lestep))<min(yfield(minh2-lestep:minh2))
                minh2=minh2+lestep;
            end
            if minh2+2*lestep<length(yfield)
                maxh2=minh2+lestep;
            end
            while maxh2<length(yfield)-2*lestep & max(yfield(maxh2:maxh2+lestep))>max(yfield(maxh2-lestep:maxh2))
                maxh2=maxh2+lestep;
            end
            
            [maxval,sech]=max(yfield(maxh-lestep:maxh));
            sech=sech+maxh-lestep;
            [minval,firsth2]=min(yfield(minh2-lestep:minh2));
            firsth2=firsth2+minh2-lestep;
            [maxval2,sech2]=max(yfield(maxh2:maxh2+lestep));
            sech2=sech2+maxh2;
            
            if  strcmp(type,'peak') % the waveform is upside-down if the recording is intracellular
                maxval=maxval*-1;
                minval=minval*-1;
                maxval2=maxval2*-1;
            end
            
            time=[-stepback:stepforward]*si;
            fieldidx=[-stepback:stepforward]+firsth2;
            if strcmp(type,'peak') % the waveform is upside-down if the recording is intracellular
                fieldy=-1*yfield(fieldidx);
            else
                fieldy=yfield(fieldidx);
            end
%             fieldtoshow=yfieldtoshow(fieldidx);
            intracell=icy(fieldidx);
%             intracell=filtfilt(bic,aic,intracell);
            relativepeakhs=[sech,firsth2,sech2]-firsth2+stepback;
            relativepeakhs(relativepeakhs<1)=1;
            relativepeakhs(relativepeakhs>stepback+stepforward)=stepback+stepforward;
            NEXT=NEXT+1;
            amplitudes=abs(diff(fieldy(relativepeakhs)));
            FieldData(NEXT).time=time';
            FieldData(NEXT).si=mode(diff(time));
            FieldData(NEXT).fieldtodetect=fieldy';
%             FieldData(NEXT).fieldfilt=fieldtoshow';
            FieldData(NEXT).ic=intracell';
            FieldData(NEXT).amplitudes=amplitudes;
            FieldData(NEXT).meanamplitude=nanmean(amplitudes);
            FieldData(NEXT).maxamplitude=nanmax(amplitudes);
            FieldData(NEXT).peakvalue1=fieldy(relativepeakhs(1));
            FieldData(NEXT).troughvalue=fieldy(relativepeakhs(2));
            FieldData(NEXT).peakvalue2=fieldy(relativepeakhs(3));
            FieldData(NEXT).relativepeakhs=relativepeakhs;
            FieldData(NEXT).prevpeak=sech;
            FieldData(NEXT).troughh=firsth2;
            FieldData(NEXT).nextpeak=sech2;
            FieldData(NEXT).fieldsweepnum=fieldsweepnum;
            FieldData(NEXT).sweepnum=sweepnum;
            FieldData(NEXT).troughtime=sweeptime(firsth2);
            
            eventdataa=apdata;
            neededevents=find([eventdataa.maxtime]>=FieldData(NEXT).troughtime-timebefore & [eventdataa.maxtime]<=FieldData(NEXT).troughtime+timeafter);
            FieldData(NEXT).relativeAPtimes=[eventdataa(neededevents).maxtime]-FieldData(NEXT).troughtime;
            if isempty(FieldData(NEXT).relativeAPtimes)
                FieldData(NEXT).apnum=0;
            else
                FieldData(NEXT).apnum=length(FieldData(NEXT).relativeAPtimes);
            end
            
            eventdataa=aapdata;
            neededevents=find([eventdataa.maxtime]>=FieldData(NEXT).troughtime-timebefore & [eventdataa.maxtime]<=FieldData(NEXT).troughtime+timeafter);
            FieldData(NEXT).relativeaxonalAPtimes=[eventdataa(neededevents).maxtime]-FieldData(NEXT).troughtime;
            if isempty(FieldData(NEXT).relativeaxonalAPtimes)
                FieldData(NEXT).axonalapnum=0;
            else
                FieldData(NEXT).axonalapnum=length(FieldData(NEXT).relativeaxonalAPtimes);
            end
            
            if isfield(aapdata,'axonalAP_sporadic')
                eventdataa=aapdata([aapdata.axonalAP_sporadic]==1);
                neededevents=find([eventdataa.maxtime]>=FieldData(NEXT).troughtime-timebefore & [eventdataa.maxtime]<=FieldData(NEXT).troughtime+timeafter);
                FieldData(NEXT).relativeaxonalAPtimes_sporadic=[eventdataa(neededevents).maxtime]-FieldData(NEXT).troughtime;
                if isempty(FieldData(NEXT).relativeaxonalAPtimes_sporadic)
                    FieldData(NEXT).axonalapnum_sporadic=0;
                else
                    FieldData(NEXT).axonalapnum=length(FieldData(NEXT).relativeaxonalAPtimes_sporadic);
                end
                
                eventdataa=aapdata([aapdata.axonalAP_persistent]==1);
                neededevents=find([eventdataa.maxtime]>=FieldData(NEXT).troughtime-timebefore & [eventdataa.maxtime]<=FieldData(NEXT).troughtime+timeafter);
                FieldData(NEXT).relativeaxonalAPtimes_persistent=[eventdataa(neededevents).maxtime]-FieldData(NEXT).troughtime;
                if isempty(FieldData(NEXT).relativeaxonalAPtimes_persistent)
                    FieldData(NEXT).axonalapnum_persistent=0;
                else
                    FieldData(NEXT).axonalapnum=length(FieldData(NEXT).relativeaxonalAPtimes_sporadic);
                end
            end
            
            eventdataa=epdata;
            neededevents=find([eventdataa.maxtime]>=FieldData(NEXT).troughtime-timebefore & [eventdataa.maxtime]<=FieldData(NEXT).troughtime+timeafter);
            FieldData(NEXT).relativeeptimes=[eventdataa(neededevents).maxtime]-FieldData(NEXT).troughtime;
            if isempty(FieldData(NEXT).relativeeptimes)
                FieldData(NEXT).epnum=0;
            else
                FieldData(NEXT).epnum=length(FieldData(NEXT).relativeeptimes);
            end
            
            eventdataa=ipdata;
            neededevents=find([eventdataa.maxtime]>=FieldData(NEXT).troughtime-timebefore & [eventdataa.maxtime]<=FieldData(NEXT).troughtime+timeafter);
            FieldData(NEXT).relativeiptimes=[eventdataa(neededevents).maxtime]-FieldData(NEXT).troughtime;
            if isempty(FieldData(NEXT).relativeiptimes)
                FieldData(NEXT).ipnum=0;
            else
                FieldData(NEXT).ipnum=length(FieldData(NEXT).relativeiptimes);
            end
            FieldData(NEXT).eventnum=FieldData(NEXT).ipnum+FieldData(NEXT).epnum+FieldData(NEXT).apnum;
            
            %             figure(3)
            %             clf
            %             subplot(3,1,1)
            %             plot(time,fieldy,'k-')
            %             hold on;
            %             plot(time(relativepeakhs),fieldy(relativepeakhs),'ro')
            %             axis tight
            %             subplot(3,1,2)
            %             plot(time,fieldtoshow,'k-')
            %             hold on;
            %             plot(time(relativepeakhs),fieldtoshow(relativepeakhs),'ro')
            %             axis tight
            %             subplot(3,1,3)
            %             plot(time,intracell,'k-')
            %             axis tight
            %             pause
        end
        %         if NEXT>0
        %         figure(2)
        %         clf
        %         needed=[FieldData.fieldsweepnum]==fieldsweepnum;
        %          subplot(3,1,1)
        %             plot(sweeptime,yfieldorig,'k-')
        %             hold on;
        %             plot(sweeptime([FieldData(needed).troughh]),yfieldorig([FieldData(needed).troughh]),'ro')
        %             axis tight
        %             subplot(3,1,2)
        %             plot(sweeptime,yfieldtoshow,'k-')
        %             hold on;
        %             plot(sweeptime([FieldData(needed).troughh]),yfieldtoshow([FieldData(needed).troughh]),'ro')
        %             axis tight
        %             axis tight
        %             subplot(3,1,3)
        %             plot(sweeptime,icy,'k-')
        %             axis tight
        % %             pause
        %         end
    end
    %     end
end

if isfield(handles.data,'BrainStateData')
    statestodo=[unique({handles.data.BrainStateData.name}),'All'];
else
    statestodo={'All'};
end

for i=1:length(FieldData)
    FieldData(i).medicV=median(FieldData(i).ic);
    if isfield(handles.data,'BrainStateData')
        idx=find(FieldData(i).troughtime>[handles.data.BrainStateData.starttime] & FieldData(i).troughtime<[handles.data.BrainStateData.endtime]);
        if ~isempty(idx)
            FieldData(i).brainstatename=handles.data.BrainStateData(idx).name;
        else
            FieldData(i).brainstatename='none';
        end
    end
end
%% EZT LECSEKKOLNI!!!

for statei=1:length(statestodo)
    statename=statestodo{statei};
    %%
%     FieldData=handles.data.FieldData;

    if ~strcmp(statename,'All')
        idxnow=strcmp({FieldData.brainstatename},statename);
    elseif length(statestodo)==1
        idxnow=ones(size(FieldData));
    else
        idxnow=zeros(size(FieldData));
    end
    %%
    FieldDatanow=FieldData(find(idxnow));
    
    %%
    %     FieldData([FieldData.eventnum]==0)=[];
    timestep=.05;
%     bins=[-timebefore:timestep:timeafter];
%     bins_corr=[-timebefore_corr:timestep:timeafter_corr];
    timewindowforfieldamplitude=10;
    [FieldData(:).medianamplitude] = deal(NaN);
    for i = 1:length(FieldDatanow)
        needed=FieldDatanow(i).troughtime-timewindowforfieldamplitude/2<=[FieldDatanow.troughtime] & FieldDatanow(i).troughtime+timewindowforfieldamplitude/2>=[FieldDatanow.troughtime];
        FieldDatanow(i).medianamplitude=median([FieldDatanow(needed).maxamplitude]);
    end
    FieldData=[FieldData,FieldDatanow];
    FieldData(find(idxnow))=[];
    disp('field analysis finished')
end
%%
set(handles.pushbutton6,'String','Detect troughs')
disp('lol')
handles.data.FieldData=FieldData;
guidata(hObject,handles)


function plotthedata(handles,testthefilter)
%%
if nargin<2
    testthefilter=false;
end
% testfilter(handles)
filters=get(handles.popupmenu1,'String');
filter=filters{get(handles.popupmenu1,'Value')};
highpass_stop=str2double(get(handles.edit2,'String'));
highpass_pass=str2double(get(handles.edit3,'String'));
lowpass_stop=str2double(get(handles.edit4,'String'));
lowpass_pass=str2double(get(handles.edit5,'String'));
starttime = get(handles.slider1,'Value');
endtime = get(handles.slider3,'Value');
startsweep = find([handles.data.fielddata.realtime]<=starttime,1,'last');
tapertime=str2double(get(handles.edit6,'String'));
if isempty(startsweep)
    startsweep = 1;
end
endsweep = find([handles.data.fielddata.realtime]>=endtime,1,'first');
if isempty(endsweep)
    endsweep = length(handles.data.fielddata);
end
maxpointnum=str2double(get(handles.edit1,'String'));
toplot = struct;
hossz=length([handles.data.fielddata(startsweep:endsweep).y]);
ratio=round(hossz/maxpointnum);
if isinf(ratio)
    ratio=1;
end
i=0;
for sweepnum=startsweep:endsweep
    i=i+1;
    si=handles.data.fielddata(sweepnum).si;
    realtime=handles.data.fielddata(sweepnum).realtime;
    y=handles.data.fielddata(sweepnum).y;
    time=([1:length(y)]-1)*si+realtime;
    
    if ~strcmp(filter,'Raw') & testthefilter
        filteridx=find([handles.data.filter.si]==si);
        if tapertime > 0
            taperlength = min(round(tapertime/si),length(y));
            y=[y(taperlength:-1:1),y,y(end:-1:length(y)-taperlength+1)];
        end
        for filti=1:length(handles.data.filter(filteridx).d)
            y=handles.data.filter(filteridx).d{filti}.filtfilt(y);
        end
        if tapertime > 0
            y=[y(taperlength+1:end-taperlength)];
        end
    elseif ~strcmp(filter,'Raw') & isfield(handles.data.fielddata(sweepnum),'y_filt')
        y=handles.data.fielddata(sweepnum).y_filt;
    end
    
    
    
    %     end
    if maxpointnum>0 & ~isinf(maxpointnum) %downsampling
        y=downsample(y,ratio);
        time=downsample(time,ratio);
    end
    toplot(i).y=y;
    toplot(i).time=time;
end



axes(handles.axes1)
cla
hold on
for sweep=1:length(toplot)
    plot(toplot(sweep).time,toplot(sweep).y,'k-')
    if isfield(handles.data,'FieldData')
        sweeptime=toplot(sweep).time(1);
        fieldsweepnum=find([handles.data.fielddata.realtime]==sweeptime);
        FieldData=handles.data.FieldData([handles.data.FieldData.fieldsweepnum]==fieldsweepnum);
        plot(toplot(sweep).time(round([FieldData.troughh]/ratio)),toplot(sweep).y(round([FieldData.troughh]/ratio)),'ro');
    end
end
ylim([nanmin([toplot.y]),nanmax([toplot.y])])
xlim([starttime,endtime])
aAPidxs=find([handles.data.eventdata.axonalAP] & [handles.data.eventdata.maxtime]<=endtime & [handles.data.eventdata.maxtime]>=starttime);
plot([handles.data.eventdata(aAPidxs).maxtime],ones(size(aAPidxs))*nanmax([toplot.y]),'b^')
% if isfield(handles.data,'FieldData')
%     timetotroughs=[];
%     for aAPi=1:length(aAPidxs)
%         %%
%         maxtime=handles.data.eventdata(aAPidxs(aAPi)).maxtime;
%         [~,minidx]=min(abs([handles.data.FieldData.troughtime]-maxtime));
%         timetotrough=maxtime-handles.data.FieldData(minidx).troughtime;
%         timetotroughs=[timetotroughs,timetotrough];
%         
% %         nearesttrough
%     end
%     plot([handles.data.eventdata(aAPidxs).maxtime],timetotroughs,'b^')
% end
axis tight
handles.axes1=gca;
hObject=findall(gcf,'Name','aE_ic_field_analysis_GUI');
guidata(hObject,handles);

function [filters,filterdata]=testfilter(handles)

sis=unique([handles.data.fielddata.si]);
filters=get(handles.popupmenu1,'String');
filter=filters{get(handles.popupmenu1,'Value')};
highpass_stop=str2double(get(handles.edit2,'String'));
highpass_pass=str2double(get(handles.edit3,'String'));
lowpass_stop=str2double(get(handles.edit4,'String'));
lowpass_pass=str2double(get(handles.edit5,'String'));

filterdata.filter=filter;
filterdata.highpass_stop=highpass_stop;
filterdata.highpass_pass=highpass_pass;
filterdata.lowpass_stop=lowpass_stop;
filterdata.lowpass_pass=lowpass_pass;

filters=struct;
for sii=1:length(sis)
    si=sis(sii);
    if ~strcmp(filter,'Raw') & ((highpass_stop>0 & highpass_pass > 0) | (lowpass_stop>0 & lowpass_pass > 0))
        if any(strcmp({'equiripple','kaiserwin'},filter))%FIR filters
            filtertype='fir';
        elseif any(strcmp(filter,{'butter','cheby1','cheby2','ellip'}))%IIF filters
            filtertype='iir';
        end
        d={};
        if highpass_stop == 0 | highpass_pass == 0% lowpass filter
            d{1} = designfilt(['lowpass',filtertype],'PassbandFrequency',lowpass_pass,'StopbandFrequency',lowpass_stop,'SampleRate',1/si,'DesignMethod',filter);
        elseif lowpass_stop == 0 | lowpass_pass == 0
            d{1} = designfilt(['highpass',filtertype],'PassbandFrequency',highpass_pass,'StopbandFrequency',highpass_stop,'SampleRate',1/si,'DesignMethod',filter);
        else %bandpass filter
%             d{} = designfilt(['bandpass',filtertype],'StopbandFrequency1',highpass_stop, 'PassbandFrequency1',highpass_pass, 'PassbandFrequency2',lowpass_pass, 'StopbandFrequency2',lowpass_stop,'SampleRate',1/si,'DesignMethod',filter);
            d{1} = designfilt(['lowpass',filtertype],'PassbandFrequency',lowpass_pass,'StopbandFrequency',lowpass_stop,'SampleRate',1/si,'DesignMethod',filter);
            d{2} = designfilt(['highpass',filtertype],'PassbandFrequency',highpass_pass,'StopbandFrequency',highpass_stop,'SampleRate',1/si,'DesignMethod',filter);
        end
        for i=1:length(d)
            [h(:,i),w(:,i)]=freqz(d{i},1000000,1/si);
        end
        if size(h,2)>1
            h=h(:,1).*h(:,2);
        end
        figure(30+sii)
        clf
        semilogx(w,20*log10(abs(h.^2)))
        ylim([-60 10])
        xlabel('Frequency(Hz)')
        ylabel('Amplitude response (dB)')
        set(gca, 'XTick', [10.^-1 10.^0 10.^1 10.^2 10.^3 10.^4 10.^5 10.^6])
%         fvtool(d)
        filters(sii).si=si;
        filters(sii).d=d;
%         info=[{'filter properties'};{d{1}.info}; {['order: ',num2str(d{1}.filtord)]}];
%         set(handles.text8,'String',info)
        
    end
end

function checkfiltershouldbedesigned(handles)
filters=get(handles.popupmenu1,'String');
filter=filters{get(handles.popupmenu1,'Value')};
highpass_stop=str2double(get(handles.edit2,'String'));
highpass_pass=str2double(get(handles.edit3,'String'));
lowpass_stop=str2double(get(handles.edit4,'String'));
lowpass_pass=str2double(get(handles.edit5,'String'));

filterdata.filter=filter;
filterdata.highpass_stop=highpass_stop;
filterdata.highpass_pass=highpass_pass;
filterdata.lowpass_stop=lowpass_stop;
filterdata.lowpass_pass=lowpass_pass;

if ~isfield(handles.data,'filterdata') | ~isequal(filterdata,handles.data.filterdata)
    set(handles.pushbutton3,'Backgroundcolor','r');
else
    set(handles.pushbutton3,'Backgroundcolor',[.94,.94,.94]);
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% plotting 
aAPtypes=get(handles.popupmenu4,'String');
aAPtype=aAPtypes{get(handles.popupmenu4,'Value')};
starttime = get(handles.slider1,'Value');
endtime = get(handles.slider3,'Value');

valtozok_fieldplot.timebefore=.5;%1;%.5;
valtozok_fieldplot.timestep=.025;
valtozok_fieldplot.timeafter=.5;%;1;%.5;
% valtozok_fieldplot.ICylim
% valtozok_fieldplot.Fieldylim
valtozok_fieldplot.highlightaAP=1;
valtozok_fieldplot.highlightaAP_timeback=.002;
valtozok_fieldplot.highlightaAP_timeforward=.02;

if isfield(handles.data,'BrainStateData')
    statestodo=[unique({handles.data.BrainStateData.name}),'All'];
else
    statestodo={'All'};
end
for statei=1:length(statestodo)
    statename=statestodo{statei};
    %%
    FieldData=handles.data.FieldData;
    FieldData=FieldData([FieldData.troughtime]>starttime & [FieldData.troughtime]<endtime);
    if ~strcmp(statename,'All')
        FieldData=FieldData(strcmp({FieldData.brainstatename},statename));
    end
    valtozok_fieldplot.figurenum=10+statei;
    if ~strcmp(aAPtype,'all')
        valtozok_fieldplot.aAPtype=aAPtype;
    end
    persistent_generatefigures_2017_fieldanal(FieldData,valtozok_fieldplot)
    subplot(4,1,1);
    title(statename)
end

function handles=loadthedata(handles)
sourcenames=get(handles.popupmenu2,'String');
sourceidx=get(handles.popupmenu2,'Value');
source=sourcenames{sourceidx};
if strcmp(source,'field')
    load([handles.data.dirs.rawexporteddir,handles.data.xlsdata(handles.data.fieldxlsnum).ID],'rawdata');
    handles.data.fielddata=rawdata;
elseif strcmp(source,'breathing')
    load([handles.data.dirs.breathingdir,handles.data.xlsdata(handles.data.fieldxlsnum).ID],'rawdata');
    handles.data.fielddata=rawdata;
elseif strcmp(source,'intracell')
    load([handles.data.dirs.bridgeddir,handles.data.xlsdata(handles.data.icxlsnum).ID],'bridgeddata');
    handles.data.fielddata=bridgeddata;
end
if isfield(handles.data,'FieldData')
    handles.data=rmfield(handles.data,'FieldData');
end
if isfield(handles.data,'filterdata')
    handles.data=rmfield(handles.data,'filterdata');
end
set(handles.slider1,'Max',handles.data.endtime);
set(handles.slider1,'Min',handles.data.starttime);
set(handles.slider1,'Value',handles.data.starttime);
set(handles.slider2,'Max',handles.data.endtime);
set(handles.slider2,'Min',handles.data.starttime);
set(handles.slider2,'Value',mean([handles.data.endtime,handles.data.starttime]));
set(handles.slider3,'Max',handles.data.endtime);
set(handles.slider3,'Min',handles.data.starttime);
set(handles.slider3,'Value',handles.data.endtime);
set(handles.popupmenu1,'String',{'Raw','equiripple','kaiserwin','butter','cheby1','cheby2','ellip'})
axes(handles.axes1)
xlim([handles.data.starttime,handles.data.endtime])


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
handles=loadthedata(handles);
guidata(hObject, handles);

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


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


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
