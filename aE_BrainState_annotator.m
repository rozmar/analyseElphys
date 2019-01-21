function varargout = aE_BrainState_annotator(varargin)
% AE_BRAINSTATE_ANNOTATOR MATLAB code for aE_BrainState_annotator.fig
%      AE_BRAINSTATE_ANNOTATOR, by itself, creates a new AE_BRAINSTATE_ANNOTATOR or raises the existing
%      singleton*.
%
%      H = AE_BRAINSTATE_ANNOTATOR returns the handle to a new AE_BRAINSTATE_ANNOTATOR or the handle to
%      the existing singleton*.
%
%      AE_BRAINSTATE_ANNOTATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AE_BRAINSTATE_ANNOTATOR.M with the given input arguments.
%
%      AE_BRAINSTATE_ANNOTATOR('Property','Value',...) creates a new AE_BRAINSTATE_ANNOTATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before aE_BrainState_annotator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to aE_BrainState_annotator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help aE_BrainState_annotator

% Last Modified by GUIDE v2.5 08-Oct-2018 19:15:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @aE_BrainState_annotator_OpeningFcn, ...
                   'gui_OutputFcn',  @aE_BrainState_annotator_OutputFcn, ...
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


% --- Executes just before aE_BrainState_annotator is made visible.
function aE_BrainState_annotator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to aE_BrainState_annotator (see VARARGIN)

% Choose default command line output for aE_BrainState_annotator

handles.data=varargin{1};
handles.data.BrainStateProps=varargin{2};
set(handles.popupmenu1,'String',handles.data.BrainStateProps.Statevalues(1:end-1));
set(handles.popupmenu2,'String',[handles.data.movementnames,{'pupilsize','duration'}]);
set(handles.popupmenu3,'String',{'<','>','='});
set(handles.edit1,'String','0')
set(handles.listbox2,'String',[handles.data.movementnames,{'pupilsize'}]);
if ~isfield(handles.data,'BrainStateRules')
    handles.data.BrainStateRules=struct;
    set(handles.listbox1,'String','none');
elseif isempty(fieldnames(handles.data.BrainStateRules))
    set(handles.listbox1,'String','none');
else
    set(handles.listbox1,'String',{handles.data.BrainStateRules.string2display});
end
set(handles.listbox2,'Max',20,'Min',1);


%BEFEJEZNI%!!! - új sampling ratere mindent..
%%
si=.2;
BrainStateTime=floor(handles.data.starttime):si:ceil(handles.data.endtime);
BrainStateVariables=struct;
BrainStateVariables.BrainStateTime=BrainStateTime;
varnames=handles.data.movementnames;
if isfield(handles.data,'pupildata')
    varnames=[varnames,{'pupilsize'}];
end
starttime=handles.data.starttime;
endtime=handles.data.endtime;
for vari=1:length(varnames)
    varname=varnames{vari};
    disp(['resampling ',varname,' ...'])
    if strcmp(varname,'pupilsize')
        time=handles.data.pupildata.time;
        y=handles.data.pupildata.diameter;
    else
        time=handles.data.movementdata.time;
        y=handles.data.movementdata.(varname);
    end
    BrainStateVariables.(varname)=zeros(size(BrainStateTime));
    for i=1:length(BrainStateTime)
        timenow=BrainStateTime(i);
        idx=time>timenow-si/2 & time<=timenow+si/2;
        BrainStateVariables.(varname)(i)=nanmedian(y(idx));
    end
    disp(['done ',varname,' ...'])
end
handles.data.BrainStateVariables=BrainStateVariables;
%%
handles.output = struct;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes aE_BrainState_annotator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = aE_BrainState_annotator_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
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
BrainStateRules=handles.data.BrainStateRules;
if isempty(fieldnames(BrainStateRules))
    NEXT=1;
else
    NEXT=length(BrainStateRules)+1;
end
states=handles.data.BrainStateProps.Statevalues(1:end-1);
BrainStateRules(NEXT).selectedstate=states{get(handles.popupmenu1,'Value')};
variables=get(handles.popupmenu2,'String');
BrainStateRules(NEXT).variable=variables{get(handles.popupmenu2,'Value')};
functions= get(handles.popupmenu3,'String');
BrainStateRules(NEXT).function_now= functions{get(handles.popupmenu3,'Value')};
BrainStateRules(NEXT).value=str2num(get(handles.edit1,'String'));
BrainStateRules(NEXT).string2display=[BrainStateRules(NEXT).selectedstate,': ', BrainStateRules(NEXT).variable, ' ',BrainStateRules(NEXT).function_now,' ', num2str(BrainStateRules(NEXT).value)];

[~,idx]=sort({BrainStateRules.string2display});
BrainStateRules=BrainStateRules(idx);
handles.data.BrainStateRules=BrainStateRules;
guidata(hObject, handles)
applyrules(handles)
% plotandupdate(handles)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(fieldnames(handles.data.BrainStateRules))
    selectedrule=get(handles.listbox1,'Value');
    handles.data.BrainStateRules(selectedrule)=[];
end
set(handles.listbox1,'Value',1);
applyrules(handles);
% plotandupdate(handles);

% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2
plotandupdate(handles)

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

function plotandupdate(handles)
selectedstuff=get(handles.listbox2,'Value');
labels=get(handles.listbox2,'String');
labels=labels(selectedstuff);
axes(handles.axes1)
cla
hold all
minvals=[];
maxvals=[];
for i=1:length(labels)
    labelnow=labels{i};
    x=handles.data.BrainStateVariables.BrainStateTime;
    y=handles.data.BrainStateVariables.(labelnow);
    plot(x,y)
    minvals=[minvals,nanmin(y)];
    maxvals=[maxvals,nanmax(y)];
end
minval=min(minvals);
maxval=max(maxvals);
BrainStateProps=handles.data.BrainStateProps;
BrainStateData=handles.data.BrainStateData;
if ~isempty(fieldnames(BrainStateData))
    ylimits=[minval,maxval];
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
handles.axes1=gca;
hObject=findall(gcf,'Name','aE_BrainState_annotator');
guidata(hObject,handles);


function applyrules(handles)
BrainStateRules=handles.data.BrainStateRules;
BrainStateData=struct;
if isempty(fieldnames(BrainStateRules))
    set(handles.listbox1,'String','None')
else
    set(handles.listbox1,'String',{BrainStateRules.string2display})
    [stagestodetect,stagenums,stageidxs]=unique({BrainStateRules.selectedstate});
    for stagei=1:length(stagestodetect)
        stageidxsnow=find(stageidxs==stagei);
        doduration=false;
        durationidx=[];
        stagevals=true(size(handles.data.BrainStateVariables.BrainStateTime));
        for rulei=1:length(stageidxsnow)
            if strcmp(BrainStateRules(stageidxsnow(rulei)).variable,'duration')
                doduration=true;
                durationidx=[durationidx,stageidxsnow(rulei)];
            else
%                 x=handles.data.BrainStateVariables.BrainStateTime;
                y=handles.data.BrainStateVariables.(BrainStateRules(stageidxsnow(rulei)).variable);
                value=BrainStateRules(stageidxsnow(rulei)).value;
                if strcmp(BrainStateRules(stageidxsnow(rulei)).function_now,'<')
                    stagevals=stagevals & y <= value;
                elseif strcmp(BrainStateRules(stageidxsnow(rulei)).function_now,'>')
                    stagevals=stagevals & y >= value;
                elseif strcmp(BrainStateRules(stageidxsnow(rulei)).function_now,'=')
                    stagevals=stagevals & y == value;
                end
            end
        end
        if doduration==true
            for idxi=1:length(durationidx)
                si = mode(diff(handles.data.BrainStateVariables.BrainStateTime));
                minlength = round(BrainStateRules(durationidx(idxi)).value/si);
                stagevals=bwareaopen(stagevals,minlength);
            end
        end
        [L,num]=bwlabel(stagevals);
        if num > 0
            for i=1:num
                if isempty(fieldnames(BrainStateData))
                    NEXT=1;
                else
                    NEXT=length(BrainStateData)+1;
                end
                startidx=find(L==i,1,'first');
                endidx=find(L==i,1,'last');
                BrainStateData(NEXT).name=BrainStateRules(stageidxsnow(rulei)).selectedstate;
                BrainStateData(NEXT).starttime=handles.data.BrainStateVariables.BrainStateTime(startidx);
                BrainStateData(NEXT).endtime=handles.data.BrainStateVariables.BrainStateTime(endidx);
                
            end
        end
    end
    handles.data.BrainStateData=BrainStateData;
end
hObject=findall(gcf,'Name','aE_BrainState_annotator');
output=struct;
output.BrainStateData=handles.data.BrainStateData;
output.BrainStateVariables=handles.data.BrainStateVariables;
handles.output=output;
guidata(hObject,handles);
plotandupdate(handles);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles_orig=handles;
BrainStateData=handles_orig.data.BrainStateData;
BrainStateRules=handles_orig.data.BrainStateRules;
%%

h=findobj('Name','aE_InspectTraces');
handles = guidata(h);
selectedsamplenum=handles.data.selectedsamplenum;
handles.data.samples(selectedsamplenum).BrainStateData=BrainStateData;
handles.data.samples(selectedsamplenum).BrainStateRules=BrainStateRules;
guidata(h,handles);

save([handles.data.dirs.brainstatedir,handles.data.IDs{handles.data.samples(selectedsamplenum).selectedID}],'BrainStateData','BrainStateRules');
%%
handles=handles_orig;
disp('lol')