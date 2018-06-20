function varargout = aE_videoanalyzer_ROIselector(varargin)
% aE_VIDEOANALYZER_ROISELECTOR MATLAB code for aE_videoanalyzer_ROIselector.fig
% USAGE: aE_videoanalyzer_ROIselector(dirs, xlsdata); after starting
% analysElphys_main.m
%      AE_VIDEOANALYZER_ROISELECTOR, by itself, creates a new AE_VIDEOANALYZER_ROISELECTOR or raises the existing
%      singleton*.
%
%      H = AE_VIDEOANALYZER_ROISELECTOR returns the handle to a new AE_VIDEOANALYZER_ROISELECTOR or the handle to
%      the existing singleton*.
%
%      AE_VIDEOANALYZER_ROISELECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AE_VIDEOANALYZER_ROISELECTOR.M with the given input arguments.
%
%      AE_VIDEOANALYZER_ROISELECTOR('Property','Value',...) creates a new AE_VIDEOANALYZER_ROISELECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before aE_videoanalyzer_ROIselector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to aE_videoanalyzer_ROIselector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help aE_videoanalyzer_ROIselector

% Last Modified by GUIDE v2.5 19-Jun-2018 17:29:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @aE_videoanalyzer_ROIselector_OpeningFcn, ...
                   'gui_OutputFcn',  @aE_videoanalyzer_ROIselector_OutputFcn, ...
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


% --- Executes just before aE_videoanalyzer_ROIselector is made visible.
function aE_videoanalyzer_ROIselector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to aE_videoanalyzer_ROIselector (see VARARGIN)

% Choose default command line output for aE_videoanalyzer_ROIselector
%%

% timernamesstodel={'videoplay'};
% timersrunning=timerfind;
% timerstodel=zeros(size(timersrunning));
% for i=1:length(timersrunning)
%     a=timersrunning(i).TimerFcn;
%     if any(strcmp(func2str(a),timernamesstodel))
%         timerstodel(i)=1;
%     else
%         timerstodel(i)=0;
%     end
% end
% if ~isempty(timersrunning)
% stop(timersrunning(find(timerstodel)));
% delete(timersrunning(find(timerstodel)));
% end

handles.t = timer('BusyMode', 'drop', 'ExecutionMode', 'fixedRate', 'Period', .05);%0.05
set(handles.t, 'TimerFcn', {@videoplay, hObject});

data.ROIS={'Snout','Nose','Mounth','Body','Ear','Hands','Whisker','Eye','Background'};
data.ROIcolors={'red','green','blue','cyan','magenta','yellow','white','black','blue'};
data.ROIcolors_short={'r','g','b','c','m','y','w','k','b'};

data.dirs=varargin{1};
data.xlsdata=varargin{2};

rawvideodir=[data.dirs.videodir,'movement/'];
data.dirs.rawvideodir=rawvideodir;
files=dir(rawvideodir);
files([files.isdir])=[];
data.files=files;

set(handles.popupmenu2,'String',[{'Not selected'},{files.name}]);
set(handles.popupmenu3,'String',[{'Not selected'}]);
set(handles.popupmenu1,'String',[data.ROIS]);

set(handles.listbox1,'String',[0:10]);
set(handles.listbox1,'Min',0,'Max',2);
handles.data=data;


%%

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes aE_videoanalyzer_ROIselector wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = aE_videoanalyzer_ROIselector_OutputFcn(hObject, eventdata, handles) 
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
handles=updateGUI(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.pushbutton1,'string'),'Play')
    start(handles.t);
    set(handles.pushbutton1,'string','Stop')
    
else
    stop(handles.t);
    set(handles.pushbutton1,'string','Play')
end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
plotpca(handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.InspectTraces
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
%%
if get(handles.popupmenu2,'Value')-1>0
 handles=loadthevideo(handles);
 handles=updateGUI(handles);
 guidata(hObject,handles);
end





% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgrupdateGUI(handles)oundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%
selectedfile=get(handles.popupmenu2,'Value')-1;
selectedROI=get(handles.popupmenu1,'Value');
if selectedROI>0
    h=imfreehand;
    handles.data.ROIdata(selectedROI).contour=h.getPosition;
    handles.data.ROIdata(selectedROI).mask=h.createMask;
    h.delete;
    ROIdata=handles.data.ROIdata;
    save([handles.data.dirs.videodir,'ROIs/',handles.data.files(selectedfile).name],'ROIdata');
    %this is not done at all
end
handles=updateGUI(handles);
guidata(hObject,handles);


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3

handles=updateGUI(handles);
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


function handles=loadthevideo(handles)
%%
selectedfile=get(handles.popupmenu2,'Value')-1;
handles.data.VIDEOdata=load([handles.data.dirs.rawvideodir,handles.data.files(selectedfile).name]);
set(handles.popupmenu3,'Value',1);
set(handles.popupmenu3,'String',1:length(handles.data.VIDEOdata.rawvideo));
a=dir([handles.data.dirs.videodir,'ROIs/',handles.data.files(selectedfile).name]);
if isempty(a);
    ROIdata=struct;
    for i=1:length(handles.data.ROIS)
        ROIdata(i).ROIname=handles.data.ROIS{i};
        ROIdata(i).ROIcolor=handles.data.ROIcolors{i};
        ROIdata(i).ROIcolor_short=handles.data.ROIcolors_short{i};
        %     data.ROIdata(i).handle=NaN;
        ROIdata(i).properties=struct;
    end
    
else
    load([handles.data.dirs.videodir,'ROIs/',handles.data.files(selectedfile).name]);
end
handles.data.ROIdata=ROIdata;
axes(handles.axes1)
pause(3)
cla;
handles.axes1=gca;
if isfield(handles.data,'imagedata')
    handles.data=rmfield(handles.data,'imagedata');
end
if isfield(handles,'ROIs')
    handles=rmfield(handles,'ROIs');
end


function handles=updateGUI(handles)
%%
selectedvideo=get(handles.popupmenu3,'Value');
pointintime=get(handles.slider1,'Value');
timescale=handles.data.VIDEOdata.videodata(selectedvideo).time;
frameidx=round(pointintime*length(timescale));
if frameidx==0
    frameidx=1;
end
axes(handles.axes1);
if isfield(handles.data,'imagedata');
    set(handles.data.imagedata,'CData',handles.data.VIDEOdata.rawvideo(selectedvideo).vid(:,:,frameidx));
else
    handles.data.imagedata=imagesc(handles.data.VIDEOdata.rawvideo(selectedvideo).vid(:,:,frameidx));
	colormap gray
    hold on;
end
if isfield(handles.data.ROIdata,'mask')
    if ~isfield(handles,'ROIs')
        handles.ROIs=struct;
    end
    for ROIi=1:length(handles.data.ROIdata)
        if ~isempty(handles.data.ROIdata(ROIi).contour)
            if isfield(handles.ROIs,['contour',num2str(ROIi)])
                handles.ROIs.(['contour',num2str(ROIi)]).XData=handles.data.ROIdata(ROIi).contour(:,1);
                handles.ROIs.(['contour',num2str(ROIi)]).YData=handles.data.ROIdata(ROIi).contour(:,2);
            else
                handles.ROIs.(['contour',num2str(ROIi)])=plot(handles.data.ROIdata(ROIi).contour(:,1),handles.data.ROIdata(ROIi).contour(:,2),[handles.data.ROIdata(ROIi).ROIcolor_short,'-'],'LineWidth',2);
                %%
                x=median(handles.data.ROIdata(ROIi).contour(:,1));
                y=median(handles.data.ROIdata(ROIi).contour(:,2));
                handles.ROIs.(['text',num2str(ROIi)])=text(x,y,handles.data.ROIdata(ROIi).ROIname,'Color',handles.data.ROIdata(ROIi).ROIcolor);
            end
        end
    end
end
handles.axes1=gca;

function videoplay(hObject, eventdata, parent_GUI)
handles = guidata(parent_GUI);
selectedvideo=get(handles.popupmenu3,'Value');
displayperiod=get(handles.t,'Period');
videolength=length(handles.data.VIDEOdata.videodata(selectedvideo).time);
videoSI=mode(diff(handles.data.VIDEOdata.videodata(selectedvideo).time));
stepsize=ceil(displayperiod/videoSI);
pointintime=get(handles.slider1,'Value');

timescale=handles.data.VIDEOdata.videodata(selectedvideo).time;

frameidx=round(pointintime*length(timescale));

if frameidx<videolength-stepsize;
    set(handles.slider1,'Value',pointintime+stepsize/videolength);
    handles=updateGUI(handles);
    guidata(parent_GUI, handles);
end

function plotpca(handles)
if isfield(handles.data.ROIdata,'PCA')
    selectedROI=get(handles.popupmenu1,'Value');
    selectedvideo=get(handles.popupmenu3,'Value');
    selectedPC=get(handles.listbox1,'Value')-1;
    axes(handles.axes2);
    cla
    if selectedPC==0
         plot(handles.data.ROIdata(selectedROI).PCA(selectedvideo).bestPC)
%         disp('a')
    else
        plot(handles.data.ROIdata(selectedROI).PCA(selectedvideo).PC(:,selectedPC))
    end
    handles.axes2=gca;
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
plotpca(handles)

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
