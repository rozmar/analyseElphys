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

% Last Modified by GUIDE v2.5 07-Feb-2019 17:48:52

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
handles.data.eventdata=varargin{3};
handles.data.eventdata_orig=varargin{3};
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
set(handles.popupmenu6,'String',{'=','~=','<','>'});
set(handles.popupmenu5,'String',fieldnames(handles.data.eventdata));
set(handles.popupmenu4,'String',[1:4]);
set(handles.popupmenu3,'String',fieldnames(handles.data.eventdata));
set(handles.listbox1,'Max',inf,'Min',0);
set(handles.popupmenu1,'String',{'t','v'});
set(handles.popupmenu2,'String',{'v','dv','ddv','dddv'});
handles = updateGUI(handles);
%%
% Update handles structure
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
            if strcmp(handles.data.rules(i).comparefunc,'=')
                kell=kell& [handles.data.eventdata_orig.(handles.data.rules(i).selectedeventdatafield)] == handles.data.rules(i).value;
            elseif strcmp(handles.data.rules(i).comparefunc,'~=')
                kell=kell& [handles.data.eventdata_orig.(handles.data.rules(i).selectedeventdatafield)] ~= handles.data.rules(i).value;
            elseif strcmp(handles.data.rules(i).comparefunc,'>')
                kell=kell& [handles.data.eventdata_orig.(handles.data.rules(i).selectedeventdatafield)] > handles.data.rules(i).value;
            elseif strcmp(handles.data.rules(i).comparefunc,'<')
                kell=kell& [handles.data.eventdata_orig.(handles.data.rules(i).selectedeventdatafield)] < handles.data.rules(i).value;
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
else
    handles.data.eventdata=handles.data.eventdata_orig;
end

function handles = updateDATA(handles)
%%

% valtozok.timeback=.001;
% valtozok.timeforward=.001;
% valtozok.movingn=3;

valtozok.timeback=str2num(get(handles.edit2,'String'))/1000;
valtozok.timeforward=str2num(get(handles.edit3,'String'))/1000;
valtozok.movingn=round(str2num(get(handles.edit4,'String')));
if valtozok.movingn<1
    valtozok.movingn=1;
end

valtozok.depolratewindow=.001;
valtozok.baselineVwindow=.005;
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


APwaves=struct;
eventdata=handles.data.eventdata;
stimdata=handles.data.stimdata;
bridgeddata=handles.data.bridgeddata;


% if ~isempty(eventdata) & ~isempty(fieldnames(eventdata))
apidxes=neededevents;
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
        stepback=round(valtozok.timeback/si);
        stepforward=round(valtozok.timeforward/si);
    end
    onseth=eventdata(api).onseth;
    maxh=eventdata(api).maxh;
    threshh=maxh-stepbackforthresh;
    while dyfiltered(threshh)<max(dyfiltered(max(1,threshh-5):threshh)) & threshh>stepback+3
        threshh=threshh-1;
    end
    while ~((dyfiltered(threshh)<50) & mean(dyfiltered(max(1,threshh-round(.0001/si)):threshh))<50)  & threshh>stepback+3 %
        threshh=threshh-1;
    end
    threshv=yfiltered(threshh);
    eventdata(api).threshv=threshv;
    eventdata(api).APamplitude=eventdata(api).maxval-eventdata(api).threshv;
    depolrateV=yfiltered(max(threshh-round(valtozok.depolratewindow/si),1):threshh-round(.0001/si));
    depolratet=(1:length(depolrateV))*si;
    p=polyfit(depolratet',depolrateV,1);
    eventdata(api).depolrate=p(1);
    eventdata(api).baselineV=mean(yfiltered(max(threshh-round(valtozok.baselineVwindow/si),1):threshh));
    
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
    
    %             v=y([-stepback:stepforward]+centerh);
    if valtozok.movingn>1
        dv=diff(moving(v,valtozok.movingn))'/si;
        ddv=diff(moving(dv,valtozok.movingn))'/si;
        dddv=diff(moving(ddv,valtozok.movingn))'/si;
    else
        dv=diff(v/si);
        ddv=diff(dv/si);
        dddv=diff(ddv/si);
    end
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
    APwaves(apii).threshv=eventdata(api).threshv;
    APwaves(apii).depolrate=eventdata(api).depolrate;
    APwaves(apii).baselineV=eventdata(api).baselineV;
    
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
    %                 end_def
end
handles.data.eventwaves=APwaves;

function handles = plotthedata(handles)
if isfield(handles.data,'eventwaves') & ~isempty(fieldnames(handles.data.eventwaves))
    for ploti=1:4
        Xvaltoplot=get(handles.popupmenu1,'String');
        Xvaltoplot=Xvaltoplot{get(handles.popupmenu1,'Val')};
        Yvaltoplot=get(handles.popupmenu2,'String');
        Yvaltoplot=Yvaltoplot{get(handles.popupmenu2,'Val')};
        
        Xvaltoplot=handles.data.plotdata(ploti).X;
        Yvaltoplot=handles.data.plotdata(ploti).Y;
        if any(strfind(Yvaltoplot,'d'))
            Xvaltoplot=[Xvaltoplot,Yvaltoplot];
        end
        axes(handles.(['axes',num2str(ploti)]));
        cla
        sis=unique([handles.data.eventwaves.si]);
        for i=1:length(sis)
            needed=find([handles.data.eventwaves.si]==sis(i));
            hold on
            plot([handles.data.eventwaves(needed).(Xvaltoplot)],[handles.data.eventwaves(needed).(Yvaltoplot)]);
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
        % disp('lol')
    end
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
comparefunc=get(handles.popupmenu6,'String');
comparefunc=comparefunc{get(handles.popupmenu6,'Value')};

if isempty(fieldnames(handles.data.rules))
    next=1;
else
    next=length(handles.data.rules)+1;
end

handles.data.rules(next).selectedeventdatafield=selectedeventdatafield;
handles.data.rules(next).comparefunc=comparefunc;

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
