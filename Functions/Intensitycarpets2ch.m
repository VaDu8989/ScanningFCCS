function varargout = Intensitycarpets2ch(varargin)
% INTENSITYCARPETS2CH MATLAB code for Intensitycarpets2ch.fig
%      INTENSITYCARPETS2CH, by itself, creates a new INTENSITYCARPETS2CH or raises the existing
%      singleton*.
%
%      H = INTENSITYCARPETS2CH returns the handle to a new INTENSITYCARPETS2CH or the handle to
%      the existing singleton*.
%
%      INTENSITYCARPETS2CH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTENSITYCARPETS2CH.M with the given input arguments.
%
%      INTENSITYCARPETS2CH('Property','Value',...) creates a new INTENSITYCARPETS2CH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Intensitycarpets2ch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Intensitycarpets2ch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Intensitycarpets2ch

% Last Modified by GUIDE v2.5 02-Jan-2017 09:29:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Intensitycarpets2ch_OpeningFcn, ...
                   'gui_OutputFcn',  @Intensitycarpets2ch_OutputFcn, ...
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


% --- Executes just before Intensitycarpets2ch is made visible.
function Intensitycarpets2ch_OpeningFcn(hObject, eventdata, handles, varargin)

global carpetsegments carpetCh1 carpetCh2 windowlength redmap greenmap
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Intensitycarpets2ch (see VARARGIN)

% Choose default command line output for Intensitycarpets2ch
handles.output = hObject;

load('redmap.mat');
load('greenmap.mat')

% Update handles structure
guidata(hObject, handles);
set(handles.slider1,'Max',carpetsegments);
set(handles.slider1,'Min',1);
set(handles.slider1,'Value',1);
set(handles.slider1,'SliderStep',[1/(carpetsegments-1) 1/(carpetsegments-1)]);
set(handles.slider2,'Max',carpetsegments);
set(handles.slider2,'Min',1);
set(handles.slider2,'Value',1);
set(handles.slider2,'SliderStep',[1/(carpetsegments-1) 1/(carpetsegments-1)]);


windowlength=length(carpetCh1(:,1))/carpetsegments;
%=========================================================================%
% Here: Plot only first carpetsegment. Set using slider!!!!!!!
get(handles.slider1,'Value')
carpetwindowCh1=carpetCh1((get(handles.slider1,'Value')-1)*windowlength+1:get(handles.slider1,'Value')*windowlength,:);
carpetwindowCh2=carpetCh2((get(handles.slider2,'Value')-1)*windowlength+1:get(handles.slider2,'Value')*windowlength,:);
%=========================================================================%



imagesc(carpetwindowCh1,'Parent',handles.axes1);
colormap(handles.axes1,greenmap)
imagesc(carpetwindowCh2,'Parent',handles.axes2);
colormap(handles.axes2,redmap)

% UIWAIT makes Intensitycarpets2ch wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Intensitycarpets2ch_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
global carpetCh1 windowlength factorgreen
load('greenmap.mat')
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

carpetwindowCh1=carpetCh1((uint32(get(hObject,'Value'))-1)*windowlength+1:uint32(get(hObject,'Value'))*windowlength,:);
get(hObject,'Value')
imagesc(carpetwindowCh1,'Parent',handles.axes1);
colormap(handles.axes1,greenmap)
ax=handles.axes1;
ax.CLim(2)=ax.CLim(2)/get(handles.slider3, 'Value');
factorgreen=get(handles.slider3, 'Value');
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
global carpetCh2 windowlength factorred
load('redmap.mat');
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

carpetwindowCh2=carpetCh2((uint32(get(hObject,'Value'))-1)*windowlength+1:uint32(get(hObject,'Value'))*windowlength,:);
imagesc(carpetwindowCh2,'Parent',handles.axes2);
colormap(handles.axes2,redmap)
ax=handles.axes2;
ax.CLim(2)=ax.CLim(2)/get(handles.slider6, 'Value');
factorred=get(handles.slider6, 'Value');


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


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min',0)
set(hObject, 'Max',10)
set(hObject, 'Value',2)
set(hObject, 'SliderStep',[.1 .1])


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
global gooncarpets
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gooncarpets=1;
close Intensitycarpets2ch


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global carpetCh1 carpetCh2 windowlength
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('greenmap.mat');

carpetwindowCh1=carpetCh1((uint32(get(handles.slider1,'Value'))-1)*windowlength+1:uint32(get(handles.slider1,'Value'))*windowlength,:);
%=========================================================================%



imagesc(carpetwindowCh1,'Parent',handles.axes1);
colormap(handles.axes1,greenmap)
ax=handles.axes1;
ax.CLim(2)=ax.CLim(2)/get(handles.slider3, 'Value');
axes(handles.axes1);
windowout=roipoly;
carpetwindowCh1=carpetwindowCh1.*(1-windowout);
carpetCh1((uint32(get(handles.slider1,'Value'))-1)*windowlength+1:uint32(get(handles.slider1,'Value'))*windowlength,:)=carpetwindowCh1;
imagesc(carpetwindowCh1,'Parent',handles.axes1);
colormap(handles.axes1,greenmap)
ax=handles.axes1;
ax.CLim(2)=ax.CLim(2)/get(handles.slider3, 'Value');



% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
global carpetCh1 carpetCh2 windowlength
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('redmap.mat');

carpetwindowCh2=carpetCh2((uint32(get(handles.slider2,'Value'))-1)*windowlength+1:uint32(get(handles.slider2,'Value'))*windowlength,:);
%=========================================================================%



imagesc(carpetwindowCh2,'Parent',handles.axes2);
colormap(handles.axes2,redmap)
ax=handles.axes2;
ax.CLim(2)=ax.CLim(2)/get(handles.slider6, 'Value');
axes(handles.axes2);
windowout=roipoly;
carpetwindowCh2=carpetwindowCh2.*(1-windowout);
carpetCh2((uint32(get(handles.slider2,'Value'))-1)*windowlength+1:uint32(get(handles.slider2,'Value'))*windowlength,:)=carpetwindowCh2;
imagesc(carpetwindowCh2,'Parent',handles.axes2);
colormap(handles.axes2,redmap)
ax=handles.axes2;
ax.CLim(2)=ax.CLim(2)/get(handles.slider6, 'Value');
% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min',0)
set(hObject, 'Max',10)
set(hObject, 'Value',2)
set(hObject, 'SliderStep',[.1 .1])

% --- Executes on button press in pushbutton3.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
