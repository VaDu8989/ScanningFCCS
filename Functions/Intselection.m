function varargout = Intselection(varargin)
% INTSELECTION MATLAB code for Intselection.fig
%      INTSELECTION, by itself, creates a new INTSELECTION or raises the existing
%      singleton*.
%
%      H = INTSELECTION returns the handle to a new INTSELECTION or the handle to
%      the existing singleton*.
%
%      INTSELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTSELECTION.M with the given input arguments.
%
%      INTSELECTION('Property','Value',...) creates a new INTSELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Intselection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Intselection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Intselection

% Last Modified by GUIDE v2.5 13-Jun-2016 17:30:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Intselection_OpeningFcn, ...
                   'gui_OutputFcn',  @Intselection_OutputFcn, ...
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

% --- Executes just before Intselection is made visible.
function Intselection_OpeningFcn(hObject, eventdata, handles, varargin)
    global segmentsincl Ifull Ifull2 timelinebinned numberofsegments segmentlength correctionfitseg correctionfitseg2 correctionfitseg2a polygrad
    %curveincl correlationcurves corfit meancorrcurve ItraceII correlationcurves2
    
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to Intselection (see VARARGIN)
    
    % Choose default command line output for Intselection
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);
    set(handles.slider3,'Max',numberofsegments);
    set(handles.slider3,'SliderStep',[1/(numberofsegments-1) 1/(numberofsegments-1)]);
    %Isegments2=Isegments;
    Ifull2=Ifull;

    for kl=1:numberofsegments
        if segmentsincl(kl)==0
            Ifull2((kl-1)*segmentlength+1:kl*segmentlength)=NaN;
        end
    end
    Ifullfit=Ifull2(isnan(Ifull2)==0);
    timelinebinnedfit=timelinebinned(isnan(Ifull2)==0);
    correctionfitseg=expfit(Ifullfit',timelinebinnedfit',timelinebinned');
    correctionfitseg2a = polyfit(timelinebinnedfit,Ifullfit,polygrad+3);
    correctionfitseg2=(polyval(correctionfitseg2a,timelinebinned))';
    % figure; plot(timelinebinned, Ifull2)
    % hold on;   plot(timelinebinned, polyval(correctionfitseg2,timelinebinned),'r-')
    %meancorrcurve=nanmean(correlationcurves2(:,2:end),2);
    % rather: exp.Fit mit updated segments!!!!!!
    %f=fit(timeline,Ifull2,'exp2');

    set(handles.checkbox2,'Value',segmentsincl(get(handles.slider3,'Value'),1));
    %plot(handles.axes2,ItraceII(:,1),ItraceII(:,get(handles.slider3,'Value')+1))
%   size(timelinebinned)
%   size(Ifull)
%   timelinebinned;
    plot(handles.axes2,timelinebinned,Ifull)
%   plot(handles.axes2,f,timeline,Ifull2); %need binning
    hold (handles.axes2,'on')
    plot(handles.axes2,timelinebinned,correctionfitseg)
    plot(handles.axes2,timelinebinned,correctionfitseg2, 'g-')
%   plot(handles.axes2,[ ItraceII(1,1) ItraceII(size(ItraceII,1),1) ],[mean(ItraceII(:,get(handles.slider3,'Value')+1),1) mean(ItraceII(:,get(handles.slider3,'Value')+1),1)],'r-')
    hold (handles.axes2,'off')
    Isegment=Ifull((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
    timesegment=timelinebinned((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
    correctionfitsegment=correctionfitseg((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
    correctionfitsegment2=correctionfitseg2((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
    %semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(handles.slider3,'Value')+1),'b+')
    plot(handles.axes1,timesegment,Isegment)
    hold (handles.axes1,'on')
    plot(handles.axes1,timesegment,correctionfitsegment);
    plot(handles.axes1,timesegment,correctionfitsegment2,'g');
    hold (handles.axes1,'off')
    %semilogx(handles.axes1,correlationcurves(:,1),meancorrcurve(:,1),'r-')
    %semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(handles.slider3,'Value')+1),'k-','LineWidth',2)
    %axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 
    %hold (handles.axes1,'off')

% UIWAIT makes Intselection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Intselection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
global segmentsincl Ifull Ifull2 timelinebinned numberofsegments segmentlength correctionfitseg correctionfitseg2 correctionfitseg2a polygrad
%curveincl correlationcurves corfit meancorrcurve ItraceII correlationcurves2

% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


if get(hObject,'Value')<numberofsegments+1 % plus 1?????  
    Ifull2=Ifull;
    for kl=1:numberofsegments
        if segmentsincl(kl)==0
            Ifull2((kl-1)*segmentlength+1:kl*segmentlength)=NaN;
        end
    end
    size(Ifull2);
    size(timelinebinned);
    Ifullfit=Ifull2(isnan(Ifull2)==0);
    timelinebinnedfit=timelinebinned(isnan(Ifull2)==0);
    correctionfitseg=expfit(Ifullfit',timelinebinnedfit',timelinebinned');
    correctionfitseg2a = polyfit(timelinebinnedfit,Ifullfit,polygrad+3);
    correctionfitseg2=(polyval(correctionfitseg2a,timelinebinned))';
    %meancorrcurve=nanmean(correlationcurves2(:,2:end),2);
    % rather: exp.Fit mit updated segments !!!!!
    get(hObject,'Value')
    set(handles.checkbox2,'Value',segmentsincl(get(hObject,'Value')));
    %get(hObject,'Min')
    %get(hObject,'Max')
    %get(hObject,'Value')+1;
    plot(handles.axes2,timelinebinned,Ifull)
    hold (handles.axes2,'on')
    plot(handles.axes2,timelinebinned,correctionfitseg)
    plot(handles.axes2,timelinebinned,correctionfitseg2, 'g-')
    %plot(handles.axes2,[ ItraceII(1,1) ItraceII(size(ItraceII,1),1) ],[mean(ItraceII(:,get(hObject,'Value')+1),1) mean(ItraceII(:,get(hObject,'Value')+1),1)],'r-')
    hold (handles.axes2,'off')
    %semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(hObject,'Value')+1),'b+')
    %hold (handles.axes1,'on')
    Isegment=Ifull((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
    timesegment=timelinebinned((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
    correctionfitsegment=correctionfitseg((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
    correctionfitsegment2=correctionfitseg2((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
    %semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(handles.slider3,'Value')+1),'b+')
    size(timelinebinned);
    size(Ifull);
    size(timesegment);
    size(Isegment);
    plot(handles.axes1,timesegment,Isegment)
    hold (handles.axes1,'on')
    plot(handles.axes1,timesegment,correctionfitsegment);
    plot(handles.axes1,timesegment,correctionfitsegment2,'g');
    hold (handles.axes1,'off')
    %semilogx(handles.axes1,correlationcurves(:,1),meancorrcurve(:,1),'r-')
    %semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(hObject,'Value')+1),'k-','LineWidth',2)
    %axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(hObject,'Value')+1))) max(max(correlationcurves(1:10,get(hObject,'Value')+1)))]) 
    %hold (handles.axes1,'off')
end

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)

% modify this part still!!!!!!

global segmentsincl Ifull Ifull2 timelinebinned numberofsegments segmentlength correctionfitseg correctionfitseg2 correctionfitseg2a polygrad

% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
segmentsincl(get(handles.slider3,'value'))=get(hObject,'Value');
Ifull2=Ifull;

for kl=1:numberofsegments
    if segmentsincl(kl)==0
         Ifull2((kl-1)*segmentlength+1:kl*segmentlength)=NaN;
    end
end

Ifullfit=Ifull2(isnan(Ifull2)==0);
timelinebinnedfit=timelinebinned(isnan(Ifull2)==0);
correctionfitseg=expfit(Ifullfit',timelinebinnedfit',timelinebinned');
correctionfitseg2a = polyfit(timelinebinnedfit,Ifullfit,polygrad+3);
correctionfitseg2=(polyval(correctionfitseg2a,timelinebinned))';
%meancorrcurve=nanmean(correlationcurves2(:,2:end),2);
set(hObject,'Value',segmentsincl(get(handles.slider3,'Value')));
plot(handles.axes2,timelinebinned,Ifull)
%plot(handles.axes2,ItraceII(:,1),ItraceII(:,get(handles.slider3,'Value')+1))
hold (handles.axes2,'on')
plot(handles.axes2,timelinebinned,correctionfitseg)
plot(handles.axes2,timelinebinned,correctionfitseg2, 'g-')
%plot(handles.axes2,[ ItraceII(1,1) ItraceII(size(ItraceII,1),1) ],[mean(ItraceII(:,get(handles.slider3,'Value')+1),1) mean(ItraceII(:,get(handles.slider3,'Value')+1),1)],'r-')
hold (handles.axes2,'off')
Isegment=Ifull((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
timesegment=timelinebinned((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
correctionfitsegment=correctionfitseg((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
correctionfitsegment2=correctionfitseg2((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
 %semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(handles.slider3,'Value')+1),'b+')
plot(handles.axes1,timesegment,Isegment)
hold (handles.axes1,'on')
plot(handles.axes1,timesegment,correctionfitsegment);
plot(handles.axes1,timesegment,correctionfitsegment2,'g');
hold (handles.axes1,'off')
% semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(handles.slider3,'Value')+1),'b+')
%hold (handles.axes1,'on')
%semilogx(handles.axes1,correlationcurves(:,1),meancorrcurve(:,1),'r-')
%semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(handles.slider3,'Value')+1),'k-','LineWidth',2)
%axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 
%hold (handles.axes1,'off')


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global goon
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
goon=1;
close Intselection
