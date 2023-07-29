

function varargout = Intselection2colors_new(varargin)
% INTSELECTION2COLORS_NEW MATLAB code for Intselection2colors_new.fig
%      INTSELECTION2COLORS_NEW, by itself, creates a new INTSELECTION2COLORS_NEW or raises the existing
%      singleton*.
%
%      H = INTSELECTION2COLORS_NEW returns the handle to a new INTSELECTION2COLORS_NEW or the handle to
%      the existing singleton*.
%
%      INTSELECTION2COLORS_NEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTSELECTION2COLORS_NEW.M with the given input arguments.
%
%      INTSELECTION2COLORS_NEW('Property','Value',...) creates a new INTSELECTION2COLORS_NEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Intselection2colors_new_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Intselection2colors_new_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Intselection2colors_new

% Last Modified by GUIDE v2.5 22-Oct-2022 19:07:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Intselection2colors_new_OpeningFcn, ...
                   'gui_OutputFcn',  @Intselection2colors_new_OutputFcn, ...
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

% --- Executes just before Intselection2colors_new is made visible.
function Intselection2colors_new_OpeningFcn(hObject, eventdata, handles, varargin)
global segmentsincl1 correctionfitseg1a correctionfitseg2a correctionfitseg2aa polygrad correctionfitseg1aa segmentsincl2 Ifull1 Ifull1_2 Ifull2 Ifull2_2 timeline1binned timeline2binned numberofsegments segmentlength correctionfitseg1 correctionfitseg2
%curveincl correlationcurves corfit meancorrcurve ItraceII correlationcurves2

% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Intselection2colors (see VARARGIN)

% Choose default command line output for Intselection2colors
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(handles.slider3,'Max',numberofsegments);
set(handles.slider3,'SliderStep',[1/(numberofsegments-1) 1/(numberofsegments-1)]);
%Isegments2=Isegments;
Ifull1_2=Ifull1;
Ifull2_2=Ifull2;

for kl=1:numberofsegments
if segmentsincl1(kl)==0
    Ifull1_2((kl-1)*segmentlength+1:kl*segmentlength)=NaN;
end
if segmentsincl2(kl)==0
    Ifull2_2((kl-1)*segmentlength+1:kl*segmentlength)=NaN;
end
end
Ifull1fit=Ifull1_2(isnan(Ifull1_2)==0);
Ifull2fit=Ifull2_2(isnan(Ifull2_2)==0);
timeline1binnedfit=timeline1binned(isnan(Ifull1_2)==0);
timeline2binnedfit=timeline2binned(isnan(Ifull2_2)==0);
[correctionfitseg1,correctionfitseg2]=expfit2(Ifull1fit',timeline1binnedfit',timeline1binned',Ifull2fit',timeline2binnedfit',timeline2binned');

% figure; plot (correctionfitseg2)
%     correctionfitseg1=expfit(Ifull1fit',timeline1binnedfit',timeline1binned');
    correctionfitseg1aa = polyfit(timeline1binnedfit,Ifull1fit,polygrad+3);
    correctionfitseg1a=(polyval(correctionfitseg1aa,timeline1binned))';

%       correctionfitseg2=expfit(Ifull2fit',timeline2binnedfit',timeline2binned');
    correctionfitseg2aa = polyfit(timeline2binnedfit,Ifull2fit,polygrad+3);
    correctionfitseg2a=(polyval(correctionfitseg2aa,timeline2binned))';


set(handles.checkbox2,'Value',segmentsincl1(get(handles.slider3,'Value'),1));
set(handles.checkbox4,'Value',segmentsincl2(get(handles.slider3,'Value'),1));
    %plot(handles.axes2,ItraceII(:,1),ItraceII(:,get(handles.slider3,'Value')+1))
%     size(timelinebinned)
%     size(Ifull)
%     timelinebinned;
    plot(handles.axes2,timeline1binned,Ifull1)
%    plot(handles.axes2,f,timeline,Ifull2); %need binning
     hold (handles.axes2,'on')
      plot(handles.axes2,timeline1binned,correctionfitseg1)
         plot(handles.axes2,timeline1binned,correctionfitseg1a, 'g-')
%     plot(handles.axes2,[ ItraceII(1,1) ItraceII(size(ItraceII,1),1) ],[mean(ItraceII(:,get(handles.slider3,'Value')+1),1) mean(ItraceII(:,get(handles.slider3,'Value')+1),1)],'r-')
      hold (handles.axes2,'off')
      
      plot(handles.axes5,timeline2binned,Ifull2)
      hold (handles.axes5,'on')
      plot(handles.axes5,timeline2binned,correctionfitseg2)
         plot(handles.axes5,timeline2binned,correctionfitseg2a, 'g-')
      hold (handles.axes5,'off')

     Isegment1=Ifull1((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
     timesegment1=timeline1binned((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
     correctionfitsegment1=correctionfitseg1((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
          correctionfitsegment1a=correctionfitseg1a((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);

     Isegment2=Ifull2((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
     timesegment2=timeline2binned((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
     correctionfitsegment2=correctionfitseg2((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
          correctionfitsegment2a=correctionfitseg2a((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);

     %semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(handles.slider3,'Value')+1),'b+')
     plot(handles.axes1,timesegment1,Isegment1)
    hold (handles.axes1,'on')
    plot(handles.axes1,timesegment1,correctionfitsegment1);
       plot(handles.axes1,timesegment1,correctionfitsegment1a,'g');
    ylim
    hold (handles.axes1,'off')
    
    plot(handles.axes4,timesegment2,Isegment2)
    hold (handles.axes4,'on')
    plot(handles.axes4,timesegment2,correctionfitsegment2);
       plot(handles.axes4,timesegment2,correctionfitsegment2a,'g');
    hold (handles.axes4,'off')



%semilogx(handles.axes1,correlationcurves(:,1),meancorrcurve(:,1),'r-')
%semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(handles.slider3,'Value')+1),'k-','LineWidth',2)
%axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 
%hold (handles.axes1,'off')

% UIWAIT makes Intselection2colors_new wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Intselection2colors_new_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
global segmentsincl1 correctionfitseg1a correctionfitseg2a correctionfitseg2aa polygrad correctionfitseg1aa segmentsincl2 Ifull1 Ifull1_2 Ifull2 Ifull2_2 timeline1binned timeline2binned numberofsegments segmentlength correctionfitseg1 correctionfitseg2
%curveincl correlationcurves corfit meancorrcurve ItraceII correlationcurves2

% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


if get(hObject,'Value')<numberofsegments+1 % plus 1?????
    
    Ifull1_2=Ifull1;
    Ifull2_2=Ifull2;

for kl=1:numberofsegments
if segmentsincl1(kl)==0
    Ifull1_2((kl-1)*segmentlength+1:kl*segmentlength)=NaN;
end
if segmentsincl2(kl)==0
    Ifull2_2((kl-1)*segmentlength+1:kl*segmentlength)=NaN;
end
end
size(Ifull1_2);
size(timeline1binned);
Ifull1fit=Ifull1_2(isnan(Ifull1_2)==0);
Ifull2fit=Ifull2_2(isnan(Ifull2_2)==0);
timeline1binnedfit=timeline1binned(isnan(Ifull1_2)==0);
timeline2binnedfit=timeline2binned(isnan(Ifull2_2)==0);
[correctionfitseg1,correctionfitseg2]=expfit2(Ifull1fit',timeline1binnedfit',timeline1binned',Ifull2fit',timeline2binnedfit',timeline2binned');

%     correctionfitseg1=expfit(Ifull1fit',timeline1binnedfit',timeline1binned');

    correctionfitseg1aa = polyfit(timeline1binnedfit,Ifull1fit,polygrad+3);
    correctionfitseg1a=(polyval(correctionfitseg1aa,timeline1binned))';

%       correctionfitseg2=expfit(Ifull2fit',timeline2binnedfit',timeline2binned');
    correctionfitseg2aa = polyfit(timeline2binnedfit,Ifull2fit,polygrad+3);
    correctionfitseg2a=(polyval(correctionfitseg2aa,timeline2binned))';


size(correctionfitseg1)
%meancorrcurve=nanmean(correlationcurves2(:,2:end),2);
% rather: exp.Fit mit updated segments !!!!!
get(hObject,'Value');
set(handles.checkbox2,'Value',segmentsincl1(get(hObject,'Value')));
set(handles.checkbox4,'Value',segmentsincl2(get(hObject,'Value')));
%get(hObject,'Min')
%get(hObject,'Max')
     %get(hObject,'Value')+1;
    plot(handles.axes2,timeline1binned,Ifull1)
    hold (handles.axes2,'on')
    plot(handles.axes2,timeline1binned,correctionfitseg1)
        plot(handles.axes2,timeline1binned,correctionfitseg1a, 'g')
    %plot(handles.axes2,[ ItraceII(1,1) ItraceII(size(ItraceII,1),1) ],[mean(ItraceII(:,get(hObject,'Value')+1),1) mean(ItraceII(:,get(hObject,'Value')+1),1)],'r-')
    hold (handles.axes2,'off')
    
    plot(handles.axes5,timeline2binned,Ifull2)
    hold (handles.axes5,'on')
    plot(handles.axes5,timeline2binned,correctionfitseg2)
     plot(handles.axes5,timeline2binned,correctionfitseg2a, 'g-')
    hold (handles.axes5,'off')

     
     %semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(hObject,'Value')+1),'b+')
%hold (handles.axes1,'on')

Isegment1=Ifull1((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
timesegment1=timeline1binned((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
correctionfitsegment1=correctionfitseg1((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
correctionfitsegment1a=correctionfitseg1a((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);

Isegment2=Ifull2((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
timesegment2=timeline2binned((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
correctionfitsegment2=correctionfitseg2((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
    correctionfitsegment2a=correctionfitseg2a((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);

%semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(handles.slider3,'Value')+1),'b+')
     size(timeline1binned);
     size(Ifull1);
     size(timesegment1);
     size(Isegment1);
     plot(handles.axes1,timesegment1,Isegment1)
     hold (handles.axes1,'on')
    plot(handles.axes1,timesegment1,correctionfitsegment1);
       plot(handles.axes1,timesegment1,correctionfitsegment1a, 'g');
    hold (handles.axes1,'off')
    
    plot(handles.axes4,timesegment2,Isegment2)
     hold (handles.axes4,'on')
    plot(handles.axes4,timesegment2,correctionfitsegment2);
      plot(handles.axes4,timesegment2,correctionfitsegment2a, 'g');
    hold (handles.axes4,'off')

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

global segmentsincl1 correctionfitseg1a correctionfitseg2a correctionfitseg2aa polygrad correctionfitseg1aa segmentsincl2 Ifull1 Ifull1_2 Ifull2 Ifull2_2 timeline1binned timeline2binned numberofsegments segmentlength correctionfitseg1 correctionfitseg2

% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



segmentsincl1(get(handles.slider3,'value'))=get(hObject,'Value');
Ifull1_2=Ifull1;
Ifull2_2=Ifull2;

for kl=1:numberofsegments
if segmentsincl1(kl)==0
     Ifull1_2((kl-1)*segmentlength+1:kl*segmentlength)=NaN;
end
end
Ifull1fit=Ifull1_2(isnan(Ifull1_2)==0);
Ifull2fit=Ifull2_2(isnan(Ifull2_2)==0);
timeline1binnedfit=timeline1binned(isnan(Ifull1_2)==0);
timeline2binnedfit=timeline2binned(isnan(Ifull2_2)==0);
[correctionfitseg1,correctionfitseg2]=expfit2(Ifull1fit',timeline1binnedfit',timeline1binned',Ifull2fit',timeline2binnedfit',timeline2binned');

%  correctionfitseg1=expfit(Ifull1fit',timeline1binnedfit',timeline1binned');

    correctionfitseg1aa = polyfit(timeline1binnedfit,Ifull1fit,polygrad+3);
    correctionfitseg1a=(polyval(correctionfitseg1aa,timeline1binned))';

%       correctionfitseg2=expfit(Ifull2fit',timeline2binnedfit',timeline2binned');
    correctionfitseg2aa = polyfit(timeline2binnedfit,Ifull2fit,polygrad+3);
    correctionfitseg2a=(polyval(correctionfitseg2aa,timeline2binned))';

%meancorrcurve=nanmean(correlationcurves2(:,2:end),2);
  
set(hObject,'Value',segmentsincl1(get(handles.slider3,'Value')));
plot(handles.axes2,timeline1binned,Ifull1)
    %plot(handles.axes2,ItraceII(:,1),ItraceII(:,get(handles.slider3,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,timeline1binned,correctionfitseg1)
       plot(handles.axes2,timeline1binned,correctionfitseg1a, 'g')
    %plot(handles.axes2,[ ItraceII(1,1) ItraceII(size(ItraceII,1),1) ],[mean(ItraceII(:,get(handles.slider3,'Value')+1),1) mean(ItraceII(:,get(handles.slider3,'Value')+1),1)],'r-')
     hold (handles.axes2,'off')
     
     
     Isegment1=Ifull1((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
timesegment1=timeline1binned((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
correctionfitsegment1=correctionfitseg1((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
correctionfitsegment1a=correctionfitseg1a((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);

     %semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(handles.slider3,'Value')+1),'b+')
     plot(handles.axes1,timesegment1,Isegment1)
     hold (handles.axes1,'on')
plot(handles.axes1,timesegment1,correctionfitsegment1);
plot(handles.axes1,timesegment1,correctionfitsegment1a, 'g');
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
close Intselection2colors_new


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4

global correctionfitseg1a correctionfitseg2a correctionfitseg2aa polygrad correctionfitseg1aa segmentsincl1 segmentsincl2 Ifull1 Ifull1_2 Ifull2 Ifull2_2 timeline1binned timeline2binned numberofsegments segmentlength correctionfitseg1 correctionfitseg2

% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



segmentsincl2(get(handles.slider3,'value'))=get(hObject,'Value');
Ifull1_2=Ifull1;
Ifull2_2=Ifull2;

for kl=1:numberofsegments
if segmentsincl2(kl)==0
     Ifull2_2((kl-1)*segmentlength+1:kl*segmentlength)=NaN;
end
end

Ifull1fit=Ifull1_2(isnan(Ifull1_2)==0);
Ifull2fit=Ifull2_2(isnan(Ifull2_2)==0);
timeline1binnedfit=timeline1binned(isnan(Ifull1_2)==0);
timeline2binnedfit=timeline2binned(isnan(Ifull2_2)==0);
[correctionfitseg1,correctionfitseg2]=expfit2(Ifull1fit',timeline1binnedfit',timeline1binned',Ifull2fit',timeline2binnedfit',timeline2binned');

%  correctionfitseg1=expfit(Ifull1fit',timeline1binnedfit',timeline1binned');

    correctionfitseg1aa = polyfit(timeline1binnedfit,Ifull1fit,polygrad+3);
    correctionfitseg1a=(polyval(correctionfitseg1aa,timeline1binned))';

%       correctionfitseg2=expfit(Ifull2fit',timeline2binnedfit',timeline2binned');
    correctionfitseg2aa = polyfit(timeline2binnedfit,Ifull2fit,polygrad+3);
    correctionfitseg2a=(polyval(correctionfitseg2aa,timeline2binned))';
    
    %meancorrcurve=nanmean(correlationcurves2(:,2:end),2);
  
set(hObject,'Value',segmentsincl2(get(handles.slider3,'Value')));
     
     plot(handles.axes5,timeline2binned,Ifull2)
    %plot(handles.axes2,ItraceII(:,1),ItraceII(:,get(handles.slider3,'Value')+1))
    hold (handles.axes5,'on')
    plot(handles.axes5,timeline2binned,correctionfitseg2)
        plot(handles.axes5,timeline2binned,correctionfitseg2a, 'g')
    %plot(handles.axes2,[ ItraceII(1,1) ItraceII(size(ItraceII,1),1) ],[mean(ItraceII(:,get(handles.slider3,'Value')+1),1) mean(ItraceII(:,get(handles.slider3,'Value')+1),1)],'r-')
     hold (handles.axes5,'off')    
    
Isegment2=Ifull2((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
timesegment2=timeline2binned((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
correctionfitsegment2=correctionfitseg2((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
correctionfitsegment2a=correctionfitseg2a((get(handles.slider3,'Value')-1)*segmentlength+1:get(handles.slider3,'Value')*segmentlength);
    
%semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(handles.slider3,'Value')+1),'b+')
plot(handles.axes4,timesegment2,Isegment2)
     hold (handles.axes4,'on')
plot(handles.axes4,timesegment2,correctionfitsegment2);
plot(handles.axes4,timesegment2,correctionfitsegment2a, 'g');
hold (handles.axes4,'off')

     
    % semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(handles.slider3,'Value')+1),'b+')
%hold (handles.axes1,'on')



%semilogx(handles.axes1,correlationcurves(:,1),meancorrcurve(:,1),'r-')
%semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(handles.slider3,'Value')+1),'k-','LineWidth',2)
%axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 
%hold (handles.axes1,'off')
