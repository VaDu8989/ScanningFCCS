

function varargout = Corrselection(varargin)
% CORRSELECTION MATLAB code for Corrselection.fig
%      CORRSELECTION, by itself, creates a new CORRSELECTION or raises the existing
%      singleton*.
%
%      H = CORRSELECTION returns the handle to a new CORRSELECTION or the handle to
%      the existing singleton*.
%
%      CORRSELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORRSELECTION.M with the given input arguments.
%
%      CORRSELECTION('Property','Value',...) creates a new CORRSELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Corrselection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Corrselection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Corrselection

% Last Modified by GUIDE v2.5 31-Jan-2021 14:44:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Corrselection_OpeningFcn, ...
                   'gui_OutputFcn',  @Corrselection_OutputFcn, ...
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

% --- Executes just before Corrselection is made visible.
function Corrselection_OpeningFcn(hObject, eventdata, handles, varargin)
global curveincl curveincl_poly correlationcurves correlationcurves_poly  meancorrcurve ItraceII correlationcurves2 S sigmascurves polygrad fitting_poly_GUI
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Corrselection (see VARARGIN)

% Choose default command line output for Corrselection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(handles.text3,'String', ['You have selected ' num2str(sum(curveincl,1)) ' out of ' num2str(size(curveincl,1)) ' curves'])

set(handles.checkbox3,'String', ['Subtract ' num2str(polygrad) ' grade polynomial'])
set(handles.slider3,'Max',size(correlationcurves,2)-1);
set(handles.slider3,'SliderStep',[1/(size(correlationcurves,2)-2) 1/(size(correlationcurves,2)-2)]);
correlationcurves2=correlationcurves;

for kl=2:size(correlationcurves,2)

if curveincl_poly(kl-1)==1
    correlationcurves2(:,kl)=correlationcurves_poly(:,kl);
end
if curveincl(kl-1)==0
    correlationcurves2(:,kl)=NaN;
end

end
meancorrcurve=nanmean(correlationcurves2(:,2:end),2);



lsfitfunc=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
x0=[1,0.01 S];
fixed=[false false true];
lb=[0 0.0001 S-1];
ub=[10000 10 S+1];
[Nsegfinal,taudsegfinal,Sfitsegfinal,CIsegfinal,fitcurve,residualssegfinal]  = autocorrfit2Ddiff(correlationcurves(:,1),meancorrcurve,lb,ub,abs(meancorrcurve./mean(sigmascurves,2)),lsfitfunc,x0,fixed);
%this fit is approximated. Sigmas are not the final ones


set(handles.checkbox2,'Value',curveincl(get(handles.slider3,'Value'),1));
set(handles.checkbox3,'Value',curveincl_poly(get(handles.slider3,'Value'),1));

    plot(handles.axes2,ItraceII(:,1),ItraceII(:,get(handles.slider3,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceII(1,1) ItraceII(size(ItraceII,1),1) ],[mean(ItraceII(:,get(handles.slider3,'Value')+1),1) mean(ItraceII(:,get(handles.slider3,'Value')+1),1)],'r-')
    if curveincl_poly(get(handles.slider3,'Value'),1)==1
            plot(handles.axes2,ItraceII(:,1),polyval(fitting_poly_GUI(:,get(handles.slider3,'Value')), ItraceII(:,1)), 'k-')

    end
    
     hold (handles.axes2,'off')

     if curveincl_poly(get(handles.slider3,'Value'),1)==1
     semilogx(handles.axes1,correlationcurves(:,1),correlationcurves_poly(:,get(handles.slider3,'Value')+1),'b+')
     else
              semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(handles.slider3,'Value')+1),'b+')

    end
hold (handles.axes1,'on')



semilogx(handles.axes1,correlationcurves(:,1),meancorrcurve(:,1),'r-')
semilogx(handles.axes1,correlationcurves(:,1),fitcurve,'k-','LineWidth',1)
 if curveincl_poly(get(handles.slider3,'Value'),1)==1
       
       axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves_poly(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves_poly(1:10,get(handles.slider3,'Value')+1)))]) 

     else
     
axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 

 end
 hold (handles.axes1,'off')

% UIWAIT makes Corrselection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Corrselection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
global curveincl_poly fitting_poly_GUI correlationcurves_poly curveincl correlationcurves corfit meancorrcurve ItraceII correlationcurves2 S sigmascurves
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


if get(hObject,'Value')<size(correlationcurves,2)
    
    correlationcurves2=correlationcurves;

for kl=2:size(correlationcurves,2)

if curveincl_poly(kl-1)==1
    correlationcurves2(:,kl)=correlationcurves_poly(:,kl);
end
if curveincl(kl-1)==0
    correlationcurves2(:,kl)=NaN;
end
end
meancorrcurve=nanmean(correlationcurves2(:,2:end),2);




lsfitfunc=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
x0=[1,0.01 S];
fixed=[false false true];
lb=[0 0.0001 S-1];
ub=[10000 1000 S+1];
% abs(meancorrcurve./mean(sigmascurves,2))
[Nsegfinal,taudsegfinal,Sfitsegfinal,CIsegfinal,fitcurve,residualssegfinal]  = autocorrfit2Ddiff(correlationcurves(:,1),meancorrcurve,lb,ub,abs(meancorrcurve./mean(sigmascurves,2)),lsfitfunc,x0,fixed);

set(handles.checkbox2,'Value',curveincl(get(hObject,'Value'),1));
set(handles.checkbox3,'Value',curveincl_poly(get(hObject,'Value'),1));

get(hObject,'Min')
get(hObject,'Max')
     get(hObject,'Value')+1;
    plot(handles.axes2,ItraceII(:,1),ItraceII(:,get(hObject,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceII(1,1) ItraceII(size(ItraceII,1),1) ],[mean(ItraceII(:,get(hObject,'Value')+1),1) mean(ItraceII(:,get(hObject,'Value')+1),1)],'r-')
    if curveincl_poly(get(hObject,'Value'),1)==1
            plot(handles.axes2,ItraceII(:,1),polyval(fitting_poly_GUI(:,get(hObject,'Value')), ItraceII(:,1)), 'k-')

    end
    
    hold (handles.axes2,'off')
get(handles.slider3,'Value')
     if curveincl_poly(get(handles.slider3,'Value'),1)==1
     semilogx(handles.axes1,correlationcurves(:,1),correlationcurves_poly(:,get(hObject,'Value')+1),'b+')
     else
     
     semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(hObject,'Value')+1),'b+')

     end
     
     
     hold (handles.axes1,'on')



semilogx(handles.axes1,correlationcurves(:,1),meancorrcurve(:,1),'r-')
semilogx(handles.axes1,correlationcurves(:,1),fitcurve,'k-','LineWidth',2)
   if curveincl_poly(get(handles.slider3,'Value'),1)==1
       
       axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves_poly(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves_poly(1:10,get(handles.slider3,'Value')+1)))]) 

     else
     
axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 

     end
hold (handles.axes1,'off')
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

global curveincl_poly fitting_poly_GUI correlationcurves_poly curveincl correlationcurves correlationcurves2 meancorrcurve corfit ItraceII S sigmascurves

% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



curveincl(get(handles.slider3,'value'),1)=get(hObject,'Value');
correlationcurves2=correlationcurves;

for kl=2:size(correlationcurves,2)

if curveincl_poly(kl-1)==1
    correlationcurves2(:,kl)=correlationcurves_poly(:,kl);
end

if curveincl(kl-1)==0
    correlationcurves2(:,kl)=NaN;
end
end
meancorrcurve=nanmean(correlationcurves2(:,2:end),2);
  


lsfitfunc=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
x0=[1,0.01 S];
fixed=[false false true];
lb=[0 0.0001 S-1];
ub=[10000 1000 S+1];
% abs(meancorrcurve./mean(sigmascurves,2))
[Nsegfinal,taudsegfinal,Sfitsegfinal,CIsegfinal,fitcurve,residualssegfinal]  = autocorrfit2Ddiff(correlationcurves(:,1),meancorrcurve,lb,ub,abs(meancorrcurve./mean(sigmascurves,2)),lsfitfunc,x0,fixed);

set(hObject,'Value',curveincl(get(handles.slider3,'Value'),1));
set(handles.checkbox3,'Value',curveincl_poly(get(handles.slider3,'Value'),1));

    plot(handles.axes2,ItraceII(:,1),ItraceII(:,get(handles.slider3,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceII(1,1) ItraceII(size(ItraceII,1),1) ],[mean(ItraceII(:,get(handles.slider3,'Value')+1),1) mean(ItraceII(:,get(handles.slider3,'Value')+1),1)],'r-')
       if curveincl_poly(get(handles.slider3,'Value'),1)==1
            plot(handles.axes2,ItraceII(:,1),polyval(fitting_poly_GUI(:,get(handles.slider3,'Value')), ItraceII(:,1)), 'k-')

    end
    
    
    hold (handles.axes2,'off')
     if curveincl_poly(get(handles.slider3,'Value'),1)==1
     semilogx(handles.axes1,correlationcurves(:,1),correlationcurves_poly(:,get(handles.slider3,'Value')+1),'b+')
     else
     semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(handles.slider3,'Value')+1),'b+')
     end
     hold (handles.axes1,'on')



semilogx(handles.axes1,correlationcurves(:,1),meancorrcurve(:,1),'r-')
semilogx(handles.axes1,correlationcurves(:,1),fitcurve,'k-','LineWidth',2)
 if curveincl_poly(get(handles.slider3,'Value'),1)==1
       
       axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves_poly(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves_poly(1:10,get(handles.slider3,'Value')+1)))]) 

     else
     
axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 

 end
 hold (handles.axes1,'off')
set(handles.text3,'String', ['You have selected ' num2str(sum(curveincl,1)) ' out of ' num2str(size(curveincl,1)) ' curves'])



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global goon
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
goon=1;
close Corrselection


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
global curveincl_poly fitting_poly_GUI correlationcurves_poly curveincl correlationcurves correlationcurves2 meancorrcurve corfit ItraceII S sigmascurves
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3

curveincl_poly(get(handles.slider3,'value'),1)=get(hObject,'Value');
correlationcurves2=correlationcurves;

for kl=2:size(correlationcurves,2)

if curveincl_poly(kl-1)==1
    correlationcurves2(:,kl)=correlationcurves_poly(:,kl);
end

if curveincl(kl-1)==0
    correlationcurves2(:,kl)=NaN;
end
end
meancorrcurve=nanmean(correlationcurves2(:,2:end),2);
  


lsfitfunc=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
x0=[1,0.01 S];
fixed=[false false true];
lb=[0 0.0001 S-1];
ub=[10000 1000 S+1];
% abs(meancorrcurve./mean(sigmascurves,2))
[Nsegfinal,taudsegfinal,Sfitsegfinal,CIsegfinal,fitcurve,residualssegfinal]  = autocorrfit2Ddiff(correlationcurves(:,1),meancorrcurve,lb,ub,abs(meancorrcurve./mean(sigmascurves,2)),lsfitfunc,x0,fixed);

set(handles.checkbox2,'Value',curveincl(get(handles.slider3,'Value'),1));
set(hObject,'Value',curveincl_poly(get(handles.slider3,'Value'),1));

    plot(handles.axes2,ItraceII(:,1),ItraceII(:,get(handles.slider3,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceII(1,1) ItraceII(size(ItraceII,1),1) ],[mean(ItraceII(:,get(handles.slider3,'Value')+1),1) mean(ItraceII(:,get(handles.slider3,'Value')+1),1)],'r-')
       if curveincl_poly(get(handles.slider3,'Value'),1)==1
            plot(handles.axes2,ItraceII(:,1),polyval(fitting_poly_GUI(:,get(handles.slider3,'Value')), ItraceII(:,1)), 'k-')

       end
    
    
    hold (handles.axes2,'off')
     if curveincl_poly(get(handles.slider3,'Value'),1)==1
     semilogx(handles.axes1,correlationcurves(:,1),correlationcurves_poly(:,get(handles.slider3,'Value')+1),'b+')
     else
     semilogx(handles.axes1,correlationcurves(:,1),correlationcurves(:,get(handles.slider3,'Value')+1),'b+')
     end
     hold (handles.axes1,'on')



semilogx(handles.axes1,correlationcurves(:,1),meancorrcurve(:,1),'r-')
semilogx(handles.axes1,correlationcurves(:,1),fitcurve,'k-','LineWidth',2)
 if curveincl_poly(get(handles.slider3,'Value'),1)==1
       
       axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves_poly(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves_poly(1:10,get(handles.slider3,'Value')+1)))]) 

     else
     
axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 

 end
 hold (handles.axes1,'off')
