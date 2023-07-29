

function varargout = Corrselection2Channelsindividuellpreview_new(varargin)
% CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW_NEW MATLAB code for Corrselection2Channelsindividuellpreview_new.fig
%      CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW_NEW, by itself, creates a new CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW_NEW or raises the existing
%      singleton*.
%
%      H = CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW_NEW returns the handle to a new CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW_NEW or the handle to
%      the existing singleton*.
%
%      CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW_NEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW_NEW.M with the given input arguments.
%
%      CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW_NEW('Property','Value',...) creates a new CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW_NEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Corrselection2Channelsindividuellpreview_new_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Corrselection2Channelsindividuellpreview_new_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Corrselection2Channelsindividuellpreview_new

% Last Modified by GUIDE v2.5 25-Oct-2022 20:28:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Corrselection2Channelsindividuellpreview_new_OpeningFcn, ...
                   'gui_OutputFcn',  @Corrselection2Channelsindividuellpreview_new_OutputFcn, ...
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

% --- Executes just before Corrselection2Channelsindividuellpreview_new is made visible.
function Corrselection2Channelsindividuellpreview_new_OpeningFcn(hObject, eventdata, handles, varargin)
global x0 y0 fixed curveincl1 curveincl2 curveincl1_poly curveincl2_poly correlationcurvesCh1 correlationcurvesCh2 
global correlationcurvesCh1_poly correlationcurvesCh2_poly correlationcurvesChCC correlationcurvesChCC_poly   fitting_poly_GUI2 polygrad fitting_poly_GUI1

global sigmas1curves sigmas2curves sigmascccurves     meancorrcurve meancorrcurveCC ItraceIICh1 ItraceIICh2  correlationcurves2Ch1 correlationcurves2Ch2 correlationcurves2ChCC sigmas1curves2 sigmas2curves2 sigmascccurves2    


% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Corrselection2Channelsindividuellpreview_new (see VARARGIN)

% Choose default command line output for Corrselection2Channelsindividuellpreview_new
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.text5,'String', ['You have selected ' num2str(sum(curveincl1,1)) ' out of ' num2str(size(curveincl1,1)) ' curves'])
set(handles.text6,'String', [num2str(sum(curveincl2,1)) ' out of ' num2str(size(curveincl2,1)) ' curves. CC curves are ' num2str(sum(curveincl1+curveincl2==2)) ])

set(handles.checkbox4,'String', ['Subtract ' num2str(polygrad) ' grade polynomial'])
set(handles.checkbox5,'String', ['Subtract ' num2str(polygrad) ' grade polynomial'])

set(handles.slider3,'Max',size(correlationcurvesCh1,2)-1);
set(handles.slider3,'SliderStep',[1/(size(correlationcurvesCh1,2)-2) 1/(size(correlationcurvesCh1,2)-2)]);
correlationcurves2Ch1=correlationcurvesCh1;
correlationcurves2Ch2=correlationcurvesCh2;
correlationcurves2ChCC=correlationcurvesChCC;

sigmas1curves2=sigmas1curves;
sigmas2curves2=sigmas2curves;
sigmascccurves2=sigmascccurves;


%approximation for now: if one channes gets the polynomial correction
% the CC is automatically taken as it were from the 2 channels both
% corrected
for kl=2:size(correlationcurvesCh1,2)
    
if curveincl1_poly(kl-1)==1
    correlationcurves2Ch1(:,kl)=correlationcurvesCh1_poly(:,kl);
    correlationcurves2ChCC(:,kl)=correlationcurvesChCC_poly(:,kl); %or will this be calculated?
     
end

if curveincl2_poly(kl-1)==1
    correlationcurves2Ch2(:,kl)=correlationcurvesCh2_poly(:,kl);
    correlationcurves2ChCC(:,kl)=correlationcurvesChCC_poly(:,kl); %or will this be calculated?

end
    
    if curveincl1(kl-1)==0
    correlationcurves2Ch1(:,kl)=NaN;
    correlationcurves2ChCC(:,kl)=NaN;
    sigmas1curves2(:,kl-1)=NaN;
    sigmascccurves2(:,kl-1)=NaN;
    end

    
    
    
if curveincl2(kl-1)==0
correlationcurves2Ch2(:,kl)=NaN;
correlationcurves2ChCC(:,kl)=NaN;
sigmas2curves2(:,kl-1)=NaN;
sigmascccurves2(:,kl-1)=NaN;
end
end



meancorrcurve(:,1)=nanmean(correlationcurves2Ch1(:,2:end),2);
meancorrcurve(:,2)=nanmean(correlationcurves2Ch2(:,2:end),2);
meancorrcurveCC=nanmean(correlationcurves2ChCC(:,2:end),2);


sigmas1curvesmean=nanmean(sigmas1curves2,2);
sigmas2curvesmean=nanmean(sigmas2curves2,2);
sigmascccurvesmean=nanmean(sigmascccurves2,2);


set(handles.checkbox2,'Value',curveincl1(get(handles.slider3,'Value'),1));
set(handles.checkbox3,'Value',curveincl2(get(handles.slider3,'Value'),1));
set(handles.checkbox4,'Value',curveincl1_poly(get(handles.slider3,'Value'),1));
set(handles.checkbox5,'Value',curveincl2_poly(get(handles.slider3,'Value'),1));

    plot(handles.axes2,ItraceIICh1(:,1),ItraceIICh1(:,get(handles.slider3,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceIICh1(1,1) ItraceIICh1(size(ItraceIICh1,1),1) ],[mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1)],'r-')
      ylim(handles.axes2,[min(ItraceIICh1(:,get(handles.slider3,'Value')+1)) max(ItraceIICh1(:,get(handles.slider3,'Value')+1))])
     if curveincl1_poly(get(handles.slider3,'Value'),1)==1
              plot(handles.axes2,ItraceIICh1(:,1),polyval(fitting_poly_GUI1(:,get(handles.slider3,'Value')), ItraceIICh1(:,1)), 'k-')
 
  end
     hold (handles.axes2,'off')
%  
       plot(handles.axes3,ItraceIICh2(:,1),ItraceIICh2(:,get(handles.slider3,'Value')+1))
      hold (handles.axes3,'on')
      plot(handles.axes3,[ ItraceIICh2(1,1) ItraceIICh2(size(ItraceIICh2,1),1) ],[mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1)],'r-')
      ylim(handles.axes3,[min(ItraceIICh2(:,get(handles.slider3,'Value')+1)) max(ItraceIICh2(:,get(handles.slider3,'Value')+1))])
         if curveincl2_poly(get(handles.slider3,'Value'),1)==1
              plot(handles.axes3,ItraceIICh2(:,1),polyval(fitting_poly_GUI2(:,get(handles.slider3,'Value')), ItraceIICh2(:,1)), 'k-')
 
     end
     hold (handles.axes3,'off')

     if curveincl1_poly(get(handles.slider3,'Value'),1)==1
 semilogx(handles.axes1,correlationcurvesCh1_poly(:,1),correlationcurvesCh1_poly(:,get(handles.slider3,'Value')+1),'g+')

     else
     semilogx(handles.axes1,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+')
     end
          hold (handles.axes1,'on')

     
          if curveincl2_poly(get(handles.slider3,'Value'),1)==1
    semilogx(handles.axes1,correlationcurvesCh2_poly(:,1),correlationcurvesCh2_poly(:,get(handles.slider3,'Value')+1),'r+')

          else     
    semilogx(handles.axes1,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+')
          end
          
          
         if  curveincl1_poly(get(handles.slider3,'Value'),1) || curveincl2_poly(get(handles.slider3,'Value'),1)==1
 semilogx(handles.axes1,correlationcurvesChCC_poly(:,1),correlationcurvesChCC_poly(:,get(handles.slider3,'Value')+1),'b+')

         else
         semilogx(handles.axes1,correlationcurvesChCC(:,1),correlationcurvesChCC(:,get(handles.slider3,'Value')+1),'b+')
         end

 semilogx(handles.axes1,correlationcurvesCh1(:,1),meancorrcurve(:,1),'g-',correlationcurvesCh2(:,1),meancorrcurve(:,2),'r-',correlationcurvesChCC(:,1),meancorrcurveCC,'b-')
% 
% 
 xlim(handles.axes1,[min(correlationcurvesCh1(:,1)) max(correlationcurvesCh1(:,1))])
% % For now: Don't plot fits of CFs
% %semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(handles.slider3,'Value')+1),'k-','LineWidth',2)
% % For now: Don't set axis manually
% %axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 
 hold (handles.axes1,'off')
% 
% 
tcorr1segfinalpre=correlationcurvesCh1(:,1);
tcorr2segfinalpre=correlationcurvesCh2(:,1);
tcorrccsegfinalfitpre=correlationcurvesChCC(:,1);
lbseg=tcorr1segfinalpre(1);
ubseg=tcorr1segfinalpre(end);
lbccseg=tcorrccsegfinalfitpre(1);
ubccseg=tcorrccsegfinalfitpre(end);
% 
fcorr1segfinalpre=meancorrcurve(:,1);
fcorr2segfinalpre=meancorrcurve(:,2);
fcrosscorrsegfinalfitpre=meancorrcurveCC;
% 
% approximated, because they do not take into account poly correction
weights1segfinalpre=abs(fcorr1segfinalpre./sigmas1curvesmean);
weights2segfinalpre=abs(fcorr2segfinalpre./sigmas2curvesmean);
weightsccsegfinalfitpre=abs(fcrosscorrsegfinalfitpre./sigmascccurvesmean);

lsautofitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
flscrossfitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5+x(4);

%x0=[100,0.5, S]; % Set in header or GUI manually later!
%fixed=[false false true];
fixedcc=[false true true ];
%flscrossfitfunc=lsautofitfunc;
%flscrossfitfunc=@(y,tcc)1/y(1)*((1+4*y(2)*tcc./y(3)^2).^-0.5).*((1+4*y(2)*tcc./(y(3)*y(4))^2).^-0.5).*exp(-d^2./(y(3)^2+4*y(2).*tcc));
%y0=[5000,100,S,10^-5]; %[N0,tau0,S], Set in header or GUI manually later!
%y0=[2000,0.5,S];

% For now: separate fit. Later: Global fit of all 3 curves (2x auto & cc)
% [N1finalpre,taud1finalpre,Sfit1finalpre,CI1finalpre,fitcurve1finalpre,residuals1finalpre] = autocorrfit2Ddiff(tcorr1fit,fcorr1fit,lb,ub,weights1fit,lsautofitfuncpre,x0,fixed);
% [N2finalpre,taud2finalpre,Sfit2finalpre,CI2finalpre,fitcurve2finalpre,residuals2finalpre] = autocorrfit2Ddiff(tcorr2fit,fcorr2fit,lb,ub,weights2fit,lsautofitfuncpre,x0,fixed);
% %[Nccfinal,tauccfinal,Sfitccfinal,CIccfinal,fitcurveccfinal,residualsccfinal] =autocorrfit2Ddiffconst(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,flscrossfitfunc,y0,fixedcc);
% [Nccfinalpre,tauccfinalpre,Sfitccfinalpre,CIccfinalpre,fitcurveccfinalpre,residualsccfinalpre] =autocorrfit2Ddiff(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,lsautofitfuncpre,y0,fixed);

[N1segfinalpre,taud1segfinalpre,Sfit1segfinalpre,CI1segfinalpre,fitcurve1segfinalpre,residuals1segfinalpre] = autocorrfit2Ddiffpreview(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
[N2segfinalpre,taud2segfinalpre,Sfit2segfinalpre,CI2segfinalpre,fitcurve2segfinalpre,residuals2segfinalpre] = autocorrfit2Ddiffpreview(tcorr2segfinalpre,fcorr2segfinalpre,lbseg,ubseg,weights2segfinalpre,lsautofitfuncpre,x0,fixed);


try
%[Nccsegfinal,tauccsegfinal,Sfitccsegfinal,CIccsegfinal,fitcurveccsegfinal,residualsccsegfinal] =autocorrfit2Ddiffconst(tcorrccsegfinalfit,fcrosscorrsegfinalfit',lbccseg,ubccseg,weightsccsegfinalfit',flscrossfitfunc,y0,fixedcc);

[Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] =autocorrfit2Ddiffpreview(tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,lbccseg,ubccseg,weightsccsegfinalfitpre,lsautofitfuncpre,y0,fixed);
catch
[Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] = autocorrfit2Ddiffpreview(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
%    size(fitcurveccsegfinalpre)
 fitcurveccsegfinalpre  =zeros(size(fitcurveccsegfinalpre,1)+1,size(fitcurveccsegfinalpre,1));

end
% h=figure('OuterPosition',[scrsz(3) 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Full Curve');
% positionvector1=[0.1 0.35 0.8 0.55];
% positionvector2=[0.1 0.1 0.8 0.15];
% subplot('Position',positionvector1),semilogx(tcorr1fit,fitcurve1finalpre,'-g',tcorr2fit,fitcurve2finalpre,'-r',tcorrccfit,fitcurveccfinalpre,'-b')
% hold on
% subplot('Position',positionvector1),semilogx(tcorr1fit,fcorr1fit,'gs',tcorr2fit,fcorr2fit,'rd',tcorrccfit,fcrosscorrfit,'bx')
% xlabel('time')
% ylabel('autocorrelation')

%hh=figure('OuterPosition',[4*scrsz(3)/3 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Segment Averaged Curve');
%positionvector1=[0.1 0.35 0.8 0.55];
%positionvector2=[0.1 0.1 0.8 0.15];
semilogx(handles.axes5,tcorr2segfinalpre,fitcurve2segfinalpre,'-k',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-k')
hold (handles.axes5,'on')
semilogx(handles.axes5,tcorr2segfinalpre,fcorr2segfinalpre,'r-',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'b-')

   if curveincl2_poly(get(handles.slider3,'Value'),1)==0
semilogx(handles.axes5,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+')
hold (handles.axes5,'off')
          else     
semilogx(handles.axes5,correlationcurvesCh2(:,1),correlationcurvesCh2_poly (:,get(handles.slider3,'Value')+1),'r+')
hold (handles.axes5,'off')
   end

semilogx(handles.axes4,tcorr1segfinalpre,fitcurve1segfinalpre,'-k',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-k')
hold (handles.axes4,'on')
semilogx(handles.axes4,tcorr1segfinalpre,fcorr1segfinalpre,'g-',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'b-')

 if curveincl1_poly(get(handles.slider3,'Value'),1)==1
semilogx(handles.axes4,correlationcurvesCh1(:,1),correlationcurvesCh1_poly(:,get(handles.slider3,'Value')+1),'g+')
hold (handles.axes4,'off')
          else    

semilogx(handles.axes4,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+')
 hold (handles.axes4,'off')

 end

% UIWAIT makes Corrselection2Channelsindividuellpreview_new wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Corrselection2Channelsindividuellpreview_new_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
global x0 y0 fixed curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesChCC  sigmas1curves sigmas2curves sigmascccurves  meancorrcurve meancorrcurveCC ItraceIICh1 ItraceIICh2  correlationcurves2Ch1 correlationcurves2Ch2 correlationcurves2ChCC sigmas1curves2 sigmas2curves2 sigmascccurves2
global correlationcurvesCh1_poly correlationcurvesCh2_poly  correlationcurvesChCC_poly   fitting_poly_GUI2  fitting_poly_GUI1
global  curveincl1_poly curveincl2_poly   

% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


if get(hObject,'Value')<size(correlationcurvesCh1,2)
    
correlationcurves2Ch1=correlationcurvesCh1;
correlationcurves2Ch2=correlationcurvesCh2;
correlationcurves2ChCC=correlationcurvesChCC;

sigmas1curves2=sigmas1curves;
sigmas2curves2=sigmas2curves;
sigmascccurves2=sigmascccurves;


%approximation for now: if one channes gets the polynomial correction
% the CC is automatically taken as it were from the 2 channels both
% corrected
for kl=2:size(correlationcurvesCh1,2)
    
if curveincl1_poly(kl-1)==1
    correlationcurves2Ch1(:,kl)=correlationcurvesCh1_poly(:,kl);
    correlationcurves2ChCC(:,kl)=correlationcurvesChCC_poly(:,kl); %or will this be calculated?
     
end

if curveincl2_poly(kl-1)==1
    correlationcurves2Ch2(:,kl)=correlationcurvesCh2_poly(:,kl);
    correlationcurves2ChCC(:,kl)=correlationcurvesChCC_poly(:,kl); %or will this be calculated?

end
    
    if curveincl1(kl-1)==0
    correlationcurves2Ch1(:,kl)=NaN;
    correlationcurves2ChCC(:,kl)=NaN;
    sigmas1curves2(:,kl-1)=NaN;
    sigmascccurves2(:,kl-1)=NaN;
    end

    
    
    
if curveincl2(kl-1)==0
correlationcurves2Ch2(:,kl)=NaN;
correlationcurves2ChCC(:,kl)=NaN;
sigmas2curves2(:,kl-1)=NaN;
sigmascccurves2(:,kl-1)=NaN;
end
end



meancorrcurve(:,1)=nanmean(correlationcurves2Ch1(:,2:end),2);
meancorrcurve(:,2)=nanmean(correlationcurves2Ch2(:,2:end),2);
meancorrcurveCC=nanmean(correlationcurves2ChCC(:,2:end),2);


sigmas1curvesmean=nanmean(sigmas1curves2,2);
sigmas2curvesmean=nanmean(sigmas2curves2,2);
sigmascccurvesmean=nanmean(sigmascccurves2,2);


set(handles.checkbox2,'Value',curveincl1(get(handles.slider3,'Value'),1));
set(handles.checkbox3,'Value',curveincl2(get(handles.slider3,'Value'),1));
set(handles.checkbox4,'Value',curveincl1_poly(get(handles.slider3,'Value'),1));
set(handles.checkbox5,'Value',curveincl2_poly(get(handles.slider3,'Value'),1));

    plot(handles.axes2,ItraceIICh1(:,1),ItraceIICh1(:,get(handles.slider3,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceIICh1(1,1) ItraceIICh1(size(ItraceIICh1,1),1) ],[mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1)],'r-')
      ylim(handles.axes2,[min(ItraceIICh1(:,get(handles.slider3,'Value')+1)) max(ItraceIICh1(:,get(handles.slider3,'Value')+1))])
     if curveincl1_poly(get(handles.slider3,'Value'),1)==1
              plot(handles.axes2,ItraceIICh1(:,1),polyval(fitting_poly_GUI1(:,get(handles.slider3,'Value')), ItraceIICh1(:,1)), 'k-')
 
  end
     hold (handles.axes2,'off')
%  
       plot(handles.axes3,ItraceIICh2(:,1),ItraceIICh2(:,get(handles.slider3,'Value')+1))
      hold (handles.axes3,'on')
      plot(handles.axes3,[ ItraceIICh2(1,1) ItraceIICh2(size(ItraceIICh2,1),1) ],[mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1)],'r-')
      ylim(handles.axes3,[min(ItraceIICh2(:,get(handles.slider3,'Value')+1)) max(ItraceIICh2(:,get(handles.slider3,'Value')+1))])
         if curveincl2_poly(get(handles.slider3,'Value'),1)==1
              plot(handles.axes3,ItraceIICh2(:,1),polyval(fitting_poly_GUI2(:,get(handles.slider3,'Value')), ItraceIICh2(:,1)), 'k-')
 
     end
     hold (handles.axes3,'off')

     if curveincl1_poly(get(handles.slider3,'Value'),1)==1
 semilogx(handles.axes1,correlationcurvesCh1_poly(:,1),correlationcurvesCh1_poly(:,get(handles.slider3,'Value')+1),'g+')

     else
     semilogx(handles.axes1,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+')
     end
          hold (handles.axes1,'on')

     
          if curveincl2_poly(get(handles.slider3,'Value'),1)==1
    semilogx(handles.axes1,correlationcurvesCh2_poly(:,1),correlationcurvesCh2_poly(:,get(handles.slider3,'Value')+1),'r+')

          else     
    semilogx(handles.axes1,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+')
          end
          
          
         if  curveincl1_poly(get(handles.slider3,'Value'),1) || curveincl2_poly(get(handles.slider3,'Value'),1)==1
 semilogx(handles.axes1,correlationcurvesChCC_poly(:,1),correlationcurvesChCC_poly(:,get(handles.slider3,'Value')+1),'b+')

         else
         semilogx(handles.axes1,correlationcurvesChCC(:,1),correlationcurvesChCC(:,get(handles.slider3,'Value')+1),'b+')
         end

 semilogx(handles.axes1,correlationcurvesCh1(:,1),meancorrcurve(:,1),'g-',correlationcurvesCh2(:,1),meancorrcurve(:,2),'r-',correlationcurvesChCC(:,1),meancorrcurveCC,'b-')
% 
% 
 xlim(handles.axes1,[min(correlationcurvesCh1(:,1)) max(correlationcurvesCh1(:,1))])
% % For now: Don't plot fits of CFs
% %semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(handles.slider3,'Value')+1),'k-','LineWidth',2)
% % For now: Don't set axis manually
% %axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 
 hold (handles.axes1,'off')
% 
% 
tcorr1segfinalpre=correlationcurvesCh1(:,1);
tcorr2segfinalpre=correlationcurvesCh2(:,1);
tcorrccsegfinalfitpre=correlationcurvesChCC(:,1);
lbseg=tcorr1segfinalpre(1);
ubseg=tcorr1segfinalpre(end);
lbccseg=tcorrccsegfinalfitpre(1);
ubccseg=tcorrccsegfinalfitpre(end);
% 
fcorr1segfinalpre=meancorrcurve(:,1);
fcorr2segfinalpre=meancorrcurve(:,2);
fcrosscorrsegfinalfitpre=meancorrcurveCC;
% 
% approximated, because they do not take into account poly correction
weights1segfinalpre=abs(fcorr1segfinalpre./sigmas1curvesmean);
weights2segfinalpre=abs(fcorr2segfinalpre./sigmas2curvesmean);
weightsccsegfinalfitpre=abs(fcrosscorrsegfinalfitpre./sigmascccurvesmean);

lsautofitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
flscrossfitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5+x(4);

%x0=[100,0.5, S]; % Set in header or GUI manually later!
%fixed=[false false true];
fixedcc=[false true true ];
%flscrossfitfunc=lsautofitfunc;
%flscrossfitfunc=@(y,tcc)1/y(1)*((1+4*y(2)*tcc./y(3)^2).^-0.5).*((1+4*y(2)*tcc./(y(3)*y(4))^2).^-0.5).*exp(-d^2./(y(3)^2+4*y(2).*tcc));
%y0=[5000,100,S,10^-5]; %[N0,tau0,S], Set in header or GUI manually later!
%y0=[2000,0.5,S];

% For now: separate fit. Later: Global fit of all 3 curves (2x auto & cc)
% [N1finalpre,taud1finalpre,Sfit1finalpre,CI1finalpre,fitcurve1finalpre,residuals1finalpre] = autocorrfit2Ddiff(tcorr1fit,fcorr1fit,lb,ub,weights1fit,lsautofitfuncpre,x0,fixed);
% [N2finalpre,taud2finalpre,Sfit2finalpre,CI2finalpre,fitcurve2finalpre,residuals2finalpre] = autocorrfit2Ddiff(tcorr2fit,fcorr2fit,lb,ub,weights2fit,lsautofitfuncpre,x0,fixed);
% %[Nccfinal,tauccfinal,Sfitccfinal,CIccfinal,fitcurveccfinal,residualsccfinal] =autocorrfit2Ddiffconst(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,flscrossfitfunc,y0,fixedcc);
% [Nccfinalpre,tauccfinalpre,Sfitccfinalpre,CIccfinalpre,fitcurveccfinalpre,residualsccfinalpre] =autocorrfit2Ddiff(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,lsautofitfuncpre,y0,fixed);

[N1segfinalpre,taud1segfinalpre,Sfit1segfinalpre,CI1segfinalpre,fitcurve1segfinalpre,residuals1segfinalpre] = autocorrfit2Ddiff(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
[N2segfinalpre,taud2segfinalpre,Sfit2segfinalpre,CI2segfinalpre,fitcurve2segfinalpre,residuals2segfinalpre] = autocorrfit2Ddiff(tcorr2segfinalpre,fcorr2segfinalpre,lbseg,ubseg,weights2segfinalpre,lsautofitfuncpre,x0,fixed);

try

%[Nccsegfinal,tauccsegfinal,Sfitccsegfinal,CIccsegfinal,fitcurveccsegfinal,residualsccsegfinal] =autocorrfit2Ddiffconst(tcorrccsegfinalfit,fcrosscorrsegfinalfit',lbccseg,ubccseg,weightsccsegfinalfit',flscrossfitfunc,y0,fixedcc);

[Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] =autocorrfit2Ddiff(tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,lbccseg,ubccseg,weightsccsegfinalfitpre,lsautofitfuncpre,y0,fixed);
catch
[Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] = autocorrfit2Ddiff(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
%    size(fitcurveccsegfinalpre)
 fitcurveccsegfinalpre  =zeros(size(fitcurveccsegfinalpre,1)+1,size(fitcurveccsegfinalpre,1));

end
% h=figure('OuterPosition',[scrsz(3) 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Full Curve');
% positionvector1=[0.1 0.35 0.8 0.55];
% positionvector2=[0.1 0.1 0.8 0.15];
% subplot('Position',positionvector1),semilogx(tcorr1fit,fitcurve1finalpre,'-g',tcorr2fit,fitcurve2finalpre,'-r',tcorrccfit,fitcurveccfinalpre,'-b')
% hold on
% subplot('Position',positionvector1),semilogx(tcorr1fit,fcorr1fit,'gs',tcorr2fit,fcorr2fit,'rd',tcorrccfit,fcrosscorrfit,'bx')
% xlabel('time')
% ylabel('autocorrelation')

%hh=figure('OuterPosition',[4*scrsz(3)/3 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Segment Averaged Curve');
%positionvector1=[0.1 0.35 0.8 0.55];
%positionvector2=[0.1 0.1 0.8 0.15];
semilogx(handles.axes5,tcorr2segfinalpre,fitcurve2segfinalpre,'-k',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-k')
hold (handles.axes5,'on')
semilogx(handles.axes5,tcorr2segfinalpre,fcorr2segfinalpre,'r-',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'b-')

   if curveincl2_poly(get(handles.slider3,'Value'),1)==0
semilogx(handles.axes5,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+')
hold (handles.axes5,'off')
          else     
semilogx(handles.axes5,correlationcurvesCh2(:,1),correlationcurvesCh2_poly (:,get(handles.slider3,'Value')+1),'r+')
hold (handles.axes5,'off')
   end

semilogx(handles.axes4,tcorr1segfinalpre,fitcurve1segfinalpre,'-k',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-k')
hold (handles.axes4,'on')
semilogx(handles.axes4,tcorr1segfinalpre,fcorr1segfinalpre,'g-',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'b-')


 if curveincl1_poly(get(handles.slider3,'Value'),1)==1
semilogx(handles.axes4,correlationcurvesCh1(:,1),correlationcurvesCh1_poly(:,get(handles.slider3,'Value')+1),'g+')
hold (handles.axes4,'off')
          else    

semilogx(handles.axes4,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+')
 hold (handles.axes4,'off')

 end

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

global x0 y0 fixed curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesChCC  sigmas1curves sigmas2curves sigmascccurves  meancorrcurve meancorrcurveCC ItraceIICh1 ItraceIICh2  correlationcurves2Ch1 correlationcurves2Ch2 correlationcurves2ChCC sigmas1curves2 sigmas2curves2 sigmascccurves2
global correlationcurvesCh1_poly correlationcurvesCh2_poly  correlationcurvesChCC_poly   fitting_poly_GUI2  fitting_poly_GUI1
global  curveincl1_poly curveincl2_poly   
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



curveincl1(get(handles.slider3,'value'),1)=get(hObject,'Value');
set(handles.text5,'String', ['You have selected ' num2str(sum(curveincl1,1)) ' out of ' num2str(size(curveincl1,1)) ' curves'])
set(handles.text6,'String', [num2str(sum(curveincl2,1)) ' out of ' num2str(size(curveincl2,1)) ' curves. CC curves are ' num2str(sum(curveincl1+curveincl2==2)) ])

correlationcurves2Ch1=correlationcurvesCh1;
correlationcurves2Ch2=correlationcurvesCh2;
correlationcurves2ChCC=correlationcurvesChCC;

sigmas1curves2=sigmas1curves;
sigmas2curves2=sigmas2curves;
sigmascccurves2=sigmascccurves;


%approximation for now: if one channes gets the polynomial correction
% the CC is automatically taken as it were from the 2 channels both
% corrected
for kl=2:size(correlationcurvesCh1,2)
    
if curveincl1_poly(kl-1)==1
    correlationcurves2Ch1(:,kl)=correlationcurvesCh1_poly(:,kl);
    correlationcurves2ChCC(:,kl)=correlationcurvesChCC_poly(:,kl); %or will this be calculated?
     
end

if curveincl2_poly(kl-1)==1
    correlationcurves2Ch2(:,kl)=correlationcurvesCh2_poly(:,kl);
    correlationcurves2ChCC(:,kl)=correlationcurvesChCC_poly(:,kl); %or will this be calculated?

end
    
    if curveincl1(kl-1)==0
    correlationcurves2Ch1(:,kl)=NaN;
    correlationcurves2ChCC(:,kl)=NaN;
    sigmas1curves2(:,kl-1)=NaN;
    sigmascccurves2(:,kl-1)=NaN;
    end

    
    
    
if curveincl2(kl-1)==0
correlationcurves2Ch2(:,kl)=NaN;
correlationcurves2ChCC(:,kl)=NaN;
sigmas2curves2(:,kl-1)=NaN;
sigmascccurves2(:,kl-1)=NaN;
end
end



meancorrcurve(:,1)=nanmean(correlationcurves2Ch1(:,2:end),2);
meancorrcurve(:,2)=nanmean(correlationcurves2Ch2(:,2:end),2);
meancorrcurveCC=nanmean(correlationcurves2ChCC(:,2:end),2);


sigmas1curvesmean=nanmean(sigmas1curves2,2);
sigmas2curvesmean=nanmean(sigmas2curves2,2);
sigmascccurvesmean=nanmean(sigmascccurves2,2);


set(handles.checkbox2,'Value',curveincl1(get(handles.slider3,'Value'),1));
set(handles.checkbox3,'Value',curveincl2(get(handles.slider3,'Value'),1));
set(handles.checkbox4,'Value',curveincl1_poly(get(handles.slider3,'Value'),1));
set(handles.checkbox5,'Value',curveincl2_poly(get(handles.slider3,'Value'),1));

    plot(handles.axes2,ItraceIICh1(:,1),ItraceIICh1(:,get(handles.slider3,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceIICh1(1,1) ItraceIICh1(size(ItraceIICh1,1),1) ],[mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1)],'r-')
      ylim(handles.axes2,[min(ItraceIICh1(:,get(handles.slider3,'Value')+1)) max(ItraceIICh1(:,get(handles.slider3,'Value')+1))])
     if curveincl1_poly(get(handles.slider3,'Value'),1)==1
              plot(handles.axes2,ItraceIICh1(:,1),polyval(fitting_poly_GUI1(:,get(handles.slider3,'Value')), ItraceIICh1(:,1)), 'k-')
 
  end
     hold (handles.axes2,'off')
%  
       plot(handles.axes3,ItraceIICh2(:,1),ItraceIICh2(:,get(handles.slider3,'Value')+1))
      hold (handles.axes3,'on')
      plot(handles.axes3,[ ItraceIICh2(1,1) ItraceIICh2(size(ItraceIICh2,1),1) ],[mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1)],'r-')
      ylim(handles.axes3,[min(ItraceIICh2(:,get(handles.slider3,'Value')+1)) max(ItraceIICh2(:,get(handles.slider3,'Value')+1))])
         if curveincl2_poly(get(handles.slider3,'Value'),1)==1
              plot(handles.axes3,ItraceIICh2(:,1),polyval(fitting_poly_GUI2(:,get(handles.slider3,'Value')), ItraceIICh2(:,1)), 'k-')
 
     end
     hold (handles.axes3,'off')

     if curveincl1_poly(get(handles.slider3,'Value'),1)==1
 semilogx(handles.axes1,correlationcurvesCh1_poly(:,1),correlationcurvesCh1_poly(:,get(handles.slider3,'Value')+1),'g+')

     else
     semilogx(handles.axes1,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+')
     end
          hold (handles.axes1,'on')

     
          if curveincl2_poly(get(handles.slider3,'Value'),1)==1
    semilogx(handles.axes1,correlationcurvesCh2_poly(:,1),correlationcurvesCh2_poly(:,get(handles.slider3,'Value')+1),'r+')

          else     
    semilogx(handles.axes1,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+')
          end
          
          
         if  curveincl1_poly(get(handles.slider3,'Value'),1) || curveincl2_poly(get(handles.slider3,'Value'),1)==1
 semilogx(handles.axes1,correlationcurvesChCC_poly(:,1),correlationcurvesChCC_poly(:,get(handles.slider3,'Value')+1),'b+')

         else
         semilogx(handles.axes1,correlationcurvesChCC(:,1),correlationcurvesChCC(:,get(handles.slider3,'Value')+1),'b+')
         end

 semilogx(handles.axes1,correlationcurvesCh1(:,1),meancorrcurve(:,1),'g-',correlationcurvesCh2(:,1),meancorrcurve(:,2),'r-',correlationcurvesChCC(:,1),meancorrcurveCC,'b-')
% 
% 
 xlim(handles.axes1,[min(correlationcurvesCh1(:,1)) max(correlationcurvesCh1(:,1))])
% % For now: Don't plot fits of CFs
% %semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(handles.slider3,'Value')+1),'k-','LineWidth',2)
% % For now: Don't set axis manually
% %axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 
 hold (handles.axes1,'off')
% 
% 
tcorr1segfinalpre=correlationcurvesCh1(:,1);
tcorr2segfinalpre=correlationcurvesCh2(:,1);
tcorrccsegfinalfitpre=correlationcurvesChCC(:,1);
lbseg=tcorr1segfinalpre(1);
ubseg=tcorr1segfinalpre(end);
lbccseg=tcorrccsegfinalfitpre(1);
ubccseg=tcorrccsegfinalfitpre(end);
% 
fcorr1segfinalpre=meancorrcurve(:,1);
fcorr2segfinalpre=meancorrcurve(:,2);
fcrosscorrsegfinalfitpre=meancorrcurveCC;
% 
% approximated, because they do not take into account poly correction
weights1segfinalpre=abs(fcorr1segfinalpre./sigmas1curvesmean);
weights2segfinalpre=abs(fcorr2segfinalpre./sigmas2curvesmean);
weightsccsegfinalfitpre=abs(fcrosscorrsegfinalfitpre./sigmascccurvesmean);

lsautofitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
flscrossfitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5+x(4);

%x0=[100,0.5, S]; % Set in header or GUI manually later!
%fixed=[false false true];
fixedcc=[false true true ];
%flscrossfitfunc=lsautofitfunc;
%flscrossfitfunc=@(y,tcc)1/y(1)*((1+4*y(2)*tcc./y(3)^2).^-0.5).*((1+4*y(2)*tcc./(y(3)*y(4))^2).^-0.5).*exp(-d^2./(y(3)^2+4*y(2).*tcc));
%y0=[5000,100,S,10^-5]; %[N0,tau0,S], Set in header or GUI manually later!
%y0=[2000,0.5,S];

% For now: separate fit. Later: Global fit of all 3 curves (2x auto & cc)
% [N1finalpre,taud1finalpre,Sfit1finalpre,CI1finalpre,fitcurve1finalpre,residuals1finalpre] = autocorrfit2Ddiff(tcorr1fit,fcorr1fit,lb,ub,weights1fit,lsautofitfuncpre,x0,fixed);
% [N2finalpre,taud2finalpre,Sfit2finalpre,CI2finalpre,fitcurve2finalpre,residuals2finalpre] = autocorrfit2Ddiff(tcorr2fit,fcorr2fit,lb,ub,weights2fit,lsautofitfuncpre,x0,fixed);
% %[Nccfinal,tauccfinal,Sfitccfinal,CIccfinal,fitcurveccfinal,residualsccfinal] =autocorrfit2Ddiffconst(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,flscrossfitfunc,y0,fixedcc);
% [Nccfinalpre,tauccfinalpre,Sfitccfinalpre,CIccfinalpre,fitcurveccfinalpre,residualsccfinalpre] =autocorrfit2Ddiff(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,lsautofitfuncpre,y0,fixed);

[N1segfinalpre,taud1segfinalpre,Sfit1segfinalpre,CI1segfinalpre,fitcurve1segfinalpre,residuals1segfinalpre] = autocorrfit2Ddiff(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
[N2segfinalpre,taud2segfinalpre,Sfit2segfinalpre,CI2segfinalpre,fitcurve2segfinalpre,residuals2segfinalpre] = autocorrfit2Ddiff(tcorr2segfinalpre,fcorr2segfinalpre,lbseg,ubseg,weights2segfinalpre,lsautofitfuncpre,x0,fixed);

try

%[Nccsegfinal,tauccsegfinal,Sfitccsegfinal,CIccsegfinal,fitcurveccsegfinal,residualsccsegfinal] =autocorrfit2Ddiffconst(tcorrccsegfinalfit,fcrosscorrsegfinalfit',lbccseg,ubccseg,weightsccsegfinalfit',flscrossfitfunc,y0,fixedcc);

[Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] =autocorrfit2Ddiff(tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,lbccseg,ubccseg,weightsccsegfinalfitpre,lsautofitfuncpre,y0,fixed);
catch
[Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] = autocorrfit2Ddiff(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
%    size(fitcurveccsegfinalpre)
 fitcurveccsegfinalpre  =zeros(size(fitcurveccsegfinalpre,1)+1,size(fitcurveccsegfinalpre,1));

end
% h=figure('OuterPosition',[scrsz(3) 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Full Curve');
% positionvector1=[0.1 0.35 0.8 0.55];
% positionvector2=[0.1 0.1 0.8 0.15];
% subplot('Position',positionvector1),semilogx(tcorr1fit,fitcurve1finalpre,'-g',tcorr2fit,fitcurve2finalpre,'-r',tcorrccfit,fitcurveccfinalpre,'-b')
% hold on
% subplot('Position',positionvector1),semilogx(tcorr1fit,fcorr1fit,'gs',tcorr2fit,fcorr2fit,'rd',tcorrccfit,fcrosscorrfit,'bx')
% xlabel('time')
% ylabel('autocorrelation')

%hh=figure('OuterPosition',[4*scrsz(3)/3 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Segment Averaged Curve');
%positionvector1=[0.1 0.35 0.8 0.55];
%positionvector2=[0.1 0.1 0.8 0.15];
semilogx(handles.axes5,tcorr2segfinalpre,fitcurve2segfinalpre,'-k',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-k')
hold (handles.axes5,'on')
semilogx(handles.axes5,tcorr2segfinalpre,fcorr2segfinalpre,'r-',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'b-')

   if curveincl2_poly(get(handles.slider3,'Value'),1)==0
semilogx(handles.axes5,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+')
hold (handles.axes5,'off')
          else     
semilogx(handles.axes5,correlationcurvesCh2(:,1),correlationcurvesCh2_poly (:,get(handles.slider3,'Value')+1),'r+')
hold (handles.axes5,'off')
   end

semilogx(handles.axes4,tcorr1segfinalpre,fitcurve1segfinalpre,'-k',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-k')
hold (handles.axes4,'on')
semilogx(handles.axes4,tcorr1segfinalpre,fcorr1segfinalpre,'g-',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'b-')

 
 if curveincl1_poly(get(handles.slider3,'Value'),1)==1
semilogx(handles.axes4,correlationcurvesCh1(:,1),correlationcurvesCh1_poly(:,get(handles.slider3,'Value')+1),'g+')
hold (handles.axes4,'off')
          else    

semilogx(handles.axes4,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+')
 hold (handles.axes4,'off')

 end





% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global goon
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
goon=1;
close Corrselection2Channelsindividuellpreview_new


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
global x0 y0 fixed curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesChCC  sigmas1curves sigmas2curves sigmascccurves  meancorrcurve meancorrcurveCC ItraceIICh1 ItraceIICh2  correlationcurves2Ch1 correlationcurves2Ch2 correlationcurves2ChCC sigmas1curves2 sigmas2curves2 sigmascccurves2
global correlationcurvesCh1_poly correlationcurvesCh2_poly  correlationcurvesChCC_poly   fitting_poly_GUI2  fitting_poly_GUI1
global  curveincl1_poly curveincl2_poly   
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



curveincl2(get(handles.slider3,'value'),1)=get(hObject,'Value');
set(handles.text5,'String', ['You have selected ' num2str(sum(curveincl1,1)) ' out of ' num2str(size(curveincl1,1)) ' curves'])
set(handles.text6,'String', [num2str(sum(curveincl2,1)) ' out of ' num2str(size(curveincl2,1)) ' curves. CC curves are ' num2str(sum(curveincl1+curveincl2==2)) ])

correlationcurves2Ch1=correlationcurvesCh1;
correlationcurves2Ch2=correlationcurvesCh2;
correlationcurves2ChCC=correlationcurvesChCC;

sigmas1curves2=sigmas1curves;
sigmas2curves2=sigmas2curves;
sigmascccurves2=sigmascccurves;


%approximation for now: if one channes gets the polynomial correction
% the CC is automatically taken as it were from the 2 channels both
% corrected
for kl=2:size(correlationcurvesCh1,2)
    
if curveincl1_poly(kl-1)==1
    correlationcurves2Ch1(:,kl)=correlationcurvesCh1_poly(:,kl);
    correlationcurves2ChCC(:,kl)=correlationcurvesChCC_poly(:,kl); %or will this be calculated?
     
end

if curveincl2_poly(kl-1)==1
    correlationcurves2Ch2(:,kl)=correlationcurvesCh2_poly(:,kl);
    correlationcurves2ChCC(:,kl)=correlationcurvesChCC_poly(:,kl); %or will this be calculated?

end
    
    if curveincl1(kl-1)==0
    correlationcurves2Ch1(:,kl)=NaN;
    correlationcurves2ChCC(:,kl)=NaN;
    sigmas1curves2(:,kl-1)=NaN;
    sigmascccurves2(:,kl-1)=NaN;
    end

    
    
    
if curveincl2(kl-1)==0
correlationcurves2Ch2(:,kl)=NaN;
correlationcurves2ChCC(:,kl)=NaN;
sigmas2curves2(:,kl-1)=NaN;
sigmascccurves2(:,kl-1)=NaN;
end
end



meancorrcurve(:,1)=nanmean(correlationcurves2Ch1(:,2:end),2);
meancorrcurve(:,2)=nanmean(correlationcurves2Ch2(:,2:end),2);
meancorrcurveCC=nanmean(correlationcurves2ChCC(:,2:end),2);


sigmas1curvesmean=nanmean(sigmas1curves2,2);
sigmas2curvesmean=nanmean(sigmas2curves2,2);
sigmascccurvesmean=nanmean(sigmascccurves2,2);


set(handles.checkbox2,'Value',curveincl1(get(handles.slider3,'Value'),1));
set(handles.checkbox3,'Value',curveincl2(get(handles.slider3,'Value'),1));
set(handles.checkbox4,'Value',curveincl1_poly(get(handles.slider3,'Value'),1));
set(handles.checkbox5,'Value',curveincl2_poly(get(handles.slider3,'Value'),1));

    plot(handles.axes2,ItraceIICh1(:,1),ItraceIICh1(:,get(handles.slider3,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceIICh1(1,1) ItraceIICh1(size(ItraceIICh1,1),1) ],[mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1)],'r-')
      ylim(handles.axes2,[min(ItraceIICh1(:,get(handles.slider3,'Value')+1)) max(ItraceIICh1(:,get(handles.slider3,'Value')+1))])
     if curveincl1_poly(get(handles.slider3,'Value'),1)==1
              plot(handles.axes2,ItraceIICh1(:,1),polyval(fitting_poly_GUI1(:,get(handles.slider3,'Value')), ItraceIICh1(:,1)), 'k-')
 
  end
     hold (handles.axes2,'off')
%  
       plot(handles.axes3,ItraceIICh2(:,1),ItraceIICh2(:,get(handles.slider3,'Value')+1))
      hold (handles.axes3,'on')
      plot(handles.axes3,[ ItraceIICh2(1,1) ItraceIICh2(size(ItraceIICh2,1),1) ],[mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1)],'r-')
      ylim(handles.axes3,[min(ItraceIICh2(:,get(handles.slider3,'Value')+1)) max(ItraceIICh2(:,get(handles.slider3,'Value')+1))])
         if curveincl2_poly(get(handles.slider3,'Value'),1)==1
              plot(handles.axes3,ItraceIICh2(:,1),polyval(fitting_poly_GUI2(:,get(handles.slider3,'Value')), ItraceIICh2(:,1)), 'k-')
 
     end
     hold (handles.axes3,'off')

     if curveincl1_poly(get(handles.slider3,'Value'),1)==1
 semilogx(handles.axes1,correlationcurvesCh1_poly(:,1),correlationcurvesCh1_poly(:,get(handles.slider3,'Value')+1),'g+')

     else
     semilogx(handles.axes1,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+')
     end
          hold (handles.axes1,'on')

     
          if curveincl2_poly(get(handles.slider3,'Value'),1)==1
    semilogx(handles.axes1,correlationcurvesCh2_poly(:,1),correlationcurvesCh2_poly(:,get(handles.slider3,'Value')+1),'r+')

          else     
    semilogx(handles.axes1,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+')
          end
          
          
         if  curveincl1_poly(get(handles.slider3,'Value'),1) || curveincl2_poly(get(handles.slider3,'Value'),1)==1
 semilogx(handles.axes1,correlationcurvesChCC_poly(:,1),correlationcurvesChCC_poly(:,get(handles.slider3,'Value')+1),'b+')

         else
         semilogx(handles.axes1,correlationcurvesChCC(:,1),correlationcurvesChCC(:,get(handles.slider3,'Value')+1),'b+')
         end

 semilogx(handles.axes1,correlationcurvesCh1(:,1),meancorrcurve(:,1),'g-',correlationcurvesCh2(:,1),meancorrcurve(:,2),'r-',correlationcurvesChCC(:,1),meancorrcurveCC,'b-')
% 
% 
 xlim(handles.axes1,[min(correlationcurvesCh1(:,1)) max(correlationcurvesCh1(:,1))])
% % For now: Don't plot fits of CFs
% %semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(handles.slider3,'Value')+1),'k-','LineWidth',2)
% % For now: Don't set axis manually
% %axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 
 hold (handles.axes1,'off')
% 
% 
tcorr1segfinalpre=correlationcurvesCh1(:,1);
tcorr2segfinalpre=correlationcurvesCh2(:,1);
tcorrccsegfinalfitpre=correlationcurvesChCC(:,1);
lbseg=tcorr1segfinalpre(1);
ubseg=tcorr1segfinalpre(end);
lbccseg=tcorrccsegfinalfitpre(1);
ubccseg=tcorrccsegfinalfitpre(end);
% 
fcorr1segfinalpre=meancorrcurve(:,1);
fcorr2segfinalpre=meancorrcurve(:,2);
fcrosscorrsegfinalfitpre=meancorrcurveCC;
% 
% approximated, because they do not take into account poly correction
weights1segfinalpre=abs(fcorr1segfinalpre./sigmas1curvesmean);
weights2segfinalpre=abs(fcorr2segfinalpre./sigmas2curvesmean);
weightsccsegfinalfitpre=abs(fcrosscorrsegfinalfitpre./sigmascccurvesmean);

lsautofitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
flscrossfitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5+x(4);

%x0=[100,0.5, S]; % Set in header or GUI manually later!
%fixed=[false false true];
fixedcc=[false true true ];
%flscrossfitfunc=lsautofitfunc;
%flscrossfitfunc=@(y,tcc)1/y(1)*((1+4*y(2)*tcc./y(3)^2).^-0.5).*((1+4*y(2)*tcc./(y(3)*y(4))^2).^-0.5).*exp(-d^2./(y(3)^2+4*y(2).*tcc));
%y0=[5000,100,S,10^-5]; %[N0,tau0,S], Set in header or GUI manually later!
%y0=[2000,0.5,S];

% For now: separate fit. Later: Global fit of all 3 curves (2x auto & cc)
% [N1finalpre,taud1finalpre,Sfit1finalpre,CI1finalpre,fitcurve1finalpre,residuals1finalpre] = autocorrfit2Ddiff(tcorr1fit,fcorr1fit,lb,ub,weights1fit,lsautofitfuncpre,x0,fixed);
% [N2finalpre,taud2finalpre,Sfit2finalpre,CI2finalpre,fitcurve2finalpre,residuals2finalpre] = autocorrfit2Ddiff(tcorr2fit,fcorr2fit,lb,ub,weights2fit,lsautofitfuncpre,x0,fixed);
% %[Nccfinal,tauccfinal,Sfitccfinal,CIccfinal,fitcurveccfinal,residualsccfinal] =autocorrfit2Ddiffconst(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,flscrossfitfunc,y0,fixedcc);
% [Nccfinalpre,tauccfinalpre,Sfitccfinalpre,CIccfinalpre,fitcurveccfinalpre,residualsccfinalpre] =autocorrfit2Ddiff(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,lsautofitfuncpre,y0,fixed);
try
[N1segfinalpre,taud1segfinalpre,Sfit1segfinalpre,CI1segfinalpre,fitcurve1segfinalpre,residuals1segfinalpre] = autocorrfit2Ddiff(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
[N2segfinalpre,taud2segfinalpre,Sfit2segfinalpre,CI2segfinalpre,fitcurve2segfinalpre,residuals2segfinalpre] = autocorrfit2Ddiff(tcorr2segfinalpre,fcorr2segfinalpre,lbseg,ubseg,weights2segfinalpre,lsautofitfuncpre,x0,fixed);
end

try

%[Nccsegfinal,tauccsegfinal,Sfitccsegfinal,CIccsegfinal,fitcurveccsegfinal,residualsccsegfinal] =autocorrfit2Ddiffconst(tcorrccsegfinalfit,fcrosscorrsegfinalfit',lbccseg,ubccseg,weightsccsegfinalfit',flscrossfitfunc,y0,fixedcc);

[Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] =autocorrfit2Ddiff(tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,lbccseg,ubccseg,weightsccsegfinalfitpre,lsautofitfuncpre,y0,fixed);
catch
[Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] = autocorrfit2Ddiff(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
%    size(fitcurveccsegfinalpre)
 fitcurveccsegfinalpre  =zeros(size(fitcurveccsegfinalpre,1)+1,size(fitcurveccsegfinalpre,1));

end
% h=figure('OuterPosition',[scrsz(3) 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Full Curve');
% positionvector1=[0.1 0.35 0.8 0.55];
% positionvector2=[0.1 0.1 0.8 0.15];
% subplot('Position',positionvector1),semilogx(tcorr1fit,fitcurve1finalpre,'-g',tcorr2fit,fitcurve2finalpre,'-r',tcorrccfit,fitcurveccfinalpre,'-b')
% hold on
% subplot('Position',positionvector1),semilogx(tcorr1fit,fcorr1fit,'gs',tcorr2fit,fcorr2fit,'rd',tcorrccfit,fcrosscorrfit,'bx')
% xlabel('time')
% ylabel('autocorrelation')

%hh=figure('OuterPosition',[4*scrsz(3)/3 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Segment Averaged Curve');
%positionvector1=[0.1 0.35 0.8 0.55];
%positionvector2=[0.1 0.1 0.8 0.15];
try
semilogx(handles.axes5,tcorr2segfinalpre,fitcurve2segfinalpre,'-k',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-k')
end
hold (handles.axes5,'on')
semilogx(handles.axes5,tcorr2segfinalpre,fcorr2segfinalpre,'r-',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'b-')

   if curveincl2_poly(get(handles.slider3,'Value'),1)==0
semilogx(handles.axes5,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+')
hold (handles.axes5,'off')
          else     
semilogx(handles.axes5,correlationcurvesCh2(:,1),correlationcurvesCh2_poly (:,get(handles.slider3,'Value')+1),'r+')
hold (handles.axes5,'off')
   end
try
semilogx(handles.axes4,tcorr1segfinalpre,fitcurve1segfinalpre,'-k',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-k')
end
hold (handles.axes4,'on')
semilogx(handles.axes4,tcorr1segfinalpre,fcorr1segfinalpre,'g-',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'b-')

 
 if curveincl1_poly(get(handles.slider3,'Value'),1)==1
semilogx(handles.axes4,correlationcurvesCh1(:,1),correlationcurvesCh1_poly(:,get(handles.slider3,'Value')+1),'g+')
hold (handles.axes4,'off')
          else    

semilogx(handles.axes4,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+')
 hold (handles.axes4,'off')

 end




% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global x0 y0 fixed curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesChCC  sigmas1curves sigmas2curves sigmascccurves  meancorrcurve meancorrcurveCC ItraceIICh1 ItraceIICh2  correlationcurves2Ch1 correlationcurves2Ch2 correlationcurves2ChCC sigmas1curves2 sigmas2curves2 sigmascccurves2
global correlationcurvesCh1_poly correlationcurvesCh2_poly  correlationcurvesChCC_poly   fitting_poly_GUI2  fitting_poly_GUI1
global  curveincl1_poly curveincl2_poly   
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



curveincl1_poly(get(handles.slider3,'value'),1)=get(hObject,'Value');
set(handles.text5,'String', ['You have selected ' num2str(sum(curveincl1,1)) ' out of ' num2str(size(curveincl1,1)) ' curves'])
set(handles.text6,'String', [num2str(sum(curveincl2,1)) ' out of ' num2str(size(curveincl2,1)) ' curves. CC curves are ' num2str(sum(curveincl1+curveincl2==2)) ])

correlationcurves2Ch1=correlationcurvesCh1;
correlationcurves2Ch2=correlationcurvesCh2;
correlationcurves2ChCC=correlationcurvesChCC;

sigmas1curves2=sigmas1curves;
sigmas2curves2=sigmas2curves;
sigmascccurves2=sigmascccurves;


%approximation for now: if one channes gets the polynomial correction
% the CC is automatically taken as it were from the 2 channels both
% corrected
for kl=2:size(correlationcurvesCh1,2)
    
if curveincl1_poly(kl-1)==1
    correlationcurves2Ch1(:,kl)=correlationcurvesCh1_poly(:,kl);
    correlationcurves2ChCC(:,kl)=correlationcurvesChCC_poly(:,kl); %or will this be calculated?
     
end

if curveincl2_poly(kl-1)==1
    correlationcurves2Ch2(:,kl)=correlationcurvesCh2_poly(:,kl);
    correlationcurves2ChCC(:,kl)=correlationcurvesChCC_poly(:,kl); %or will this be calculated?

end
    
    if curveincl1(kl-1)==0
    correlationcurves2Ch1(:,kl)=NaN;
    correlationcurves2ChCC(:,kl)=NaN;
    sigmas1curves2(:,kl-1)=NaN;
    sigmascccurves2(:,kl-1)=NaN;
    end

    
    
    
if curveincl2(kl-1)==0
correlationcurves2Ch2(:,kl)=NaN;
correlationcurves2ChCC(:,kl)=NaN;
sigmas2curves2(:,kl-1)=NaN;
sigmascccurves2(:,kl-1)=NaN;
end
end



meancorrcurve(:,1)=nanmean(correlationcurves2Ch1(:,2:end),2);
meancorrcurve(:,2)=nanmean(correlationcurves2Ch2(:,2:end),2);
meancorrcurveCC=nanmean(correlationcurves2ChCC(:,2:end),2);


sigmas1curvesmean=nanmean(sigmas1curves2,2);
sigmas2curvesmean=nanmean(sigmas2curves2,2);
sigmascccurvesmean=nanmean(sigmascccurves2,2);


set(handles.checkbox2,'Value',curveincl1(get(handles.slider3,'Value'),1));
set(handles.checkbox3,'Value',curveincl2(get(handles.slider3,'Value'),1));
set(handles.checkbox4,'Value',curveincl1_poly(get(handles.slider3,'Value'),1));
set(handles.checkbox5,'Value',curveincl2_poly(get(handles.slider3,'Value'),1));

    plot(handles.axes2,ItraceIICh1(:,1),ItraceIICh1(:,get(handles.slider3,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceIICh1(1,1) ItraceIICh1(size(ItraceIICh1,1),1) ],[mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1)],'r-')
      ylim(handles.axes2,[min(ItraceIICh1(:,get(handles.slider3,'Value')+1)) max(ItraceIICh1(:,get(handles.slider3,'Value')+1))])
     if curveincl1_poly(get(handles.slider3,'Value'),1)==1
              plot(handles.axes2,ItraceIICh1(:,1),polyval(fitting_poly_GUI1(:,get(handles.slider3,'Value')), ItraceIICh1(:,1)), 'k-')
 
  end
     hold (handles.axes2,'off')
%  
       plot(handles.axes3,ItraceIICh2(:,1),ItraceIICh2(:,get(handles.slider3,'Value')+1))
      hold (handles.axes3,'on')
      plot(handles.axes3,[ ItraceIICh2(1,1) ItraceIICh2(size(ItraceIICh2,1),1) ],[mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1)],'r-')
      ylim(handles.axes3,[min(ItraceIICh2(:,get(handles.slider3,'Value')+1)) max(ItraceIICh2(:,get(handles.slider3,'Value')+1))])
         if curveincl2_poly(get(handles.slider3,'Value'),1)==1
              plot(handles.axes3,ItraceIICh2(:,1),polyval(fitting_poly_GUI2(:,get(handles.slider3,'Value')), ItraceIICh2(:,1)), 'k-')
 
     end
     hold (handles.axes3,'off')

     if curveincl1_poly(get(handles.slider3,'Value'),1)==1
 semilogx(handles.axes1,correlationcurvesCh1_poly(:,1),correlationcurvesCh1_poly(:,get(handles.slider3,'Value')+1),'g+')

     else
     semilogx(handles.axes1,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+')
     end
          hold (handles.axes1,'on')

     
          if curveincl2_poly(get(handles.slider3,'Value'),1)==1
    semilogx(handles.axes1,correlationcurvesCh2_poly(:,1),correlationcurvesCh2_poly(:,get(handles.slider3,'Value')+1),'r+')

          else     
    semilogx(handles.axes1,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+')
          end
          
          
         if  curveincl1_poly(get(handles.slider3,'Value'),1) || curveincl2_poly(get(handles.slider3,'Value'),1)==1
 semilogx(handles.axes1,correlationcurvesChCC_poly(:,1),correlationcurvesChCC_poly(:,get(handles.slider3,'Value')+1),'b+')

         else
         semilogx(handles.axes1,correlationcurvesChCC(:,1),correlationcurvesChCC(:,get(handles.slider3,'Value')+1),'b+')
         end

 semilogx(handles.axes1,correlationcurvesCh1(:,1),meancorrcurve(:,1),'g-',correlationcurvesCh2(:,1),meancorrcurve(:,2),'r-',correlationcurvesChCC(:,1),meancorrcurveCC,'b-')
% 
% 
 xlim(handles.axes1,[min(correlationcurvesCh1(:,1)) max(correlationcurvesCh1(:,1))])
% % For now: Don't plot fits of CFs
% %semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(handles.slider3,'Value')+1),'k-','LineWidth',2)
% % For now: Don't set axis manually
% %axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 
 hold (handles.axes1,'off')
% 
% 
tcorr1segfinalpre=correlationcurvesCh1(:,1);
tcorr2segfinalpre=correlationcurvesCh2(:,1);
tcorrccsegfinalfitpre=correlationcurvesChCC(:,1);
lbseg=tcorr1segfinalpre(1);
ubseg=tcorr1segfinalpre(end);
lbccseg=tcorrccsegfinalfitpre(1);
ubccseg=tcorrccsegfinalfitpre(end);
% 
fcorr1segfinalpre=meancorrcurve(:,1);
fcorr2segfinalpre=meancorrcurve(:,2);
fcrosscorrsegfinalfitpre=meancorrcurveCC;
% 
% approximated, because they do not take into account poly correction
weights1segfinalpre=abs(fcorr1segfinalpre./sigmas1curvesmean);
weights2segfinalpre=abs(fcorr2segfinalpre./sigmas2curvesmean);
weightsccsegfinalfitpre=abs(fcrosscorrsegfinalfitpre./sigmascccurvesmean);

lsautofitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
flscrossfitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5+x(4);

%x0=[100,0.5, S]; % Set in header or GUI manually later!
%fixed=[false false true];
fixedcc=[false true true ];
%flscrossfitfunc=lsautofitfunc;
%flscrossfitfunc=@(y,tcc)1/y(1)*((1+4*y(2)*tcc./y(3)^2).^-0.5).*((1+4*y(2)*tcc./(y(3)*y(4))^2).^-0.5).*exp(-d^2./(y(3)^2+4*y(2).*tcc));
%y0=[5000,100,S,10^-5]; %[N0,tau0,S], Set in header or GUI manually later!
%y0=[2000,0.5,S];

% For now: separate fit. Later: Global fit of all 3 curves (2x auto & cc)
% [N1finalpre,taud1finalpre,Sfit1finalpre,CI1finalpre,fitcurve1finalpre,residuals1finalpre] = autocorrfit2Ddiff(tcorr1fit,fcorr1fit,lb,ub,weights1fit,lsautofitfuncpre,x0,fixed);
% [N2finalpre,taud2finalpre,Sfit2finalpre,CI2finalpre,fitcurve2finalpre,residuals2finalpre] = autocorrfit2Ddiff(tcorr2fit,fcorr2fit,lb,ub,weights2fit,lsautofitfuncpre,x0,fixed);
% %[Nccfinal,tauccfinal,Sfitccfinal,CIccfinal,fitcurveccfinal,residualsccfinal] =autocorrfit2Ddiffconst(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,flscrossfitfunc,y0,fixedcc);
% [Nccfinalpre,tauccfinalpre,Sfitccfinalpre,CIccfinalpre,fitcurveccfinalpre,residualsccfinalpre] =autocorrfit2Ddiff(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,lsautofitfuncpre,y0,fixed);
try
[N1segfinalpre,taud1segfinalpre,Sfit1segfinalpre,CI1segfinalpre,fitcurve1segfinalpre,residuals1segfinalpre] = autocorrfit2Ddiff(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
[N2segfinalpre,taud2segfinalpre,Sfit2segfinalpre,CI2segfinalpre,fitcurve2segfinalpre,residuals2segfinalpre] = autocorrfit2Ddiff(tcorr2segfinalpre,fcorr2segfinalpre,lbseg,ubseg,weights2segfinalpre,lsautofitfuncpre,x0,fixed);
end
try

%[Nccsegfinal,tauccsegfinal,Sfitccsegfinal,CIccsegfinal,fitcurveccsegfinal,residualsccsegfinal] =autocorrfit2Ddiffconst(tcorrccsegfinalfit,fcrosscorrsegfinalfit',lbccseg,ubccseg,weightsccsegfinalfit',flscrossfitfunc,y0,fixedcc);

[Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] =autocorrfit2Ddiff(tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,lbccseg,ubccseg,weightsccsegfinalfitpre,lsautofitfuncpre,y0,fixed);
catch
[Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] = autocorrfit2Ddiff(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
%    size(fitcurveccsegfinalpre)
 fitcurveccsegfinalpre  =zeros(size(fitcurveccsegfinalpre,1)+1,size(fitcurveccsegfinalpre,1));

end
% h=figure('OuterPosition',[scrsz(3) 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Full Curve');
% positionvector1=[0.1 0.35 0.8 0.55];
% positionvector2=[0.1 0.1 0.8 0.15];
% subplot('Position',positionvector1),semilogx(tcorr1fit,fitcurve1finalpre,'-g',tcorr2fit,fitcurve2finalpre,'-r',tcorrccfit,fitcurveccfinalpre,'-b')
% hold on
% subplot('Position',positionvector1),semilogx(tcorr1fit,fcorr1fit,'gs',tcorr2fit,fcorr2fit,'rd',tcorrccfit,fcrosscorrfit,'bx')
% xlabel('time')
% ylabel('autocorrelation')

%hh=figure('OuterPosition',[4*scrsz(3)/3 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Segment Averaged Curve');
%positionvector1=[0.1 0.35 0.8 0.55];
%positionvector2=[0.1 0.1 0.8 0.15];
try
semilogx(handles.axes5,tcorr2segfinalpre,fitcurve2segfinalpre,'-k',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-k')
end
hold (handles.axes5,'on')
semilogx(handles.axes5,tcorr2segfinalpre,fcorr2segfinalpre,'r-',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'b-')

   if curveincl2_poly(get(handles.slider3,'Value'),1)==0
semilogx(handles.axes5,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+')
hold (handles.axes5,'off')
          else     
semilogx(handles.axes5,correlationcurvesCh2(:,1),correlationcurvesCh2_poly (:,get(handles.slider3,'Value')+1),'r+')
hold (handles.axes5,'off')
   end
try
semilogx(handles.axes4,tcorr1segfinalpre,fitcurve1segfinalpre,'-k',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-k')
end
hold (handles.axes4,'on')
semilogx(handles.axes4,tcorr1segfinalpre,fcorr1segfinalpre,'g-',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'b-')

 
 if curveincl1_poly(get(handles.slider3,'Value'),1)==0
semilogx(handles.axes4,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+')
hold (handles.axes4,'off')
          else    

semilogx(handles.axes4,correlationcurvesCh1(:,1),correlationcurvesCh1_poly(:,get(handles.slider3,'Value')+1),'g+')
 hold (handles.axes4,'off')

 end

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global x0 y0 fixed curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesChCC  sigmas1curves sigmas2curves sigmascccurves  meancorrcurve meancorrcurveCC ItraceIICh1 ItraceIICh2  correlationcurves2Ch1 correlationcurves2Ch2 correlationcurves2ChCC sigmas1curves2 sigmas2curves2 sigmascccurves2
global correlationcurvesCh1_poly correlationcurvesCh2_poly  correlationcurvesChCC_poly   fitting_poly_GUI2  fitting_poly_GUI1
global  curveincl1_poly curveincl2_poly   
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



curveincl2_poly(get(handles.slider3,'value'),1)=get(hObject,'Value');
set(handles.text5,'String', ['You have selected ' num2str(sum(curveincl1,1)) ' out of ' num2str(size(curveincl1,1)) ' curves'])
set(handles.text6,'String', [num2str(sum(curveincl2,1)) ' out of ' num2str(size(curveincl2,1)) ' curves. CC curves are ' num2str(sum(curveincl1+curveincl2==2)) ])

correlationcurves2Ch1=correlationcurvesCh1;
correlationcurves2Ch2=correlationcurvesCh2;
correlationcurves2ChCC=correlationcurvesChCC;

sigmas1curves2=sigmas1curves;
sigmas2curves2=sigmas2curves;
sigmascccurves2=sigmascccurves;


%approximation for now: if one channes gets the polynomial correction
% the CC is automatically taken as it were from the 2 channels both
% corrected
for kl=2:size(correlationcurvesCh1,2)
    
if curveincl1_poly(kl-1)==1
    correlationcurves2Ch1(:,kl)=correlationcurvesCh1_poly(:,kl);
    correlationcurves2ChCC(:,kl)=correlationcurvesChCC_poly(:,kl); %or will this be calculated?
     
end

if curveincl2_poly(kl-1)==1
    correlationcurves2Ch2(:,kl)=correlationcurvesCh2_poly(:,kl);
    correlationcurves2ChCC(:,kl)=correlationcurvesChCC_poly(:,kl); %or will this be calculated?

end
    
    if curveincl1(kl-1)==0
    correlationcurves2Ch1(:,kl)=NaN;
    correlationcurves2ChCC(:,kl)=NaN;
    sigmas1curves2(:,kl-1)=NaN;
    sigmascccurves2(:,kl-1)=NaN;
    end

    
    
    
if curveincl2(kl-1)==0
correlationcurves2Ch2(:,kl)=NaN;
correlationcurves2ChCC(:,kl)=NaN;
sigmas2curves2(:,kl-1)=NaN;
sigmascccurves2(:,kl-1)=NaN;
end
end



meancorrcurve(:,1)=nanmean(correlationcurves2Ch1(:,2:end),2);
meancorrcurve(:,2)=nanmean(correlationcurves2Ch2(:,2:end),2);
meancorrcurveCC=nanmean(correlationcurves2ChCC(:,2:end),2);


sigmas1curvesmean=nanmean(sigmas1curves2,2);
sigmas2curvesmean=nanmean(sigmas2curves2,2);
sigmascccurvesmean=nanmean(sigmascccurves2,2);


set(handles.checkbox2,'Value',curveincl1(get(handles.slider3,'Value'),1));
set(handles.checkbox3,'Value',curveincl2(get(handles.slider3,'Value'),1));
set(handles.checkbox4,'Value',curveincl1_poly(get(handles.slider3,'Value'),1));
set(handles.checkbox5,'Value',curveincl2_poly(get(handles.slider3,'Value'),1));

    plot(handles.axes2,ItraceIICh1(:,1),ItraceIICh1(:,get(handles.slider3,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceIICh1(1,1) ItraceIICh1(size(ItraceIICh1,1),1) ],[mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1)],'r-')
      ylim(handles.axes2,[min(ItraceIICh1(:,get(handles.slider3,'Value')+1)) max(ItraceIICh1(:,get(handles.slider3,'Value')+1))])
     if curveincl1_poly(get(handles.slider3,'Value'),1)==1
              plot(handles.axes2,ItraceIICh1(:,1),polyval(fitting_poly_GUI1(:,get(handles.slider3,'Value')), ItraceIICh1(:,1)), 'k-')
 
  end
     hold (handles.axes2,'off')
%  
       plot(handles.axes3,ItraceIICh2(:,1),ItraceIICh2(:,get(handles.slider3,'Value')+1))
      hold (handles.axes3,'on')
      plot(handles.axes3,[ ItraceIICh2(1,1) ItraceIICh2(size(ItraceIICh2,1),1) ],[mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1)],'r-')
      ylim(handles.axes3,[min(ItraceIICh2(:,get(handles.slider3,'Value')+1)) max(ItraceIICh2(:,get(handles.slider3,'Value')+1))])
         if curveincl2_poly(get(handles.slider3,'Value'),1)==1
              plot(handles.axes3,ItraceIICh2(:,1),polyval(fitting_poly_GUI2(:,get(handles.slider3,'Value')), ItraceIICh2(:,1)), 'k-')
 
     end
     hold (handles.axes3,'off')

     if curveincl1_poly(get(handles.slider3,'Value'),1)==1
 semilogx(handles.axes1,correlationcurvesCh1_poly(:,1),correlationcurvesCh1_poly(:,get(handles.slider3,'Value')+1),'g+')

     else
     semilogx(handles.axes1,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+')
     end
          hold (handles.axes1,'on')

     
          if curveincl2_poly(get(handles.slider3,'Value'),1)==1
    semilogx(handles.axes1,correlationcurvesCh2_poly(:,1),correlationcurvesCh2_poly(:,get(handles.slider3,'Value')+1),'r+')

          else     
    semilogx(handles.axes1,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+')
          end
          
          
         if  curveincl1_poly(get(handles.slider3,'Value'),1) || curveincl2_poly(get(handles.slider3,'Value'),1)==1
 semilogx(handles.axes1,correlationcurvesChCC_poly(:,1),correlationcurvesChCC_poly(:,get(handles.slider3,'Value')+1),'b+')

         else
         semilogx(handles.axes1,correlationcurvesChCC(:,1),correlationcurvesChCC(:,get(handles.slider3,'Value')+1),'b+')
         end

 semilogx(handles.axes1,correlationcurvesCh1(:,1),meancorrcurve(:,1),'g-',correlationcurvesCh2(:,1),meancorrcurve(:,2),'r-',correlationcurvesChCC(:,1),meancorrcurveCC,'b-')
% 
% 
 xlim(handles.axes1,[min(correlationcurvesCh1(:,1)) max(correlationcurvesCh1(:,1))])
% % For now: Don't plot fits of CFs
% %semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(handles.slider3,'Value')+1),'k-','LineWidth',2)
% % For now: Don't set axis manually
% %axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 
 hold (handles.axes1,'off')
% 
% 
tcorr1segfinalpre=correlationcurvesCh1(:,1);
tcorr2segfinalpre=correlationcurvesCh2(:,1);
tcorrccsegfinalfitpre=correlationcurvesChCC(:,1);
lbseg=tcorr1segfinalpre(1);
ubseg=tcorr1segfinalpre(end);
lbccseg=tcorrccsegfinalfitpre(1);
ubccseg=tcorrccsegfinalfitpre(end);
% 
fcorr1segfinalpre=meancorrcurve(:,1);
fcorr2segfinalpre=meancorrcurve(:,2);
fcrosscorrsegfinalfitpre=meancorrcurveCC;
% 
% approximated, because they do not take into account poly correction
weights1segfinalpre=abs(fcorr1segfinalpre./sigmas1curvesmean);
weights2segfinalpre=abs(fcorr2segfinalpre./sigmas2curvesmean);
weightsccsegfinalfitpre=abs(fcrosscorrsegfinalfitpre./sigmascccurvesmean);

lsautofitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
flscrossfitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5+x(4);

%x0=[100,0.5, S]; % Set in header or GUI manually later!
%fixed=[false false true];
fixedcc=[false true true ];
%flscrossfitfunc=lsautofitfunc;
%flscrossfitfunc=@(y,tcc)1/y(1)*((1+4*y(2)*tcc./y(3)^2).^-0.5).*((1+4*y(2)*tcc./(y(3)*y(4))^2).^-0.5).*exp(-d^2./(y(3)^2+4*y(2).*tcc));
%y0=[5000,100,S,10^-5]; %[N0,tau0,S], Set in header or GUI manually later!
%y0=[2000,0.5,S];

% For now: separate fit. Later: Global fit of all 3 curves (2x auto & cc)
% [N1finalpre,taud1finalpre,Sfit1finalpre,CI1finalpre,fitcurve1finalpre,residuals1finalpre] = autocorrfit2Ddiff(tcorr1fit,fcorr1fit,lb,ub,weights1fit,lsautofitfuncpre,x0,fixed);
% [N2finalpre,taud2finalpre,Sfit2finalpre,CI2finalpre,fitcurve2finalpre,residuals2finalpre] = autocorrfit2Ddiff(tcorr2fit,fcorr2fit,lb,ub,weights2fit,lsautofitfuncpre,x0,fixed);
% %[Nccfinal,tauccfinal,Sfitccfinal,CIccfinal,fitcurveccfinal,residualsccfinal] =autocorrfit2Ddiffconst(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,flscrossfitfunc,y0,fixedcc);
% [Nccfinalpre,tauccfinalpre,Sfitccfinalpre,CIccfinalpre,fitcurveccfinalpre,residualsccfinalpre] =autocorrfit2Ddiff(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,lsautofitfuncpre,y0,fixed);
try
[N1segfinalpre,taud1segfinalpre,Sfit1segfinalpre,CI1segfinalpre,fitcurve1segfinalpre,residuals1segfinalpre] = autocorrfit2Ddiff(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
[N2segfinalpre,taud2segfinalpre,Sfit2segfinalpre,CI2segfinalpre,fitcurve2segfinalpre,residuals2segfinalpre] = autocorrfit2Ddiff(tcorr2segfinalpre,fcorr2segfinalpre,lbseg,ubseg,weights2segfinalpre,lsautofitfuncpre,x0,fixed);
end
try

%[Nccsegfinal,tauccsegfinal,Sfitccsegfinal,CIccsegfinal,fitcurveccsegfinal,residualsccsegfinal] =autocorrfit2Ddiffconst(tcorrccsegfinalfit,fcrosscorrsegfinalfit',lbccseg,ubccseg,weightsccsegfinalfit',flscrossfitfunc,y0,fixedcc);

[Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] =autocorrfit2Ddiff(tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,lbccseg,ubccseg,weightsccsegfinalfitpre,lsautofitfuncpre,y0,fixed);
catch
[Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] = autocorrfit2Ddiff(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
%    size(fitcurveccsegfinalpre)
 fitcurveccsegfinalpre  =zeros(size(fitcurveccsegfinalpre,1)+1,size(fitcurveccsegfinalpre,1));

end
% h=figure('OuterPosition',[scrsz(3) 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Full Curve');
% positionvector1=[0.1 0.35 0.8 0.55];
% positionvector2=[0.1 0.1 0.8 0.15];
% subplot('Position',positionvector1),semilogx(tcorr1fit,fitcurve1finalpre,'-g',tcorr2fit,fitcurve2finalpre,'-r',tcorrccfit,fitcurveccfinalpre,'-b')
% hold on
% subplot('Position',positionvector1),semilogx(tcorr1fit,fcorr1fit,'gs',tcorr2fit,fcorr2fit,'rd',tcorrccfit,fcrosscorrfit,'bx')
% xlabel('time')
% ylabel('autocorrelation')

%hh=figure('OuterPosition',[4*scrsz(3)/3 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Segment Averaged Curve');
%positionvector1=[0.1 0.35 0.8 0.55];
%positionvector2=[0.1 0.1 0.8 0.15];
try
semilogx(handles.axes5,tcorr2segfinalpre,fitcurve2segfinalpre,'-k',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-k')
end
hold (handles.axes5,'on')
semilogx(handles.axes5,tcorr2segfinalpre,fcorr2segfinalpre,'r-',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'b-')

   if curveincl2_poly(get(handles.slider3,'Value'),1)==0
semilogx(handles.axes5,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+')
hold (handles.axes5,'off')
          else     
semilogx(handles.axes5,correlationcurvesCh2(:,1),correlationcurvesCh2_poly (:,get(handles.slider3,'Value')+1),'r+')
hold (handles.axes5,'off')
   end
try
semilogx(handles.axes4,tcorr1segfinalpre,fitcurve1segfinalpre,'-k',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-k')
end
hold (handles.axes4,'on')
semilogx(handles.axes4,tcorr1segfinalpre,fcorr1segfinalpre,'g-',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'b-')


 if curveincl1_poly(get(handles.slider3,'Value'),1)==1
semilogx(handles.axes4,correlationcurvesCh1(:,1),correlationcurvesCh1_poly(:,get(handles.slider3,'Value')+1),'g+')
hold (handles.axes4,'off')
          else    

semilogx(handles.axes4,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+')
 hold (handles.axes4,'off')

 end

% Hint: get(hObject,'Value') returns toggle state of checkbox5
