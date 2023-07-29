function [N,taud,Sfit,CI,fitcurve,residuals] = autocorrfit2Ddiff(tcorr,fcorr,lb,ub,weights,fitfunc,x0,fixed)
[fcorrfitparameters,residuals,J,COVB,MSE] = nlinfitsome(fixed,tcorr,fcorr,fitfunc,x0,'Weights',weights);
%[fcorrfitparameters,fcorrfitresiduals] = lsqcurvefit(lsfitfunc,x0,tcorrfit,fcorrfit,lb,ub,'Weights',weights);
N=fcorrfitparameters(1);
taud=fcorrfitparameters(2);
Sfit=fcorrfitparameters(3);
% Compute confidence intervals of coefficients
%CI=0;
fitparametersnotfixed=fcorrfitparameters(fixed==false);
CInotfixed = nlparci(fitparametersnotfixed,residuals,'covar',COVB);
CI=zeros(length(fcorrfitparameters),2);
k=1;
for i=1:length(fcorrfitparameters)
    if fixed(i)==false
        CI(i,:)=CInotfixed(k,:);
        k=k+1;
    else
        CI(i,:)=[0 0];
    end
end    
fitcurve=1/N*((1+tcorr./taud).^-0.5).*(1+tcorr./(taud*Sfit^2)).^-0.5;

