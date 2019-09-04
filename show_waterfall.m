function show_waterfall(Fr,chemshift,cs_range,Nspectra,xshift,yshift,pyr_index,lac_index,bic_index,ala_index)
% show_waterfall(Fr,chemshift,cs_range,Nspectra,xshift,yshift,pyr_index,lac_index,bic_index,ala_index) 
%gets real part of spectra Fr and plots the waterplot plot.
% 
% Input parameters:
%
% Fr: real part of spectrum
% chemshift: vector of chemical shifts.
% cs_range: vector range of chemical shifts to be displayed. [start end]
% Nspectra: Number of spectra to show.
% xshift: shift in chemical shift on the x axis.
% ysfhit: shift on psectra in y direction. 
%
% optional input parameters:
% pyr_index....: the index of peaks to make colored. 
%
% No output arguments.   
% 
% Created by Mehrdad pourfathi on 7/10/2013.
%
% Updated by MP on 1/26/2015.
%
delta_cs = chemshift(1)-chemshift(2);
xshift = round(xshift/delta_cs);
cs_win = intersect(find(chemshift>cs_range(1)),find(chemshift<cs_range(2)));   % C1 injection
% norm= max((max(Fr(pyr_index,:))));
norm = max((max(Fr)));
yshift = yshift * norm;
for ii = 1:length(Nspectra)-1
    x = cs_win + (ii-1) * xshift;
    signal = sum(Fr(:,[Nspectra(ii);Nspectra(ii+1)-1]),2);
    if nargin > 6
        plot(chemshift(x),signal(cs_win)+(ii-1)*yshift,'k'); hold on;
        pyr = plot(chemshift(pyr_index + (ii-1) * xshift),signal(pyr_index)+(ii-1)*yshift,'b'); hold on;
        lac = plot(chemshift(lac_index + (ii-1) * xshift),signal(lac_index)+(ii-1)*yshift,'r'); hold on;
        ala = plot(chemshift(ala_index + (ii-1) * xshift),signal(ala_index)+(ii-1)*yshift,'g'); hold on;
        bic = plot(chemshift(bic_index + (ii-1) * xshift),signal(bic_index)+(ii-1)*yshift,'m'); hold on;
        legend([pyr,lac,ala,bic],'Pyruvate','Lactate','Alanine','Bicarbonate')
    else
        fill(chemshift(x),signal(cs_win)+(ii-1)*yshift,'w','LineStyle','none'); hold on;
        plot(chemshift(x),signal(cs_win)+(ii-1)*yshift,'k','LineWidth',1.5); hold on;
    end
    set(gca, 'XDir', 'reverse');
end
hold off;  