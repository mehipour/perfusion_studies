function [win, cs, pyr_index,lac_index, ala_index, bic_index, hyd_index co2_index]...
    = select_pyr_C1_peaks_ppm(F,ppm_ref,t_plot,path,NP,chemshift_limit)

%[win, cs, pyr_index,lac_index, ala_index, bic_index, hyd_index co2_index]...
%    = select_pyr_C1_peaks_ppm(F,ppm_ref,t_plot,path,NP,chemshift_limit)
%
% select the indices for the peaks after injection of 1-13C pyruvate. 
%
% input parameters:
%
% F is the complex fid after phase correction and base removal. 
% t_plot the time point at which the spectra is looked at to determine
% the maximum peak (which is pyruvate).
% ppm_ref is the reference ppm of the largest peak. (ex: 171 ppm for c1
% pyruvate.)
% path: fid file path. 
% NP = number of points. 
% chemshift_limit : range of chemical shifts to show.
%
% output parameters:
%
% win: indicies asociates with the total ppm range. 
% cs: vector of chemical shifts.
% indices: .... vectors of each chemical shifts of each peak. 
%
% created by Mehrdad Pourfathi on 07/03/2013
%
% Update 1 by Mehrdad Pourfathi on 08/06/2013
% range of chemical shift range was made as an input to the function.
fid_path = path(1:end-4);   % remove .fid from the end of the file path. 

sw = readprocpar(fid_path, 'sw'); sw = sw(2);
% tr = readprocpar(fid_path, 'tr'); tr = t(2);
sfrq = readprocpar(fid_path, 'sfrq'); sfrq = sfrq(2);

[p,Npyr] = max(abs(F(:,t_plot)));
cs = index2ppm(ppm_ref,Npyr,NP,path); 
win = intersect(find(cs>chemshift_limit(1)),find(cs<chemshift_limit(2)));
pyr_index = intersect(find(cs>168.6),find(cs<173.6));
% lac_index = intersect(find(cs>182.4),find(cs<183.5));
lac_index = intersect(find(cs>182.8),find(cs<184));
ala_index = intersect(find(cs>176.3),find(cs<177));
% ala_index = intersect(find(cs>172),find(cs<177.9));
bic_index = intersect(find(cs>160),find(cs<162));
hyd_index = intersect(find(cs>178.9),find(cs<180.7));
co2_index = intersect(find(cs>122),find(cs<128));

