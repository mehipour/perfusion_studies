function [Mpyr,Mlac,Mala,Mbic,Mhyd,tim,tim_idx,tstart] = integrate_peaks_and_show_kinetics_C1_pyr(path,Fr,tr,Navg,Nmove,tend,pyr_index,lac_index,ala_index,bic_index,hyd_index,met_ratio,baseline_pnts,find_tstart,save_peakfits)

%
% function [Mpyr,Mlac,Mala,Mbic,Mhyd,tim,t_ind,tstart]...
% = integrate_peaks_and_show_kinetics_C1_pyr...
% (path,Fr,tr,Navg,Nmove,tend,pyr_index,lac_index,ala_index,bic_index,hyd_index,met_ratio,baseline_pnts,find_tstart,save_peakfits)
% Takes fid file path, real part of the spectra and the indicies for
% chemical shifts of 1-c-pyruvate and its metabolites finds the integrals
% under the peaks. 
% 
% if the file is already analyzed (peakfit.txt) file exists then it loads
% the contents of the peakfit.txt file. 
%
%% Created by Mehrdad Pourfathi on 07/30/2013.
%
% Update 1 by Mehrdad Pourfathi on 02/17/2014
% 1.added the find_tstart flag to make the finding of the injection time
% optional.
% 2. also optionally saves the peakfit file. 
%
% Update 2 by Mehrdad Pourfathi 04/08/2014
% modified how the code finds the starting time of the injection.
% added tr to in the input argument of the function.
%
% Update 3 by Mehrdad Pourfathi 10/07/2014
% Bugs with writing the text file during integration was fixed.


% if peakfit file does not exists find integrals, else load them. 
if (~exist([path '/peakfit.txt'],'file'))   % If the peakfit.txt file does not exists
  
    M = floor((size(Fr,2)-Navg)/Nmove)+1;
    Fr_avg = zeros(size(Fr,1),M);
    tim = zeros(M,1);
    tim_idx = zeros(M,1);
    jj = 0;
    for ii = 1:Nmove:size(Fr,2)
        jj = jj+1; 
        tim_idx(jj) = jj; % create time index vector
        if jj-1 % create ime vector
            tim(tim_idx(jj)) = Nmove + tim(tim_idx(jj-1));
        end
        idx = Nmove*(jj-1)+1:Nmove*jj;
        Fr_avg(:,jj) = mean(Fr(:,idx),2);
    end

    Mpyr = sum(Fr_avg(pyr_index,:)); % no baseline correction
    Mlac = remove_peak_baseline(Fr_avg,lac_index,[baseline_pnts(1) baseline_pnts(2)]);
    Mala = remove_peak_baseline(Fr_avg,ala_index,[baseline_pnts(3) baseline_pnts(4)]);
    Mbic = remove_peak_baseline(Fr_avg,bic_index,[baseline_pnts(5) baseline_pnts(6)]);
    Mhyd = remove_peak_baseline(Fr_avg,hyd_index,[baseline_pnts(7) baseline_pnts(8)]);
    tim = [1:size(Mpyr,2)]-1;
    
    Mlac = Mlac/max(Mpyr);
    Mala = Mala/max(Mpyr);
    Mbic = Mbic/max(Mpyr);
    Mhyd = Mhyd/max(Mpyr);
    Mpyr = Mpyr/max(Mpyr);
    
    tim = tim*tr;
    y = [tim; Mpyr; Mlac; Mala; Mbic; Mhyd];
%     save_to_txtfile(path,'C1ala',y)  
    if save_peakfits
        save_to_txtfile(path,'C1',y)  
    end
else
%    loads peaks if peakfit.txt file exists
    peaks = rearrange_peaks_from_txtfile(path);
    tim_idx = round(peaks(:,1)/tr);
    tim = peaks(:,1);
    Mpyr = peaks(:,2);
    Mlac = peaks(:,3);
    Mala = peaks(:,4);
    Mbic = peaks(:,5);
    Mhyd = peaks(:,6);
end    
% Normalize Peaks
Mlac = Mlac/sum(Mpyr);
Mbic = Mbic/sum(Mpyr);
Mala = Mala/sum(Mpyr);
Mhyd = Mhyd/sum(Mpyr);
Mpyr = Mpyr/sum(Mpyr);

if find_tstart
    tstart_idx = find_injection_start_time(Mpyr);
else
    tstart_idx = 1;
end
% tstart = 0;
% tstart_idx =0;
% tim_idx1 = tstart_idx + 1:length(Mala);
tim = tim - tim(tstart_idx);
tim_idx = intersect(tim_idx(find(tim<tend)),tim_idx(find(tim>=0)))+1;
tim = tim(tim_idx);
tstart = tstart_idx*tr;
 
pyr_ratio = met_ratio(1);
bic_ratio = met_ratio(4);

plot(tim,Mpyr(tim_idx)*pyr_ratio,'LineWidth',1.5); hold on;
plot(tim,Mlac(tim_idx),'r','LineWidth',1.5); hold on;
plot(tim,Mala(tim_idx),'g','LineWidth',1.5); hold on;
plot(tim,bic_ratio*Mbic(tim_idx),'m','LineWidth',1.5); hold on;
xlabel('Time (seconds)');
ylabel('Relative Signal Intensity (A.U.)');
% ylim([0 max(Mpyr/sum(Mpyr))*pyr_ratio]); grid on;
hold off;   
legend(strcat('Pyruvate (x ',num2str(pyr_ratio),')')...
    ,'Lacate','Alanine',...
    strcat('Bicarbonate (x ',num2str(bic_ratio),')'));