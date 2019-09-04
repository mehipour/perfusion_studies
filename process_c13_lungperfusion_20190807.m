% function process_c13_lungperfusion_xxxxxxxx
%
% this analysis tool was used for Ascorbic Acid Project to analyze C13
% spectra and fit the pyruvate, lactate, alanine and bicarbonate peaks to
% the kinetic model for one injection (or several differet files)
%
% File created: 1/28/2013 by Mehrdad Pourfathi
%
% Update 1: 2/21/2013 by Mehrdad Pourfathi
% % averages the carbon spectra for demonstration..
%
% Update 2: 2/22/2013 by Mehrdad Pourfathi
% % 
%
% Update 3: 5/10/2013 by Mehrdad Pourfathi
% Added baseline correction, modified chemical shift labeling function to
% use global variables. Added waterplot capability. 
%
% Update 4: 5/13/2013 by Mehrdad Pourfathi
% Added first order phase correction. 
% Added condition to check for C1 or C2 labeled carbon on the pyruvate and
% switch limits for baseline adjustment accordingly.
%
% Update 5: 5/15/2013 by Mehrdad Pourfathi
% Changed first order phase correction to occur outside the loop. 
% Also change the mechanism for first order correction. (Corrected the
% frequency spacing)
%
% Update 6: 5/17/2013 by Mehrdad Pourfathi
% Added two-site exchange model for pyruvate lactate interconversion to fit
% to data and acquire apparente forward and backward rate constants.
%
% Update 7: 5/22/2013 by Mehrdad Pourfathi
% Completed the two-site exchange fit model for pyruvate and lactate. The
% cost function is for simultanous fittin of pyruvate and lactate data and
% uses weighted mean square error of the data and the simulated results.
% The lactate MSE is weighted higher to penalize its error more than
% pyruvate since the data is much smaller in values. 
% however the fitting does not work perfectly.
%
% Update 8: 6/03/2013 by Mehrdad Pourfathi
% Added color to the peaks when showing the spectra.
%
% Update 9: 6/04/2013 by Mehrdad Pourfathi
% Added section to save time course of the integral of the peaks in a  text
% file. Also added the the pyruvate hydrate peak.
%
% Update 10: 6/07/2013 by Mehrdad Pourfathi
% Fixed a bug associated with coloring peaks differently in the waterfall
% plot. Also added legents to the colored peaks. 
% Updated the visualization of the graphs and fitting parameters of the
% kinetic model. 
% Also adds the varian function path automatically.
% A major bug in the T1 estimation of the pyruvate and lactate species for
% the fit was corrected. 
%
% Update 11: 6/08/2013 by Mehrdad Pourfathi
% Made the code modular; parts of the code where turned into functions that
% are executed in order.
%
% Update 12: 6/11/2013 by Mehrdad Pourfathi
% Added non-linear regression to fit kinetic data to the kinetic model as
% oppposed to running an MSE or LSQ curve fitting. The benefit of uinsg
% non-linear fitting is that non-linear fitting toolbox in matlab provides
% covariance matrix of the estimated parameters which can be used to
% estimate the confidence interval of the fitting parameters. 
% Also added a function to estimate the local SNRs of the peaks to weight
% the MSEs of each peak inversly to the its local SNR to penalize deviation
% of the fit from smaller peaks more. 
%
% Update 13: 6/13/2013 by Mehrdad Pourfathi
% Added the effect of RF pulse on the fit, assuming that it is constant
% over the fit. (Not working though)
%
% Update 14: 6/27/2013 by Mehrdad Pourfathi
% Added the save_phasing_info function to save ph0 and ph1 for each file as
% a text file inside the '.fid' folder.
%
% Update 15: 6/30/2013 by Mehrdad Pourfathi
% Normalize spectra to maximum value when plotting.
%
% Update 16: 7/10/2013 by Mehrdad Pourfathi
% Added the function "find_injection_start_time()" in the
% integrate_peaks_and_show_kinetics() function.
% Modified the show_waterfall function and created it as a function.
% all data saved to text files are now using one single function named
% "save_to_txtfile(path,'ph',y)".
%
% Update 17: 7/30/2013 by Mehrdad Pourfathi
% -Renamed the function integrate_peaks_and_show_kinetics to
% integrate_peaks_and_show_kinetics_C1 and made as a separate function
% it was also modified modified to load peakfit.txt file if available
% -Also updated the function correct_phase_and_baseline and use that instaed
% of the internal funciton. 
% -Additionally the new change includes reading the files by going through
% the directory and selecting "13C" .fid folders.
% -The figures are automatically saved in png format (300 dpi)
%
% Update 18: 8/02/2013 by Mehrdad Pourfathi
% -It is not a function anymore.
% -The kinetic fit is a separate funciton now.
%
% Update 19: 8/08/2013 by Mehrdad Pourfathi
% -Autophasing is added!!
% -kinetics (Mpyr,Mlac,etc.) are not arrayed over the files anymore. 
% -Net metabolite fractions are not calculated anymore. 
% -All global variables are dismissed
%
% Update 20: 2/17/2014 by Mehrdad Pourfathi
% -integrate_peaks_and_show_kinetics_C1 function was updated to make the
% finding of the injection time optional. 
% - optionally saves the peakfit file if fitting is not used. 
%
% Update 21: 2/21/2014 by Mehrdad Pourfathi
% - Saves the matlab figure files as well as the png files. 
%
% Update 22: 3/22/2014 by Mehrdad Pourfathi
% - disables latex interpreters for presenting the file name in figure
% titles
%
% Update 23: 4/9/2014 by Mehrdad Pourfathi
% - used save_figure function to save figures.
% - (major) started using integrate_peaks_and_show_kinetics_C1_pyr
% function.
%
% Update 24 (MAJOR) 8/7/2019 by Mehrdad Pourfathi
% - Combined Steve's function fit spectra with this code
% - has a flag to give user option to fit data.
% -cleaned up the code, and deleted excess fitting code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize variables and paths.

% clear all

% where the varian function is.
%addpath('/Users/mehipour/Documents/MATLAB/Varian matlab Files/C13 Data Processing');
%addpath('/Users/mehipour/Documents/MATLAB/Varian matlab Files/Perfused Lung Studies');

close all; format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load Directory
data_list_perfusion_studies_20190807();

% Processing parameters
lb = 30;   % not used for the fit
fit_spectra = 1;

% Output parameters
show_in_subplot = 0;
save_png = 1;
save_eps = 0;
save_fig = 1;

% Waterplot visualization parameters
waterfall_plot = 1;
show_color = 0;
Nsubplotx = 1;
Nsubploty = 1;
Navg = 1;
Nmove = 1;
tend = 250;
N_interval_waterplot = 5;
xshift = 0.1;   % shift in chemical shifts
yshift = -3e-2; % shift in y axis.

% Timecourse visualization parameters
pyr_scale = 0.03; 
bic_scale = 2;
find_tstart = 0;
time_range = 30:120;    % time to average and plot.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main function


for kk = 9:9

    file = dir(path(kk).name);
    ii = 0;
    % for ii = numel(cohort):-1:1
    for jj = 1:numel(file)   
        % searches for 13C files. 
        if or(~isempty(strfind(lower(file(jj).name),'13c')), ~isempty(strfind(lower(file(jj).name),'c13')))
            ii = ii + 1;
            disp(file(jj).name)
            filepath = [path(kk).name '/' file(jj).name];
            fig_path = filepath;
        else
            continue
        end

        %% Load Spectra from raw-data files and apply line broadening. 
        [F,NP,NB,sw,sfrq,tr] = read_and_broaden_fid(filepath,lb);

        % Select peaks and convert indicies to ppm (chemical shift)
        [win, chemshift, pyr_index,lac_index, ala_index, bic_index, hyd_index]...
            = select_pyr_C1_peaks_ppm(F,171,time_range,filepath,NP,[150,210]);

        %% Correct phase and baseline. 
        ph0 = pi; ph1 = 0;
        [Fr_avg,Fc,Fr,Fi_avg] = correct_phase_and_baseline...
        (F,NP,NB,win,time_range,filepath,ph0,ph1,[0,1,0,1]);
        Fr_avg = Fr_avg/max(Fr_avg);

        %% Fit spectra
        if fit_spectra
            if (~exist([filepath '/peakfit.txt'],'file'))   % If the peakfit.txt file does not exists
               perfused_lung_13c_fit_20190807(filepath);
            end
            mets = rearrange_peaks_from_txtfile(filepath);
        end

        %%
        figure(1);
        if show_in_subplot
            subplot(Nsubploty,Nsubplotx,ii);
            fig_path = filepath;
        end
        plot(chemshift(win),Fr_avg(win),'LineWidth',2); hold off;
        if show_color
            plot(chemshift(win),Fr_avg(win),'k'); hold on;
            pyr = plot(chemshift(pyr_index),Fr_avg(pyr_index),'b'); hold on;
            lac = plot(chemshift(lac_index),Fr_avg(lac_index),'r'); hold on;
            ala = plot(chemshift(ala_index),Fr_avg(ala_index),'g'); hold on;
            bic = plot(chemshift(bic_index),Fr_avg(bic_index),'m'); hold off;
            legend([pyr,lac,ala,bic],'Pyruvate','Lactate','Alanine','Bicarbonate')
        end
        ylim([-0.05*pyr_scale 1*pyr_scale]);
        set(gca, 'XDir', 'reverse');
        xlabel('Chemical Shift (ppm)'); ylabel('Signal Intensity (A.U.)'); 
        title(file(jj).name,'interpreter','none');

        save_figure('/spectrum',fig_path,save_fig,save_png,save_eps);


        %% Integrate under the pyruvate, lactate, alanine and bicarbonate peaks.    
        % if data spectra was fitted, it will plot the fits
        figure(2)
        baseline_pnts = [100 100 100 100 100 100 100 100];
        [Mpyr,Mlac,Mala,Mbic,Mhyd,tim,tim_idx,tstart]...
            = integrate_peaks_and_show_kinetics_C1_pyr(filepath...
            ,Fr,tr,Navg,Nmove,tend,pyr_index,lac_index,ala_index,bic_index,...
            hyd_index,[pyr_scale,1,1,bic_scale],baseline_pnts,find_tstart,0);
        title(file(jj).name,'interpreter','none');
        save_figure('/kinetics1',fig_path,save_fig,save_png,save_eps);

        %% plot waterplot
        if waterfall_plot
            figure(3);
            Nspectra = [1:N_interval_waterplot:tend]+tstart;
            Nspectra = round(Nspectra/tr);
            Nspectra = 1:N_interval_waterplot:tend;
            if show_in_subplot
               subplot(Nsubploty,Nsubplotx,ii);
               fig_path = filepath;
            end
            cs_win = [160 185];     % range of the chemical shift to be displayed.
            if show_color
                show_waterfall(Fr,chemshift,cs_win,Nspectra,xshift,yshift,pyr_index,lac_index,bic_index,ala_index)
            else
                show_waterfall(Fr,chemshift,cs_win,Nspectra,xshift,yshift)
            end    
            xlabel('Chemical Shift (ppm)'); ylabel('Normalized Spectra (A.U.)'); 
            title(file(jj).name,'interpreter','none');
            save_figure('/waterfall',fig_path,save_fig,save_png,save_eps);
        end

    end
end

