function [c1,F1,Fr] = analyzeP31()
%
% this analysis tool was used for Ascorbic Acud Project to plot the
% phosphorous spectra taken over time for a given perfused lung. 
%
% File created: 2/20/2013 by Mehrdad Pourfathi
%
% Update 1: 1 2/21/2013 by Mehrdad Pourfathi
% 


%% Load FIle
clear all

% Experiment parameters.
% data_directory = '/Users/mehipour/Documents/MATLAB/HUP-B data/20130128/';
% data_directory = '/Users/mehipour/Documents/MATLAB/HUP-B data/20130129/';
% data_directory = '/Users/mehipour/Documents/MATLAB/HUP-B data/20130130/';
% data_directory = '/Users/mehipour/Documents/MATLAB/HUP-B data/ASA/';
% data_directory = '/Users/mehipour/Documents/MATLAB/HUP-B data/20130218/';
data_directory = '/Volumes/Macintosh HD 1/Data/proccessig folder 1/processing folder 2/';
Nsubplotx = 1;
Nsubploty = 1;


% Load files
% files = { ...     
%     struct('filename', 'C13-perfusedlung_1+4mMpyr_012813', 'name', 'result/C13-perfusedlung_1+4mMpyr_012813', 'lb', 5, 'autophase', 1, 'ph0', 1.36*pi, 'ph1', 0, 'tstart', 18,'title','Without ASA') ... 
%     struct('filename', 'C13-perfusedlung+2mM-ASA_2+4mMpyr_012813', 'name', 'result/C13-perfusedlung+2mM-ASA_2+4mMpyr_012813', 'lb', 5, 'autophase', 1, 'ph0', 1.7*pi, 'ph1', 0, 'tstart', 11,'title','2mM ASA') ...  
%     struct('filename', 'C13-perfusedlung_1+2mMASA_4mM_012913', 'name', 'result/C13-perfusedlung_1+2mMASA_4mM_012913', 'lb', 5, 'autophase', 1, 'ph0', 1.71*pi, 'ph1', 0, 'tstart', 18,'title','2mM ASA') ... 
%     struct('filename', 'C13-perfusedlung_2_4mM_012913', 'name', 'result/C13-perfusedlung_2_4mM_012913', 'lb', 5, 'autophase', 1, 'ph0', 1.21*pi, 'ph1', 0, 'tstart', 18,'title','Without ASA') ... 
%     struct('filename', 'C13-perfusedlung_2+2mMASA_4mM_012913', 'name', 'result/C13-perfusedlung_2+2mMASA_4mM_012913', 'lb', 5, 'autophase', 1, 'ph0', 1.21*pi, 'ph1', 0, 'tstart', 18,'title','2mM ASA') ... 
%     struct('filename', 'C13-perfusedlung_1+2mMASA_4mM_f_013013', 'name', 'result/C13-perfusedlung_1+2mMASA_4mM_f_013013', 'lb', 5, 'autophase', 1, 'ph0', 1.21*pi, 'ph1', 0, 'tstart', 18,'title','2mM ASA') ... 
%     struct('filename', 'C13-perfusedlung_2+2mMASA_4mM_f_102min_013013', 'name', 'result/C13-perfusedlung_2+2mMASA_4mM_f_102min_013013', 'lb', 5, 'autophase', 1, 'ph0', 1.21*pi, 'ph1', 0, 'tstart', 18,'title','2mM ASA') ... 
%     struct('filename', 'C13-perfusdlung_normal1+2mMASA_46min_020713_1', 'name', 'result/C13-perfusdlung_normal1+2mMASA_46min_020713_1', 'lb', 5, 'autophase', 1, 'ph0', 1.7*pi, 'ph1', 0, 'tstart', 18,'title','2mM ASA') ... 
%     struct('filename', 'C13-perfusdlung_normal1+2mMASA_46min_020713_1', 'name', 'result/C13-perfusdlung_normal1+2mMASA_46min_020713_1', 'lb', 5, 'autophase', 1, 'ph0', 1.7*pi, 'ph1', 0, 'tstart', 18,'title','2mM ASA') ... 
%     struct('filename', 'C13-perfusdlung_normal2+2mMASA_140min_020713_2', 'name', 'result/C13-perfusdlung_normal2+2mMASA_140min_020713_2', 'lb', 5, 'autophase', 1, 'ph0', 1.71*pi, 'ph1', 0, 'tstart', 18,'title','2mM ASA') ... 
%     struct('filename', 'C13-perfusdlung_normal1+2mMASAcoinjection_020813', 'name', 'result/C13-perfusdlung_normal1+2mMASAcoinjection_020813', 'lb', 5, 'autophase', 1, 'ph0', 1.7*pi, 'ph1', 0, 'tstart', 18,'title','Without ASA') ...   
%     struct('filename', '31p-perfusdlung1_normal1_before_021813', 'name', 'result/31p-perfusdlung1_normal1_before_021813', 'lb', 5, 'autophase', 1, 'ph0', pi+3*pi/10, 'ph1', 0, 'tstart', 18,'title','Averaged ^{31}P Spectum 60 mins after Perfusion') ...   
%     struct('filename', 'P31-perfusedlung1_Normal2_4mM_021813', 'name', 'result/P31-perfusedlung1_Normal2_4mM_021813', 'lb', 5, 'autophase', 1, 'ph0', pi+3*pi/10, 'ph1', 0, 'tstart', 18,'title','Averaged ^{31}P Spectrum 110 mins after Perfusion') ...   

% };

files = { ...
%     struct('filename', '20140304_p31_alanine_lung_control_after', 'name', 'result/20130509_31p-perfusedlung_1_2mM_GSH', 'lb', 100, 'autophase', 1, 'ph0', 4*pi/3+pi/10, 'ph1', 0, 'tstart', 18,'title','^{31}^P spectrum') ...   
%     struct('filename', '20140304_p31_alanine_lung_control', 'name', 'result/20130509_31p-perfusedlung_1_2mM_GSH', 'lb', 100, 'autophase', 1, 'ph0', 4*pi/3-pi/10, 'ph1', 0, 'tstart', 18,'title','^{31}^P spectrum') ...   
%     struct('filename', '20131223_31p-perfusedlung_control_1Cala', 'name', 'result/20130509_31p-perfusedlung_1_2mM_GSH', 'lb', 100, 'autophase', 1, 'ph0', 2*pi/3, 'ph1', 0, 'tstart', 18,'title','^{31}^P spectrum') ...   
    struct('filename', '20140522_perfusedmouselung_metcancer20_1024_after', 'name', 'result/20130509_31p-perfusedlung_1_2mM_GSH', 'lb', 120, 'autophase', 1, 'ph0',7*pi/8, 'ph1', 0, 'tstart', 18,'title','^{31}^P spectrum') ...   


    %     struct('filename', '20131009_1_31p-perfusedlung_TV_25_severe', 'name', 'result/20130509_31p-perfusedlung_1_2mM_GSH', 'lb', 100, 'autophase', 1, 'ph0', pi/3, 'ph1', 0, 'tstart', 18,'title','^{31}^P spectrum') ...   
%     struct('filename', '20130925_2_31p-perfusedlung_TV_25_mild', 'name', 'result/20130509_31p-perfusedlung_1_2mM_GSH', 'lb', 100, 'autophase', 1, 'ph0', pi/3, 'ph1', 0, 'tstart', 18,'title','^{31}^P spectrum') ...   

    %     struct('filename', '20130411_31p-perfusedlung_2_0p5mM_ASA_after', 'name', 'result/20130509_31p-perfusedlung_1_2mM_GSH_after', 'lb', 20, 'autophase', 1, 'ph0', -pi/10, 'ph1', 0, 'tstart', 18,'title','^{31}^P spectumr of Lung 1 after 75 mins with 2 mM GSH') ...   
%     struct('filename', '20130509_31p-perfusedlung_2_2mM_GSH', 'name', 'result/20130509_31p-perfusedlung_2_2mM_GSH', 'lb', 5, 'autophase', 1, 'ph0', pi+3*pi/10, 'ph1', 0, 'tstart', 18,'title','^{31}^P spectumr of Lung 2 after 15 mins with 2 mM GSH') ...   
%     struct('filename', '20130509_31p-perfusedlung_2_2mM_GSH_after', 'name', 'result/20130509_31p-perfusedlung_2_2mM_GSH_after', 'lb', 5, 'autophase', 1, 'ph0', pi+3*pi/10, 'ph1', 0, 'tstart', 18,'title','^{31}^P spectumr of Lung 2 after 80 mins with 2 mM GSH') ...   
%     struct('filename', '01', 'name', 'result/20130509_31p-perfusedlung_2_2mM_GSH_after', 'lb', 5, 'autophase', 1, 'ph0', pi+3*pi/10, 'ph1', 0, 'tstart', 18,'title','^{31}^P spectumr of Lung 2 after 80 mins with 2 mM GSH') ...   

    };



for ii = 1:length(files)    

    % Load Spectra from raw-data files. 
    [RE,IM,NP,NB,NT,HDR]=varianloadfid([data_directory files{ii}.filename], 1,1);
    
    % Read spectral parameters from proc file. 
    sw = readprocpar([data_directory files{ii}.filename], 'sw'); sw = sw(2);
    tr = readprocpar([data_directory files{ii}.filename], 'tr'); tr = tr(2);
    sfrq = readprocpar([data_directory files{ii}.filename], 'sfrq'); sfrq = sfrq(2);

    
    n = [0:NP-1]';
    for kk = 1:NB
        C(:,kk) = squeeze(RE(:,kk))+ 1i*squeeze(IM(:,kk));    % read raw FID
        C(:,kk) = C(:,kk).*exp(-n/sw*files{ii}.lb);           % Exponential Line Broadening
        F(:,kk) = fftshift(fft(C(:,kk)));                     % Apply Fourier Transform
    end
    
    %% Apply zero-order phase correction. 
    Fr = real(F*exp(1i*files{ii}.ph0)); 
    Fi = imag(F*exp(1i*files{ii}.ph0));
    
    % translate sample number to ppm values. 
    % chemical shifts of gamma-ATP, alfa-ATP and beta-ATP peaks in lung
    % extracts are -3 ppm, -7.9 ppm and -18.4 ppm respectively. [31P
    % spectroscopy of isolated perfusated rat lungs, Y. Hayashi et al.,
    % Journal of Applied Physiology, 1993]
    [p,Npyr] = max(Fr);
    cs = index2ppm(NP,sw,4.7,Npyr,sfrq);    % reference would be alfa-ATP peak
%     cs = index2ppm(NP,sw,3.665,Npyr,sfrq);    % reference would be alfa-ATP peak

    win = intersect(find(cs>-25),find(cs<15)); 

    figure(1); hold on;
    subplot(Nsubploty,Nsubplotx,ii);
    plot(cs(win),Fr(win)/max(Fr(win)) );
    c1 = cs(win);
    F1 = Fr(win)/max(Fr(win));
%     plot(cs(win),zeros(size(win))+(ii-1)*1,'k');
    title(files{ii}.title);
    xlabel('ppm'); ylabel('A.U.');
    set(gca, 'XDir', 'reverse');
    hold off
    ylim([-0.2 2]); box on;
        hold off;
%     figure(4);
%     subplot(2,1,ii);
%     plot(cs(win),Fr(win)/max(Fr(win)));
%     title(files{ii}.title);
%     xlabel('ppm'); ylabel('A.U.');
%     set(gca, 'XDir', 'reverse');
    
end

%% Functions
function a = readprocpar(f, param)
% Before the end of file, read a line, if the first character in that line
% is parameter read and return the value on the next line!    
% g = fopen([f '.fid\procpar']);
    g = fopen([f '.fid/procpar']);
    while(~feof(g))
        s = fgets(g);
        if (strcmp(strtok(s),param))
            a = strread(fgets(g));
        end
    end
    
function cs = index2ppm(NP,sw,ref_ppm,ref_index,f0)
% creates chemical shift vector given reference ref_ppm (typically pyruvate
% at 171 ppm) and center frequency f0 in MHz and sw in Hz.
    sw_ppm = (sw/f0);
    delta_cs = sw_ppm/NP;
    cs1 = [ref_index+1:-1:1]' * delta_cs;
    cs2 = [-1:-1:-NP+ref_index]' * delta_cs; 
    cs = [cs1;0;cs2] + ref_ppm;
        