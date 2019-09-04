function [F,NP,NB,sw1,sfrq1,tr1] = read_and_broaden_fid(path,lb)
% opens the fid file at path and applies line broadening by lb Hz and
% returns the complex spectra as an output. 
%
% input parameters:
%
% path: path of the fid file.
% lb  : line broadening applied in Hz. 
% 
% ouptut parameters: 
%
% F: complex spectra.
% NP: number of points per spectrum
% NB: number of spectra.
% 
% optional outputs:
%
% sw1 = spectral bandwidth
% sfrq1 = center frequency
% tr1 = repetition time
%
% Created by Mehrdad Pourfathi on 07/03/2013.
% 
% Update 1: by Mehrdad Pourfathi on 07/10/2013
% Gives sw,sfrq and tr as outputs.
%

    fid_path = path(1:end-4);   % remove .fid from the end of the file path. 
    [RE,IM,NP,NB,~,HDR]=varianloadfid(fid_path,1,1);
    % Read spectral parameters from proc file. 
    sw = readprocpar(fid_path, 'sw'); sw = sw(2);
    tr = readprocpar(fid_path, 'tr'); tr = tr(2);
    sfrq = readprocpar(fid_path, 'sfrq'); sfrq = sfrq(2);
    % Fourier Transform Loop.
    n = [0:NP-1]';
    
    for kk = 1:NB
        C(:,kk) = squeeze(RE(:,kk))+ 1i*squeeze(IM(:,kk));   % read raw FID
%         C(:,kk) = C(:,kk) - mean(C(36000:end,kk));
        C(:,kk) = C(:,kk).*exp(-n/sw*lb);                    % Exponential Line Broadening
        A(:,kk) = fftshift(fft(C(:,kk)));                    % Apply Fourier Transform  
    end
    F = A;
    if nargout == 4
        sw1 = sw;
    elseif nargout == 5
        sw1 = sw;
        sfrq1 = sfrq;
    elseif nargout == 6;
        tr1 = tr;
%         tr1 = 3;
        sw1 = sw;
        sfrq1 = sfrq;
    end
  
        
        