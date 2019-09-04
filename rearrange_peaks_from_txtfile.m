function peaks = rearrange_peaks_from_txtfile(path) 
% Reads the saved executed kinetic model from the path and rearrage the
% peaks as followed:
% time pyruvate lactate alanine bicarbonate pyruvate-hydrate
% 
% Created by Mehrdad Pourfathi on 07/08/2013
%
addpath('/Users/mehipour/Documents/MATLAB/Varian matlab Files/');
kinetic_file_path = [path '/peakfit.txt'];
data = importdata(kinetic_file_path);

% rearrange the coloumns in the peakfit.txt file. 
    for kk = 1:length(data.textdata)
        if strcmp(data.textdata(kk),'time')
            peaks(:,1) = data.data(:,kk);
        elseif strcmp(data.textdata(kk),'pyruvate')
            peaks(:,2) = data.data(:,kk);
        elseif strcmp(data.textdata(kk),'lactate')
            peaks(:,3) = data.data(:,kk);
        elseif strcmp(data.textdata(kk),'alanine')
            peaks(:,4) = data.data(:,kk);
        elseif strcmp(data.textdata(kk),'bicarb')
            peaks(:,5) = data.data(:,kk);
        elseif strcmp(data.textdata(kk),'hydrate')
            peaks(:,6) = data.data(:,kk);
        end
    end