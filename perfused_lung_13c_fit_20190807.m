function perfused_lung_13c_fit_20190807(fidfile)
global RE;
global IM;
global line;
global current_line;
global current_spec;
global c;
global g;
global f;
global xscale;
% global subdir_slash;
global sn;

    fitting_options = optimset('TolFun', 1e-2, 'MaxFunEvals' , 3000, 'MaxIter', 3000, 'Display', 'on');

    clear_peakfits = 0;
    dirname = fidfile(1:(end-4));
    
    % chemical shift values
    column_values = { ' pyruvate ', '  lactate ', '   bicarb ', '  alanine ', '  hydrate ' };
    pyruvate_chemical_shift = 170.88;
    chemical_shift_estimates = [0 12.09 -10.11 5.70 8.285]; % relative to pyruvate
    
    % read acquisition parameters
    sfrq = readprocpar(dirname, 'sfrq'); sfrq = sfrq(2);
    sw = readprocpar(dirname, 'sw');    sw = sw(2);
    tr = readprocpar(dirname, 'tr');    tr = tr(2);
%     if (tr == -1)
%         at = readprocpar(dirname, 'at');
%         d1 = readprocpar(dirname, 'd1');
%         pw = readprocpar(dirname, 'pw');
%         tr = at + d1 + pw * 1E-6;
%     end
    % read spectra
    [RE,IM,NP,NB,NT,HDR] = varianloadfid(dirname, 1, 1);
    if (~make_complex_spectra(2, sw)) % initial line broadening estimate = 2 Hz, gets refined later
        % no pyruvate found, make an empty file and skip it
        fclose(fopen([fidfile '/peakfit.txt'], 'w'));
    end
    window_complex_spectra(pyruvate_chemical_shift, chemical_shift_estimates, sw, sfrq);
    [ph0 ph1] = autophase(atan2(-imag(sum(g)), real(sum(g))));
    % peak pick- peaks are maxima that are > 10x noise and it goes down by at least 10% before going up again
    gp = g * 0;
    line = [];
    for m = size(c, 1):-1:1
        if (sn(m) < 12)
            continue;
        end
        gp = gp + real(f(m,:));
        peaks = gp > 8 * std(gp(1:50));
        for q = 1:200
            peaks((1+q):(length(gp)-q)) = peaks((1+q):(length(gp)-q)) & ...
                (gp((1+q):(length(gp)-q)) > gp(1:(length(gp)-2*q))) & ...
                (gp((1+q):(length(gp)-q)) > gp((2*q+1):length(gp)));
        end
        peaksn = gp(peaks) / std(gp(1:50));
        peakindices = find(peaks == 1);
        if (length(peakindices) == 0)
            continue;
        end
        % adjust each line frequency to be that of the closest picked line
        [max_val max_idx] = max(g);
        for k = 1:length(chemical_shift_estimates)
            % find the peak that's closest
            [min_val pi] = min(abs(xscale(peaks) - pyruvate_chemical_shift - chemical_shift_estimates(k)));
            if (min_val > 0.5)
                continue;
            end
            found = cell2mat(arrayfun(@(x)(strcmp(x.name, column_values(k))), line, 'uniformoutput', false));
            if (sum(found) == 0)
                a.center = xscale(peakindices(pi));
                a.name = column_values{k};
                a.ispyruvate = (chemical_shift_estimates(k) == 0);
                a.width = 0.2;
                line = [line a];
            end
            if ((sum(found) > 0) && (peaksn(pi) < 20))
                line(found).center = xscale(peakindices(pi));
                for lo = peakindices(pi):-1:2
                    if (g(lo - 1) < g(peakindices(pi)) / 2)
                        break;
                    end
                end
                for hi = peakindices(pi):(length(g) - 1)
                    if (g(hi + 1) < g(peakindices(pi)) / 2)
                        break;
                    end
                end
                line(found).width = (xscale(hi) - xscale(lo)) / 2;
            end
        end
        plot(xscale, gp);
        set(gca, 'Xdir', 'reverse');
        hold on;
        for k = 1:length(line)
            plot([line(k).center line(k).center], [0 max(gp)], 'r');
            text(line(k).center, max(gp) * (1 - k/10), ['\leftarrow ' line(k).name],'HorizontalAlignment','left')
        end
        hold off;
        drawnow;
    end
    avgwidth = 0;
    for k = 1:length(line)
        [minval x_index] = min(abs(xscale - line(k).center));
        line(k).amplitude = real(g(x_index));
        avgwidth = avgwidth + line(k).width / length(line);
    end
    for k = 1:length(line)
        line(k).width = avgwidth;
    end

    % redo fft with line broadening equal to 1/2 of the linewidth we just found
    max_idx = make_complex_spectra(avgwidth * sfrq / 2, sw);
    window_complex_spectra(pyruvate_chemical_shift, chemical_shift_estimates, sw, sfrq);
    [ph0 ph1] = autophase(atan2(-imag(sum(g)), real(sum(g))));
    figure(1);
    plot(xscale,real(g),xscale,linefiteval);
    set(gca,'XDir','reverse');
    set(get(gca,'XLabel'),'String','chemical shift (ppm)');
    for k = 1:length(line)
        [minval x_index] = min(abs(xscale - line(k).center));
        text(xscale(x_index), line(k).amplitude, ['\leftarrow ' line(k).name],'HorizontalAlignment','left')
    end
    drawnow;
    % fit each peak independently to get amplitudes
    fitampl = zeros(size(c, 1),length(column_values)+1);
    for k = 1:size(c, 1)
        current_spec = squeeze(f(k, :));
        oldline = line;
        for q = 1:length(line)
            [min_val min_idx] = min(abs(xscale - line(q).center));
            x(q) = real(current_spec(min_idx));
            line(q).amplitude = x(q);
        end
        for q = 1:length(line)
            current_line = q;
            if (sn(k) > 12 && line(q).ispyruvate)
                [min_val min_idx] = max(current_spec);
                line(q).center = xscale(min_idx);
                xp = fminsearch(@linefit_fixed, [real(min_val)]);
                % xp = fminsearch(@linefit_movable, [line(q).amplitude line(q).center]);
                % line(q).center = xp(2);
            else
                xp = fminsearch(@linefit_fixed, [line(q).amplitude]);
            end
            line(q).amplitude = xp(1);
            fit_ampl(k, 1+q) = xp(1);
        end
        fit_ampl(k,1) = (k-1) * tr;
        yp = linefiteval;
        figure(2);
        plot(xscale, real(current_spec),xscale,yp);
        set(gca, 'XDir', 'reverse');
        drawnow;
        %pause;
        fx = fopen([fidfile '/peakfit.txt'], 'w');
        fprintf(fx, '     time ');
        for q = 1:length(column_values)
            fprintf(fx, '%s', column_values{q});
        end
        fprintf(fx, '\r\n');
        for q = 1:size(fit_ampl, 1)
            for p = 1:(length(column_values) + 1)
                printval = fit_ampl(q, 1);
                if (p > 1)
                    printval = 0.0;
                    for pp = 1:length(line)
                        if (strcmp(line(pp).name, column_values{p - 1}))
                            printval = fit_ampl(q, pp + 1);
                        end
                    end
                end
                fprintf(fx, '%9.5f ', printval);
            end
            fprintf(fx, '\r\n');
        end
        fclose(fx);
        line = oldline;
    end
    
function spectra_ok = make_complex_spectra(lb, sw)
global RE;
global IM;
global c;
global g;
global f;
global sn;
    c = (RE + 1i * IM)';
    baseline_subtraction_points = floor(size(c, 2) / 10);
    c = c - sum(sum(c(:,(size(c, 2) - baseline_subtraction_points + 1):size(c, 2)))) / (baseline_subtraction_points * size(c, 1));
    found_pyruvate = 0;
    f = zeros(size(c));
    g = zeros(1, size(c, 2));
    spectra_ok = 0;
    max_val = 0;
    for k = 1:size(c, 1)
        c(k,:) = c(k,:) .* exp(-lb*(0:(size(c, 2)-1))/sw);
        f(k,:) = fftshift(fft(c(k,:)));
        fr = real(squeeze(f(k,:)));
        fi = imag(squeeze(f(k,:)));
        fn = fr .* fr + fi .* fi;
        peakheight = sqrt(max(fn));
        snr_calculation_points = floor(length(fn) / 100);
        sn(k) = sqrt(peakheight^2 / sum(fn(1:snr_calculation_points)) * snr_calculation_points);
        if (sn(k) > 30) % only include spectra with S/N>30 in the average
            max_val = max(max(fn), max_val);
            g = g + squeeze(f(k,:));
            spectra_ok = 1;
        end
    end
    if (spectra_ok)
        f = f / sqrt(max_val);
        g = g / max(vecnorm(g));
    end
  
function window_complex_spectra(pyruvate_chemical_shift, chemical_shift_estimates, sw, sfrq)
global xscale;
global g;
global f;
    % find maximum of averaged spectra
    [max_val max_idx] = max(vecnorm(g));
    % set x scale and window
    xscale = ((-max_idx):(length(g) - max_idx-1)) / length(g) * sw / sfrq + pyruvate_chemical_shift;
    window = (max_idx + min(chemical_shift_estimates) * sfrq / sw * length(g)):(max_idx + max(chemical_shift_estimates) * sfrq / sw * length(g));
    window = (floor(max(0,window(1)-0.20*(window(end)-window(1))))):(ceil(min(length(g),window(end)+0.20*(window(end)-window(1)))));
    xscale = xscale(window);
    g = g(window);
    f = f(:, window);

function [ph0 ph1] = autophase(ph0)
global c;
global f;
global g;
global xscale;
    done = 0;
    dph = .001;
    while (~done)
        gw0 = g * exp(1i * ph0);
        gwp = g * exp(1i * (ph0 + dph));
        gwn = g * exp(1i * (ph0 - dph));
        sngw0 = sum(real(gw0)<0);
        sngwp = sum(real(gwp)<0);
        sngwn = sum(real(gwn)<0);
        if (sngw0==0 || sngwp==0 || sngwn==0)
            n0 = min(real(gw0));
            np = min(real(gwp));
            nn = min(real(gwn));
        else
            n0 = sum((real(gw0) < 0) .* real(gw0));
            np = sum((real(gwp) < 0) .* real(gwp));
            nn = sum((real(gwn) < 0) .* real(gwn));
        end
        if (n0 > np && n0 > nn)
            done = 1;
        elseif (np > n0)
            ph0 = ph0 + dph;
        else
            ph0 = ph0 - dph;
        end
        figure(1); plot(real(g*exp(1i*ph0))); pause(0.5);
    end
    ph1 = 0;
    % apply phase
    spectrum_baseline_points = floor(length(g) / 20);
    g = g .* exp(1i * ph0);
    g = g - real(sum(g(1:spectrum_baseline_points) + g((end-spectrum_baseline_points + 1):end)))/(2*spectrum_baseline_points);
    for k = 1:size(c, 1)
        f(k,:) = f(k,:) .* exp(1i * ph0);
        f(k,:) = f(k,:) - real(sum(sum(f(k,1:spectrum_baseline_points) + f(k,(end-spectrum_baseline_points + 1):end))))/(2*spectrum_baseline_points);
    end

function a = widthfit(xp)
global line;
global xscale;
global current_spec;
    for k = 1:length(line)
        line(k).width = xp;
    end
    yp = linefiteval;
    err = real(current_spec) - yp;
    a = sum(err.*err);

  
function yp = linefiteval()
global line;
global xscale;
    yp = zeros(1,length(xscale));
    for p=1:length(line)
        for k = 1:length(xscale)
            %yp(k) = yp(k) + line(p).amplitude/(1+(xscale(k)-line(p).center)^2/line(p).width^2);
            yp(k) = yp(k) + line(p).amplitude*exp(-(xscale(k)-line(p).center)^2/line(p).width^2);
        end
    end

function a = linefit_fixed(xp)
global line;
global current_line;
global current_spec;
    store = line(current_line).amplitude;
    line(current_line).amplitude = xp(1);
    yp = linefiteval;
    err = yp - current_spec;
    a = sum(err.*err);
    line(current_line).amplitude = store;
  
function a = linefit_movable(xp)
global line;
global current_line;
global current_spec;
    store_amplitude = line(current_line).amplitude;
    store_center = line(current_line).center;
    line(current_line).amplitude = xp(1);
    line(current_line).center = xp(2);
    yp = linefiteval;
    err = yp - current_spec;
    a = sum(err.*err);
    line(current_line).amplitude = store_amplitude;
    line(current_line).center = store_center;

function a=vecnorm(c)
for j=1:length(c)
    a(j)=norm(c(j));
end