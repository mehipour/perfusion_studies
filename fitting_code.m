function fitting_code(Mpyr1,Mlac1,Fr_avg1,chemshift1,pyr_index1,lac_index1)
global Mpyr Mlac counter TR pyr_ratio ii Fr_avg chemshift pyr_index lac_index

ii = 1;
Mpyr = Mpyr1;
Mlac = Mlac1;
Fr_avg = Fr_avg1;
chemshift = chemshift1;
counter = 0;
pyr_index = pyr_index1;
lac_index = lac_index1;
TR = 1;
pyr_ratio = 0.05;
c1_ode_nonlin;

[p(:,ii),l(:,ii)] = deal(x(:,1),x(:,2));
figure(3); hold on;
plot(t,p(:,ii)/sum(p(:,ii))*pyr_ratio,'--',t,l(:,ii)/sum(p(:,ii)),'r--');
[p(:,ii),l(:,ii)] = deal(x(:,1),x(:,2));
figure(3); hold on;
plot(t,p(:,ii)/sum(p(:,ii))*pyr_ratio,'--',t,l(:,ii)/sum(p(:,ii)),'r--');
hold off;hold off;


% convert figures to PNG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions

function [t x] = c1_ode_lsq
% Fit kinetic model to C1-pyr injection spectra to estimate apparent
% forward and backward rate constants using LSQ curve fitting. 
global kopt Nmetabolites Tspan tstart
    tic
    find_snr_peaks()
    Tspan = 0:400;
    Nmetabolites = 2;
    [p0,l0] = deal(0,0);                   % Initial magnetizations
    % k0 = [1,0.5,1];                  % initial kinetic parameters.
    % k0 = [0.0035,0.03,1];
%     k0 = [0.008040713450315, 0.066500372252357, 0.702747014971904,9.392859384899749,5.285189191197075];
%     k0 = [0.008040713450315, 0.066500372252357, 0.702747014971904,9.392859384899749,5.285189191197075,55];
    k0 = [0.008040713450315, 0.066500372252357, 0.702747014971904,9.392859384899749,5.285189191197075,55,10];

    kopt = fminsearch(@solve_c1_kinetic_model_lsq,k0); % solve for optimum k values.
    disp(kopt');
    kopt1 = ones(Nmetabolites,1)*kopt; % optimum value of k with adjusted dimensions
    [t x] = ode45(@c1_kinetic_model,Tspan,[p0,l0],1,kopt1);
    toc
    % Injection occurs for 120 seconds starting at ti for dur seconds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve equations
function J = solve_c1_kinetic_model_lsq(y)
% Solve differential equations and define cost function. 
global t ii Mpyr Mlac Tspan pyr_ratio bic_ratio tstart tend counter snrp snrl theta_rf
global T1_pyr_fit T1_lac_fit k_pyr2lac_fit k_lac2pyr_fit alpha_fit cost_func T1_pyr_inj
    [p0,l0] = deal(0,0);                     % Initial magnetizations
    y1 = ones(2,1)*y; % inital values for [kp,kl,alpha] with adjusted dimensions
    [t x] = ode45(@c1_kinetic_model,Tspan,[p0 l0],1,y1);
    % x = interp1(t,x1,[0:299]);
    J = sum(((x(:,1)/sum(x(:,1))-Mpyr(1:200,ii))/snrp).^2)+...
        10*sum(((x(:,2)/sum(x(:,1))-Mlac(1:200,ii))/snrl).^2);
    %     -0.01*(y(1,1)+y(1,2)+y(1,3)+y(1,4)+y(1,5));   % cost function (simultaneous fitting)
    figure(10);
    subplot(224);
    plot(pyr_ratio*Mpyr(:,ii)); hold on;
    plot(pyr_ratio*x(:,1)/sum(x(:,1)),'--');
    plot(Mlac(:,ii),'r');  
    plot(x(:,2)/sum(x(:,1)),'r--'); hold off;
    legend(strcat('Pyruvate (x ',num2str(pyr_ratio),')'),...
        strcat('Pyruvate Fit (x ',num2str(pyr_ratio),')'),...
        'Lacate/Pyruvate','Lacate/Pyruvate Fit');
    title('Fitting'); ylabel('A.U.'); xlabel('Time (sec)');
    counter = counter +1;
    % create fitting parameter vectors for visualization. 
    k_pyr2lac_fit(counter) = y1(1,1);
    k_lac2pyr_fit(counter) = y1(1,2);
    alpha_fit(counter) = y1(1,3);
    % time_start(counter) = y(1,6);
    T1_pyr_fit(counter) = y1(1,4);
    T1_lac_fit(counter) = y1(1,5);
    T1_pyr_inj(counter) = y(1,6);
    theta_rf(counter) = y(1,7);
    cost_func(counter) = J;
    subplot(241)
    plot(k_pyr2lac_fit,'b'); title('Rate Constants'); hold on;
    plot(k_lac2pyr_fit,'r'); xlabel('iteration'); ylabel('s^{-1}'); hold off;
    legend('k+','k^-');
    subplot(242)
    plot(T1_pyr_fit,'b'); title('T_1 Relaxation Times'); hold on;
    plot(T1_lac_fit,'r'); xlabel('iteration'); ylabel('seconds'); hold off;
    legend('T_{1-pyr}','T_{1-lac}');
    subplot(243)
    [haxes,hline1,hline2] = plotyy(1:length(alpha_fit),alpha_fit,1:length(alpha_fit),T1_pyr_inj);
%     title('Injection Parameters');
    axes(haxes(1))
    ylabel('Injection Constant (A.U.)');
    axes(haxes(2));
    ylabel('T_1 of Injected Pyruvate');
    xlabel('iteration');
    subplot(244)
    plot(theta_rf,'b'); title('RF pulse flip angle'); hold on;
    xlabel('iteration'); ylabel('\theta (degrees)'); hold off;
    subplot(223);
    semilogy(cost_func); title('Cost Function of Non-linear Fit');
    ylabel('A.U.'); xlabel('iteration');
    drawnow();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function [t x] = c1_ode_nonlin
% Fit kinetic model to C1-pyr injection spectra to estimate apparent
% forward and backward rate constants using non-linear fit.    
global kopt Nmetabolites Tspan tstart Mpyr Mlac data t snrp snrl
    tic
    find_snr_peaks()   % find snr of the peaks to weight MSE of peaks
    Tspan = 0:199;
    Nmetabolites = 2;
    [p0,l0] = deal(0,0);                   % Initial magnetizations
    % k0 = [1,0.5,1];                  % initial kinetic parameters.
    % k0 = [0.0035,0.03,1];
    k0 = [0.008040713450315, 0.066500372252357, 0.702747014971904,9.392859384899749,5.285189191197075,68,10];
%     k0 = [0.008040713450315, 0.066500372252357, 0.702747014971904,9.392859384899749,5.285189191197075,68];
%     k0 = [0.008040713450315, 0.066500372252357, 0.702747014971904,9.392859384899749,5.285189191197075 55];
    t = 1:200;
    time = [t;t+max(t)+1];
    data = [Mpyr(1:length(t))/snrp;Mlac(1:length(t))/snrl];
    [kopt,R,J,COVB,MSE] = nlinfit(time,data,@solve_c1_kinetic_model_nl,k0); % solve for optimum k values.
    disp(kopt');
    disp(MSE)
    xe = nlparci(kopt,R,'jacobian',J);
    xe = (xe(:,2)-xe(:,1))/2;
    disp(xe);
    kopt1 = ones(Nmetabolites,1)*kopt;   % optimum value of k with adjusted dimensions
    [t x] = ode45(@c1_kinetic_model,Tspan,[p0,l0],1,kopt1);
    toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Solve equations
function output = solve_c1_kinetic_model_nl(y,time)
% Solve differential equations and define cost function. 
global t ii Mpyr Mlac Tspan pyr_ratio bic_ratio tstart tend counter theta_rf
global T1_pyr_fit T1_lac_fit k_pyr2lac_fit k_lac2pyr_fit alpha_fit cost_func T1_pyr_inj
global snrp snrl
    [p0,l0] = deal(0,0);                     % Initial magnetizations
    y1 = ones(2,1)*y; % inital values for [kp,kl,alpha] with adjusted dimensions
    [t x] = ode45(@c1_kinetic_model,Tspan,[p0 l0],1,y1);
    output = [x(:,1)/snrp;x(:,2)/snrl]/sum(x(:,1)); 
%     x = interp1(t,x1,[0:299]);
    J = sum(((x(:,1)/sum(x(:,1))-Mpyr(1:200,ii))/snrp).^2)+...
        10*sum(((x(:,2)/sum(x(:,1))-Mlac(1:200,ii))/snrl).^2);
%         -0.01*(y(1,1)+y(1,2)+y(1,3)+y(1,4)+y(1,5));   % cost function (simultaneous fitting)
    figure(11);
    subplot(224);
    plot(pyr_ratio*Mpyr(:,ii)); hold on;
    plot(pyr_ratio*x(:,1)/sum(x(:,1)),'--');
    plot(Mlac(:,ii),'r');  
    plot(x(:,2)/sum(x(:,1)),'r--'); hold off;
    legend(strcat('Pyruvate (x ',num2str(pyr_ratio),')'),...
        strcat('Pyruvate Fit (x ',num2str(pyr_ratio),')'),...
        'Lacate/Pyruvate','Lacate/Pyruvate Fit');
    title('Fitting'); ylabel('A.U.'); xlabel('Time (sec)');
    counter = counter +1;
    % create fitting parameter vectors for visualization. 
    k_pyr2lac_fit(counter) = y1(1,1);
    k_lac2pyr_fit(counter) = y1(1,2);
    alpha_fit(counter) = y1(1,3);
    % time_start(counter) = y(1,6);
    T1_pyr_fit(counter) = y1(1,4);
    T1_lac_fit(counter) = y1(1,5);
    T1_pyr_inj(counter) = y1(1,6);
    theta_rf(counter) = y1(1,7);
    cost_func(counter) = J;
    subplot(241)
    plot(k_pyr2lac_fit,'b'); title('Rate Constants'); hold on;
    plot(k_lac2pyr_fit,'r'); xlabel('iteration'); ylabel('s^{-1}'); hold off;
    legend('k+','k^-');
    subplot(242)
    plot(T1_pyr_fit,'b'); title('T_1 Relaxation Times'); hold on;
    plot(T1_lac_fit,'r'); xlabel('iteration'); ylabel('seconds'); hold off;
    legend('T_{1-pyr}','T_{1-lac}');
    subplot(243)
    [haxes,hline1,hline2] = plotyy(1:length(alpha_fit),alpha_fit,1:length(alpha_fit),T1_pyr_inj);
    title('Injection Parameters');
    axes(haxes(1))
    ylabel('Injection Constant (A.U.)');
    axes(haxes(2));
    ylabel('T_1 of Injected Pyruvate');
    xlabel('iteration');
    subplot(244)
    plot(theta_rf,'b'); title('RF pulse flip angle'); hold on;
    xlabel('iteration'); ylabel('\theta (degrees)'); hold off;
    subplot(223);
    semilogy(cost_func); title('Cost Function of Non-linear Fit');
    ylabel('A.U.'); xlabel('iteration');
    drawnow();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function dxdt = c1_kinetic_model(t,x,k)
% Define the differental equations sytem.
global tstart tend TR
    [T1,ti,Flow,dur]= deal(55,tstart,10,130);    % constant kinetic parameters.
    % [ti,Flow,dur]= deal(0,10,120);    % constant kinetic parameters.
    % disp(k);
    dxdt = zeros(2,1);
%     injection = k(1,3)*Flow*exp(-(t-ti)/k(1,6)).*((heaviside(t-ti)-heaviside(t-dur-ti)).*(t-ti));
%     injection = k(1,3)*Flow*exp(-(t-ti)/k(1,6)).*(heaviside(t-ti)-heaviside(t-dur-ti));
%     dxdt(1,1) = -k(1,1)*x(1) + k(1,2)*x(2) - x(1)/k(1,4) + injection - x(1)/TR*(1-cos(pi*k(1,7)/180));
%     dxdt(2,1) =  k(2,1)*x(1) - k(2,2)*x(2) - x(2)/k(2,5) - x(2)/TR*(1-cos(pi*k(1,7)/180));
    injection = k(1,3)*Flow*exp(-(t-ti)/T1).*((heaviside(t-ti)-heaviside(t-dur-ti)).*(t-ti));
    dxdt(1,1) = -k(1,1)*x(1) + k(1,2)*x(2) - x(1)/k(1,4) + injection;
    dxdt(2,1) =  k(2,1)*x(1) - k(2,2)*x(2) - x(2)/k(2,5);
%     dxdt(1,1) = -k(1,1)*x(1) + k(1,2)*x(2) - x(1)/k(1,4) + injection - x(1)/TR*(1-cos(pi*10/180));
%     dxdt(2,1) =  k(2,1)*x(1) - k(2,2)*x(2) - x(2)/k(2,5) - x(2)/TR*(1-cos(pi*10/180));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function find_snr_peaks
global chemshift Fr_avg noise_index snrp snrl pyr_index lac_index
    noise_index = intersect(find(chemshift>192),find(chemshift<195));
    snrp = sum(Fr_avg(pyr_index))/std(Fr_avg(noise_index));
    snrl = sum(Fr_avg(lac_index))/std(Fr_avg(noise_index));  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function FitT1(y,x0,tstart,tend)
% Fit the signal y beginning from time point tstart to tend to a single 
% exponential to estimate T1 relaxation time. 
    tstart = 120; tend = 200; 
    t = (tstart:tend)-tstart;
    x = fminsearch(@(x)sum((x(1)*exp(-(t-tstart)'/x(2))+x(3)-y(tstart:end)/max(y)).^2./(y(tstart:end)/max(y))),x0);
    x = fminsearch(@(x)sum((x(1)*exp(-(t-tstart)'/x(2))+x(3)-y(tstart:end)/max(y)).^2),x0);
    plot(y,'.'); hold on;
    plot([tstart:tend],x,'r');   