% CFANVC
%
% Code for the work:
% Paulo A. C. Lopes, José A. B. Gerald
% Careful Feedback Active Noise and Vibration Control Algorithm Robust to Large Secondary Path Changes
% European Control Journal, to be published.
% 
% Paulo Alexandre Crisóstomo Lopes, 26 June 2023
% Instituto de Engenharia de Sistemas e Computadores - Investigação e Desenvolvimento
% Instituto Superior Técnico
% Universidade de Lisboa

% simulation parameters
Nsim = 10;                  % number of simulations
fs = 2000;                  % sampling frequency
simulation_time = 10*fs;    % simulation samples
change_at = simulation_time/2;
change_time = 1*fs;         % 1s fast
qn_steady = 0*1/400/fs;     % state noise
qn_change = 2/change_time;  % relative state noise at change

stability_margin = 0.1;     % distance of model poles and zeros from the unit circle
f = 200*(1:3)';             % primary noise sinusoids frequencies
amplitudes = [0.5, 1.2, 0.3]';      % primary noise sinusoids amplitudes
phases = [56, 170, -23]'*pi/180;    % primary noise sinusoids phases
frequency_noise = 0;        % rms Hz
Nx = 6;                     % model order (size-1)
on_id = 100;                % system identification start
on = 1000;                  % ANC start
qv0 = 0.01;                 % background noise power

% algorithm parameters
N = 15;                     % model order (size-1)
L = 64;                     % the equations are 2*L, from i = -L+1 to L
M = 32*(N+1);               % algorithm memory
P = M;                      % size of qv estimation blocks
R = N;                      % number of samples used to obtain uM
Lu = 10;                    % saturation of the antinoise signal
wd = 10000;                 % cost weigth of past control signal changes
wu = 0;                     % cost weigth of control effort
deltax = 1e-3;              % prior on x variance
delta = 1e-9;               % actual value added Rxx to calc Sxx (<deltax)
alpha = 0;                  % controls auxiliary noise power
eta = 1.01;                 % max grouth rate of u (per Q samples)
Q = 8;                      % controller update period

% LOGS
log_e = nan*zeros(simulation_time, Nsim);
log_u = nan*zeros(simulation_time, Nsim);
log_xi = nan*zeros(simulation_time, Nsim);
log_xi0 = nan*zeros(simulation_time, Nsim);
log_qv = nan*zeros(simulation_time, Nsim);

% warning('off', 'MATLAB:singularMatrix');
% warning('off', 'MATLAB:nearlySingularMatrix');

for n_sim = 1:Nsim
    tic
    rng(343989 + n_sim);
   
    % simulation intializations
    frequency = f' + frequency_noise*randn(simulation_time+L, length(f));
    phase = 2*pi*cumsum(frequency)/fs + phases';
    d0 = sin(phase)*amplitudes;
    d = d0 + sqrt(qv0)*randn(simulation_time+L, 1);  % primary noise signal

    [a,b] = generate_plant(Nx, stability_margin);
    
    uv = zeros(max(M+N+1, 2*L), 1);        % anti-noise buffer
    e1v = zeros(Nx+1, 1);    % residual noise minus background noise buffer

    % good values for ah and bh
    n = 0:N-length(a)+1;
    Pz = exp(2*pi*f/fs*n*1i);
    Pz = [real(Pz); imag(Pz)];
    Hz = exp(-2*pi*f*n/fs*1i)*Pz.';
    Hz = [real(Hz); imag(Hz)];
    Az = Hz\[ones(length(f),1); zeros(length(f),1)];
    pz = - Pz.'*Az; pz(1) = pz(1) + 1;
    pz = pz/pz(1);

    % algorithm initialization
    u = 0;                   % anti-noise signal
    ev = zeros(M+N,1);       % residual noise buffer
    wuv = wu*ones(L+N,1);
    wev = ones(L,1);         % cost weigth of each residual noise sample
    bc = zeros(1,N);
    ac = zeros(1,N);
    udv = zeros(L+N,1);

    for k = 1:simulation_time
        
        % simulation
        qn = qn_steady + qn_change*(abs(k-change_at)<=change_time/2);
        a(2:end) = a(2:end) + rms(a(2:end))*sqrt(qn)*randn(Nx,1);
        b = b + rms(b)*sqrt(qn)*randn(Nx+1,1);
        if qn > 0
            [a,b] = adjust_plant(a,b,stability_margin);
        end
%         if k == change_at
%             b = - b;
%         end
        
        log_u(k,n_sim) = u;    % logs u(n) and not u(n+1)
              
        uv = [u; uv(1:end-1)]; % simulation and algorithm
        e1v = [0; e1v(1:end-1)];
        e1 = b'*uv(1:Nx+1) - a'*e1v;
        e1v(1) = e1;
        e = e1 + d(k);
        
        % algorithm
        ev = [e; ev(1:end-1)];
        if k >= on_id

            if k>=on

                if mod(k, Q) == 0

                    % Calculate the parameters
                    Em = ev((1:M)'+(1:N));
                    Um = uv((0:M-1)'+(1:N+1));
                    H = [-Em, Um];
                    % eb = H x
                    Rxx = H'*H;
                    Deltax = deltax*eye(2*N+1);
                    x = (Rxx+Deltax)\H'*ev(1:M);
                    ab = x(1:N); bb=x(N+1:end);

                    qvM = max(mean(reshape(abs(ev(1:M)-H*x).^2, P, [])));
        
                    Delta = delta*eye(2*N+1);
                    Sxx = qvM*inv(Rxx + Delta);
                    
                    % good values for ah and bh
                    %a0 = conv(a,pz); a0 = a0(2:end); a0 = [a0; zeros(N-length(a0),1)];
                    %b0 = conv(b,pz); b0 = [b0; zeros(N+1-length(b0),1)];
            
                    Saa = Sxx(1:N, 1:N);
                    Saax = [zeros(N+1,1), [zeros(1,N); Saa]];
                    Sbb = Sxx(N+1:end, N+1:end);
                    Sab = Sxx(1:N, N+1:end);
                    Sba = Sxx(N+1:end, 1:N);
                    Sbax = [zeros(N+1,1), Sba];
                    Qbb = bb*bb'+ Sbb;
                    Qba = bb*[1;ab]' + Sbax;
                    Qaa = [1;ab]*[1;ab]' + Saax;
                    Qab = Qba';

                    Imx = diag([1./wev; zeros(N,1)]);
                    Raax = conv2(Imx, flip(flip(Qaa),2), 'valid');
                    RaaI = inv(Raax);
                    Raa = conv2(RaaI, Qaa);
                    Rbb = conv2(RaaI, Qbb) + diag(wuv);
                    Rba = conv2(RaaI, Qba);

                    wv = [zeros(L,1); ones(N,1)];
                    Wm = diag(wv);
                    Rmd = decomposition(Rbb + wd*Wm);
                    Bc = Rmd\Rba;
                    Ac = Rmd\Wm*wd;
                    bc0 = bc; ac0 = ac;
                    bc = Bc(L, L+1:end);
                    ac = - Ac(L, L+1:end);

                    edv = zeros(L+N, 1);
                    edv(L+1:end) = ev(1:N);
                    u0v = [zeros(L,1); uv(1:N)];
                    udv0 = udv;
                    udv = Bc*edv + Ac*u0v;
                    xi = edv'*Raa*edv + udv'*Rbb*udv - 2*udv'*Rba*edv;
                    xi0 = edv'*Raa*edv + u0v'*Rbb*u0v - 2*u0v'*Rba*edv;
                    qr = alpha*xi/trace(Rbb);
                    uM = eta*max(abs(uv(1:R)));
                end
                % u = udv0(L-mod(k,Q)-Q);
                u = bc0*ev(1:N) - ac0*uv(1:N);
                u = u + sqrt(qr)*randn;
                u = min(Lu, max(-Lu,u));
                u = min(uM, max(-uM,u));

                log_qv(k,n_sim) = qvM;
                log_xi(k,n_sim) = xi;
                log_xi0(k,n_sim) = xi0;
            else
                u = randn;
            end
        else
            u = 0;
        end

        log_e(k,n_sim) = e;

    end

    qe = mean(log_e(end-simulation_time/10:end,n_sim).^2);
    fprintf(1, 'sim: %d, residual noise power: %f\n', n_sim, qe);
    toc
end

warning('on', 'MATLAB:singularMatrix');
warning('on', 'MATLAB:nearlySingularMatrix');

figure(1)
plot_xy_p3((0:size(log_e,1)-1)/fs, 10*log10(sline(log_e.^2,min(round(simulation_time/200),200))));
set(gca, 'YLim', [-20,30]);
xlabel('time (s)');
ylabel('Noise (dB)');
grid on;
%title('Residual noise versus time percentile plot');

ix = (1:min(Nsim,10))+0;
figure(2)
plot((0:size(log_e,1)-1)/fs, 10*log10(sline(log_e(:,ix).^2,min(round(simulation_time/200),200))));
set(gca, 'YLim', [-20,30]);
xlabel('time (s)');
ylabel('Noise (dB)');
grid on;
%legend('simulation 1', 'simulation 2', 'simulation 3');
%title('Residual noise versus time plot of 3 simulations');
%legend('1','2','3','4','5','6','7','8','9','10');

figure(3);
histogram(10*log10(mean(log_e(end-simulation_time/10:end,:).^2)),[-25:0.5:0,inf]);
xlabel('Noise Power (dB)');
ylabel(['Frequency (out of ', num2str(Nsim),')']);
grid on;
%title('Final residual noise power histogram');

figure(4);
plot_xy_p3((0:size(log_qv,1)-1)/fs, 10*log10(sline(log_qv,100)));
grid on;
title('Model identification error signal');

ix = 1:min(Nsim,10);
figure(5)
plot((0:size(log_xi0,1)-1)/fs, 10*log10(sline(log_xi0(:,ix)/sum(wev), min(round(simulation_time/200),200))));
hold on;
plot((0:size(log_xi0,1)-1)/fs, 10*log10(sline(log_xi(:,ix)/sum(wev),min(round(simulation_time/200),200))));
hold off;
xlabel('time (s)');
ylabel('Noise (dB)');
grid on;
title('xi and xi0');
legend('xi0', 'xi');
