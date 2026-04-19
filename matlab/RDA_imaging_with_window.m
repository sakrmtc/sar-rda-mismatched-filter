%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              side view
%                RDA
%              point target simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
% clear;
% close all;
% clc;

%%
% --------------------------------------------------------------------
% define parameters
% --------------------------------------------------------------------
R_nc = 20e3;            % Slope distance from center of view
Vr = 150;               % Radar effective speed
Tr = 2.5e-6;            % Transmit pulse duration
Kr = 15e12;             % Distance modulation frequency
f0 = 5.3e9;             % Radar operating frequency
BW_dop = 80;            % Doppler bandwidth
Fr = 60e6;              % distance sampling rate
Fa = 200;               % Azimuth sampling rate
Naz = 1024;          	% Number of distance lines (ie data matrix, number of rows) - modified to 1024 here.
Nrg = 320;              % The number of sample points for the distance line (ie, the data matrix, the number of columns)
sita_r_c = (0*pi)/180;	% Beam oblique angle, 0 degrees, converted to radians here
c = 3e8;                % speed of light

R0 = R_nc*cos(sita_r_c);	% The nearest slope distance corresponding to R_nc, denoted as R0
Nr = Tr*Fr;             % The number of sample points of the chirp signal


BW_range = Kr*Tr;       % range bandwidth
lamda = c/f0;           % wavelength
fnc = 2*Vr*sin(sita_r_c)/lamda;     % Doppler center frequency, calculated according to Equation (4.33).
La_real = 0.886*2*Vr*cos(sita_r_c)/BW_dop; %Azimuth antenna length, according to equation (4.36)
beta_bw = 0.886*lamda/La_real;           % Radar 3dB beam
La = beta_bw*R0;        % Synthetic Aperture Length
a_sr = Fr / BW_range;   % distance oversampling factor
a_sa = Fa / BW_dop;     % Azimuth oversampling factor

Mamb = round(fnc/Fa);   % Doppler blur

NFFT_r = Nrg;           % Range FFT length
NFFT_a = Naz;           % Azimuth FFT length

% --------------------------------------------------------------------
% Set the location of the simulation point target
% Take the distance direction as the positive direction of the x-axis
% Take the azimuth direction as the positive direction of the y-axis
% -------------------------------------------------------------------- 
delta_R0 = 0;       % The time when the beam center of target 1 crosses is defined as the azimuth time zero.
delta_R1 = 120; 	% Azimuth distance difference between target 1 and target 2, 120m
delta_R2 = 50;      % Distance to distance difference between target 2 and target 3, 50m

% Goal 1
x1 = R0;            % distance to target 1
y1 = delta_R0 + x1*tan(sita_r_c);	% Azimuth distance to target 1

% Goal 2
x2 = x1;            % Target 2 and target 1 have the same distance to distance
y2 = y1 + delta_R1; % Azimuth distance to target 2
% Goal 3
x3 = x2 + delta_R2;                 % There is a distance difference between target 3 and target 2, which is delta_R2
y3 = y2 + delta_R2*tan(sita_r_c);  	% Azimuth distance to target 3
% Define the following arrays for easy handling
x_range = [x1,x2,x3];
y_azimuth = [y1,y2,y3];

% Calculate the beam center crossing time of each of the three targets
nc_1 = (y1-x1*tan(sita_r_c))/Vr;    % The moment when the beam center of target 1 crosses.
nc_2 = (y2-x2*tan(sita_r_c))/Vr;    % The beam center crossing moment of target 2.
nc_3 = (y3-x3*tan(sita_r_c))/Vr;    % The beam center crossing moment of target 3.
nc_target = [nc_1,nc_2,nc_3];       %This array is defined for easy processing.

%% 
% --------------------------------------------------------------------
% Range (azimuth) versus time, frequency dependent definitions
% --------------------------------------------------------------------
% distance
tr = 2*x1/c + ( -Nrg/2 : (Nrg/2-1) )/Fr;                % distance timeline
fr = ( -NFFT_r/2 : NFFT_r/2-1 )*( Fr/NFFT_r );          % distance frequency axis
% position
ta = ( -Naz/2: Naz/2-1 )/Fa;                            % Azimuth Timeline
fa = fnc + ( -NFFT_a/2 : NFFT_a/2-1 )*( Fa/NFFT_a );	% Azimuth frequency axis

% Generate a distance (azimuth) time (frequency) matrix
tr_mtx = ones(Naz,1)*tr;    % Distance time axis matrix, size: Naz*Nrg
ta_mtx = ta.'*ones(1,Nrg);  % Azimuth time axis matrix, size: Naz*Nrg

%% 
% --------------------------------------------------------------------
% Generate point target raw data
% --------------------------------------------------------------------
s_echo = zeros(Naz,Nrg);    % Used to store the generated echo data

A0 = 1;                     % The target echo amplitude is set to 1.
for k = 1:1                 % Generate raw echo data for k targets
    R_n = sqrt( (x_range(k).*ones(Naz,Nrg)).^2 + (Vr.*ta_mtx-y_azimuth(k).*ones(Naz,Nrg)).^2 );% Instantaneous slope distance of target k
    w_range = ((abs(tr_mtx-2.*R_n./c)) <= ((Tr/2).*ones(Naz,Nrg)));     % distance envelope, the distance window
    % =====================================================================    
    % The azimuth envelope, that is, the two-way pattern contribution factor of the antenna.
    %{
    % way 1
    % sinc square function, calculated according to equation (4.31)   
    sita = atan( Vr.*(ta_mtx-nc_target(k).*ones(Naz,Nrg))/x_range(k) );
    w_azimuth = (sinc(0.886.*sita./beta_bw)).^2;    
    % Use the beam center crossing time corresponding to each target instead of nc in the previous parameter.
    %}
    %
    % 方式2
    % Use the synthetic aperture length to directly construct a rectangular window (in fact, this is just to limit the data range, no real windowing)
    w_azimuth = (abs(ta - nc_target(k)) <= (La/2)/Vr);    % row vector
    w_azimuth = w_azimuth.'*ones(1,Nrg);    % Generate a matrix of Naz*Nrg
    %}
    % =====================================================================     
    s_k = A0.*w_range.*w_azimuth.*exp(-(1j*4*pi*f0).*R_n./c).*exp((1j*pi*Kr).*(tr_mtx-2.*R_n./c).^2);
    % The above formula is the generated echo signal of a certain point target (target k).
    % After several cycles, the echo signals of several point targets are generated and added.
    if k == 1
        s_1 = s_k;          % Echo signal of target 1
    end
    if k == 2   
        s_2 = s_k;          % Echo signal of target 2
    end
    if k == 3
        s_3 = s_k;          % Echo signal of target 3
    end
    s_echo = s_echo + s_k;  % Sum of all point target echo signals  
end
% s_echo is the original data we need, the point target echo signal。

% drawing
% Figure 1 - Raw data
% figure;
% subplot(2,2,1);
% imagesc(real(s_echo));
% title('(a) Real part');
% xlabel('Distance Time Domain (Sampling Points)');
% ylabel('Azimuth time domain (sample points)');
% text(300,-60,'Figure 1. Raw data');       % Text description for Figure 1
% text function: add a text box at the specified coordinate position of the image

% subplot(2,2,2);
% imagesc(imag(s_echo));
% title('(b) Imaginary part');
% xlabel('Distance time domain (sample point）');
% ylabel('Azimuth time domain (sample points）');

% subplot(2,2,3);
% imagesc(abs(s_echo));
% title('(c) Amplitude');
% xlabel('Distance Time Domain (Sampling Points)');
% ylabel('Azimuth time domain (sample points)');

% subplot(2,2,4);
% imagesc(angle(s_echo));
% title('(d) Phase');
% xlabel('Distance Time Domain (Sampling Points)');
% ylabel('Azimuth time domain (sample points)');
% colormap(gray);

% 
% figure;
% subplot(2,2,1);
% imagesc(abs(fft(s_echo,[],1)));
% title('RD spectral magnitude');
% subplot(2,2,2);
% imagesc(angle(fft(s_echo,[],1)));
% title('RD Spectral Phase');
% subplot(2,2,3);
% imagesc(abs(fft2(s_echo)));
% title('2D spectral magnitude');
% subplot(2,2,4);
% imagesc(angle(fft2(s_echo)));
% title('2D spectral phase');
% colormap(gray);

%%
% --------------------------------------------------------------------
% range compression
% --------------------------------------------------------------------
S_range = fft(s_echo,NFFT_r,2);     % Perform a distance-to-Fourier transform with zero frequencies at both ends

%
% drawing
% Figure 2 - Range Frequency Domain, Azimuth Time Domain, Spectrum (Not Range Compressed)
% figure;
% subplot(1,2,1);
% imagesc(real(S_range));
% title('(a) Real part');
% xlabel('Distance frequency domain (sample points)');
% ylabel('Azimuth time domain (sample points)');
% text(280,-60,'Figure 2. Distance frequency domain');       %Text description for Figure 2
% text(340,-10,'uncompressed');       
% 
% subplot(1,2,2);
% imagesc(abs(S_range));
% title('(b) Amplitude');
% xlabel('Distance frequency domain (sample points)');
% ylabel('Azimuth time domain (sample points)');
%}

%　Generating Range Matched Filters
% ====================================================
% Use method 2
% Time domain replica pulse, zero-fill at the end, fft, and then take the complex conjugate.
t_ref = ( -Nr/2 : (Nr/2-1) )/Fr;    % Distance timeline used to generate distance MF
t_ref_mtx = ones(Naz,1)*t_ref;      % matrix form
w_ref = kaiser(Nr,2.5);             % Distance direction, build Kaiser window, this is a column vector.
%w_ref = hann(Nr)
%w_ref = taylorwin(Nr);
w_ref = ones(Naz,1)*(w_ref.');      % In matrix form, each row is equally windowed.

%s_ref = exp((1j*pi*Kr).*((t_ref_mtx).^2)); % Duplicate (transmit) pulse, unwindowed.
 s_ref = w_ref.*exp((1j*pi*Kr).*((t_ref_mtx).^2)); % Duplicate (transmit) pulse, windowed.

s_ref = [s_ref,zeros(Naz,Nrg-Nr)];      % For copy pulses, the trailing end is filled with zeros.
 
S_ref = fft(s_ref,NFFT_r,2);            %The distance Fourier transform of the replicated pulse, with zero frequency at both ends.
H_range = conj(S_ref);                  % Distance matched filter with zero frequency at both ends.
S_range_c = S_range.*H_range;           % Multiply by a matched filter, with zero frequencies at both ends.   
s_rc = ifft(S_range_c,[],2);            % Complete the distance compression and go back to the two-dimensional time domain.
% The length of s_rc is: Naz*Nrg. The disposal area was not removed.

% Perform the operation of removing the waste area on s_rc
%The length of the disposal area is：2*（Nr-1）
% The length we intercepted: (Nrg-Nr+1), denoted as N_rg.
N_rg = Nrg-Nr+1;                        % length of full convolution
s_rc_c = zeros(Naz,N_rg);               % Used to store data after removing the trash area
s_rc_c = s_rc(:,1:N_rg);                % Take the first N_rg columns.
% ====================================================

%
% drawing
% Figure 3 - Range frequency domain, azimuth time domain, spectrum (range compressed)
% figure;
% subplot(1,2,1);
% imagesc(real(S_range_c));
% title('(a) Real part');
% xlabel('Distance frequency domain (sample points)');
% ylabel('Azimuth time domain (sample points)');
% text(280,-60,'Figure 3. Distance frequency domain');       % Text description for Figure 3
% text(340,-10,'compressed');       
% 
% subplot(1,2,2);
% imagesc(abs(S_range_c));
% title('(b) Amplitude');
% xlabel('Distance frequency domain (sample points)');
% ylabel('Azimuth time domain (sample points)');
%}
%
% drawing
% Figure 4 - 2D Time Domain (Complete Distance Compression)
% figure;
% subplot(1,2,1);
% imagesc(real(s_rc_c));  %　This and the following, all directly use the result after removing the waste area
% title('(a) Real part');
% xlabel('Distance Time Domain (Sampling Points)');
% ylabel('Azimuth time domain (sample points)');
% text(150,-60,'Figure 4. Two-dimensional time domain');       % Text description for Figure 4
% text(172,-10,'complete compression');       

% subplot(1,2,2);
% imagesc(abs(s_rc_c));
% title('(b) Amplitude');
% xlabel('Distance Time Domain (Sampling Points)');
% ylabel('Azimuth time domain (sample points)');
%}

%%
% --------------------------------------------------------------------
% Transform to range Doppler domain for range migration correction
% --------------------------------------------------------------------
s_rc_c = s_rc_c.*exp(-1j*2*pi*fnc.*(ta.'*ones(1,N_rg)));    % data movement
S_rd = fft(s_rc_c,NFFT_a,1);            %Azimuth Fourier Transform, to Range Doppler Domain
% ====================================================================
% Set Azimuth Frequency Axis - Key Points! ! !
fa = fnc + fftshift(-NFFT_a/2:NFFT_a/2-1)/NFFT_a*Fa;    %The azimuth frequency axis is set as such.
% =====================================================================
% The following one is an improvement, each closest slope distance (R0) varies with the distance gate.


tr_RCMC = 2*x1/c + ( -N_rg/2 : (N_rg/2-1) )/Fr;   % Timeline under the new distance line length.

R0_RCMC = (c/2).*tr_RCMC;       % The R0 that varies with the distance line, denoted as R0_RCMC, is used to calculate RCM and Ka.
delta_Rrd_fn = lamda^2.*((fa.').^2)*(R0_RCMC)/(8*Vr^2);

num_range = c/(2*Fr);   % A distance sampling unit, the corresponding length
delta_Rrd_fn_num = delta_Rrd_fn./num_range; % For each azimuth frequency, the number of distance sampling units corresponding to its RCM

R = 8;  % sinc interpolation kernel length
S_rd_rcmc = zeros(NFFT_a,N_rg); % Used to store the value after RCMC
for p = 1 : NFFT_a
    for q = 1 : N_rg   % At this time, the length of the distance direction is (Nrg-Nr+1)=N_rg        
        delta_Rrd_fn_p = delta_Rrd_fn_num(p,q);
        Rrd_fn_p = q + delta_Rrd_fn_p;
        Rrd_fn_p_zheng = ceil(Rrd_fn_p);        % ceil, rounded up.
        ii = ( Rrd_fn_p-(Rrd_fn_p_zheng-R/2):-1:Rrd_fn_p-(Rrd_fn_p_zheng+R/2-1)  );        
        rcmc_sinc = sinc(ii);
        rcmc_sinc = rcmc_sinc/sum(rcmc_sinc);   % Normalization of the interpolation kernel
        % ii are variables of the sinc interpolation process;
        % g(x)=sum(h(ii)*g_d(x-ii)) = sum(h(ii)*g_d(ll));
               
        % Since S_rd has only integer point values, and the range is limited. Therefore, its value overflow boundary problem should be considered in the interpolation.
        % Here I take the idea of ​​circular shift to solve the value overflow problem.
        if (Rrd_fn_p_zheng-R/2) > N_rg    % Full right overflow
            ll = (Rrd_fn_p_zheng-R/2-N_rg:1:Rrd_fn_p_zheng+R/2-1-N_rg);
        else
            if (Rrd_fn_p_zheng+R/2-1) > N_rg    % partial right overflow
                ll_1 = (Rrd_fn_p_zheng-R/2:1:N_rg);
                ll_2 = (1:1:Rrd_fn_p_zheng+R/2-1-N_rg);
                ll = [ll_1,ll_2];
            else
                if (Rrd_fn_p_zheng+R/2-1) < 1    % Full left overflow (unlikely, but still something to consider)
                    ll = (Rrd_fn_p_zheng-R/2+N_rg:1:Rrd_fn_p_zheng+R/2-1+N_rg);
                else
                    if (Rrd_fn_p_zheng-R/2) < 1       % partial left overflow
                        ll_1 = (Rrd_fn_p_zheng-R/2+N_rg:1:N_rg);
                        ll_2 = (1:1:Rrd_fn_p_zheng+R/2-1);
                        ll = [ll_1,ll_2];
                    else
                        ll = (Rrd_fn_p_zheng-R/2:1:Rrd_fn_p_zheng+R/2-1);
                    end                    
                end
            end
        end   
        rcmc_S_rd = S_rd(p,ll);
        S_rd_rcmc(p,q) = sum( rcmc_sinc.*rcmc_S_rd );
    end
end
% S_rd_rcmc is the range Doppler domain spectrum after RCMC.

%drawing
% Figure 5 - Range Doppler Domain (not RCMC)
% figure;
% subplot(1,2,1);
% imagesc(real(S_rd));
% title('（a) Real part');
% xlabel('Distance Time Domain (Sampling Points)');
% ylabel('Azimuth frequency domain (sample points)');
% text(150,-60,'Figure 5. Range Doppler Domain');       % Text description for Figure 5
% text(172,-10,'Not RCMC');       
% subplot(1,2,2);
% imagesc(abs(S_rd));
% title('(b) Amplitude');
% xlabel('Distance time domain (sample point）');
% ylabel('Azimuth frequency domain (sample points)');

% drawing
% Figure 6 - Range Doppler domain, results after RCMC
% figure;
% subplot(1,2,1);
% imagesc(real(S_rd_rcmc));
% title('(a) Real part');
% xlabel('Distance Time Domain (Sampling Points)');
% ylabel('Azimuth frequency domain (sample points)');
% text(150,-60,'Figure 6. Range Doppler Domain');       % Text description for Figure 6
% text(172,-10,'RCMC');       
% 
% subplot(1,2,2);
% imagesc(abs(S_rd_rcmc));
% title('(b) Amplitude');
% xlabel('Distance Time Domain (Sampling Points)');
% ylabel('Azimuth frequency domain (sample points)');
%}

%%
% --------------------------------------------------------------------
% Azimuth compression
% --------------------------------------------------------------------
fa_azimuth_MF = fa;         % The azimuth frequency axis is the same as the frequency axis used in RCMC.
Ka = 2*Vr^2*(cos(sita_r_c))^3./(lamda.* R0_RCMC);  	% The azimuth modulation frequency changes with the nearest slope distance R0.
Ka_1 = 1./Ka;                                       % For the convenience of calculation, take the reciprocal first.
Haz = exp( -1j*pi.*(((fa_azimuth_MF).').^2*Ka_1) );	% Azimuth matched filter
% It should be noted here that the zero frequency of the generated MF is neither at the ends nor at the center.
%Consider what the frequency axis looks like and where the discontinuities are. Note the composition of fa.
% The frequency axis here corresponds to the azimuth spectrum in the range Doppler domain.

S_rd_c = S_rd_rcmc.*Haz;            % Multiply by the matched filter
s_ac = ifft(S_rd_c,[],1);       	% Complete the azimuth compression, change to the image domain. Finish.

% drawing
% Figure 7 - Imaging results
% figure;
% imagesc(abs(s_ac));
% title('point target imaging');
% xlabel('Distance time domain (sample point）');
% ylabel('Azimuth time domain (sample points)');     

%%
% Next, by calling the function, the respective slices of the three point targets are obtained and upsampled
% At the same time, make distance slices and azimuth slices for the center of the point target
% Calculate the corresponding indicators: PSLR, ISLR, IRW

NN = 16;
% Get the slice magnification of each point target separately; row slice, column slice; and corresponding indicators
% Target 1, point target center at ((Naz/2+1), 86)
target_1 = target_analysis( s_ac((Naz/2+1)-NN:(Naz/2+1)+NN,86-NN:86+NN),Fr,Fa,Vr);

% Target 2, point target center at （（Naz/2+1+delta_R1/Vr*Fa）,86）
tg_2_delatx = (Naz/2+1 + delta_R1/Vr*Fa);
% target_2 = target_analysis( s_ac(tg_2_delatx-NN:tg_2_delatx+NN,86-NN:86+NN),Fr,Fa,Vr);

% Goal 3
tg_3_delatx = tg_2_delatx + delta_R2*tan(sita_r_c)/Vr*Fa;
tg_3_delaty = 2*delta_R2/c*Fr;
% target_3 = target_analysis( s_ac(tg_3_delatx-NN:tg_3_delatx+NN,86+tg_3_delaty-NN:86+tg_3_delaty+NN),Fr,Fa,Vr);


