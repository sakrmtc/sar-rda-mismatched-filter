%% Optimization of Mismatched Filter using via QCQP Approach
clear all
clc
close all
Tr = 2.5e-6;            % Transmit pulse duration
Kr = 6e12;             % range modulation frequency
Fs = 60e6;              % range sampling rate
Nr = Tr*Fs;       % The number of sample points of the chirp signal
t_ref = ( -Nr/2 : (Nr/2-1) )/Fs;
s = exp((1j*pi*Kr).*((t_ref).^2)); %generate the lfm signal

s_ref_flip = flip(s);
NN = size(s, 2);    % No. of sampled signal vector
p = 1;
KK = NN;
TT = KK + NN - 1;

% Generate Toeplitz matrix
Matrix_A = zeros(KK, TT);
for j = 1:KK
    Matrix_A(j, j:length(s) + j - 1) = s_ref_flip;
end
build_T = Matrix_A';

% Generate the diagonal matrix for main lobe and side lobe control
Totallenth = TT + 1;
AcutallLenth = Totallenth - 9;
period = 2;

z = zeros(1, AcutallLenth);
for c = 1:AcutallLenth
    r = rem(c, period);
    if (r == 0)
        z(c) = 1;
    else
        z(c) = 1;
    end
end
breakDown = AcutallLenth / 2;
M = [z(1:breakDown) 0 0 0 0 0 0 0 0 0 z(breakDown + 1:end)];
F = diag(M);
delta = F * build_T;  % The output product of the Toeplitz matrix and diagonal matrix

% Ideal matched filter
matched = xcorr(s);
[ee, oo] = max(abs(matched));
b = zeros(1, 299);
b(146:154) = matched(146:154);
BB = b.^2;

% Plot the matched filter output
figure()  % Create a new figure for the matched filter output
plot(abs(matched))
mag2 = abs(matched);
max2 = max(mag2);
max2_n = mag2 / max2;
max2_n_dB = mag2db(max2_n);
figure()  % Create a new figure for the normalized matched filter output in dB
plot(max2_n_dB)
ylim([-60 0])
xlabel('time(in samples)')
ylabel('filter output')
legend('Matched filter')
grid on
hold on

%% Begin CVX model

cvx_begin 
   variable q(150, 1) complex 
   variable a
   minimize a
   subject to
       S_T' * q == S_T' * S_T
       q' * q <= Alpha * S_T' * S_T
        
   for i = 1:299
       if ((1 <= i) < 146) && ((154 < i) <= 299)
           (delta(i, :) * q)' * (delta(i, :) * q) < a; 
       elseif ((146 <= i) <= 154)
           (delta(i, :) * q)' * (delta(i, :) * q) <= BB(i);
       end
   end
cvx_end

% Plot the mismatched filter output
q_T = q.'; % Transpose from column to row
DD = xcorr(s, q_T);
figure()  % Create a new figure for the mismatched filter output
plot(abs(DD))
grid on
mag1 = abs(DD);
max1 = max(mag1);
max1_n = mag1 / max1;
max1_n_dB = mag2db(max1_n);
figure()  % Create a new figure for the normalized mismatched filter output in dB
plot(max1_n_dB, 'r')
ylim([-60 0])
grid on
legend('Matched filter', 'Mismatched filter')
xlabel('time(in samples)')
ylabel('filter output')
