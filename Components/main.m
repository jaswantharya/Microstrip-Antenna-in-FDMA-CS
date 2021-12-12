%Clear all history
clear,clc,close all
% Parameters
CF1 = 2.43e9; % First carrier
CF2 = 2.63e9; % Second carrier

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1) Transmitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose input message
%first signal
[FileName,PathName] = uigetfile('*.wav*');
[sig1,Fsa1] = audioread([PathName FileName]);
Fs1=Fsa1*1e4;
S1 = size(sig1);
%Make sure signal is mono not (stereo or more)
if S1(2)>1
    sig1 = sum(sig1,2);
end

%second signal
[FileName,PathName] = uigetfile('*.wav*');
[sig2,Fsa2] = audioread([PathName FileName]);
Fs2=Fsa2*1e4;
S2 = size(sig2);
%Make sure signal is mono not (stereo or more)
if S2(2)>1
    sig2 = sum(sig2,2); 
end

%% Plotting frequency spectrum
%preforming fourier transform
df = Fs1/S1(1); 
ff1 = (-Fs1/2:df:Fs1/2-df)';
df = Fs2/S2(1); 
ff2 = (-Fs2/2:df:Fs2/2-df)';
FFT1=abs(fftshift(fft(sig1)));
FFT2=abs(fftshift(fft(sig2)));

%Plotting
figure('Position',[85 85 1151 588])
subplot(2,2,1),plot(ff1,FFT1); title('Spectrum of Sig1');
xlabel('Frequency'); ylabel('Amplitude');
subplot(2,2,2),plot(ff2,FFT2); title('Spectrum of Sig2');
xlabel('Frequency'); ylabel('Amplitude');

%Frequency change
F1=CF1*(1e9)/CF1;
F2=CF2*(1.6e9)/CF2;
mar = 0.3*abs(F1-F2); % Margin for filters

%% Upsampling signals
% Both signals are upsampled to make nyquist condition met 
r1 = fix(5*max(F1,F2)/(Fs1));
sig1 = interp(sig1,r1);
Fs1 = r1*Fs1; 
S1 = size(sig1);

r2 = fix(5*max(F1,F2)/(Fs2));
sig2 = interp(sig2,r2); 
Fs2 = r2*Fs2;
S2 = size(sig2);

%% Make the two signals of same length
if S1(1)>S2(1)
    sig2=[sig2;zeros(S1(1)-S2(1),1)];
    S = S1;
elseif S2(1)>S1(1)
    sig1=[sig1;zeros(S2(1)-S1(1),1)];
    S = S2;
else
    S = S2; 
end
%clear S1 S2
%% Creating two carrier signals
n = (0:S(1)-1)';
C1 = 2*cos(2*pi*F1*n*(1/Fs1));
C2 = 2*cos(2*pi*F2*n*(1/Fs2));

%Modulation of carrier signals
C1m = C1.*sig1;
C2m = C2.*sig2;

%preforming fourier transform
df = Fs1/S(1);
ff1 = (-Fs1/2:df:Fs1/2-df)';

df = Fs2/S(1);
ff2 = (-Fs2/2:df:Fs2/2-df)';

FFT1=abs(fftshift(fft(C1m)));
FFT2=abs(fftshift(fft(C2m)));

%Plotting
subplot(2,2,3),plot(ff1,FFT1); 
title('Spectrum of modulated Sig1');
xlabel('Frequency'); 
ylabel('Amplitude');

subplot(2,2,4),plot(ff2,FFT2); 
title('Spectrum of modulated Sig2');
xlabel('Frequency'); 
ylabel('Amplitude');

%% Adding two signals in one message
MSG = C1m + C2m; % Sending through air
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2) Receiver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spectrum of received messsage
%preforming fourier transform
df = Fs2/S(1); ff = (-Fs2/2:df:Fs2/2-df)';
FFTmsg=abs(fftshift(fft(MSG)));
%Plotting
figure(),plot(ff,FFTmsg); 
title('Spectrum of received Signal & filters');
xlabel('Frequency'); ylabel('Amplitude');
hold on; axis manual

%% Band pass filters
% 1st filter
FPra=[min(ff),F1-mar,F1-0.5*mar,F1+0.5*mar,F1+mar,max(ff)]; %Freq parameters
APra=max(FFTmsg)*[0,0,1,1,0,0]; %Amplitude prameters
plot([-flip(FPra),FPra],[APra,APra],'r--') %Plotting filter response
FPra=(2/Fs1)*FPra; % Normalize frequency parameters
d=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',FPra(2),FPra(3),FPra(4),FPra(5),60,1,60);
BFilt1 = design(d,'equiripple');

% 2nd filter
FPra=[min(ff),F2-mar,F2-0.5*mar,F2+0.5*mar,F2+mar,max(ff)]; %Freq parameters
APra=max(FFTmsg)*[0,0,1,1,0,0]; %Amplitude prameters
plot([-flip(FPra),FPra],[APra,APra],'g--') %Plotting filter response
FPra=(2/Fs2)*FPra; % Normalize frequency parameters
d=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',FPra(2),FPra(3),FPra(4),FPra(5),60,1,60);
BFilt2 = design(d,'equiripple'); hold off

% Apply filter
Ch1 = filter(BFilt1,MSG); %Channel 1
Ch2 = filter(BFilt2,MSG); %Channel 2

%preforming fourier transform
df = Fs1/S(1); ff1 = (-Fs1/2:df:Fs1/2-df)';
df = Fs2/S(1); ff2 = (-Fs2/2:df:Fs2/2-df)';
FFT1=abs(fftshift(fft(Ch1)));
FFT2=abs(fftshift(fft(Ch2)));

%Plotting
figure('Position',[85 85 1151 588])
subplot(3,2,1),plot(ff1,FFT1); title('Spectrum of Channel 1');
xlabel('Frequency'); ylabel('Amplitude');
subplot(3,2,2),plot(ff2,FFT2); title('Spectrum of Channel 2');
xlabel('Frequency'); ylabel('Amplitude');

%% Demodulation
Ch1 = Ch1.*C1; Ch2 = Ch2.*C2;
%preforming fourier transform
df = Fs1/S(1); ff1 = (-Fs1/2:df:Fs1/2-df)';
df = Fs2/S(1); ff2 = (-Fs2/2:df:Fs2/2-df)';
FFT1=abs(fftshift(fft(Ch1)));
FFT2=abs(fftshift(fft(Ch2)));
%Plotting
subplot(3,2,3),plot(ff1,FFT1); title('Demodulated Channel 1 & LBF');
xlabel('Frequency'); ylabel('Amplitude');
subplot(3,2,4),plot(ff2,FFT2); title('Demodulated Channel 2 & LBF');
xlabel('Frequency'); ylabel('Amplitude');

%% Low pass filter
F = min(F1,F2);
FPra=[0,0.8*F,F]; %Freq parameters
APra=max(FFT1)*[1,1,0]; %Amplitude prameters
subplot(3,2,3),hold on
plot([-flip(FPra),FPra],[flip(APra),APra],'r--') %Plotting filter response
APra=max(FFT2)*[1,1,0]; %Amplitude prameters
subplot(3,2,4),hold on
plot([-flip(FPra),FPra],[flip(APra),APra],'r--') %Plotting filter response

d=fdesign.lowpass('Fp,Fst,Ap,Ast',(2/Fs1)*FPra(2),(2/Fs1)*FPra(3),1,60);
LFilt = design(d,'equiripple'); hold off
% Apply filter
Ch1 = filter(LFilt,Ch1); %Channel 1

d=fdesign.lowpass('Fp,Fst,Ap,Ast',(2/Fs2)*FPra(2),(2/Fs2)*FPra(3),1,60);
LFilt = design(d,'equiripple'); hold off
% Apply filter
Ch2 = filter(LFilt,Ch2); %Channel 2

%preforming fourier transform
FFT1=abs(fftshift(fft(Ch1)));
FFT2=abs(fftshift(fft(Ch2)));
%Plotting
subplot(3,2,5),plot(ff1,FFT1); title('Filtered Channel 1');
xlabel('Frequency'); ylabel('Amplitude');
subplot(3,2,6),plot(ff2,FFT2); title('Filtered Channel 2');
xlabel('Frequency'); ylabel('Amplitude');

%% Downsampling
Ch1 = downsample(Ch1,r1); Fs1=Fs1/r1;
Ch2 = downsample(Ch2,r2); Fs2=Fs2/r2;
f1 = 'output1.wav';
audiowrite(f1,Ch1,Fs1/1e4);
f2 = 'output2.wav';
audiowrite(f2,Ch2,Fs2/1e4);
%% Play messages
sound(Ch1,Fs1/1e4);
pause(5)
sound(Ch2,Fs2/1e4);
fvtool(BFilt1,BFilt2,LFilt); %preformance of filters