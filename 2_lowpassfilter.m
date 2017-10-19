clear;clc;

N   = 200;        % FIR filter order
Fp  = 5.5e3;      % 5.5 kHz passband-edge frequency
Fs  = 16e3;       % 16 kHz sampling frequency
Rp  = 0.00057565; % Corresponds to 0.01 dB peak-to-peak ripple
Rst = 1e-4;       % Corresponds to 80 dB stopband attenuation
    
eqnum = firceqrip(N,Fp/(Fs/2),[Rp Rst],'passedge'); % eqnum = vec of coeffs
%fvtool(eqnum,'Fs',Fs,'Color','White') % Visualize filter

x=audioread('mythisisatest.wav');
y=filter(eqnum,1,x);

sound(y,Fs);
% sound(x,Fs)
