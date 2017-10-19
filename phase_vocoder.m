close all; clear;clc;
[y,fs] = audioread('si552.wav');
output = pvoc(y, 0.75, 1024);
audiowrite('prolonged speech.wav',output,fs)
sound(output,fs)

function x = istft(d, ftsize, w, h)     %inverse short-time Fourier transform
%   d: interpolated array at new rate
%   ftsize: size of FFT
%   w: frame size of STFT
%   h: window offset
%   x: output rate

hop = round(h);
s = size(d);
cols = s(2);
xlen = ftsize + cols * (hop);
x = zeros(1,xlen);

% force window be odd
    if rem(w, 2) == 0       
        w = w + 1;
    end
    
% create hann window
    win = zeros(1, ftsize);
    halff = ftsize/2;   % midpoint of win
    halflen = (w-1)/2;
    acthalflen = min(halff, halflen);
    
% create symmetric hann window
    halfwin = 0.5 * ( 1 + cos( pi * (0:halflen)/halflen));
    win((halff+1):(halff+acthalflen)) = halfwin(1:acthalflen);
    win((halff+1):-1:(halff-acthalflen+2)) = halfwin(1:acthalflen);
    
% calculate new rate from inverse STFT with hann windows
    for b = 0:hop:(hop*(cols-1))
        ft = d(:,1+b/hop)';
        ft = [ft, conj(ft([((ftsize/2)):-1:2]))];
        px = real(ifft(ft));
        x((b+1):(b+ftsize)) = x((b+1):(b+ftsize))+px.*win;
    end;

end 


function d = stft(x, f, w, h)       %short-time Fourier transform
%   x: input
%   f: FFT size in samples (usually same as window size)
%   w: window size in samples
%   h: window shift in samples
%   d: matrix of short-time Fourier transforms

hop=round(h);
s = length(x);

% force window be odd
    if rem(w, 2) == 0   
        w = w + 1;
    end
    
% define hann window as symmetric half cosine
    halflen = (w-1)/2;
    halff = f/2;   % midpoint of win
    acthalflen = min(halff, halflen);
    halfwin = 0.5 * ( 1 + cos( pi * (0:halflen)/halflen));
    win = zeros(1, f);
    win((halff+1):(halff+acthalflen)) = halfwin(1:acthalflen);
    win((halff+1):-1:(halff-acthalflen+2)) = halfwin(1:acthalflen);
    
c = 1;
d = zeros((1+f/2),1+fix((s-f)/hop));

% perform STFT using FFTs; store results in d(1+f/2,fix((sif)/h)
    for b = 0:hop:(s-f)
        u = win.*x((b+1):(b+f));
        t = fft(u);
        d(:,c) = t(1:(1+f/2))';
        c = c+1;
    end;

end 


function c = pvsample(b, t, hop)        %interpolate an STFT array
%   b: matrix of STFT of input signal (rate=1)
%   t: new set of sampling times according to new rate (rate=r)
%   hop: window shift of original STFT; defaults to half the window length
%   c: interpolated STFT array based on rate r

if nargin < 3
    hop = 0;
end

[rows,cols] = size(b);
N = 2*(rows-1);

if hop == 0
    hop = N/2;
end

c = zeros(rows, length(t));
ph = angle(b(:,1));
b = [b,zeros(rows,1)];
ocol = 1;
for tt = t

% Interpolate magnitude from the two columns
    bcols = b(:,floor(tt)+[1 2]);
    tf = tt - floor(tt);
    bmag = (1-tf)*abs(bcols(:,1)) + tf*(abs(bcols(:,2)));    

    if (ocol <= 5)
        fprintf('tt:%6.2f, ocol:%d, cols: %d %d, tf:%8.2f \n',tt,ocol,floor(tt)+[1,2],tf);
    end
    
%phase advance
    dp = angle(bcols(:,2)) - angle(bcols(:,1));
    
    c(:,ocol) = bmag .* exp(j*ph);
    ph = ph + dp;
    ocol = ocol+1;
end

end 


function y = pvoc(x, r, n)      %time-scale a signal to r times faster
%   x: input waveform
%   r: rate of speedup (1<=r<=4) or slowdown (0.25<=r<=1)
%   n: size of FFT for transforms (default of 1024)
%   y: speeded up or slowed down waveform

if nargin < 3
    n = 1024;
end

hop = n/4;
scf = 2/3;
X = scf * stft(x', n, n, hop);

[rows, cols] = size(X);
t = 0:r:(cols-2);

X2 = pvsample(X, t, hop);

y = istft(X2, n, n, hop)';
 
end 


function y = pvoc(x, r, n)      %time-scale a signal to r times faster
%   x: input waveform
%   r: rate of speedup (1<=r<=4) or slowdown (0.25<=r<=1)
%   n: size of FFT for transforms (default of 1024)
%   y: speeded up or slowed down waveform

if nargin < 3
    n = 1024;
end

hop = n/4;
scf = 2/3;
X = scf * stft(x', n, n, hop);

[rows, cols] = size(X);
t = 0:r:(cols-2);

X2 = pvsample(X, t, hop);

y = istft(X2, n, n, hop)';
 
end 



function c = pvsample(b, t, hop)        %interpolate an STFT array
%   b: matrix of STFT of input signal (rate=1)
%   t: new set of sampling times according to new rate (rate=r)
%   hop: window shift of original STFT; defaults to half the window length
%   c: interpolated STFT array based on rate r

if nargin < 3
    hop = 0;
end

[rows,cols] = size(b);
N = 2*(rows-1);

if hop == 0
    hop = N/2;
end

c = zeros(rows, length(t));
ph = angle(b(:,1));
b = [b,zeros(rows,1)];
ocol = 1;
for tt = t

% Interpolate magnitude from the two columns
    bcols = b(:,floor(tt)+[1 2]);
    tf = tt - floor(tt);
    bmag = (1-tf)*abs(bcols(:,1)) + tf*(abs(bcols(:,2)));    

    if (ocol <= 5)
        fprintf('tt:%6.2f, ocol:%d, cols: %d %d, tf:%8.2f \n',tt,ocol,floor(tt)+[1,2],tf);
    end
    
%phase advance
    dp = angle(bcols(:,2)) - angle(bcols(:,1));
    
    c(:,ocol) = bmag .* exp(j*ph);
    ph = ph + dp;
    ocol = ocol+1;
end

end 


function d = stft(x, f, w, h)       %short-time Fourier transform
%   x: input
%   f: FFT size in samples (usually same as window size)
%   w: window size in samples
%   h: window shift in samples
%   d: matrix of short-time Fourier transforms

hop=round(h);
s = length(x);

% force window be odd
    if rem(w, 2) == 0   
        w = w + 1;
    end
    
% define hann window as symmetric half cosine
    halflen = (w-1)/2;
    halff = f/2;   % midpoint of win
    acthalflen = min(halff, halflen);
    halfwin = 0.5 * ( 1 + cos( pi * (0:halflen)/halflen));
    win = zeros(1, f);
    win((halff+1):(halff+acthalflen)) = halfwin(1:acthalflen);
    win((halff+1):-1:(halff-acthalflen+2)) = halfwin(1:acthalflen);
    
c = 1;
d = zeros((1+f/2),1+fix((s-f)/hop));

% perform STFT using FFTs; store results in d(1+f/2,fix((sif)/h)
    for b = 0:hop:(s-f)
        u = win.*x((b+1):(b+f));
        t = fft(u);
        d(:,c) = t(1:(1+f/2))';
        c = c+1;
    end;

end 

