clear;clc;

%% create a mat file containing required information
cd('C:\Users\...); % go to the desired diretory
speechfiles  = dir();  %read the filescontains the list of ASCII speech files 
N= length(speechfiles)-2; % the number of speech files;

files=zeros(1,N);
files=speechfiles(3:N+2,1);

%% lmax= 26500;  % maximum number of samples to be played for any of the files

%% create a beep signal
Fs = 16000;                         % Default Sampling Frequency (Hz)
T = 0:1/Fs:((Fs-1)/Fs);             % The interval time is 1 second
Ftone = 1000;                       % Tone Frequency
Y = 0.01*sin(2*pi*Ftone*T);          
Y=Y';                               % Create the beep signal

%% play out the sequence of files
audio=[];
for i = 1:length(files)
    cd('C:\Users\...');
    [speech,fs] = audioread(files(i).name);    % read the speech signal
    lmax(i) = length(speech);                  % maximum number of samples to be played for any of the files
    audio=[audio; speech; Y];                  % create the speech file combined with seperate speech and beep signal 
end

sound (audio,16000);          % sound the combined speech signal
