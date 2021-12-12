clear,clc,close all
Fs = 8000;
% Record your voice for 5 seconds.
recObj = audiorecorder(Fs,8,1);
disp('Get ready ...')
pause(1.8)
disp('Start speaking.')
recordblocking(recObj, 5);
disp('End of Recording.');
% Play back the recording.
play(recObj);
% Store data in double-precision array.
myRecording = getaudiodata(recObj);
% Plot the waveform.
plot(myRecording);
FileName = input('Enter saving name for the message : ','s');
audiowrite([FileName '.wav'],myRecording,Fs);