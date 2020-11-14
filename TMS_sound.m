%% Auditory Mask for TMS

%% load click
% data taken from:
% https://data.mendeley.com/datasets/jb93vbfkwr/1
% sound structure (renaiming to sounds to not interfer with soudn function)
% contains recorded clicks from different coils and intensities (10 pulses 
% per intensity), at 2 different distances
load('C:\Users\ckohl\Desktop\Current\Other\Auditory masking\sound_mat\sound.mat')
sounds=sound;
clear sound

coil_oi='Magstim 70mm Double Coil';
coil_oi_i=2;
Fs=192000;

% %plot sound wave at different intensities
% figure
% subplot(10,1,1)
% hold on
% for pulse=1:length(sounds.stimulator_output(coil_oi_i,:))
%     if pulse/10==round(pulse/10) & pulse ~=length(sounds.stimulator_output(coil_oi_i,:))
%         subplot(10,1,pulse/10+1)
%         hold on
%     end
%     plot(sounds.time,sounds.microphone_025{coil_oi_i,1}(:,pulse))
%     title(num2str(sounds.stimulator_output(coil_oi_i,pulse)))
%     xlim([sounds.time(1),sounds.time(end)])
%     ylim([-4 4])
% end

% % play click sequence (10 pulses each at 10 intensities)
% % one coil and one distance selected
% for pulse=1:length(sounds.stimulator_output(coil_oi_i,:))
%     sound(sounds.microphone_025{coil_oi_i,1}(:,pulse),Fs)
%     pause(.5)
% end

click=sounds.microphone_025{coil_oi_i,1}(:,length(sounds.stimulator_output(coil_oi_i,:)));
sound(click,Fs)


%% Phase Randomise
% match click spectrum but eliminate temporal information 
% (Jason Ritt gave us this code) 
% to use his example
% load handel -> use y Fs that comes with that

snd_orig = click;
snd_fft = fft(snd_orig,2^14);%fourier
N = length(snd_fft); % Should be power of two
rand_phase = exp((sqrt(-1)*2*pi)*rand(N/2-1,1));
rand_phase = [1; rand_phase; 1; conj(rand_phase(end:-1:1))]; % Conjugate symmetry

snd_rand = real(ifft(rand_phase.*snd_fft)); % Back to time domain
% The "real" is to eliminate any spurious imaginary components due to round off error.

sound(snd_orig,Fs)% original sound
pause(5)
sound(snd_rand,Fs)  % Same “sound” but no temporal info

%% attempt to make it longer:
% looping doesn't work
% simply repeating doesn't work either
% long_snd_rand=repmat(snd_rand,100,1);% 
% sound(long_snd_rand,Fs) 
% so here I just sample from the snd-random lots 
seconds=10;
long_snd_rand=datasample(snd_rand,Fs*seconds);
sound(long_snd_rand,Fs)

% this result doesn't contain white noise, but why would it?
% this does not sound the same
% If I sample only as many samples as there are in snd_rand (so it's the same
% length), those don't sound the same
sound(snd_rand,Fs)
long_snd_rand=datasample(snd_rand, 16384)
sound(long_snd_rand,Fs)


% another attempt to make it longer
long_snd_rand = real(ifft(repmat(rand_phase,5,1).*repmat(snd_fft,5,1))); 
sound(snd_rand,Fs)  % Same “sound” but no temporal info
sound(long_snd_rand,Fs)

% also no
snd_fft = fft(snd_orig,2^14);%fourier
N = length(snd_fft)*10; % Should be power of two
rand_phase = exp((sqrt(-1)*2*pi)*rand(N/2-1,1));
rand_phase = [1; rand_phase; 1; conj(rand_phase(end:-1:1))]; 
long_snd_rand = real(ifft(rand_phase.*repmat(snd_fft,10,1))); 





%plot the new spectrum
figure
subplot(1,2,1)
hold on
power=[];
x=long_snd_rand;
fs=Fs;
y = fft(x);
n = length(x);          % number of samples
f = (0:n-1)*(fs/n);     % frequency range
power = abs(y).^2/n;    % power of the DFT
plot(f,pow2db(power))
xlabel('Freq')
ylabel('DB')
title('long_snd_rand')
xlim([100 50000])
subplot(1,2,2)
hold on
power=[];
x=snd_rand;
fs=Fs;
y = fft(x);
n = length(x);          % number of samples
f = (0:n-1)*(fs/n);     % frequency range
power = abs(y).^2/n;    % power of the DFT
plot(f,pow2db(power))
xlabel('Freq')
ylabel('DB')
title('snd_rand')
xlim([100 50000])


%% SO the question is still: how do we make this longer
















%% Code Graveyard
% 
% 
% 
% figure
% hold on
% for pulse=1:length(sounds.stimulator_output(coil_oi_i,:))  
%     power=[];
%     x=sounds.microphone_025{coil_oi_i,1}(:,pulse);
%     fs=Fs;
%     y = fft(x);
%     n = length(x);          % number of samples
%     f = (0:n-1)*(fs/n);     % frequency range
%     power = abs(y).^2/n;    % power of the DFT
% 
%     plot(f,pow2db(power))
%     xlim([100 50000])
% end
% 
% 
%     x=sounds.microphone_025{coil_oi_i,1}(:,pulse);
%     fs=Fs;
%     y = fft(x,1000);
%     n = length(x);  
%     n=30000;% number of samples
%     f = (0:n-1)*(fs/n);     % frequency range
%     power = abs(y).^2/n;    % power of the DFT
% 
%     plot(f,pow2db(power))
%     
%     
%     
% sound(ifft(y),Fs)
% %plot 1/3 octave power spectrum
% %https://www.mathworks.com/help/audio/ug/octave-band-and-fractional-octave-band-filters.html
% https://www.researchgate.net/post/How_do_I_generate_time_series_data_from_given_PSD_of_random_vibration_input
% 
% x=sounds.microphone_025{coil_oi_i,1}(:,pulse);
% y = fft(x,'symflag');
% amp=sqrt(2*y);
% %assign random phase
% phase=0+(360-0).*rand(size(y));
% freqdomain=amp.*exp(i*phase);
% 
% new_y=ifft(freqdomain)
% sound(new_y)
% 
% %make white noise
% % noise=wgn
% 
% figure
% power=[];
% data=sounds.microphone_025{coil_oi_i,1}(:,pulse);
% x=data;
% fs=Fs;
% y = fft(x);
% n = length(x);          % number of samples
% f = (0:n-1)*(fs/n);     % frequency range
% data_power = abs(y).^2/n;    % power of the DFT
% subplot(1,3,1)
% plot(f,pow2db(data_power))
% xlim([0 20000])
% ylim([-50 10])
%     
% mu=0;sigma=1;    
% noise= sigma *randn(1,length(sounds.microphone_025{coil_oi_i,1}(:,pulse)))+mu;
% x=noise;
% fs=Fs;
% y = fft(x);
% n = length(x);          % number of samples
% f = (0:n-1)*(fs/n);     % frequency range
% noise_power = abs(y).^2/n;    % power of the DFT
% subplot(1,3,2)
% plot(f,pow2db(noise_power))
% xlim([0 20000])
% ylim([-50 10])
% 
% comb=noise'.*data;
% x=comb;
% fs=Fs;
% y = fft(x);
% n = length(x);          % number of samples
% f = (0:n-1)*(fs/n);     % frequency range
% comb_power = abs(y).^2/n;    % power of the DFT
% subplot(1,3,3)
% plot(f,pow2db(comb_power))
% xlim([0 20000])
% ylim([-50 10])
% 
% 
% sound(data,Fs)
% sound(noise,Fs)
% sound(comb,Fs)
% 
% 
% for i=1:10
%   soundsc(comb,Fs)
%   pause(.2)
% end
% 
% 
% 
% 
% 
% Fs = 192000;                                     % Sampling Frequency
% secs = 1;
% t  = linspace(0, secs, Fs*secs+1);              % Time Vector + 1 sample
% t(end) = [];                                    % remove extra sample
% for freq=1:length(f)
%     % freq = 1000;
%     w = 2*pi*f(freq);                                  % Radian Value To Create 1kHz Tone
%     s = 100*data_power(freq)*sin(w*t);                  % Create Tone
% %     sound(s, Fs)                                    % Produce Tone As Sound
%     if sum(s)>0
%         
%         if freq==1
%             combined=s;
%         else
%             combined=combined.*s;
%         end
%     end
% end
% soundsc(combined, Fs)    
% 
% 
% 
% 
% 
% 
% 
% 
% 
% amp=10 
% fs=14400  % sampling frequency
% duration=1
% freq=[100]
% values=0:1/fs:duration;
% a=amp*sin(2*pi* freq*values)
% sound(a)
% 
% 
% amp = 0.5;
% Fs = 50 + 100 + 200 + 400 + 1000 + 2000 + 4000 + 6000 + 8000 + 10000 + 12000 + 14000 + 16000 + 18000;                         
% Ts = 1/Fs;
% T = 0:Ts:3;
% x = [];
% for freq = [50 100 200 400 1000 2000 4000 6000 8000 10000 12000 14000 16000 18000]
%     x = [ x amp * sin(2*pi*freq*T)];
% end
% sound(x,Fs)
% 
% 
% 
% 
% 
% Fs = 14400;                                     % Sampling Frequency
% secs = 1;
% t  = linspace(0, secs, Fs*secs+1);              % Time Vector + 1 sample
% t(end) = [];                                    % remove extra sample
% freq = 1000;
% w = 2*pi*freq;                                  % Radian Value To Create 1kHz Tone
% s = sin(w*t);                                   % Create Tone
% sound(s, Fs)                                    % Produce Tone As Sound
% 
% 
% sound(y,Fs)
% 
