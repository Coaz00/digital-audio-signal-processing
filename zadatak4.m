%Aleksandar Djordjevic 2019/0086

%% Ucitavanje signala
clear all;
close all;
clc;

[x, Fs] = audioread('violina_3.wav');

%% Vremenski domen i amplitudska karakteristika
t = 0:1/Fs:(length(x)-1)/Fs;

N = 2^nextpow2(length(x));
f1 = 0:Fs/N:Fs/2;
X = fft(x,N)/length(x);
amp_x = abs(X(1:N/2+1));
amp_x(2:N/2+1) = 2*amp_x(2:N/2+1);

figure(1)

subplot(2,1,1)
plot(t,x);
title('Signal u vremenskom domenu');
xlabel('t[s]');
ylabel('x(t)');

subplot(2,1,2)
plot(f1,amp_x);
title('Amplitudska frekvencijska karakteristika signala');
xlabel('f[Hz]');

%% Izdvajanje prvog pika

% FIR filtar sa hamming-ovim prozorom

n = 32;
window = hamming(n+1);
Wn =[1000 1100]/(Fs/2);
b = fir1(n,Wn,window);
a = 1;
[h, f] = freqz(b,a,N/2+1,Fs);

figure(2)

subplot(2,1,1)
stem(window);
xlabel('n[odb]');
title('Hamming-ova prozorska funkcija reda');

subplot(2,1,2)
plot(f, 20*log10(abs(h)));
xlabel('f[Hz]');
title('Amplitudska karakteristika filtra')

% Filtriranje signala

x_fil = filter(b,a,x);
X_fil = fft(x_fil,N)/length(x_fil);
amp_x_fil = abs(X_fil(1:N/2+1));
amp_x_fil(2:N/2+1) = 2*amp_x_fil(2:N/2+1);

figure(3)

plot(f1,amp_x);
hold on;
plot(f1,amp_x_fil);
title('Amplitudska frekvencijska karakteristika signala');
xlabel('f[Hz]');
legend('Pre filtriranja','Posle filtriranja');

audiowrite('violina_3_fil.wav',x_fil,Fs);

%% Zadatak4 pod b)

m = 8;

x_fil_zam = zeros(1,round(length(x_fil)/m));
count = 1;

for i =1:length(x_fil)
   if(mod(i,m) == 0)
       x_fil_zam(count) = x_fil(i);
       count = count + 1;
   end
end


X_fil_zam = fft(x_fil_zam,N)/length(x_fil_zam);
amp_x_fil_zam = abs(X_fil_zam(1:N/2+1));
amp_x_fil_zam(2:N/2+1) = 2*amp_x_fil_zam(2:N/2+1);

f2 = linspace(0,Fs/16,N/2+1);

figure(4)
plot(f1,amp_x_fil);
hold on;
plot(f2,amp_x_fil_zam);
hold off;
title('Amplitudska frekvencijska karakteristika signala');
xlabel('f[Hz]');
legend('Adekvatna ucestanost odabiranja', 'Neadekvatna ucestanost odabiranja');

audiowrite('violina_3_zam.wav',x_fil_zam,Fs);