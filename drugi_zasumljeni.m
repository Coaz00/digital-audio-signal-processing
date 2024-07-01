%Aleksandar Djordjevic 2019/0086

%% Ucitavanje signala, vremenski domen i amp karakteristika
clear all;
close all;
clc

[x, Fs] = audioread('star_wars_zasumljen2.wav');


t=0:1/Fs:(length(x) - 1)/Fs;

N = 2^nextpow2(length(x));
f1 = 0:Fs/N:Fs/2;
X = fft(x,N)/length(x);
amp_x = abs(X(1:N/2+1));
amp_x(2:N/2+1) = 2*amp_x(2:N/2+1);

figure(1)

subplot(2,1,1)
plot(t,x)
xlabel('t[s]');
ylabel('x(t)');
title('Zasumljeni signal');

subplot(2,1,2);
plot(f1,amp_x);
xlim([0 0.25*10^4])
xlabel('f[Hz]');
title('Amplitudska frekvencijska karakteristika zasumljenog signala');

%% Uklanjanje suma

% Analogni elipticni filtar(lowpass)
Wp = 2000*2*pi;
Ws = 2500*2*pi;
Rp = 2;
Rs = 40;

[n,Wn] = ellipord(Wp,Ws,Rp,Rs,'s');
[b, a] = ellip(n,Rp,Rs,Wp,'s');

[h, w] = freqs(b,a,N/2+1);

[bz,az] = bilinear(b,a,Fs);

[hz, fz] = freqz(bz,az,N/2+1,Fs);

figure(2)
plot(w/(2*pi),20*log10(abs(h)));
hold on;
plot(fz,20*log10(abs(hz)));
hold off;
xlabel('f[Hz]');
title('Amplitudska frekvencijska karakteristika filtra');
legend('Analogni','Digitalni');



% Filtriranje signala
x_bez_suma = filter(bz,az,x);
X_bez_suma = fft(x_bez_suma,N)/length(x_bez_suma);
amp_x_bez_suma = abs(X_bez_suma(1:N/2+1));
amp_x_bez_suma(2:N/2+1) = 2*amp_x_bez_suma(2:N/2+1);

figure(3);
plot(f1,amp_x);
hold on;
plot(f1,amp_x_bez_suma);
hold off;
title('AFK signala');
xlabel('f[Hz]');
legend('Pre filtriranja','Posle filtriranja')


figure(4)
plot(t,x);
hold on;
plot(t,x_bez_suma);
hold off;
title('Signal');
xlabel('t[s]');
ylabel('x(t)');
legend('Pre filtriranja','Posle filtriranja');

audiowrite('isfiltriran2.wav',x_bez_suma,Fs);
