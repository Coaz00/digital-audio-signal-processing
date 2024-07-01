%Aleksandar Djordjevic 2019/0086

clear all;
close all
clc

%% Ucitavanje signala

[x, Fs] = audioread('star_wars_zasumljen1.wav');

t = 0:1/Fs:(length(x) - 1)/Fs;

%% Vremenski domen i amplitudska karakteristika

figure(1)
plot(t,x);
xlabel('t[s]');
ylabel('x(t)');
title('Zasumljeni signal');

N = 2^nextpow2(length(x));
f1 = 0:Fs/N:Fs/2; 
X = fft(x,N)/length(x);
amp_x = abs(X(1:N/2+1));
amp_x(2:N/2+1) = 2*amp_x(2:N/2+1);

figure(2)
plot(f1,amp_x);
xlim([0 10^4]);
xlabel('f[Hz]');
ylabel('|X(j2 \pi f)|');
title('Amplitudska frekvencijska karakteristika zasumljenog signala');

%% Uklanjanje frekvencij na 4000Hz

% Digitalni elipticki filtar(bandstop)

Wp = [3600 4400]/(Fs/2);
Ws = [3500 4500]/(Fs/2);
Rp = 2;
Rs = 40;

[n, Wn] = ellipord(Wp,Ws,Rp,Rs);
[b, a] = ellip(n,Rp,Rs,Wn, 'stop');
[h, f]= freqz(b, a, N/2 + 1, Fs);

figure(3)

plot(f, 20*log10(abs(h)));
xlabel('f[Hz]');
title('Amplitudska frekvencijska karakteristika bandstop filtra');
ylim([-100 10]);

% Filtriranje signala

x_bez1 = filter(b,a,x);
X_bez1 = fft(x_bez1,N)/length(x_bez1);
amp_x_bez1 = abs(X_bez1(1:N/2+1));
amp_x_bez1(2:N/2+1) = 2*amp_x_bez1(2:N/2+1);

figure(4)

subplot(2,2,[1 3]);
plot(f1,amp_x);
hold on;
plot(f1,amp_x_bez1);
hold off;
xlabel('f[Hz]');
title('AFK posle bandstop filtra');
legend('Pre filtriranja', 'Posle filtriranja');
xlim([0 10^4])

subplot(2,2,2);
plot(t,x);
xlabel('t[s]');
ylabel('x(t)');
title('Zasumljeni signal');

subplot(2,2,4)
plot(t,x_bez1);
xlabel('t[s]');
ylabel('x_ bez1(t)');
title('Signal sa sumom na 6500Hz');

%% Uklanjanje suma na frekveniciji 6500Hz

% Digitalni elipticki filtar(lowpass)

Wp = 4500/(Fs/2);
Ws = 4600/(Fs/2);
Rp = 2;
Rs = 40;

[n, Wn] = ellipord(Wp,Ws,Rp,Rs);
[b2, a2] = ellip(n, Rp, Rs, Wp);
[h, f] = freqz(b2,a2,N/2+1,Fs);

figure(5)
plot(f, 20*log10(abs(h)));
xlabel('f[Hz]');
title('Amplitudska frekvencijska karakteristika lowpass filtra')

% Filtriranje signala

x_bez2 = filter(b2,a2,x);
X_bez2 = fft(x_bez2,N)/length(x_bez2);
amp_x_bez2 = abs(X_bez2(1:N/2+1));
amp_x_bez2(2:N/2+1) = 2*amp_x_bez2(2:N/2+1);

figure(6)

subplot(2,2,[1 3]);
plot(f1,amp_x);
hold on;
plot(f1,amp_x_bez2);
hold off;
xlabel('f[Hz]');
title('AFK posle lowpass filtra');
legend('Pre filtriranja', 'Posle filtriranja');
xlim([0 10^4])

subplot(2,2,2);
plot(t,x);
xlabel('t[s]');
ylabel('x(t)');
title('Zasumljeni signal');

subplot(2,2,4)
plot(t,x_bez2);
xlabel('t[s]');
ylabel('x_ bez1(t)');
title('Signal sa sumom na 4000Hz');

%% Potpuno uklanjanje suma

x_bez_suma = filter(b,a,x);
x_bez_suma = filter(b2,a2,x_bez_suma);

X_bez_suma = fft(x_bez_suma,N)/length(x_bez_suma);
amp_x_bez_suma = abs(X_bez_suma(1:N/2+1));
amp_x_bez_suma(2:N/2+1) = amp_x_bez_suma(2:N/2+1);

figure(7)

subplot(2,1,2);
plot(f1,amp_x);
hold on;
plot(f1,amp_x_bez_suma);
hold off;
xlabel('f[Hz]');
title('Amplitudska frekvencijska karakteristika signala');
legend('Pre filtriranja', 'Posle filtriranja');
xlim([0 10^4])

subplot(2,1,1);
plot(t,x);
hold on;
plot(t,x_bez_suma);
hold off;
title('Signal');
xlabel('t[s]');
ylabel('x(t)');
legend('Pre filtriranja', 'Posle filtriranja');

audiowrite('isfiltriran1.wav', x_bez_suma,Fs);

