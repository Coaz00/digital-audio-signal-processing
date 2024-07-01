%Aleksandar Djordjevic 2019/0086

%% Ucitavanje signala

clear all
close all
clc

[x_flauta, Fs_flauta] = audioread('flauta_3.wav');
[x_klavir, Fs_klavir] = audioread('klavir_3.wav');
[x_truba, Fs_truba] = audioread('truba_3.wav');
[x_violina, Fs_violina] = audioread('violina_3.wav');

t_flauta = 0:1/Fs_flauta:(length(x_flauta) - 1)/Fs_flauta;
t_klavir = 0:1/Fs_klavir:(length(x_klavir) - 1)/Fs_klavir;
t_truba = 0:1/Fs_truba:(length(x_truba) - 1)/Fs_truba;
t_violina = 0:1/Fs_violina:(length(x_violina) - 1)/Fs_violina;

%% Zadatak 2 pod a)

figure(1)

subplot(2,2,1)
plot(t_flauta,x_flauta);
xlabel('t[s]');
ylabel('x_{flauta}(t)');
title('Flauta');

subplot(2,2,2)
plot(t_klavir,x_klavir);
xlabel('t[s]');
ylabel('x_{klavir}(t)');
title('Klavir');

subplot(2,2,3)
plot(t_truba,x_truba);
xlabel('t[s]');
ylabel('x_{truba}(t)');
title('Truba');

subplot(2,2,4)
plot(t_violina,x_violina);
xlabel('t[s]');
ylabel('x_{violina}(t)');
title('Violina');

%% Zadatak 2 pod b)

N_flauta = 2^nextpow2(length(x_flauta));
N_klavir = 2^nextpow2(length(x_klavir));
N_truba = 2^nextpow2(length(x_truba));
N_violina = 2^nextpow2(length(x_violina));

X_flauta = fft(x_flauta,N_flauta)/length(x_flauta);
X_klavir = fft(x_klavir,N_klavir)/length(x_klavir);
X_truba = fft(x_truba,N_truba)/length(x_truba);
X_violina = fft(x_violina,N_violina)/length(x_violina);

f_flauta = 0:Fs_flauta/N_flauta:Fs_flauta/2;
f_klavir = 0:Fs_klavir/N_klavir:Fs_klavir/2;
f_truba = 0:Fs_truba/N_truba:Fs_truba/2;
f_violina = 0:Fs_violina/N_violina:Fs_violina/2;

amp_flauta = abs(X_flauta(1:N_flauta/2+1));
amp_flauta(2:N_flauta/2+1) = 2*amp_flauta(2:N_flauta/2+1);

amp_klavir = abs(X_klavir(1:N_klavir/2+1));
amp_klavir(2:N_klavir/2+1) = 2*amp_klavir(2:N_klavir/2+1);

amp_truba = abs(X_truba(1:N_truba/2+1));
amp_truba(2:N_truba/2+1) = 2*amp_truba(2:N_truba/2+1);

amp_violina = abs(X_violina(1:N_violina/2+1));
amp_violina(2:N_violina/2+1) = 2*amp_violina(2:N_violina/2+1);

figure(5)

subplot(4,1,1)
plot(f_flauta,amp_flauta);
xlabel('f[Hz]');
ylabel('|X_{flauta}(jf)|');
title('Amplitudska frekvencijska karakteristika za flautu');

subplot(4,1,2)
plot(f_klavir,amp_klavir);
xlabel('f[Hz]');
ylabel('|X_{klavir}(jf)|');
title('Amplitudska frekvencijska karakteristika za klavir');

subplot(4,1,3)
plot(f_truba,amp_truba);
xlabel('f[Hz]');
ylabel('|X_{truba}(jf)|');
title('Amplitudska frekvencijska karakteristika za trubu');

subplot(4,1,4)
plot(f_violina,amp_violina);
xlabel('f[Hz]');
ylabel('|X_{violina}(jf)|');
title('Amplitudska frekvencijska karakteristika za violinu');

ph_flauta = unwrap(angle(X_flauta(1:N_flauta/2+1)));
ph_klavir = unwrap(angle(X_klavir(1:N_klavir/2+1)));
ph_truba = unwrap(angle(X_truba(1:N_truba/2+1)));
ph_violina = unwrap(angle(X_violina(1:N_violina/2+1)));

figure(6)

subplot(4,1,1)
plot(f_flauta,ph_flauta);
xlabel('f[Hz]');
ylabel('arg(X_{flauta}(jf))');
title('Fazna frekvencijska karakteristika za flautu');

subplot(4,1,2)
plot(f_klavir,ph_klavir);
xlabel('f[Hz]');
ylabel('arg(X_{klavir}(jf))');
title('Fazna frekvencijska karakteristika za klavir');

subplot(4,1,3)
plot(f_truba,ph_truba);
xlabel('f[Hz]');
ylabel('arg(X_{truba}(jf))');
title('Fazna frekvencijska karakteristika za trubu');

subplot(4,1,4)
plot(f_violina,ph_violina);
xlabel('f[Hz]');
ylabel('arg(X_{violina}(jf))');
title('Fazna frekvencijska karakteristika za violinu');

%% Zadatak 2 pod c) i pod d)

%FLAUTA

[pik_vred,pik_frekv] = findpeaks(amp_flauta);

pik_frekv = pik_frekv(pik_vred > 0.05*max(pik_vred));
pik_vred = pik_vred(find(pik_vred > 0.05*max(pik_vred)));

pik_vred = pik_vred(3:length(pik_vred));

pik_frekv = pik_frekv*Fs_flauta/N_flauta;
pik_frekv = pik_frekv(3:length(pik_frekv));

i = 1;
count = 1;


while(i < length(pik_frekv))
    j = i + 1;
    while(j < length(pik_frekv) && abs(pik_frekv(i) - pik_frekv(j)) < 100)
       j = j + 1;
    end
    frekvencije_flauta(count) = pik_frekv(find(pik_vred ==  max(pik_vred(i:j))));
    count = count + 1;
    i = j + 1;
end

frekvencije_flauta;
ton_flauta = frekvencije_flauta(1);
frekvencije_flauta = diff(frekvencije_flauta);
prosecan_razmak_flauta = mean(frekvencije_flauta);

% KLAVIR

[pik_vred,pik_frekv] = findpeaks(amp_klavir);

pik_frekv = pik_frekv(find(pik_vred > 0.05*max(pik_vred)));
pik_vred = pik_vred(find(pik_vred > 0.05*max(pik_vred)));

pik_vred = pik_vred(3:length(pik_vred));

pik_frekv = pik_frekv*Fs_klavir/N_klavir;
pik_frekv = pik_frekv(3:length(pik_frekv));

i = 1;
count = 1;


while(i < length(pik_frekv))
    j = i + 1;
    while(j < length(pik_frekv) && abs(pik_frekv(i) - pik_frekv(j)) < 100)
       j = j + 1;
    end
    frekvencije_klavir(count) = pik_frekv(find(pik_vred ==  max(pik_vred(i:j))));
    count = count + 1;
    i = j + 1;
end

frekvencije_klavir;
ton_klavir = frekvencije_klavir(1);
frekvencije_klavir = diff(frekvencije_klavir);
prosecan_razmak_klavir = mean(frekvencije_klavir);

%TRUBA

[pik_vred,pik_frekv] = findpeaks(amp_truba);

pik_frekv = pik_frekv(find(pik_vred > 0.05*max(pik_vred)));
pik_vred = pik_vred(find(pik_vred > 0.05*max(pik_vred)));

pik_vred = pik_vred(3:length(pik_vred));

pik_frekv = pik_frekv*Fs_truba/N_truba;
pik_frekv = pik_frekv(3:length(pik_frekv));

i = 1;
count = 1;


while(i < length(pik_frekv))
    j = i + 1;
    while(j < length(pik_frekv) && abs(pik_frekv(i) - pik_frekv(j)) < 100)
       j = j + 1;
    end
    frekvencije_truba(count) = pik_frekv(find(pik_vred ==  max(pik_vred(i:j))));
    count = count + 1;
    i = j + 1;
end

frekvencije_truba;
ton_truba= frekvencije_truba(1);
frekvencije_truba = diff(frekvencije_truba);
prosecan_razmak_truba = mean(frekvencije_truba);

%VIOLINA

[pik_vred,pik_frekv] = findpeaks(amp_violina);

pik_frekv = pik_frekv(find(pik_vred > 0.05*max(pik_vred)));
pik_vred = pik_vred(find(pik_vred > 0.05*max(pik_vred)));

pik_vred = pik_vred(3:length(pik_vred));

pik_frekv = pik_frekv*Fs_violina/N_violina;
pik_frekv = pik_frekv(3:length(pik_frekv));

i = 1;
count = 1;


while(i < length(pik_frekv))
    j = i + 1;
    while(j < length(pik_frekv) && abs(pik_frekv(i) - pik_frekv(j)) < 100)
       j = j + 1;
    end
    frekvencije_violina(count) = pik_frekv(find(pik_vred ==  max(pik_vred(i:j))));
    count = count + 1;
    i = j + 1;
end

frekvencije_violina;
ton_violina = frekvencije_violina(1);
frekvencije_violina = diff(frekvencije_violina);
prosecan_razmak_violina = mean(frekvencije_violina);

disp(ton_flauta);
disp(ton_klavir);
disp(ton_truba);
disp(ton_violina);

disp(prosecan_razmak_truba);