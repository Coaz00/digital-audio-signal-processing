%Aleksandar Djordjevic 2019/0086

%% Zadatak 1 - ucitavanje signala

clear all;
close all;
clc

godina_upisa = 2019;
broj_indeksa = 86;

P = mod(broj_indeksa,4) + 4;

n = 0:1:(P - 1);

n1 = n(1:floor(P/2));
n2 = n(floor(P/2) + 1:P);

x1 = sin(n1) + 2*cos(2*n1) + P/2;
x2 = 2.^n2;

x = [x1 x2];

y1 = (-1).^n1 + mod(n1,2);
y2 = n2 - P/2;

y = [y1 y2];

%% Zadatak 1 pod a)

figure(1)

subplot(2,1,1)
stem(n,x,'r');
xlabel('n[odb]');
ylabel('x(n))');
title('Signal x(n)');

subplot(2,1,2)
stem(n,y, 'r');
xlabel('n[odb]');
ylabel('y(n)');
title('Signal y(n)');

N = P + 2;

c = zeros(1,N);

add_zeros = zeros(1, N - length(x)); % niz nula koji treba dodati na signale

x_p2 = [x add_zeros];
y_p2 = [y add_zeros];

% racunanje ciklicne konvolucije
for k = 0:N-1
    c(k+1) = x_p2*y_p2(mod(k:-1:k-(N - 1),N)+1)';
end

n_p2 = 0:1:N-1;

c2 = cconv(x,y,N); % provera

figure(2)
stem(n_p2,c2);
hold on;
stem(n_p2,c);
hold off;
xlabel('n[odb]');
ylabel('c(n)');
title('Ciklicna konvolucija signala x i y u P + 2 tacaka');
legend('cconv','manuelno');

%% Zadatak 1 pod b)

N = length(x) + length(y) - 1;

n_l = 0:1:N-1;

cc = zeros(1,N);

add_zeros = zeros(1, N - length(x));

x_2pm1 = [x add_zeros];
y_2pm1 = [y add_zeros];

% racunanje ciklicne konvolucije u 2P - 1 tacaka
for k = 0:N-1
    cc(k+1) = x_2pm1*y_2pm1(mod(k:-1:k-(N - 1),N)+1)';
end

lc = zeros(1,N);

% racunanje linearne konvolucije
for i = 1:N
    for j = 1:i
        lc(i) = lc(i) + x_2pm1(j)*y_2pm1(i - j + 1);
    end
end

figure(3)

subplot(2,1,1);
stem(n_l,lc,'r')
title('Linearna konvolucija singala x i y');
xlabel('n[odb]');
ylabel('lc(n)');

subplot(2,1,2)
stem(n_l,cc,'r')
title('Ciklicna konvolucija singala x i y u 2P - 1 tacaka');
xlabel('n[odb]');
ylabel('cc(n)');

disp(lc == conv(x,y)); % provera jednakosti ugradjene funckije za lin konv i moje
disp(lc == cc); % provera jednakosti ciklicne i lin konv