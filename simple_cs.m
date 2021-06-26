clear all
close all
clc

%% Some Initialization of arbitary signal
n = 500; %Original Number sample
t = linspace(0, 1, n);
sampling_f=n;
x = cos(2 * 10 * pi * t) + 0.7*cos(2 * 80 * pi * t);
xt = fft(x); 
PSD = xt.* conj(xt)/n; 
f = sampling_f*(0:(length(xt)/2))/length(xt);
%% Randomly sample that signal
p = 150;  %Samples under CS
perm = round(rand(p, 1) * n);
y = x(perm); 
t_permuted=t(perm);
%% Construct standard CS form y=Theta*x
Psi = dct(eye(n, n)); 
Theta = Psi(perm, :);
Phi=perm;
for k=1:n-1
    Phi=[Phi perm];
end
%% L1 norm minimzation
cvx_begin;
variable s(n);
minimize( norm(s,1) );
subject to
Theta*s==y';
cvx_end;
%% Reconstruction
xrec = idct(s); % reconstruct full signal
xtrec = fft(xrec); 
PSDrec = xtrec.* conj(xtrec)/n; 
%% Visualization
figure
subplot(2,2,1)
plot(t,x);
title('Original Signal in T-domain');
xlabel('time(sec)');
ylabel('Amplitude(au)');
ylim([-2 2]);

subplot(2,2,2)
plot(f,PSD(1:end/2+1));
title('Original Signal in F-domain');
xlabel('frequency(\nu)');
ylabel('PSD(|F|^2)');
xlim([0 140]);

subplot(2,2,3)
plot(t,xrec);
title(join(["Reconstructed Signal in T-domain with p/n=",num2str(p/n)]));
xlabel('time(sec)');
ylabel('Amplitude(au)');
ylim([-2 2]);

subplot(2,2,4)
plot(f,PSDrec(1:end/2+1));
title('Reconstructed Signal in F-domain');
xlabel('frequency(\nu)');
ylabel('PSD(|F|^2)');
xlim([0 140]);

figure
subplot(2,2,1)
plot(y)
title('y');
xlabel('sample #');
ylabel('Amplitude(au)');
xlim([0 p]);

subplot(2,2,2)
plot(xrec)
title('x');
xlabel('sample #');
ylabel('Amplitude(au)');
xlim([0 n]);

subplot(2,2,3)
plot(s)
title('\psi');
xlabel('sample #');
ylabel('Amplitude(au)');
xlim([0 n]);

figure
heatmap(Theta, 'Colormap',gray,'GridVisible','off');
title('\Theta');

figure
heatmap(Psi, 'Colormap',gray,'GridVisible','off');
title('D');

figure
heatmap(Phi, 'Colormap',gray,'GridVisible','off');
title('\Phi');

figure
subplot(3,1,1)
plot(t,xrec)
xlabel('time(sec)');
ylabel('Amplitude(au)');
title('Reconstructed Signal')
ylim([-2 2]);
xlim([0 0.5]);

subplot(3,1,2)
plot(t,x)
xlabel('time(sec)');
ylabel('Amplitude(au)');
title('Original Signal')
ylim([-2 2]);
xlim([0 0.5]);

subplot(3,1,3)
plot(t,abs(x-xrec'))
xlabel('time(sec)');
ylabel('Amplitude(au)');
title('Absolute Error')
ylim([0 2]);
xlim([0 0.5]);






