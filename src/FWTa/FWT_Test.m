% L = 4*1024;
% fmin = 0.01;
% fmax = 0.4;
% x0 = chirp(1:L,fmin,L,fmax,[],90);
% % x = repmat(x.',[1 10]);
% x = zeros(1,L);
% x(2*1024+1) = 1;
% x(2*1024) = -1;


L = 2001;        % Length of every sequence.
Ls = 801;        % Length of the signal part.
Ls0 = 600;       % Starting point of the signal

sz = 2*round(0.4*Ls/2);
tmp = hanning(sz);
win = [tmp(1:sz/2).' ones(1,Ls-sz) tmp(sz/2+1:end).'];

t = 0:Ls-1;
f0 = 0.005; f1 = 0.1;
f = f0 + (f1 - f0)/3 * (t/t(end)).^2;

xs = sin(2*pi*f.*t) .* win;
x0 = [zeros(1,Ls0) xs zeros(1, L-Ls0-Ls)];



L = 1001;        % Length of every sequence.
Ls = 801;        % Length of the signal part.
Ls0 = 100;       % Starting point of the signal

sz = 2*round(0.4*Ls/2);
sze = 2*round(0.4*Ls/2);
tmp = hanning(sz);
tmpe = hanning(sze);
win = [tmp(1:sz/2).' ones(1,Ls-(sz+sze)/2) tmpe(sze/2+1:end).'];

t = 0:Ls-1;
f0 = 0.005; f1 = 0.05;
f = f0 + (f1 - f0)/3 * (t/t(end)).^2;

f0 = 0.005; f1 = 0.03;
a = (f1/f0)^(1/t(end));
f = f0*(a.^t);

xs = sin(2*pi*f.*t) .* win;
x0 = [zeros(1,Ls0) xs zeros(1, L-Ls0-Ls)];




% x = randn(L,10);
x = x0; %  + 0.1*randn(size(x0));


WT.Family = -1;
WT.J = 8;
WT.V = 6;
WT.s0 = 4;
WT.b0 = 1;
WT.op1 = 5*(2*sqrt(log(2))); % pi*sqrt(2/log(2));
WT.continuous = 0;
WT.convtype = 0;
WT.RealInvOutput = 0;

% WT.Family = -1;
% WT.J = 7;
% WT.V = 6;
% WT.s0 = 2;
% WT.b0 = 0.5;
% WT.op1 = pi*sqrt(2/log(2));
% WT.continuous = 0;
% WT.convtype = 0;
% WT.RealInvOutput = 0;

tic; [y, WF] = gw_FWT(x, WT); toc
tic; xrec  = gw_IFWT(y, WT, L); toc
