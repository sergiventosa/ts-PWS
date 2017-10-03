n = randn(1000,30);
s = randn(1000,1);
x = n + repmat(s,1,30);

Data.x = x;
Data.reference = s;
tspws.rm = 1;
tspws.wu = 2;
tspws.unbiased = 1;
tspws.Nmax = 10;
tspws.convergence = 1;

out = gw_ts_pws(Data, tspws);