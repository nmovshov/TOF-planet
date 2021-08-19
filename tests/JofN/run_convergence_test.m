%% How many density levels and shape levels are required for desired precision?
clear
clc
close all

%% Choose layer numbers and distribution strategies to investigate
N = [128, 256];
nx = [32, 64];
zstrat = {@zvecs.equal_dr};
toforder = [4,7];

%% The experiment MAY be sensitive to the chosen density model
% Here we use a smooth polynomial model, with no discontinuities or even steep
% slopes. The goal is to find the *minimal* resolution needed for a desired J
% precisoion. More complex density profiles will surely require higher
% resolution, and/or special treatment of discontinuous density.

rho0 = 0.199;
p=[-9.1352277e+02, -2.0616194e+03, -2.5768977e+02, -6.8877550e+01,...
    8.6817818e+03, 4.2076235e+03, -1.3579913e+04, 0.0, 3.9924165e+03];
mrot = 0.1574;

%% Create a table to store experiment results
vars = {'toforder','N','nx','zstrat','runtime','J2','J4'};
tps = {'double','double','double','string','double','double','double'};
nN = length(N);
nnx = length(nx);
nZ = length(zstrat);
no = length(toforder);
T = table('Size',[no*nN*nnx*nZ, length(vars)],...
    'VariableTypes', tps,...
    'VariableNames', vars);

%% Call tof, this will take a while with tof7!
fname = sprintf('%f.mat',now());
fprintf("Running TOF.\n")
for i=1:nZ
    for j=1:nN
        for k=1:nnx
            for o=1:no
                row = (i-1)*(nN*nnx*no) + (j-1)*(nnx*no) + (k-1)*(no) + o;
                fprintf('working on row %d of %d (N=%d, nx=%d, tof%d)...',...
                    row,height(T),N(j),nx(k),toforder(o));
                T.zstrat(row) = string(char(zstrat{i}));
                T.toforder(row) = toforder(o);
                T.N(row) = N(j);
                T.nx(row) = nx(k);
                zvec = zstrat{i}(N(j));
                dvec = polyval(p, zvec);
                tic
                if toforder(o) == 4
                    Js = tof4(zvec, dvec, mrot, 'xlevels', nx(k));
                else
                    Js = tof7(zvec, dvec, mrot, 'xlevels', nx(k));
                end
                T.runtime(row) = toc;
                T.J2(row) = Js(2);
                T.J4(row) = Js(3);
                fprintf('done.\n')
                save(fname);
            end
        end
    end
end
fprintf('All done. (%s)\n',seconds2human(sum(T.runtime)));
