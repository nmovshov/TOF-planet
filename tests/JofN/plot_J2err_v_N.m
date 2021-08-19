function [fh, ah] = plot_J2err_v_N(T, kwargs)
arguments
    T (:,:) table
    kwargs.logx (1,1) logical = false
    kwargs.logy (1,1) logical = false
end
N = unique(T.N);
nx = unique(T.nx);
tofo = unique(T.toforder);
J2_err = nan(length(N),length(tofo));

for j=1:length(tofo)
    tind = find(T.N == max(N) & T.nx == max(nx) & T.toforder == tofo(j));
    Jinf = T.J2(tind);
    for k=1:length(N)
        ind = T.N == N(k) & T.nx == max(nx) & T.toforder == tofo(j);
        J2_err(k,j) = abs(T.J2(ind) - Jinf);
    end
end
[fh, ah] = ngraf.get_canvas('proj');
for j=1:length(tofo)
    plot(N, J2_err(:,j)*1e6,'--+','DisplayName',sprintf('tof%d',tofo(j)));
end
if kwargs.logx, ah.XScale = 'log'; end
if kwargs.logy, ah.YScale = 'log'; end
xlabel(sprintf('N (nx = %d)', max(nx)))
ylabel('$10^6(J_2-J_2^c)$')
legend('Location','ne')
hline(0.04)
end
