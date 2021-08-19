function [fh, ah] = plot_runtime_v_nx(T, kwargs)
arguments
    T (:,:) table
    kwargs.logx (1,1) logical = false
    kwargs.logy (1,1) logical = false
end
N = unique(T.N);
tofo = unique(T.toforder);
[fh, ah] = ngraf.get_canvas('proj');
for k=1:length(N)
    for j=1:length(tofo)
        ind = (T.N == N(k)) & (T.toforder == tofo(j));
        lh = plot(T.nx(ind), T.runtime(ind), '--+');
        lh.DisplayName = sprintf('N=%d, tof%d',N(k),tofo(j));
    end
end
if kwargs.logx, ah.XScale = 'log'; end
if kwargs.logy, ah.YScale = 'log'; end
xlabel('nx')
ylabel('run time [sec]')
legend('Location','nw')

end
