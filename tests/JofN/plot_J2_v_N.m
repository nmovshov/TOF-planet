function [fh, ah] = plot_J2_v_N(T, kwargs)
arguments
    T (:,:) table
    kwargs.logx (1,1) logical = false
    kwargs.logy (1,1) logical = false
end
nx = unique(T.nx);
tofo = unique(T.toforder);
[fh, ah] = ngraf.get_canvas('proj');
for k=1:length(nx)
    for j=1:length(tofo)
        ind = (T.nx == nx(k)) & (T.toforder == tofo(j));
        lh = plot(T.N(ind), T.J2(ind), '--+');
        lh.DisplayName = sprintf('nx=%d, tof%d',nx(k),tofo(j));
    end
end
if kwargs.logx, ah.XScale = 'log'; end
if kwargs.logy, ah.YScale = 'log'; end
xlabel('N')
ylabel('$J_2$')
legend('Location','nw')

end
