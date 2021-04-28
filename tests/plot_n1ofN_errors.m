function fh = plot_n1ofN_errors(T_errs)
%PLOT_N1OFN_ERRORS Plot fractional diffs between my Js and WH16; Nadine style.

% Organize data to plot
n = height(T_errs)/2; % it should always be even; for tof4 and tof7
N = T_errs.N(1:n);
j2err4 = abs(T_errs.J2(1:n));
j2err7 = abs(T_errs.J2(n+1:end));
j4err4 = abs(T_errs.J4(1:n));
j4err7 = abs(T_errs.J4(n+1:end));
j6err4 = abs(T_errs.J6(1:n));
j6err7 = abs(T_errs.J6(n+1:end));
j8err4 = abs(T_errs.J8(1:n));
j8err7 = abs(T_errs.J8(n+1:end));
j10err4 = abs(T_errs.J10(1:n));
j10err7 = abs(T_errs.J10(n+1:end));
j12err4 = abs(T_errs.J12(1:n));
j12err7 = abs(T_errs.J12(n+1:end));

% plot, nadine style
subplot(2,3,1);
loglog(N, j2err4, 'b--', 'LineWidth',2);
hold
loglog(N, j2err7, 'b-^', 'LineWidth',2);
title('J2')

subplot(2,3,2);
loglog(N, j4err4, 'b--', 'LineWidth',2);
hold
loglog(N, j4err7, 'b-^', 'LineWidth',2);
title('J4')

subplot(2,3,3);
loglog(N, j6err4, 'b--', 'LineWidth',2);
hold
loglog(N, j6err7, 'b-^', 'LineWidth',2);
title('J6')

subplot(2,3,4);
loglog(N, j8err4, 'b--', 'LineWidth',2);
hold
loglog(N, j8err7, 'b-^', 'LineWidth',2);
title('J8')

subplot(2,3,5);
loglog(N, j10err4, 'b--', 'LineWidth',2);
hold
loglog(N, j10err7, 'b-^', 'LineWidth',2);
title('J10')

subplot(2,3,6);
loglog(N, j12err4, 'b--', 'LineWidth',2);
hold
loglog(N, j12err7, 'b-^', 'LineWidth',2);
title('J12')

% Style
fh = gcf;
fh.WindowState = 'max';
end
