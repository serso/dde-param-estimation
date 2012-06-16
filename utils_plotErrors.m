t = get(gco, 'XData');
x = get(gco, 'YData');

tSol = get(gco, 'XData');
xSol = get(gco, 'YData');

figure;
hold on;
grid on;

t = t(end - length(tSol) + 1:end);
x = x(end - length(tSol) + 1:end);

errors = x - xSol;

plot(t, errors, '.b');