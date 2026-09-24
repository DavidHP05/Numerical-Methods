%% =========================================================================
%  PRACTICA 1 - INTERPOLACION NUMERICA (LAGRANGE)
%  =========================================================================
%  Problema 1:
%   A. Polinomios de base de Lagrange y polinomio interpolador de Lagrange.
%   1. Aplicacion al ejemplo de clase: f(x) = sin(x) en [0, pi/4],
%      con distintos soportes.
%   2. Aplicacion a f(x) = exp(-x) + cos(4x/pi) en [0,2],
%      con distintos soportes (2, 3 y mas puntos).
%
%  En cada apartado se compara la funcion con el polinomio interpolador
%  mediante representaciones graficas (funcion vs. interpolante, y error
%  absoluto cometido).
%  =========================================================================

clear; clc; close all;

%% ------------------------------------------------------------------------
%  APARTADO 1: f(x) = sin(3*x) en [0, pi/4]  (ejemplo resuelto en clase)
%  ------------------------------------------------------------------------
f1 = @(x) sin(3*x);
a1 = 0; b1 = pi/4;

% Distintos soportes pedidos en el enunciado:
soportes1 = { ...
    [0, pi/4], ...              % 2 puntos
    [0, pi/8, pi/4], ...        % 3 puntos
    linspace(a1, b1, 11) };     % 11 puntos equidistantes (con MATLAB)

nombres1 = { ...
    '2 puntos: \{0, \pi/4\}', ...
    '3 puntos: \{0, \pi/8, \pi/4\}', ...
    '11 puntos equidistantes' };

xx1 = linspace(a1, b1, 500);   % malla fina para representar f y p_n

figure('Name', 'Apartado 1: f(x) = sin(x) en [0,pi/4]', 'NumberTitle', 'off');
for k = 1:numel(soportes1)
    xs = soportes1{k};
    ys = f1(xs);
    pk = lagrange_interp(xs, ys, xx1);

    subplot(numel(soportes1), 2, 2*k-1);
    plot(xx1, f1(xx1), 'b-', 'LineWidth', 1.5); hold on;
    plot(xx1, pk, 'r--', 'LineWidth', 1.5);
    plot(xs, ys, 'ko', 'MarkerFaceColor', 'k');
    title(['f(x) vs p_n(x) - ' nombres1{k}]);
    xlabel('x'); ylabel('y');
    legend('f(x) = sin(x)', 'p_n(x)', 'Soporte', 'Location', 'best');
    grid on;

    subplot(numel(soportes1), 2, 2*k);
    plot(xx1, abs(f1(xx1) - pk), 'm-', 'LineWidth', 1.5);
    title(['Error absoluto |f(x)-p_n(x)| - ' nombres1{k}]);
    xlabel('x'); ylabel('Error'); grid on;
end
sgtitle('Interpolacion de Lagrange para f(x) = sin(x) en [0, \pi/4]');

%% ------------------------------------------------------------------------
%  APARTADO 2: f(x) = e^{-x} + cos(4x/pi) en [0,2]
%  ------------------------------------------------------------------------
f2 = @(x) exp(-x) + cos(4*x/pi);
a2 = 0; b2 = 2;

npuntos2 = [2, 3, 4, 6, 11];   % soportes de 2, 3 y mas puntos (equidistantes)
xx2 = linspace(a2, b2, 500);

figure('Name', 'Apartado 2: f(x) = e^{-x}+cos(4x/pi) en [0,2]', 'NumberTitle', 'off');
for k = 1:numel(npuntos2)
    n = npuntos2(k);
    xs = linspace(a2, b2, n);
    ys = f2(xs);
    pk = lagrange_interp(xs, ys, xx2);

    subplot(numel(npuntos2), 2, 2*k-1);
    plot(xx2, f2(xx2), 'b-', 'LineWidth', 1.5); hold on;
    plot(xx2, pk, 'r--', 'LineWidth', 1.5);
    plot(xs, ys, 'ko', 'MarkerFaceColor', 'k');
    title(sprintf('f(x) vs p_{%d}(x)  (%d puntos)', n-1, n));
    xlabel('x'); ylabel('y');
    legend('f(x)', 'p_n(x)', 'Soporte', 'Location', 'best');
    grid on;

    subplot(numel(npuntos2), 2, 2*k);
    plot(xx2, abs(f2(xx2) - pk), 'm-', 'LineWidth', 1.5);
    title(sprintf('Error absoluto (%d puntos)', n));
    xlabel('x'); ylabel('Error'); grid on;
end
sgtitle('Interpolacion de Lagrange para f(x) = e^{-x} + cos(4x/\pi) en [0,2]');

%% ------------------------------------------------------------------------
%  Comparativa conjunta: evolucion del error maximo segun el soporte
%  ------------------------------------------------------------------------
errores_max = zeros(size(npuntos2));
for k = 1:numel(npuntos2)
    n = npuntos2(k);
    xs = linspace(a2, b2, n);
    ys = f2(xs);
    pk = lagrange_interp(xs, ys, xx2);
    errores_max(k) = max(abs(f2(xx2) - pk));
end

figure('Name', 'Error maximo vs numero de puntos', 'NumberTitle', 'off');
semilogy(npuntos2, errores_max, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlabel('Numero de puntos del soporte (n+1)');
ylabel('Error maximo en [0,2] (escala log)');
title('Evolucion del error maximo de interpolacion segun el soporte');
grid on;

%% =========================================================================
%  APARTADO A: FUNCIONES LOCALES
%  =========================================================================

function Li = lagrange_base(xs, i, x)
% LAGRANGE_BASE  Polinomio de base de Lagrange L_i(x) de grado n = length(xs)-1
%
%   xs : vector con el soporte {x0, x1, ..., xn}
%   i  : indice del polinomio base a calcular (1-indexado, i = 1..n+1)
%   x  : punto o vector de puntos donde se evalua L_i(x)
%
%   L_i(x) = prod_{j=0, j~=i} (x - x_j) / (x_i - x_j)
    n = length(xs);
    Li = ones(size(x));
    for j = 1:n
        if j ~= i
            Li = Li .* (x - xs(j)) / (xs(i) - xs(j));
        end
    end
end

function pn = lagrange_interp(xs, ys, x)
% LAGRANGE_INTERP  Polinomio interpolador de Lagrange p_n(x)
%
%   xs : vector con el soporte {x0, ..., xn}
%   ys : vector con los valores f(x0), ..., f(xn)
%   x  : punto o vector de puntos donde se evalua p_n(x)
%
%   p_n(x) = sum_i f_i * L_i(x)
    n = length(xs);
    pn = zeros(size(x));
    for i = 1:n
        pn = pn + ys(i) * lagrange_base(xs, i, x);
    end
end
