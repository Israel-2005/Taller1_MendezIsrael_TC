% Práctica 1 - Identificación Gráfica de Sistemas

clear;

clc;

close all;

%% ===================== 1. LECTURA DE DATOS =====================

table = readtable("data_motor.csv");

tM = table.time_t_;

uM = table.ex_signal_u_;

yM = table.system_response_y_;

%% ===================== 2. GRAFICA DE LOS DATOS REALES =====================

figure;

plot(tM,uM,'g','LineWidth',1.5); hold on;

plot(tM,yM,'b','LineWidth',1.8);

grid on;

xlabel('Tiempo [s]');

ylabel('Amplitud');

title('Señal de entrada y respuesta real del sistema');

legend('Entrada u(t)','Salida real y(t)','Location','best');

%% ===================== 3. METODO ANALITICO =====================

km = 1/1.5;

tau = 0.5;

theta = 0.2;

num = km;

den = [tau 1];

G1 = tf(num, den, 'InputDelay', theta);

ym1 = lsim(G1,uM,tM);

fit1 = 100 * (1 - norm(yM - ym1)/norm(yM - mean(yM)));

%% ===================== 4. DATOS BASICOS DEL SISTEMA =====================

% Promedios para evitar ruido

u_inicial = mean(uM(1:5));

u_final = mean(uM(end-5:end));

y_inicial = mean(yM(1:5));

y_final = mean(yM(end-5:end));

% Ganancia del sistema

K_base = (y_final - y_inicial)/(u_final - u_inicial);

%% ===================== 5. DETECCION DEL ESCALON =====================

% Detectar el escalón usando un umbral más robusto

delta_u_total = u_final - u_inicial;

umbral = 0.05 * abs(delta_u_total); % 5% del cambio total

indice_escalon = find(abs(uM - u_inicial) > umbral, 1, 'first');

t_escalon = tM(indice_escalon);

%% ===================== 6. PUNTOS DE LA TANGENTE =====================

% MODIFICACIÓN: Para un sistema que se asemeja más a un primer orden puro,

% la pendiente máxima está al inicio. Tomamos el 15% y 45% para capturar

% la parte más empinada de la subida y esquivar el ruido inicial de los datos.

y_15 = y_inicial + 0.15 * (y_final - y_inicial);

y_45 = y_inicial + 0.45 * (y_final - y_inicial);

indice_15 = find(yM >= y_15, 1, 'first');

indice_45 = find(yM >= y_45, 1, 'first');

t1 = tM(indice_15);

y1_tan = yM(indice_15);

t2 = tM(indice_45);

y2_tan = yM(indice_45);

% Pendiente de la tangente (m = dy/dt)

m = (y2_tan - y1_tan)/(t2 - t1);

% Intercepto de la tangente (y = mx + b -> b = y - m*x)

b = y1_tan - m*t1;

% Corte con línea base

t_base = (y_inicial - b)/m;

% Corte con línea final

t_final_tan = (y_final - b)/m;

%% ===================== 7. METODO MILLER =====================

K_miller = K_base;

% Retardo: desde que ocurre el escalón hasta que la tangente corta la base

theta_miller = t_base - t_escalon;

% Si sale negativo, lo forzamos a cero

if theta_miller < 0

theta_miller = 0;

end

% Nivel 63.2%

y_63 = y_inicial + 0.632*(y_final - y_inicial);

% Buscar el instante donde la salida alcanza el 63.2%

if y_final >= y_inicial

indice_63 = find(yM >= y_63, 1, 'first');

else

indice_63 = find(yM <= y_63, 1, 'first');

end

t_63 = tM(indice_63);

% Constante de tiempo: T_63 menos el tiempo de escalón y el retardo de Miller

tau_miller = t_63 - t_escalon - theta_miller;

% Si tau sale negativa o cero, la corregimos

if tau_miller <= 0

tau_miller = 0.001;

end

% Función de transferencia Miller

num_miller = K_miller;

den_miller = [tau_miller 1];

G_miller = tf(num_miller, den_miller, 'InputDelay', theta_miller);

% Respuesta Miller

y_miller = lsim(G_miller,uM,tM);

% Fit Miller

fit_miller = 100 * (1 - norm(yM - y_miller)/norm(yM - mean(yM)));

%% ===================== 8. METODO ZIEGLER-NICHOLS =====================

K_zn = K_base;

% Retardo Ziegler-Nichols (L)

theta_zn = t_base - t_escalon;

if theta_zn < 0

theta_zn = 0;

end

% Constante de tiempo Ziegler-Nichols (T)

tau_zn = t_final_tan - t_base;

if tau_zn <= 0

tau_zn = 0.001;

end

num_zn = K_zn;

den_zn = [tau_zn 1];

G_zn = tf(num_zn, den_zn, 'InputDelay', theta_zn);

y_zn = lsim(G_zn,uM,tM);

fit_zn = 100 * (1 - norm(yM - y_zn)/norm(yM - mean(yM)));

%% ===================== 9. GRAFICA GENERAL =====================

figure;

plot(tM, yM, 'k', 'LineWidth', 2); hold on;

plot(tM, ym1, 'b--', 'LineWidth', 1.8);

plot(tM, y_miller, 'r-.', 'LineWidth', 1.8);

plot(tM, y_zn, 'm:', 'LineWidth', 2);

plot(tM, uM, 'g', 'LineWidth', 1.2);

grid on;

xlabel('Tiempo [s]');

ylabel('Amplitud');

title('Comparación: Real vs Analítico vs Miller vs Ziegler-Nichols');

legend('Salida real','Analítico','Miller','Ziegler-Nichols','Entrada','Location','best');

%% ===================== 10. GRAFICA AUXILIAR DE LA TANGENTE =====================

figure;

plot(tM,yM,'b','LineWidth',1.8); hold on;

y_tangente = m*tM + b;

plot(tM,y_tangente,'k--','LineWidth',1.5);

yline(y_inicial,':','Línea base');

yline(y_final,':','Línea final');

xline(t_escalon,':','t escalón');

xline(t_base,':','t base');

xline(t_final_tan,':','t final tangente');

% Añadimos los puntos usados para calcular la tangente (opcional, muy útil visualmente)

plot(t1, y1_tan, 'ko', 'MarkerFaceColor', 'r');

plot(t2, y2_tan, 'ko', 'MarkerFaceColor', 'r');

grid on;

% Limitamos el eje Y para que la tangente no estire demasiado el gráfico

ylim([min(y_inicial,0)-0.2, y_final+0.2]);

xlabel('Tiempo [s]');

ylabel('Amplitud');

title('Visualización de la tangente');

legend('Salida real','Tangente','Línea base','Línea final','t escalón','t base','t final tang','Puntos 15% y 45%','Location','best');

%% ===================== 11. RESULTADOS =====================

disp('==================== RESULTADOS ====================');

disp(' ');

disp('--- METODO ANALITICO ---');

fprintf('K = %.4f\n', km);

fprintf('tau = %.4f s\n', tau);

fprintf('theta = %.4f s\n', theta);

fprintf('fit analitico = %.2f %%\n', fit1);

disp(G1);

disp(' ');

disp('--- METODO MILLER ---');

fprintf('K = %.4f\n', K_miller);

fprintf('tau = %.4f s\n', tau_miller);

fprintf('theta = %.4f s\n', theta_miller);

fprintf('fit Miller = %.2f %%\n', fit_miller);

disp(G_miller);

disp(' ');

disp('--- METODO ZIEGLER-NICHOLS ---');

fprintf('K = %.4f\n', K_zn);

fprintf('tau = %.4f s\n', tau_zn);

fprintf('theta = %.4f s\n', theta_zn);

fprintf('fit Ziegler-Nichols = %.2f %%\n', fit_zn);

disp(G_zn);