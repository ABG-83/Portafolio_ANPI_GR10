clear;
clc;
close all;

source('parte1.m');

% Unificar la legibilidad de las figuras y conservar su tamano al exportarlas.
set(0,'defaultaxesfontsize',11,'defaulttextfontsize',11, ...
      'defaultfigurecolor','w','defaultfigurepaperpositionmode','auto');


% Datos del problema
L = 1;
alpha = 0.01;

u0 = @(x) 80*sin(pi*x/L);

T_exacta = @(x,t) 20 + 80*exp(-alpha*(pi/L)^2*t).*sin(pi*x/L);

% Experimento base
M = 50;
dt = 0.01;
tf = 10;

dx = L/M;
mu = alpha*dt/dx^2;

tic;

[x_exp,t_exp,U_exp] = euler_explicito(L,alpha,u0,M,dt,tf);

T_exp = U_exp + 20;

tiempo_exp = toc;

% Persona 2: medir el esquema implicito incluyendo matriz y factorizacion de Thomas.
tic;
[x_imp,t_imp,U_imp,residuos_imp,conv_imp] = euler_implicito(L,alpha,u0,M,dt,tf);
tiempo_imp = toc;
if ~conv_imp
    error('Se detuvo el experimento base de Euler implicito.');
end
T_imp = U_imp + 20;


tic;
[x_cn,t_cn,U_cn,residuos_cn,conv_cn] = crank_nicolson(L,alpha,u0,M,dt,tf);
tiempo_cn = toc;
if ~conv_cn
    error('Se detuvo el experimento base de Crank-Nicolson.');
end
  T_cn = U_cn + 20;

% Verificar estabilidad
fprintf('\nEuler explicito\n');
fprintf('mu = %.4f\n',mu);

if mu <= 0.5
    fprintf('Condicion de estabilidad: cumple\n');
else
    fprintf('Condicion de estabilidad: no cumple\n');
end

% Error al tiempo final
T_ex_final = T_exacta(x_exp,tf);

E_f_exp = max(abs(T_exp(end,:) - T_ex_final));

% Error en toda la malla
E_malla_exp = 0;

for n = 1:length(t_exp)

    T_ex_n = T_exacta(x_exp,t_exp(n));

    error_n = max(abs(T_exp(n,:) - T_ex_n));

    E_malla_exp = max(E_malla_exp,error_n);

end

% Persona 2: validar con la solucion exacta fuera de la medicion de tiempo.
E_f_imp = max(abs(T_imp(end,:) - T_exacta(x_imp,tf)));
E_malla_imp = 0;
for n = 1:length(t_imp)
    error_n = max(abs(T_imp(n,:) - T_exacta(x_imp,t_imp(n))));
    E_malla_imp = max(E_malla_imp,error_n);
end
residuo_max_imp = max(residuos_imp);


E_f_cn = max(abs(T_cn(end,:) - T_exacta(x_cn,tf)));
E_malla_cn = 0;
for n = 1:length(t_cn)
    E_malla_cn = max(E_malla_cn,max(abs(T_cn(n,:) - T_exacta(x_cn,t_cn(n)))));
end
residuo_max_cn = max(residuos_cn);

% Resultados
fprintf('\n');
fprintf('%-18s %-14s %-14s %-12s %-14s %-14s %-14s\n','Metodo','E_f','E_malla','Tiempo(s)','Residuo max','Iter. acum.','Iter. max');
fprintf('%-18s %-14.6e %-14.6e %-12.6f %-14s %-14s %-14s\n','Euler explicito',E_f_exp,E_malla_exp,tiempo_exp,'no aplica','no aplica','no aplica');
fprintf('%-18s %-14.6e %-14.6e %-12.6f %-14.6e %-14s %-14s\n','Euler implicito',E_f_imp,E_malla_imp,tiempo_imp,residuo_max_imp,'no aplica','no aplica');
fprintf('%-18s %-14.6e %-14.6e %-12.6f %-14.6e %-14s %-14s\n','Crank-Nicolson',E_f_cn,E_malla_cn,tiempo_cn,residuo_max_cn,'no aplica','no aplica');
fprintf('Euler implicito y Crank-Nicolson: Thomas con factorizacion reutilizada; todos los residuos < 1e-10.\n');

% Reservar una franja inferior para una leyenda compartida de los cuatro perfiles.
tiempos_perfiles = [0,1,5,10];
posiciones_perfiles = [0.09,0.64,0.37,0.28; 0.59,0.64,0.37,0.28; ...
                       0.09,0.22,0.37,0.28; 0.59,0.22,0.37,0.28];
figure('position',[80,80,1200,900],'paperunits','inches', ...
       'paperposition',[0,0,12,9],'paperpositionmode','manual');
for k = 1:length(tiempos_perfiles)
    tk = tiempos_perfiles(k);
    indice = round(tk/dt) + 1;
    eje_perfil = axes('position',posiciones_perfiles(k,:));
    curva_exp = plot(x_exp,T_exp(indice,:),'o-','MarkerSize',4,'LineWidth',1.2);
    hold on;
    curva_exacta = plot(x_exp,T_exacta(x_exp,tk),'-','LineWidth',1.5);
    xlabel('x (m)');
    ylabel('T (C)');
    title(sprintf('t = %.1f s',tk));
    grid on;
end

% Ubicar una sola leyenda fuera de los datos sin reducir el ultimo panel.
leyenda_perfiles = legend([curva_exp,curva_exacta],{'Euler explicito','Solucion exacta'}, ...
    'location','southoutside','orientation','horizontal','fontsize',10,'box','off');
set(leyenda_perfiles,'units','normalized','position',[0.23,0.025,0.54,0.07]);
set(eje_perfil,'position',posiciones_perfiles(end,:));
print('perfiles_euler_explicito.png','-dpng','-r150');

% Estabilidad de Euler explicito
M_est = 50;
dx_est = L/M_est;
tf_est = 9.6;

mu_1 = 0.4;
mu_2 = 0.6;

dt_1 = mu_1*dx_est^2/alpha;
dt_2 = mu_2*dx_est^2/alpha;

% Condicion inicial perturbada
u0_pert = @(x) 80*sin(pi*x/L) + 1e-6*(-1).^round(x/dx_est);

% Caso mu = 0.4
[x_04,t_04,U_04] = euler_explicito(L,alpha,u0_pert,M_est,dt_1,tf_est);

T_04 = U_04 + 20;

% Caso mu = 0.6
[x_06,t_06,U_06] = euler_explicito(L,alpha,u0_pert,M_est,dt_2,tf_est);

T_06 = U_06 + 20;

% max |T_i^n - 20|
max_04 = max(abs(T_04 - 20),[],2);
max_06 = max(abs(T_06 - 20),[],2);

% Verificar caso inestable
indice_inestable = find(max_06 > 1e6,1);

if ~isempty(indice_inestable)

    tiempo_inestable = t_06(indice_inestable);

    fprintf('\nmu = 0.6 supera 10^6 C en t = %.6f s\n',tiempo_inestable);

end

% Mostrar el caso inestable solo hasta el primer cruce, sin alterar los resultados calculados.
fin_inestable = length(t_06);
if ~isempty(indice_inestable)
    fin_inestable = indice_inestable;
end
figure('position',[80,80,1100,700],'paperunits','inches', ...
       'paperposition',[0,0,11,7],'paperpositionmode','manual');
eje_estabilidad = axes('position',[0.10,0.25,0.85,0.64]);
curva_04 = semilogy(t_04,max_04,'-','LineWidth',1.6);
hold on;
curva_06 = semilogy(t_06(1:fin_inestable),max_06(1:fin_inestable),'-','LineWidth',1.6);
xlabel('Tiempo (s)');
ylabel('max |T_i^n - 20| (C)');
title('Estabilidad de Euler explicito');
xlim([0,tf_est]);
grid on;
set(eje_estabilidad,'yminorgrid','off','xminorgrid','off','yminortick','off');

% Limitar las marcas a seis decadas legibles y senalar el instante de cruce del umbral.
if ~isempty(indice_inestable)
    ylim([10,max(2e6,1.4*max_06(indice_inestable))]);
    set(eje_estabilidad,'ytick',10.^(1:6));
    semilogy([0,tf_est],[1e6,1e6],'k:','LineWidth',1);
    semilogy(tiempo_inestable,max_06(indice_inestable),'ko','MarkerSize',5);
    text(tiempo_inestable+0.25,2e5, ...
         sprintf('Supera 10^6 C\nt = %.3f s',tiempo_inestable),'fontsize',10);
end
leyenda_estabilidad = legend([curva_04,curva_06],{'mu = 0.4','mu = 0.6'}, ...
    'location','southoutside','orientation','horizontal','fontsize',10,'box','off');
set(leyenda_estabilidad,'units','normalized','position',[0.29,0.035,0.42,0.07]);
set(eje_estabilidad,'position',[0.10,0.25,0.85,0.64]);
print('estabilidad_euler_explicito.png','-dpng','-r150');

% Resumen de estabilidad
fprintf('\nEstabilidad de Euler explicito\n');
fprintf('mu = 0.4, dt = %.6f s: cumple mu <= 0.5\n',dt_1);
fprintf('mu = 0.6, dt = %.6f s: no cumple mu <= 0.5\n',dt_2);


% Comparar los perfiles con una sola leyenda situada debajo de los cuatro paneles.
figure('position',[80,80,1200,900],'paperunits','inches', ...
       'paperposition',[0,0,12,9],'paperpositionmode','manual');
for k = 1:length(tiempos_perfiles)
    tk = tiempos_perfiles(k);
    indice = round(tk/dt) + 1;
    eje_perfil = axes('position',posiciones_perfiles(k,:));
    curvas_euler = plot(x_exp,T_exp(indice,:),'o-',x_imp,T_imp(indice,:),'s-', ...
         x_cn,T_cn(indice,:),'^-', x_imp,T_exacta(x_imp,tk),'-','LineWidth',1.2,'MarkerSize',4);
    xlabel('x (m)');
    ylabel('T (C)');
    title(sprintf('t = %.1f s',tk));
    grid on;
end
leyenda_euler = legend(curvas_euler,{'Euler explicito','Euler implicito','Crank-Nicolson','Solucion exacta'}, ...
    'location','southoutside','orientation','horizontal','fontsize',10,'box','off');
set(leyenda_euler,'units','normalized','position',[0.14,0.025,0.72,0.07]);
set(eje_perfil,'position',posiciones_perfiles(end,:));
print('perfiles_comparacion_euler.png','-dpng','-r150');

% Representar los errores con puntos evita la base cero de las barras en escala logaritmica.
figure('position',[80,80,1400,700],'paperunits','inches', ...
       'paperposition',[0,0,14,7],'paperpositionmode','manual');
eje_errores = axes('position',[0.08,0.26,0.37,0.62]);
semilogy(1:3,[E_f_exp,E_f_imp,E_f_cn], 'o-','LineWidth',1.5,'MarkerSize',7);
hold on;
semilogy(1:3,[E_malla_exp,E_malla_imp,E_malla_cn],'x--','LineWidth',1.5,'MarkerSize',9);
set(eje_errores,'xtick',1:3,'xticklabel',{'Euler explicito','Euler implicito','Crank-Nicolson'}, ...
    'yminorgrid','off','xminorgrid','off');
xlim([0.75,3.25]);
ylabel('Error maximo (C)');
title('Errores del experimento base');
grid on;
set(eje_errores,'yminorgrid','off','xminorgrid','off');
leyenda_errores = legend('E_f','E_{malla}','location','southoutside', ...
    'orientation','horizontal','fontsize',10,'box','off');
set(leyenda_errores,'units','normalized','position',[0.14,0.055,0.25,0.06]);
set(eje_errores,'position',[0.08,0.26,0.37,0.62]);

% Usar la misma representacion para comparar tiempos sin barras invertidas.
eje_tiempos = axes('position',[0.59,0.26,0.37,0.62]);
semilogy(1:3,[tiempo_exp,tiempo_imp,tiempo_cn],'o-','LineWidth',1.5,'MarkerSize',7);
set(eje_tiempos,'xtick',1:3,'xticklabel',{'Euler explicito','Euler implicito','Crank-Nicolson'}, ...
    'yminorgrid','off','xminorgrid','off');
xlim([0.75,3.25]);
ylabel('Tiempo (s)');
title('Tiempo del procedimiento numerico');
grid on;
set(eje_tiempos,'yminorgrid','off','xminorgrid','off');
print('comparacion_errores_tiempos_euler.png','-dpng','-r150');

% El refinamiento corresponde a los tres esquemas; incorporar el de Persona 3 al integrarlo.
metodos_ref = {@euler_explicito,@euler_implicito};
nombres_ref = {'Euler explicito','Euler implicito'};
if exist('crank_nicolson')
    metodos_ref{end+1} = @crank_nicolson;
    nombres_ref{end+1} = 'Crank-Nicolson';
else
    fprintf('\nRefinamiento: falta integrar crank_nicolson de Persona 3 para comparar los tres esquemas.\n');
end

% Organizar exactamente las mallas espaciales y temporales del punto 8.
M_espacial = [25,50,100];
dt_espacial = [0.0025,0.0025,0.0025];
M_temporal = [50,50,50];
dt_temporal = [0.016,0.008,0.004];
tipos_ref = {'espacial','temporal'};
resultados_ref = zeros(2,3,length(metodos_ref),5);

% Conservar los residuos por sistema permite revisar tambien la futura integracion de Crank-Nicolson.
residuos_ref = cell(2,3,length(metodos_ref));
residuo_max_ref = NaN(2,3,length(metodos_ref));

% Guardar dx, dt, mu, E_f y tiempo para cada esquema y cada tipo de refinamiento.
for tipo = 1:2
    if tipo == 1
        valores_M = M_espacial;
        valores_dt = dt_espacial;
    else
        valores_M = M_temporal;
        valores_dt = dt_temporal;
    end
    fprintf('\nRefinamiento %s, tf = 10 s\n',tipos_ref{tipo});
    fprintf('%-18s %-6s %-12s %-12s %-12s %-14s %-12s\n', ...
            'Metodo','M','dx','dt','mu','E_f','Tiempo(s)');

    % Medir solo el procedimiento numerico y calcular despues el error final.
    for metodo = 1:length(metodos_ref)
        for caso = 1:3
            M_ref = valores_M(caso);
            dt_ref = valores_dt(caso);
            dx_ref = L/M_ref;
            mu_ref = alpha*dt_ref/dx_ref^2;
            tic;
            if metodo == 1
                [x_ref,t_ref,U_ref] = metodos_ref{metodo}(L,alpha,u0,M_ref,dt_ref,10);
            else
                [x_ref,t_ref,U_ref,res_ref,conv_ref] = metodos_ref{metodo}(L,alpha,u0,M_ref,dt_ref,10);
            end
            tiempo_ref = toc;

            % Aceptar ambos esquemas implicitos solo si todos sus sistemas cumplen la tolerancia.
            if metodo > 1
                if ~isequal(conv_ref,1)
                    error('Fallo de convergencia en el refinamiento de %s: M = %d, dt = %.6f.', ...
                          nombres_ref{metodo},M_ref,dt_ref);
                end
                N_ref = round(10/dt_ref);
                if numel(res_ref) ~= N_ref || any(~isfinite(res_ref(:))) || ...
                   any(res_ref(:) < 0) || any(res_ref(:) >= 1e-10)
                    error('Residuos no validos en el refinamiento de %s: M = %d, dt = %.6f.', ...
                          nombres_ref{metodo},M_ref,dt_ref);
                end
                residuos_ref{tipo,caso,metodo} = res_ref(:);
                residuo_max_ref(tipo,caso,metodo) = max(res_ref(:));
            end
            if abs(t_ref(end)-10) > 1e-10
                error('El experimento de refinamiento no alcanzo tf = 10 s.');
            end
            E_f_ref = max(abs(U_ref(end,:) + 20 - T_exacta(x_ref,10)));
            resultados_ref(tipo,caso,metodo,:) = reshape([dx_ref,dt_ref,mu_ref,E_f_ref,tiempo_ref],1,1,1,5);
            fprintf('%-18s %-6d %-12.6f %-12.6f %-12.6f %-14.6e %-12.6f\n', ...
                    nombres_ref{metodo},M_ref,dx_ref,dt_ref,mu_ref,E_f_ref,tiempo_ref);
        end
    end

    % Reservar espacio para la leyenda y dejar margenes alrededor de los puntos extremos.
    figure('position',[80,80,1100,700],'paperunits','inches', ...
       'paperposition',[0,0,11,7],'paperpositionmode','manual');
    eje_ref = axes('position',[0.10,0.25,0.85,0.64]);
    for metodo = 1:length(metodos_ref)
        pasos = squeeze(resultados_ref(tipo,:,metodo,tipo));
        errores = squeeze(resultados_ref(tipo,:,metodo,4));
        loglog(pasos,errores,'o-','LineWidth',1.5,'MarkerSize',6);
        hold on;
    end
    if tipo == 1
        xlabel('dx (m)');
    else
        xlabel('dt (s)');
    end
    ylabel('E_f (C)');
    title(sprintf('Refinamiento %s',tipos_ref{tipo}));
    grid on;
    set(eje_ref,'xminorgrid','off','yminorgrid','off','xtick',sort(pasos));
    xlim([min(pasos)/1.15,max(pasos)*1.15]);
    leyenda_ref = legend(nombres_ref,'location','southoutside', ...
        'orientation','horizontal','fontsize',10,'box','off');
    set(leyenda_ref,'units','normalized','position',[0.14,0.035,0.72,0.07]);
    set(eje_ref,'position',[0.10,0.25,0.85,0.64]);
    print(sprintf('refinamiento_%s.png',tipos_ref{tipo}),'-dpng','-r150');
end

% Explicar las propiedades que justifican Thomas y las tendencias del refinamiento.
fprintf('\nAnalisis de Persona 2\n');
fprintf('Euler implicito: diagonal 1+2*mu y secundarias -mu.\n');
fprintf('Crank-Nicolson: diagonal 1+mu y secundarias -mu/2.\n');
fprintf('Para mu >= 0, ambas matrices son tridiagonales, simetricas y definidas positivas.\n');
fprintf('La dominancia diagonal estricta garantiza pivotes no nulos para Thomas sin pivoteo.\n');
fprintf('La factorizacion y las sustituciones de Thomas requieren O(M) operaciones; no tiene iteraciones.\n');
fprintf('La verificacion del residuo usa A densa y cuesta O(M^2) por paso; ese costo se incluye en los tiempos.\n');
fprintf('Euler implicito resuelve (I+mu*K)*u_nuevo = u_anterior y es estable sin exigir mu <= 0.5.\n');
fprintf('Euler explicito e implicito tienen error O(dt+dx^2); Crank-Nicolson, O(dt^2+dx^2).\n');
fprintf('Al refinar solo una variable, el error de la otra puede dominar.\n');
fprintf('En Euler explicito tambien puede haber compensacion entre errores espacial y temporal.\n');
fprintf('Por eso E_f puede no disminuir monotonamente; estabilidad no equivale a precision.\n');
fprintf('Las mallas de refinamiento explicito cumplen mu <= 0.5; el mayor mu es 0.4.\n');
fprintf('Los tiempos incluyen la construccion y factorizacion; las graficas y la solucion exacta se excluyen.\n');


% Preparar un aporte al informe con los resultados medidos, sin completar tareas de Persona 3.
analisis_ref = sprintf('# Refinamiento de malla: aporte de Persona 2\n\n');
if length(metodos_ref) == 3
    analisis_ref = [analisis_ref, sprintf('Se ejecutaron los tres esquemas con tf = 10 s.\n\n')];
else
    analisis_ref = [analisis_ref, sprintf(['Comparacion parcial: falta la funcion crank_nicolson de Persona 3. ', ...
        'Sus resultados e interpretacion se incorporaran al ejecutar nuevamente parte2.m.\n\n'])];
end

% Registrar las seis mallas de cada esquema y sus controles para que el informe sea verificable.
for tipo = 1:2
    analisis_ref = [analisis_ref, sprintf('## Refinamiento %s\n\n',tipos_ref{tipo})];
    analisis_ref = [analisis_ref, sprintf('| Esquema | dx (m) | dt (s) | mu | E_f (C) | Tiempo (s) | Residuo maximo |\n')];
    analisis_ref = [analisis_ref, sprintf('| --- | --- | --- | --- | --- | --- | --- |\n')];
    for metodo = 1:length(metodos_ref)
        for caso = 1:3
            datos = reshape(resultados_ref(tipo,caso,metodo,:),1,5);
            if metodo == 1
                residuo_texto = 'no aplica';
            else
                residuo_texto = sprintf('%.6e',residuo_max_ref(tipo,caso,metodo));
            end
            analisis_ref = [analisis_ref, sprintf('| %s | %.6f | %.6f | %.6f | %.6e | %.6f | %s |\n', ...
                nombres_ref{metodo},datos(1),datos(2),datos(3),datos(4),datos(5),residuo_texto)];
        end
    end
    analisis_ref = [analisis_ref, sprintf('\n')];

    % Describir los errores reales y sus cocientes sin atribuirles un orden puro de convergencia.
    for metodo = 1:length(metodos_ref)
        errores = reshape(resultados_ref(tipo,:,metodo,4),1,3);
        razones = errores(1:2)./errores(2:3);
        analisis_ref = [analisis_ref, sprintf(['%s: los errores de las tres mallas, en el orden de la tabla, ', ...
            'son %.6e, %.6e y %.6e C; los cocientes E_anterior/E_siguiente son %.3f y %.3f. '], ...
            nombres_ref{metodo},errores(1),errores(2),errores(3),razones(1),razones(2))];
        if all(diff(errores) < 0)
            analisis_ref = [analisis_ref, sprintf('El error disminuye en los dos refinamientos. ')];
        elseif any(diff(errores) > 0)
            analisis_ref = [analisis_ref, sprintf('El error no disminuye en ambos refinamientos. ')];
        else
            analisis_ref = [analisis_ref, sprintf('Al menos un refinamiento mantiene el mismo error. ')];
        end

        % Relacionar cada tendencia con la contribucion que permanece fija en el experimento.
        if tipo == 1
            analisis_ref = [analisis_ref, sprintf(['Se reduce dx a la mitad manteniendo dt = 0.0025 s. ', ...
                'El error espacial es de segundo orden; un cociente cercano a cuatro es compatible ', ...
                'con su predominio, pero el error temporal fijo y su posible compensacion ', ...
                'impiden exigir ese cociente al error total.\n\n'])];
        elseif metodo == 1
            analisis_ref = [analisis_ref, sprintf(['Con M = 50 fijo, el error espacial permanece. ', ...
                'Para el modo senoidal, el operador espacial discreto amortigua menos que el continuo ', ...
                'y Euler explicito introduce mayor amortiguamiento temporal. Ambos errores pueden ', ...
                'compensarse; al reducir dt esa compensacion cambia, lo que explica una posible ', ...
                'subida del error aun cumpliendo mu <= 0.5.\n\n'])];
        elseif metodo == 2
            analisis_ref = [analisis_ref, sprintf(['Con M = 50 fijo, Euler implicito tiene error temporal ', ...
                'de primer orden. Para este dato senoidal, tanto el error espacial como el temporal ', ...
                'reducen el enfriamiento respecto de la solucion exacta. Al reducir dt disminuye ', ...
                'la contribucion temporal, pero permanece el error espacial; por eso los cocientes ', ...
                'del error total no tienen que ser dos.\n\n'])];
        else
            analisis_ref = [analisis_ref, sprintf(['Con M = 50 fijo, Crank-Nicolson tiene error temporal ', ...
                'de segundo orden. Un cociente cercano a cuatro solo se espera si domina ese error; ', ...
                'si los errores cambian poco al reducir dt, el error espacial fijo limita la mejora. ', ...
                'La compensacion de contribuciones tambien puede alterar los cocientes.\n\n'])];
        end
    end
end

% Resumir controles y costos usando los valores obtenidos en esta misma ejecucion.
for metodo = 2:length(metodos_ref)
    maximos = residuo_max_ref(:,:,metodo);
    analisis_ref = [analisis_ref, sprintf(['%s: todos los sistemas de las seis mallas cumplieron ', ...
        'el criterio; el mayor residuo normalizado fue %.6e, inferior a 1e-10.\n\n'], ...
        nombres_ref{metodo},max(maximos(:)))];
end
analisis_ref = [analisis_ref, sprintf(['Los tiempos de las tablas incluyen la construccion de matrices ', ...
    'y el avance numerico; excluyen errores, tablas y graficas. Euler implicito reutiliza la factorizacion ', ...
    'de Thomas y no tiene iteraciones internas. Al refinar el tiempo aumentan los pasos ', ...
    '(625, 1250 y 2500); al refinar el espacio aumentan las incognitas interiores ', ...
    '(24, 49 y 99). Las mediciones pueden fluctuar y no constituyen por si solas ', ...
    'una estimacion del orden de complejidad. ', ...
    'La factorizacion y las sustituciones de Thomas cuestan O(M), pero la verificacion ', ...
    'del residuo usa una matriz A densa y requiere O(M^2) por paso; ese costo tambien ', ...
    'esta incluido en los tiempos.\n'])];

% Mostrar el aporte y guardarlo como texto editable para que Persona 3 integre su informe.
fprintf('\n%s',analisis_ref);
archivo_ref = fopen('aporte_refinamiento_persona2.md','w');
if archivo_ref == -1
    error('No se pudo guardar el aporte de refinamiento de Persona 2.');
end
fprintf(archivo_ref,'%s',analisis_ref);
fclose(archivo_ref);



% ===== experimento final=====
tf_f = 30;
umbral = 30;

tic;
[x_f,t_f,U_f,res_f,conv_f] = crank_nicolson(L,alpha,u0,M,dt,tf_f);
tiempo_f = toc;
if ~conv_f
    error('Se detuvo el experimento final de Crank-Nicolson.');
end
T_f = U_f + 20;

% Punto 10
Tmax_f = max(T_f,[],2);
n_cruce = find(Tmax_f <= umbral,1);
if isempty(n_cruce) || n_cruce == 1
    error('El umbral no se alcanza dentro de (0, tf] o ya se cumple en t = 0.');
end
t_num = t_f(n_cruce);
Tmax_cruce = Tmax_f(n_cruce);
t_previo = t_f(n_cruce-1);
Tmax_previo = Tmax_f(n_cruce-1);
[~,i_max] = max(T_f(n_cruce,:));

% Punto 11: separar el redondeo exacto del desplazamiento entre cruces discretos.
t_star = log(80/(umbral-20))/(alpha*(pi/L)^2);
n_exacto = ceil(t_star/dt - 1e-9);
t_malla_exacto = n_exacto*dt;
error_redondeo = t_malla_exacto - t_star;
% Esta diferencia incluye el efecto del esquema sobre el cruce y su discretizacion temporal.
diferencia_cruces_discretos = t_num - t_malla_exacto;
error_total = t_num - t_star;

fprintf('\n=== Experimento final (Crank-Nicolson, M = %d, dt = %.3f s, tf = %d s) ===\n',M,dt,tf_f);
fprintf('Residuo normalizado maximo: %.3e, tiempo: %.4f s\n',max(res_f),tiempo_f);
fprintf('Primer tiempo con max T <= %d C: t = %.2f s (n = %d), max T = %.6f C, en x = %.2f m\n', ...
        umbral,t_num,n_cruce-1,Tmax_cruce,x_f(i_max));
fprintf('Paso anterior: t = %.2f s, max T = %.6f C (> %d C)\n',t_previo,Tmax_previo,umbral);
fprintf('Instante exacto t* = ln(8)/(alpha*(pi/L)^2) = %.6f s\n',t_star);
fprintf('Primer tiempo de la malla con la solucion exacta: %.2f s\n',t_malla_exacto);
fprintf('Numerico - exacto: %.6f s = redondeo del cruce exacto (%.6f s) + diferencia entre cruces discretos (%.6f s)\n', ...
        error_total,error_redondeo,diferencia_cruces_discretos);
fprintf('La diferencia entre cruces discretos no es el error puro del esquema separado del redondeo.\n');

% Punto 12: mapa de calor
figure('position',[80,80,1000,700],'paperunits','inches', ...
       'paperposition',[0,0,10,7],'paperpositionmode','manual');
imagesc(x_f,t_f,T_f);
set(gca,'YDir','normal');
xlabel('Posicion x (m)');
ylabel('Tiempo t (s)');
title('Mapa de calor de T(x,t), Crank-Nicolson');
barra = colorbar();
try
    ylabel(barra,'Temperatura (C)');
catch
    title(barra,'T (C)');
end
hold on;
plot([0,L],[t_num,t_num],'w--','LineWidth',1.5);
text(0.03*L,t_num+0.8,sprintf('t = %.2f s: max T <= 30 C',t_num),'color','w','fontsize',10);
print('mapa_calor_enfriamiento.png','-dpng','-r150');

% Punto 12: temperatura en el centro
ic = round(M/2) + 1;                           % x = L/2 (M par)
figure('position',[80,80,1100,700],'paperunits','inches', ...
       'paperposition',[0,0,11,7],'paperpositionmode','manual');
curva_c = plot(t_f,T_f(:,ic),'-','LineWidth',1.6);
hold on;
curva_ex = plot(t_f,T_exacta(x_f(ic),t_f),'k--','LineWidth',1.2);
curva_u = plot([0,tf_f],[umbral,umbral],'r:','LineWidth',1.5);
curva_t = plot(t_num,Tmax_cruce,'ro','MarkerSize',7,'MarkerFaceColor','r');
xlabel('Tiempo t (s)');
ylabel('T(L/2,t) (C)');
title('Temperatura en el centro de la barra');
grid on;
legend([curva_c,curva_ex,curva_u,curva_t], ...
    {'Crank-Nicolson','Solucion exacta','Umbral 30 C', ...
     sprintf('Primer t con max T <= 30 C (%.2f s)',t_num)}, ...
    'location','northeast','fontsize',10);
drawnow;
print('temperatura_centro.png','-dpng','-r150');

% ===== Resumen final en la ventana de comandos =====
fprintf('\n=== Tabla comparativa final (M = %d, dt = %.2f s, tf = %d s) ===\n',M,dt,tf);
fprintf('%-18s %-14s %-14s %-12s %-14s %-14s %-14s\n','Metodo','E_f','E_malla','Tiempo(s)','Residuo max','Iter. acum.','Iter. max');
fprintf('%-18s %-14.6e %-14.6e %-12.6f %-14s %-14s %-14s\n','Euler explicito',E_f_exp,E_malla_exp,tiempo_exp,'no aplica','no aplica','no aplica');
fprintf('%-18s %-14.6e %-14.6e %-12.6f %-14.6e %-14s %-14s\n','Euler implicito',E_f_imp,E_malla_imp,tiempo_imp,residuo_max_imp,'no aplica','no aplica');
fprintf('%-18s %-14.6e %-14.6e %-12.6f %-14.6e %-14s %-14s\n','Crank-Nicolson',E_f_cn,E_malla_cn,tiempo_cn,residuo_max_cn,'no aplica','no aplica');

fprintf('\n=== Resumen del analisis ===\n');
fprintf('1. Precision: con M = 50 y mu = %.2f domina el error espacial O(dx^2). En Euler explicito el error\n',mu);
fprintf('   temporal tiene signo opuesto y compensa parte del espacial: E_f = %.2e (EE), %.2e (CN), %.2e (EI).\n',E_f_exp,E_f_cn,E_f_imp);
fprintf('2. Costo: Euler explicito no resuelve sistemas; los implicitos reutilizan la factorizacion de Thomas.\n');
fprintf('   Thomas requiere O(M) por solucion; verificar el residuo con A densa requiere O(M^2) por paso.\n');
fprintf('3. Estabilidad: Euler explicito exige mu <= 1/2 (con mu = 0.6 crece la perturbacion); EI y CN son estables\n');
fprintf('   para todo mu, pero estabilidad no implica precision.\n');
fprintf('4. Residuos: maximo %.2e (EI) y %.2e (CN), muy por debajo de 1e-10; Thomas es directo (sin iteraciones).\n',residuo_max_imp,residuo_max_cn);
fprintf('5. Enfriamiento: la temperatura maxima esta en x = L/2 y decae como exp(-alpha*(pi/L)^2*t); max T <= 30 C en\n');
fprintf('   t = %.2f s (exacto %.4f s; diferencia total %.4f s).\n',t_num,t_star,error_total);
fprintf('   Redondeo del cruce exacto: %.6f s; diferencia entre cruces discretos: %.6f s.\n', ...
        error_redondeo,diferencia_cruces_discretos);
fprintf('6. Fisica: el aislamiento lateral deja salir el calor solo por los extremos, fijos a 20 C (Dirichlet);\n');
fprintf('   la barra tiende a 20 C uniforme y el perfil conserva la forma senoidal mientras decae.\n');
