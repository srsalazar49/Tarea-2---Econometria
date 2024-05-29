%% ===== Tarea 2 - Grupo =======%
clear, clc;
cd('C:\Users\HP\Documents\MAESTRÍA\1. Primer Semestre\Econometría I\Tareas\Tarea 2')
%Cargamos la base de datos
[data_num, txt_data, raw_data] = xlsread('WAGEPAN.xls');
data = readtable('WAGEPAN.xls', 'ReadVariableNames', true);
v_ind = {'educ', 'black', 'hisp', 'exper', 'expersq','married','union'};
cols = data(:, v_ind);
X = table2array(cols);
v_exp= {'lwage'};
cols_exp= data (:,v_exp);
Y = table2array(cols_exp);
periods=8;
%% ------- 1. Estimación por POLS ---------------
[beta, e]=ols(Y,X);
n=size(X,1);
[K,s2] = s2(n, beta, e);
% Errores estándar homocedasticos y sin autocorrelación
[var_beta, ee_estardar]=errores_estandar(s2,X);
%Errores estándar robustos
[var_robust, ee_robust] = errores_robustos(n, K ,X, e);
%Errores estándar agrupados
id=data.nr;
id_num=unique(id);
id_2=repelem(1:length(unique(id_num)),8)'; %Le asignamos un valor
groups=length(id_num);
e_cluster=zeros(groups,K);%Matriz de cluster
for j = 1:K
    e_cluster(:,j) = accumarray(id_2, (X(:,j))'.*e');
end
[var_cluster, ee_cluster] = errores_cluster(n, K, groups, X, e_cluster);
%% ------ 2. Estimador Between -------------------
y_mean=accumarray(id_2,Y,[],@mean);
for i = 1:size(X,2)
    x_mean(:,i) = accumarray(id_2, X(:,i), [],@mean);
end

[beta_bet, e_bet]=ols(y_mean,x_mean);
n_bet=size(x_mean,1);
K_bet=length(beta_bet);
s2_bet=(e_bet'*e_bet)/(n_bet-K_bet);
% Errores estándar homocedasticos y sin autocorrelación
[var_beta_bet, ee_estardar_bet]=errores_estandar(s2_bet,x_mean);
%Errores estándar robustos
[var_robust_bet, ee_robust_bet] = errores_robustos(n_bet, K_bet ,x_mean, e_bet);
%Errores estándar agrupados (revisar)
[var_cluster_bet, ee_cluster_bet] = errores_cluster(n_bet, K_bet, groups, x_mean, e_bet);
%% --------- 3. Estimador Within ----------------
y_with=Y-y_mean(id_2);
for i=1:size(X,2)
    x_with(:,i)=X(:,i)-x_mean(id_2,i);
end
x_with=x_with(:,4:end);
[beta_with, e_with]=ols(y_with,x_with);
n_with=size(x_with,1);
K_with=length(beta_with);
s2_with=(e_with'*e_with)/(n_with-K_with);
% Errores estándar homocedasticos y sin autocorrelación
[var_beta_with, ee_estardar_with]=errores_estandar(s2_with,x_with);
%Errores estándar robustos
[var_robust_with, ee_robust_with] = errores_robustos(n_with, K_with ,x_with, e_with);
%Errores estándar agrupados (revisar)
e_cluster_with=zeros(groups,K_with);%Matriz de cluster
for j = 1:K_with
    e_cluster_with(:,j) = accumarray(id_2, (x_with(:,j))'.*e_with');
end
[var_cluster_with, ee_cluster_with] = errores_cluster(n_with, K_with, groups, x_with, e_cluster_with);
%% ---------- 4. Estimación por efectos aleatorios ----------
sigma_eps2=(1/(n-groups-K_with)).*sum(e_with.^2);
sigma_u2=(1/(groups-K)).*(sum(e_bet.^2))-(1/periods)*sigma_eps2;
dummy_random=dummyvar(id_2);
I=eye(n,n);
omega_eps=dummy_random*dummy_random'.*sigma_u2+I.*sigma_eps2;
omega=omega_eps./sigma_eps2;
[beta_random, e_random] = fgls(Y,X,omega);

n_random=size(X,1);
K_random=length(beta_random);
s2_random=(e_random'*e_random)/(n_random-K_random);

% Errores estándar homocedasticos y sin autocorrelación
[var_beta_random, ee_estardar_random]=errores_estandar(s2_random,X);
%Errores estándar robustos
[var_robust_random, ee_robust_random] = errores_robustos(n_random, K_random ,X, e_random);
%Errores estándar agrupados (revisar)
id=data.nr;
id_num=unique(id);
id_2=repelem(1:length(unique(id_num)),8)'; %Le asignamos un valor
groups=length(id_num);
e_cluster_random=zeros(groups,K_random);%Matriz de cluster
for j = 1:K_random
    e_cluster_random(:,j) = accumarray(id_2, (X(:,j))'.*e_random');
end
[var_cluster_random, ee_cluster_random] = errores_cluster(n_random, K_random, groups, X, e_cluster_random);

%% ----------- 5. Intervalos de confianza ------------------
alpha = 0.05;                                   % Nivel de significancia
z = norminv(1 - alpha/2);                       % Valor crítico
B = 2000;                                       % Núme remuestreo

beta_union = beta(7);                              % Premio por estar sindicalizado
ee_union = ee_estardar(7);                           % Error estándar variable sindicato

% Teoría asíntótica
conf_int_asy = [beta_union - z*ee_union, beta_union + z*ee_union];     % Intervalo

% Block bootstrap
block_boots = zeros(B, 1);                           % Vector para almacenar betas
for i = 1:B                                     % Remuestreo 2000 veces
    ind_remuestreo = [];
    for j = 1:groups                                % índices para remuestreo por bloque
        start = randi([1, n - periods + 1]);          % índice de inicio  
        block = start:(start + periods - 1);          % Índice para cada obs del bloque
        ind_remuestreo = [ind_remuestreo; block'];  % Guanda en vector de índices general
    end
    ind_remuestreo = ind_remuestreo(1:n);       % Asegurar tamaño correcto
    x_block = X(ind_remuestreo, :);             % Variable resampleada según los índices
    y_block = Y(ind_remuestreo);
    b_block = inv(x_block'*x_block)*(x_block'*y_block);   % Estimación para la nueva muestra
    block_boots(i) = b_block(7);                               % Guardo el parámetro de interés
end

conf_iint_boots = prctile(block_boots, [100*alpha/2, 100*(1-alpha/2)]);  % Intervalo sgn percentiles

% Wild cluster bootstrap
wild_boots = zeros(B, 1);                        % Vector para almacenar betas
for i = 1:B                                     % Remuestreo 2000 veces  
    y_wild = Y;                                 % Parto desde el Y original                           
    for j = 1:groups                                % Resampleo para cada cluster
        cluster_ind = find(id_2 == j);            % Busco obs del cluster
        resample_factor = 2 * (rand() > 0.5) - 1;  % Factor de remuestreo aleatoreo
        y_wild(cluster_ind) = Y(cluster_ind) * resample_factor;   % Nueva muestra = Variabilidad a los datos
    end
    b_wild = inv(X'*X)*(X'*y_wild);             % Estimación para la muestra 
    wild_boots(i) = b_wild(7);                   % Guardo parámetro de interés
end

conf_int_wboots= prctile(wild_boots, [100*alpha/2, 100*(1-alpha/2)]);      % Intervalo
%% ----------- 6. Regresión por POLS -----------------------
%% 
