%% LABORATORY OF NATURAL RESOURCES MANAGEMENT - Project
clear
clc
close all

% ----------------------
%%  DATA PREPROCESSING
% ----------------------

load -ascii tudela.txt 

% calibration and validation sets
anno = tudela(1:365*20,1);
mese = tudela(1:365*20,2);
giorno = tudela(1:365*20,3);
p = tudela(1:365*20,4);
q = tudela(1:365*20,5);
temp = tudela(1:365*20,6);

anno_v = tudela((365*20+1):end,1);
mese_v = tudela((365*20+1):end,2);
giorno_v = tudela((365*20+1):end,3);
p_v = tudela((365*20+1):end,4);
q_v = tudela((365*20+1):end,5);
temp_v = tudela((365*20+1):end,6);

% normalization

% CALCOLIAMO LA MEDIA MOBILE PER L'AFFLUSSO
% mi_q = media mobile periodica 
% m_q = media mobile periodica per l'intero periodo(20 anni)
T = 365;
[mi_q, m_q] = moving_average(q, T, 7);  
[sigma2_q, s2_q] = moving_average( (q-m_q).^2, T, 7);
s_q = sqrt(s2_q);

% CALCOLIAMO LA MEDIA MOBILE PER LA PRECIPITAZIONE(ESOGENA)
[mi_p, m_p] = moving_average(p, T, 7);
[sigma2_p, s2_p] = moving_average( (p-m_p).^2, T, 7);
s_p = sqrt(s2_p);

% CALCOLIAMO LA MEDIA MOBILE PER LA TEMPERATURA(ESOGENA)
[mi_temp, m_temp] = moving_average(temp, T, 7);
[sigma2_temp, s2_temp] = moving_average( (temp-m_temp).^2, T, 7);
s_temp = sqrt(s2_temp);

x = (q-m_q) ./ s_q;
u = (p-m_p) ./ s_p;
t = (temp-m_temp) ./ s_temp;



[mi_qv, m_qv] = moving_average(q_v, T, 7);
[sigma2_qv, s2_qv] = moving_average( (q_v-m_qv).^2, T, 7);
s_qv = sqrt(s2_qv);

[mi_pv, m_pv] = moving_average(p_v, T, 7);
[sigma2_pv, s2_pv] = moving_average( (p_v-m_pv).^2, T, 7);
s_pv = sqrt(s2_pv);

[mi_tempv, m_tempv] = moving_average(temp_v, T, 7);
[sigma2_tempv, s2_tempv] = moving_average( (temp_v-m_tempv).^2, T, 7);
s_tempv = sqrt(s2_tempv);

xv = (q_v-m_qv) ./ s_qv;
uv = (p_v-m_pv) ./ s_pv;
tv = (temp_v-m_tempv) ./ s_tempv;



figure
plot(x)
xlabel('time t')
ylabel('deseasonalized inflow x_t [m^3/s]')

figure
plot(tt, x, '.')
xlabel('time t (1 year)')
ylabel('deseasonalized inflow x_t [m^3/s]')
% COMMENT : periodicity has been (partially) removed



%  -----------------
%%  AUTOCORRELATION 
%  -----------------

figure
correlogram(t, t, 20);
xlabel('k')
ylabel('r_k')
% COMMENT : it is not white and we can notice some peaks
% (with period 7) that are due to upstream hydropower
% releases (this deterministic behavior was not removed).



% ---------------------
%%  AR(1)
% ---------

% calibration
y = x(2:end);
M = x(1:end-1);
theta = M \ y;
x_ = [ x(1) ; M*theta ];
q_ar1 = x_ .* s_q + m_q;

R2_ar1 = 1 - sum((q(2:end)-q_ar1(2:end)).^2) / ...
    sum((q(2:end)-m_q(2:end)).^2);


%% 
    
Mv = xv(1:end-1);
xv_ = [xv(1) ; Mv*theta];
qv_ar1 = xv_ .* s_qv + m_qv;
ev = q_v(2:end) - qv_ar1(2:end);
figure;
correlogram(ev, ev, 20);
Qv = mean( ev.^2);
R2_ar1v = 1 - sum((q_v(2:end)-qv_ar1(2:end)).^2) / ...
    sum((q_v(2:end)-m_qv(2:end)).^2);

% 
%r2_ar = 0.929538641664889
%r2_arv = 0.915993400749009

% ---------------------
%%  AR(2)
% ---------

% calibration
y2 = x(3:end);
M2 = [x(1:end-2), x(2:end-1)];
theta2 = M2 \ y2;
x2_ = [ x(1:2) ; M2*theta2 ];
q_ar2 = x2_ .* s_q + m_q;

R2_ar2 = 1 - sum((q(3:end)-q_ar2(3:end)).^2) / ...
    sum((q(3:end)-m_q(3:end)).^2);


%% 
    
M2v = [xv(1:end-2), xv(2:end-1)];
xv2_ = [xv(1:2) ; M2v*theta2];
qv_ar2 = xv2_ .* s_qv + m_qv;
ev2 = q_v(3:end) - qv_ar2(3:end);
figure;
correlogram(ev2, ev2, 20);
Qv2 = mean( ev2.^2);
R2_ar2v = 1 - sum((q_v(3:end)-qv_ar2(3:end)).^2) / ...
    sum((q_v(3:end)-m_qv(3:end)).^2);

%%  AR(3)
% ---------

% calibration
y3 = x(4:end);
M3 = [x(1:end-3), x(2:end-2), x(3:end-1)];
theta3 = M3 \ y3;
x3_ = [ x(1:3) ; M3*theta3 ];
q_ar3 = x3_ .* s_q + m_q;

R2_ar3 = 1 - sum((q(4:end)-q_ar3(4:end)).^2) / ...
    sum((q(4:end)-m_q(4:end)).^2);


%% 
    
M3v = [xv(1:end-3), xv(2:end-2), xv(3:end-1)];
xv3_ = [xv(1:3) ; M3v*theta3];
qv_ar3 = xv3_ .* s_qv + m_qv;
ev3 = q_v(4:end) - qv_ar3(4:end);
figure;
correlogram(ev3, ev3, 20);
Qv3 = mean( ev3.^2);
R2_ar3v = 1 - sum((q_v(4:end)-qv_ar3(4:end)).^2) / ...
    sum((q_v(4:end)-m_qv(4:end)).^2);

%%  AR(4)
% ---------

% calibration
y4 = x(5:end);
M4 = [x(1:end-4), x(2:end-3), x(3:end-2), x(4:end-1)];
theta4 = M4 \ y4;
x4_ = [ x(1:4) ; M4*theta4 ];
q_ar4 = x4_ .* s_q + m_q;
R2_ar4 = 1 - sum((q(4:end)-q_ar4(4:end)).^2) / ...
    sum((q(4:end)-m_q(4:end)).^2);


%% 
    
M4v = [xv(1:end-4), xv(2:end-3), xv(3:end-2), xv(4:end-1)];
xv4_ = [xv(1:4) ; M4v*theta4];
qv_ar4 = xv4_ .* s_qv + m_qv;
ev4 = q_v(5:end) - qv_ar4(5:end);
figure;
correlogram(ev4, ev4, 20);
Qv4 = mean( ev4.^2);
R2_ar4v = 1 - sum((q_v(4:end)-qv_ar4(4:end)).^2) / ...
    sum((q_v(4:end)-m_qv(4:end)).^2);

% -------------------
%%  ARX(1,1) PROPER
% -------------------

% calibration
y = x(2:end);
M = [ x(1:end-1) u(1:end-1) ];  
theta = M \ y;
x_arx_pro_ = [ x(1) ; M*theta ];
q_arx_pro_prec = x_arx_pro_ .* s_q + m_q;

R2_arx_pro = 1 - sum((q(2:end)-q_arx_pro_prec(2:end)).^2) / ...
    sum((q(2:end)-m_q(2:end)).^2);

% validation
yv = xv(2:end);
Mv = [ xv(1:end-1) uv(1:end-1) ];
x_arx_prov = [ x(1) ; Mv*theta ];
q_arx_prov = x_arx_prov .* s_qv + m_qv;
ev = q_v(2:end) - q_arx_prov(2:end);
figure;
correlogram(ev, ev, 20)

R2_arx_prov = 1 - sum((q_v(2:end)-q_arx_prov(2:end)).^2) / ...
    sum((q_v(2:end)-m_qv(2:end)).^2);

% TEMPERATURA

% calibration
y = x(2:end);
M = [ x(1:end-1) t(1:end-1) ];  
theta = M \ y;
x_ = [ x(1) ; M*theta ];
q_arx_pro_temp = x_ .* s_q + m_q;

R2_arx_pro_temp = 1 - sum((q(2:end)-q_arx_pro_temp(2:end)).^2) / ...
    sum((q(2:end)-m_q(2:end)).^2);

% validation
yv = xv(2:end);
Mv = [ xv(1:end-1) tv(1:end-1) ];
x_v = [ x(1) ; Mv*theta ];
q_arx_prov_temp = x_v .* s_qv + m_qv;
ev = q_v(2:end) - q_arx_prov_temp(2:end);
figure;
correlogram(ev, ev, 20)

R2_arx_prov_temp = 1 - sum((q_v(2:end)-q_arx_prov_temp(2:end)).^2) / ...
    sum((q_v(2:end)-m_qv(2:end)).^2);

% con la temperatura viene peggio r2 quindi si attacca al cazzo :-)

% ---------------------
%%  ARX(1,1) IMPROPER
% ---------------------

% calibration
y = x(2:end);
M = [ x(1:end-1) u(2:end) ];
theta = M \ y;
x_ = [ x(1); M*theta ];
q_arx_imp = x_ .* s_q + m_q;

R2_arx_imp = 1 - sum((q(2:end)-q_arx_imp(2:end)).^2) / ...
    sum((q(2:end)-m_q(2:end)).^2);

% validation
yv = xv(2:end);
Mv = [ xv(1:end-1) uv(2:end) ];
x_v = [ x(1) ; Mv*theta ];
q_arx_impv = x_v .* s_qv + m_qv;
ev = q_v(2:end) - q_arx_impv(2:end);
figure;
correlogram(ev, ev, 20)

R2_arx_impv = 1 - sum((q_v(2:end)-q_arx_impv(2:end)).^2) / ...
    sum((q_v(2:end)-m_qv(2:end)).^2);

% TEMP

% calibration
y = x(2:end);
M = [ x(1:end-1) t(2:end) ];
theta = M \ y;
x_ = [ x(1); M*theta ];
q_arx_imp = x_ .* s_q + m_q;

R2_arx_imp_temp = 1 - sum((q(2:end)-q_arx_imp(2:end)).^2) / ...
    sum((q(2:end)-m_q(2:end)).^2);

% validation
yv = xv(2:end);
Mv = [ xv(1:end-1) tv(2:end) ];
x_v = [ x(1) ; Mv*theta ];
q_arx_impv = x_v .* s_qv + m_qv;
ev = q_v(2:end) - q_arx_impv(2:end);
figure;
correlogram(ev, ev, 20)

R2_arx_impv_temp = 1 - sum((q_v(2:end)-q_arx_impv(2:end)).^2) / ...
    sum((q_v(2:end)-m_qv(2:end)).^2);

% funziona peggio anche qua, spiaze

% ----------------------------
%%  ANN PROPER
% --------------

% training
X = [ x(1:end-1) u(1:end-1) ];
Y = x(2:end);
net = feedforwardnet(5);
% 1 hidden layer with 5 neurons - view(net)
net = train(net, X', Y');
Y_ = [ x(1) ; net(X')' ] ;
q_ann_pro = Y_ .* s_q + m_q ;

R2_ann_pro = 1 - sum((q(2:end)-q_ann_pro(2:end)).^2) / ...
    sum((q(2:end)-m_q(2:end)).^2);

% test
Xv = [ xv(1:end-1) uv(1:end-1) ];
Y_v = [ x(1) ; net(Xv')' ];
q_ann_prov = Y_v .* s_qv + m_qv;

ev = q_v(2:end) - q_ann_prov(2:end);
figure;
correlogram(ev, ev, 20)

R2_ann_prov = 1 - sum((q_v(2:end)-q_ann_prov(2:end)).^2) / ...
    sum((q_v(2:end)-m_qv(2:end)).^2);



% ----------------
%%  ANN IMPROPER
% ----------------

% training
X = [ x(1:end-1) u(2:end) ];
Y = x(2:end);
net = feedforwardnet(5);
% 1 hidden layer with 5 neurons - view(net)
net = train(net, X', Y');
Y_ = [ x(1) ; net(X')' ];
q_ann_imp = Y_ .* s_q + m_q;

R2_ann_imp = 1 - sum((q(2:end)-q_ann_imp(2:end)).^2) / ...
    sum((q(2:end)-m_q(2:end)).^2);

% test
Xv = [ xv(1:end-1) uv(2:end) ];
Y_v = [ x(1) ; net(Xv')' ];
q_ann_impv = Y_v .* s_qv + m_qv;

ev = q_v(2:end) - q_ann_impv(2:end);
figure;
correlogram(ev, ev, 20)
R2_ann_impv = 1 - sum((q_v(2:end)-q_ann_impv(2:end)).^2) / ...
    sum((q_v(2:end)-m_qv(2:end)).^2);



% ---------------------
%%  MODELS COMPARISON
% ---------------------

figure
bar([R2_ar1 R2_ar1v;
    R2_arx_pro R2_arx_prov;
    R2_arx_imp R2_arx_impv;
    R2_ann_pro R2_ann_prov;
    R2_ann_imp R2_ann_impv])
set(gca, 'XTickLabel',{'AR(1)', 'ARX(1,1) pro','ARX(1,1,) imp', ...
    'ANN pro', 'ANN imp'})
axis([0.5 5.5 0 1])



% -----------------------------
%%  MULTIPLE ANN CALIBRATIONS
% -----------------------------

X = [ x(1:end-1) u(2:end) ];
Y = x(2:end);
Xv = [ xv(1:end-1) uv(2:end) ];
Yv = xv(2:end);
N_runs = 30;
R2 = zeros(N_runs, 2); % col1 = train , col2 = test
for i = 1:N_runs
    net_i = feedforwardnet(5);
    net_i = train(net_i, X', Y');
    
    Y_ = [ x(1) ; net_i(X')' ];
    q_i = Y_ .* s_q + m_q ;
    R2(i, 1) = 1 - sum((q(2:end)-q_i(2:end)).^2) / ...
        sum((q(2:end)-m_q(2:end)).^2);
    
    Y_v = [ xv(1) ; net_i(Xv')' ];
    q_iv = Y_v .* s_qv + m_qv ;
    R2(i, 2) = 1 - sum((qv(2:end)-q_iv(2:end)).^2) / ...
        sum((qv(2:end)-m_qv(2:end)).^2);    
    
    % saving the best model (in test) as net_opt
    if R2(i, 2) >= max(R2(:, 2))
        net_opt = net_i ;
    end
end

figure
bar(R2)



% ----------------------------------
%%  SHALLOW & DEEP NEURAL NETWORKS
% ----------------------------------

clear
clc
close all

load -ascii 'friedman_dataset.txt'

X = friedman_dataset(1:2000, 1:15);
Y = friedman_dataset(1:2000, 16);
Xv = friedman_dataset(2001:2500, 1:15);
Yv = friedman_dataset(2001:2500, 16);

tic

N_runs = 30;
R2 = zeros(N_runs, 2); % col1 = train , col2 = test
for i = 1:N_runs
    net = feedforwardnet(20);
    net = train(net, X', Y');
    Y_ = net(X')';
    Y_v = net(Xv')';
    
    R2(i, 1) = (1-sum((Y-Y_).^2)/sum((Y-mean(Y)).^2));
    R2(i, 2) = (1-sum((Yv-Y_v).^2)/sum((Yv-mean(Yv)).^2));
end
sum(R2(:,2) >= 0.95)

N_runs = 30;
R2 = zeros(N_runs, 2); % col1 = train , col2 = test
for i = 1:N_runs
    net = feedforwardnet([5 5 5 5]);
    net = train(net, X', Y');
    Y_ = net(X')';
    Y_v = net(Xv')';
    
    R2(i, 1) = (1-sum((Y-Y_).^2)/sum((Y-mean(Y)).^2));
    R2(i, 2) = (1-sum((Yv-Y_v).^2)/sum((Yv-mean(Yv)).^2));
end
sum(R2(:,2) >= 0.95)

toc

% doc trainlm
%%  CART
% --------

X = [x(1:end-1), u(1:end-1)];
Y = x(2:end);
% training and predicting
T = fitrtree(X, Y, 'MinLeafSize', 10);
Y_pro = [ x(1); predict(T, X) ];
q_cart_pro = Y_pro .* s_q + m_q ;
R2_cart_pro = 1 - sum((q(2:end)-q_cart_pro(2:end)).^2) / ...
    sum((q(2:end)-m_q(2:end)).^2);

% validation
Xv = [xv(1:end-1), uv(1:end-1)];
Yv = xv(2:end);
Yv_pro = [ xv(1); predict(T, Xv) ];
qv_cart_pro = Yv_pro .* s_qv + m_qv;
R2_cart_prov = 1 - sum((q_v(2:end)-qv_cart_pro(2:end)).^2) / ...
    sum((q_v(2:end)-m_qv(2:end)).^2)



%%
% 'MaxNumSplits'
% 'MinLeafSize'
% 'MaxDepth' only in the latest matlab versions

% training
% MinLeafSize (default = 1) - 10; 100; 300
T1 = fitrtree(X, Y, 'MinLeafSize', 10);
view(T1, 'mode', 'graph')

Y_pro = [ x(1); predict(T1, X) ];
q_cart_pro = Y_pro .* s_q + m_q;
R2_cart_pro = 1 - sum((q(2:end)-q_cart_pro(2:end)).^2) / ...
    sum((q(2:end)-m_q(2:end)).^2)

% validation
Xv = [xv(1:end-1), uv(1:end-1)];
Yv = xv(2:end);
Y_ = predict(T1, Xv);
Yv_pro = [ xv(1); Y_ ] ;
qv_cart_pro = Yv_pro .* s_qv + m_qv ;
R2_cart_prov = 1 - sum((q_v(2:end)-qv_cart_pro(2:end)).^2) / ...
    sum((q_v(2:end)-m_qv(2:end)).^2)

% note that length(unique(Y_)) is equal to the number of leaves

% optimization of minLeafSize
% Topt = fitrtree(X, Y, 'OptimizeHyperparameters', 'auto')

% test pruning
% Tpruned = prune(T, 'Alpha', 0.2);



% ------------------------
%%  RANDOM FOREST PROPER
% ------------------------

% training
t = templateTree('MinLeafSize', 10);
RF = fitrensemble(X, Y, 'Method', 'Bag', ...
    'Learners', t, ...
    'NumLearningCycles', 20);

Y_pro = [ x(1); predict(RF, X) ];
q_RF_pro = Y_pro .* s_q + m_q ;
R2_RF_pro = 1 - sum((q(2:end)-q_RF_pro(2:end)).^2) / ...
    sum((q(2:end)-m_q(2:end)).^2)

% validation
Xv = [xv(1:end-1), uv(1:end-1)];
Yv = xv(2:end);
Yv_pro = [ xv(1); predict(RF, Xv) ];
qv_RF_pro = Yv_pro .* s_qv + m_qv;
R2_RF_prov = 1 - sum((q_v(2:end)-qv_RF_pro(2:end)).^2) / ...
    sum((q_v(2:end)-m_qv(2:end)).^2)



% --------------------------
%%  RANDOM FOREST IMPROPER
% --------------------------

X = [x(1:end-1), u(2:end)];

% training
t = templateTree('MinLeafSize', 10);
RF = fitrensemble(X, Y, 'Method', 'Bag', ...
    'Learners', t, ...
    'NumLearningCycles', 20);
Y_imp = [ x(1); predict(RF, X) ];
q_RF_imp = Y_imp .* s_q + m_q ;
R2_RF_imp = 1 - sum((q(2:end)-q_RF_imp(2:end)).^2) / ...
    sum((q(2:end)-m_q(2:end)).^2)

% validation
Xv = [xv(1:end-1), uv(2:end)];
Yv = xv(2:end);
Yv_imp = [ xv(1); predict(RF, Xv) ];
qv_RF_imp = Yv_imp .* s_qv + m_qv;
R2_RF_impv = 1 - sum((q_v(2:end)-qv_RF_imp(2:end)).^2) / ...
    sum((q_v(2:end)-m_qv(2:end)).^2)



% --------------------
%%  FRIEDMAN DATASET
% --------------------


load -ascii 'friedman_dataset.txt'

X = friedman_dataset(1:2000, 1:15);
Y = friedman_dataset(1:2000, 16);
Xv = friedman_dataset(2001:2500, 1:15);
Yv = friedman_dataset(2001:2500, 16);

tic

% CART
R2_T = zeros(1, 2);
T = fitrtree(X, Y);
Y_ = predict(T, X);
Y_v = predict(T, Xv);
R2_T(1, 1) = (1-sum((Y-Y_).^2)/sum((Y-mean(Y)).^2));
R2_T(1, 2) = (1-sum((Yv-Y_v).^2)/sum((Yv-mean(Yv)).^2));

% RANDOM FOREST
n_trees_vec = [1 5 50 100 500];
R2_RF = zeros(length(n_trees_vec), 2);
for i = 1:length(n_trees_vec)
    
    t = templateTree('NumVariablesToSample', 5);
    RF = fitrensemble(X, Y, 'Method', 'Bag', ...
        'Learners', t, ...
        'NumLearningCycles', n_trees_vec(i));
    % fitrensemble uses bagging with random predictor selections
    % at each split (random forest) by default. To use bagging
    % without the random selections, use tree learners whose
    % 'NumVariablesToSample' value is 'all'.

    Y_ = predict(RF, X);
    Y_v = predict(RF, Xv);
    R2_RF(i, 1) = (1-sum((Y-Y_).^2)/sum((Y-mean(Y)).^2));
    R2_RF(i, 2) = (1-sum((Yv-Y_v).^2)/sum((Yv-mean(Yv)).^2));
end

[R2_T; R2_RF]

toc

figure
plot(n_trees_vec, R2_RF(:,2), 'o-')

%% 

qMonth = dailyToMonthly(q, 20); % m3/s

% water demand, quantile al 40 o al 65
w = quantile(q, .40) ; % m3/s

deltaT = 60*60*24*[31 28 31 30 31 30 31 31 30 31 30 31]';
% detltaT contains the number of seconds in each month
Q = qMonth(:).*repmat(deltaT,20,1); % m3/month
W = w*ones(size(Q)).*repmat(deltaT,20,1); % m3/month
figure; plot(Q); hold on; plot(W)

% Sequent Peak Analysis
K = zeros(size(Q));
K(1) = 0; %redundant line

for t = 1:length(Q)
    K(t+1) = K(t) + W(t) - Q(t);
    if K(t+1) < 0
        K(t+1) = 0;
    end
end

figure; plot( K )
% Valore dello storage(almeno)
Kopt = max(K)

% ipotizziamo 30 metri di altezza
s = Kopt / 30;
%                                                                                    ALT CAP SUP
%Spain	ESP		Negratin		Irrigation	Granada	Freila	Guadiana Menor			75	546	14.7
%Spain	ESP		Santa Teresa Saddle Dam		Irrigation	Salamanca	Montejo	Tormes			59	560	17.4
%Spain	ESP		Itoiz Dam	Itoitz	Irrigation	Navarra	Longuida	Irati			128	586	10.4
%Spain	ESP		Garcia de Sola	Puerto Pena	Irrigation	Badajoz	Talarrubias	Guadiana			65	640	22.8
%Spain	ESP		Riano		Irrigation	Leon	Cremenes	Esla			101	664	15.6
%Spain	ESP		Contreras		Hydroelectricity	Cuenca	Minglanilla	Cabriel			129	874	6.8
 
% s = superficie = 2.057251018946008e+07 km2
% Facendo l'operazione inversa(da active capacity stimo superficie(con
% orografie simili)risultano valori di altezza analoghi cioè circa 20km2)


%%  LOAD THE PARAMETERS OF THE LAKE
% -----------------------------------

% 211500000 			  % S    (lake surface [m2] = 211.5 [km2])
% 86.6021881286577        % beta (storage-discharge relationship)
% 1.95062290363465        % alfa (storage-discharge relationship)
% -1.27                   % h0   (storage-discharge relationship)
param.nat.S = s;
param.nat.beta = lake_mod_param(2);
param.nat.alpha = lake_mod_param(3);
param.nat.h0 = lake_mod_param(4);
% regulated level-discharge relationship:
param.reg.w = 100;
param.reg.h_min = 0;
param.reg.h_max = 1.5;
param.reg.h1 = 0.3;
param.reg.h2 = 1;
param.reg.m1 = 300;
param.reg.m2 = 1000;

% -----------------------------------------------
%%  TEST REGULATED LEVEL-DISCHARGE RELATIONSHIP
% -----------------------------------------------
h_test = [ -1.5 : 0.1 : 2 ];
r_test = regulated_release( param, h_test );
figure
plot(h_test, r_test, 'o')

% -------------------------------
%%  SIMULATION OF LAKE DYNAMICS
% -------------------------------
load -ascii Verbano_inflows_74_83_nb.txt
n = Verbano_inflows_74_83_nb;
n = [ nan; n ]; % for time convention
h_init = 0.88; % initial condition

[s, h, r] = simulate_lake( n, h_init, param );

figure
plot( h )

% ---------------------------
%%  INDICATORS FOR FLOODING
% ---------------------------

% IF1 = mean annual number of days with flooding
% events (in the city of Tudela)
h_flo = ;
figure
plot( h )
hold on
plot( ones( size( h ) ) * h_flo, 'r--' )

Ny = length(h)/365
Iflo = sum( h>h_flo )/Ny

% IF2 = maximum flooded area (in the city of Locarno)
S_flo = zeros( size(h) );
idx = find( h > h_flo );
S_flo(idx) = 0.081*h(idx).^3-0.483*h(idx).^2+1.506*h(idx)-1.578;
Iflo2 = max(S_flo)

h_max = max(h);
Iflo2 = 0.081*h_max^3-0.483*h_max^2+1.506*h_max-1.578

flooded_area = @(x) max(0, 0.081*x.^3-0.483*x.^2+1.506*x-1.578);
Iflo2 = max(flooded_area(h))



