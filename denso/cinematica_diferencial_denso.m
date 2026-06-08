%% Cinematica Diferencial e Jacobiana Geometrica
% Denso VP-6242 — Parametros DH Padrao
% Robotica — Nivel Graduacao
%
% MATLAB 2026a + Robotics System Toolbox
%
% Para abrir como LiveScript: File > Save As > .mlx
%
% Conteudo:
%   0. Inicializacao
%   1. Modelo DH e rigidBodyTree
%   2. Cinematica Direta
%   3. Jacobiana Geometrica (RST + manual)
%   4. Singularidades e Manipulabilidade
%   5. Controle por Velocidade: pinv e DLS
%   6. Visualizacao 3D com show()
%   7. Simulacao: Trajetoria Circular
%   8. Validacao: geometricJacobian vs Diferencas Finitas
%   9. Exercicios Propostos

%% 0. Inicializacao

clc; clear; close all;

% Verifica disponibilidade do Robotics System Toolbox
assert(~isempty(ver('robotics')), ...
    'Robotics System Toolbox nao encontrado. Instale via Add-Ons.');

disp('Robotics System Toolbox: OK');
disp(['MATLAB versao: ' version]);

%% 1. Modelo do Robo — Parametros DH e rigidBodyTree
%
% Convencao DH PADRAO (Craig):
%
%   Junta | a (m)  | alpha (rad)  | d (m)  | theta
%   ------|--------|--------------|--------|-------
%     1   | 0.0000 | deg2rad(-90) | 0.2635 |  q1
%     2   | 0.2100 | deg2rad(  0) | 0.0000 |  q2
%     3   | 0.0750 | deg2rad(-90) | 0.0000 |  q3
%     4   | 0.0000 | deg2rad( 90) | 0.2225 |  q4
%     5   | 0.0000 | deg2rad(-90) | 0.0000 |  q5
%     6   | 0.0000 | deg2rad(  0) | 0.0365 |  q6
%
% RST usa SI (metros). setFixedTransform(joint, dhparams, 'dh')
% aceita [a, alpha, d, theta_offset] diretamente.

% Tabela DH: [a(m), alpha(rad), d(m), theta_offset(rad)]
dh = [ 0.0000, deg2rad(-90), 0.2635, 0;
       0.2100, deg2rad(  0), 0.0000, 0;
       0.0750, deg2rad(-90), 0.0000, 0;
       0.0000, deg2rad( 90), 0.2225, 0;
       0.0000, deg2rad(-90), 0.0000, 0;
       0.0000, deg2rad(  0), 0.0365, 0];

% Criar arvore cinematica
robot    = rigidBodyTree('DataFormat','column','MaxNumBodies',7);
robot.Gravity = [0 0 -9.81];
prevBody = robot.BaseName;

for i = 1:size(dh,1)
    body  = rigidBody(sprintf('link%d', i));
    joint = rigidBodyJoint(sprintf('joint%d', i), 'revolute');
    joint.PositionLimits = [-pi, pi];
    setFixedTransform(joint, dh(i,:), 'dh');
    body.Joint = joint;
    addBody(robot, body, prevBody);
    prevBody = body.Name;
end

% Frame do efetuador (TCP fixo)
tcp       = rigidBody('tcp');
tcp.Joint = rigidBodyJoint('tcp_fixed', 'fixed');
setFixedTransform(tcp.Joint, eye(4));
addBody(robot, tcp, prevBody);

disp(['Robot criado: ' num2str(robot.NumBodies) ' corpos rigidos']);
showdetails(robot);

%% 2. Cinematica Direta
%
% getTransform(robot, config, endEffectorName, baseName)
% retorna a matriz homogenea ^0T_TCP (4x4).
%
% Nota: RST trabalha em metros. Multiplicar por 1e3 para mm.

% --- Configuracao zero ---
q0  = homeConfiguration(robot);
T06 = getTransform(robot, q0, 'tcp', robot.BaseName);
pe0 = T06(1:3,4) * 1e3;
fprintf('Posicao q=0 (mm): [%.2f, %.2f, %.2f]\n', pe0(1), pe0(2), pe0(3));

% --- Configuracao arbitraria ---
q_t = deg2rad([30, -45, 60, 0, 45, -30]).';
T_t = getTransform(robot, q_t, 'tcp', robot.BaseName);
p_t = T_t(1:3,4) * 1e3;
R_t = T_t(1:3,1:3);

fprintf('\nConfiguracao q = [30, -45, 60, 0, 45, -30] graus:\n');
fprintf('Posicao (mm): [%.2f, %.2f, %.2f]\n', p_t(1), p_t(2), p_t(3));
disp('Orientacao R:'); disp(round(R_t, 4));

%% 3. Jacobiana Geometrica
%
% --- Funcao nativa do Robotics System Toolbox ---
% geometricJacobian(robot, config, 'tcp') retorna J(q) [6x6].
%   Linhas 1-3: Jv (velocidade linear)
%   Linhas 4-6: Jw (velocidade angular)
%
% --- Deducao manual (didatica) ---
% Para junta de revolucao i:
%   Jv_i = z_{i-1} x (pe - p_{i-1})
%   Jw_i = z_{i-1}
% onde z_{i-1} = 3a coluna de ^0R_{i-1}

q_j   = deg2rad([30, -45, 60, 0, 45, -30]).';
J_rst = geometricJacobian(robot, q_j, 'tcp');

fprintf('Jacobiana Geometrica J(q) [6x6] — via geometricJacobian:\n');
disp(round(J_rst, 4));
fprintf('Posto: %d   Numero de condicao: %.2e\n', rank(J_rst), cond(J_rst));

% Implementacao manual para fins didaticos
J_man = jacGeometrica(robot, q_j);
fprintf('Erro max |J_RST - J_manual| = %.2e\n', max(abs(J_rst - J_man), [], 'all'));

% Heatmap das sub-Jacobianas
figure('Position', [100 100 880 320]);
subplot(1,2,1);
imagesc(J_rst(1:3,:)); colorbar;
set(gca, 'XTick',1:6, 'XTickLabel',{'J1','J2','J3','J4','J5','J6'});
set(gca, 'YTick',1:3, 'YTickLabel',{'vx','vy','vz'});
title('Jv(q) — parte linear (3x6)');

subplot(1,2,2);
imagesc(J_rst(4:6,:)); colorbar;
set(gca, 'XTick',1:6, 'XTickLabel',{'J1','J2','J3','J4','J5','J6'});
set(gca, 'YTick',1:3, 'YTickLabel',{'wx','wy','wz'});
title('Jw(q) — parte angular (3x6)');

sgtitle('Jacobiana Geometrica — q = [30,-45,60,0,45,-30] graus');
colormap(turbo);

%% 4. Singularidades e Manipulabilidade
%
% Indice de Yoshikawa:
%   w(q) = sqrt(det(J * J^T))
%   w = 0  =>  singularidade cinematica
%
% Denso VP-6242:
%   Pulso: q5 = +/- 180 graus  (eixos J4 e J6 colineares)
%   Cotovelo: q3 ~= 0           (configuracao estendida)
%   Ombro: efetuador sobre z0   (x=0, y=0)

q_b  = deg2rad([0, -30, 45, 10, 0, 30]).';
q5v  = linspace(-180, 180, 360);
wv   = zeros(size(q5v));

for k = 1:numel(q5v)
    qs    = q_b;
    qs(5) = deg2rad(q5v(k));
    wv(k) = manipulabilidade(robot, qs);
end

figure;
plot(q5v, wv, 'Color', [0.12 0.37 0.64], 'LineWidth', 2);
xline(-180, '--r', 'Label', 'q5 = +/-180 (singular)', 'LineWidth', 1.5);
xline( 180, '--r', 'LineWidth', 1.5);
xlabel('q5 (graus)'); ylabel('w(q) [manipulabilidade]');
title('Manipulabilidade vs q5 — Singularidade de Pulso'); grid on;

qs_sing = q_b; qs_sing(5) = deg2rad(-180);
qs_ok   = q_b; qs_ok(5)   = deg2rad(90);
fprintf('w(q5=-180 graus) = %.2e  -> SINGULAR\n', manipulabilidade(robot, qs_sing));
fprintf('w(q5=  90 graus) = %.4e  -> regular\n',  manipulabilidade(robot, qs_ok));

% Mapa 2D: plano (q2, q3)
q2v = linspace(-90, 90, 50);
q3v = linspace(-90, 90, 50);
W   = zeros(50, 50);
qm  = deg2rad([0, 0, 0, 0, 90, 0]).';

for ii = 1:50
    for jj = 1:50
        qm(2)    = deg2rad(q2v(jj));
        qm(3)    = deg2rad(q3v(ii));
        W(ii,jj) = manipulabilidade(robot, qm);
    end
end

figure;
contourf(q2v, q3v, W, 30); colorbar; hold on;
contour(q2v, q3v, W, [0 0], 'r', 'LineWidth', 2);
xlabel('q2 (graus)'); ylabel('q3 (graus)');
title('Mapa de Manipulabilidade — plano (q2, q3)');
colormap(parula); grid on;

%% 5. Controle por Velocidade: Pseudo-Inversa e DLS
%
% Pseudo-inversa de Moore-Penrose:
%   qdot = J^+(q) * xdot_d,    J^+ = J^T (J J^T)^{-1}
%
% Damped Least Squares (DLS) — robusto perto de singularidades:
%   qdot = J^T (J J^T + lambda^2 I)^{-1} * xdot_d
%
% lambda > 0 limita a amplificacao de qdot.

q_near = deg2rad([0, -45, 0, 0, 0, 0]).';
q_ok2  = deg2rad([0, -30, 45, 10, 90, 30]).';
xd     = [0.01; 0; 0.005; 0; 0; 0];   % (m/s e rad/s)

lsw = logspace(-3, 0, 200);
Jn  = geometricJacobian(robot, q_near, 'tcp');
Jo  = geometricJacobian(robot, q_ok2,  'tcp');

ns = arrayfun(@(l) norm(dlsControl(Jn, xd, l)), lsw);
no = arrayfun(@(l) norm(dlsControl(Jo, xd, l)), lsw);

figure;
loglog(lsw, ns, 'r', 'LineWidth', 2, 'DisplayName', 'Perto singularidade'); hold on;
loglog(lsw, no, 'b', 'LineWidth', 2, 'DisplayName', 'Configuracao regular');
xline(0.05, '--k', 'Label', 'lambda = 0.05', 'LineWidth', 1.2);
xlabel('lambda (fator DLS)'); ylabel('||qdot|| (rad/s)');
title('Efeito do Amortecimento DLS na Norma de qdot');
legend; grid on;

%% 6. Visualizacao 3D com show()
%
% show(robot, config, 'Frames','on') renderiza links, juntas e frames.

cfgs = { zeros(6,1),                          'q = zeros (extensao)';
         deg2rad([0,-45,90,0,45,0]).',         'q = [0,-45,90,0,45,0] graus';
         deg2rad([45,-30,60,-20,70,30]).',     'q = [45,-30,60,-20,70,30] graus' };

figure('Position', [50 50 1100 400]);
for k = 1:3
    subplot(1,3,k);
    show(robot, cfgs{k,1}, 'Frames', 'on', 'PreservePlot', false);
    title(cfgs{k,2}, 'FontSize', 9);
    view(30, 20); grid on; axis equal;
    xlim([-0.65 0.65]); ylim([-0.65 0.65]); zlim([0 0.80]);
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
end
sgtitle('Denso VP-6242 — Visualizacao 3D (rigidBodyTree)', 'FontWeight', 'bold');

%% 7. Simulacao: Trajetoria Circular no Plano XY
%
% Controle por velocidade — integracao de Euler:
%   q(t+dt) = q(t) + J^+(q) * xdot_d(t) * dt
%
% Trajetoria: circulo de raio 80 mm, frequencia pi rad/s.

dt    = 0.005;
T_sim = 2.0;
t_arr = 0:dt:T_sim-dt;
N     = numel(t_arr);
Rc    = 0.08;    % raio (m)
wc    = pi;      % frequencia angular (rad/s)
lam   = 0.05;

xdot_f = @(t) [-Rc*wc*sin(wc*t); Rc*wc*cos(wc*t); 0; 0; 0; 0];

q  = deg2rad([0, -30, 60, 0, 60, 0]).';
pt = zeros(3, N);
qt = zeros(6, N);
wt = zeros(1, N);

for k = 1:N
    Tc       = getTransform(robot, q, 'tcp', robot.BaseName);
    pt(:,k)  = Tc(1:3,4);
    qt(:,k)  = q;
    wt(k)    = manipulabilidade(robot, q);
    J        = geometricJacobian(robot, q, 'tcp');
    q        = q + dlsControl(J, xdot_f(t_arr(k)), lam) * dt;
end

fprintf('Erro de fechamento: %.4f mm\n', norm(pt(:,end) - pt(:,1)) * 1e3);

figure('Position', [50 50 1200 380]);

subplot(1,3,1);
plot3(pt(1,:)*1e3, pt(2,:)*1e3, pt(3,:)*1e3, 'r-', 'LineWidth', 2);
hold on;
plot3(pt(1,1)*1e3, pt(2,1)*1e3, pt(3,1)*1e3, 'go', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
title('Trajetoria 3D do Efetuador'); grid on; view(30, 25);

subplot(1,3,2);
plot(t_arr, rad2deg(qt.'));
xlabel('Tempo (s)'); ylabel('Angulo (graus)');
legend('q1','q2','q3','q4','q5','q6', 'Location', 'best');
title('Angulos das Juntas q(t)'); grid on;

subplot(1,3,3);
plot(t_arr, wt, 'b', 'LineWidth', 2);
yline(mean(wt), '--k', sprintf('Media = %.2e', mean(wt)));
xlabel('Tempo (s)'); ylabel('w(q)');
title('Manipulabilidade ao Longo da Trajetoria'); grid on;

sgtitle('Simulacao: Trajetoria Circular — DLS lambda = 0.05', 'FontWeight', 'bold');

%% 8. Validacao: geometricJacobian vs Diferencas Finitas
%
% Jacobiana numerica por diferencas finitas centrais (parte Jv, 3x6).
% Erro esperado da ordem de eps^(2/3) ~ 1e-10.

q_val  = deg2rad([15, -30, 45, -10, 60, 20]).';
J_rst2 = geometricJacobian(robot, q_val, 'tcp');
J_num  = jacNumerica(robot, q_val);

erro = max(abs(J_rst2(1:3,:) - J_num), [], 'all');
fprintf('Erro max |Jv_RST - Jv_numerica| = %.2e\n', erro);
assert(erro < 1e-5, 'Erro acima do tolerado — verificar implementacao!');
disp('Validacao FK/Jacobiana: OK');

sv = svd(J_rst2);
fprintf('Valores singulares:'); fprintf('  %.4f', sv); fprintf('\n');
fprintf('Numero de condicao: %.2e\n', sv(1)/sv(end));

qs = deg2rad([0, -30, 45, 10, -180, 30]).';
fprintf('w(q5=-180 graus) = %.2e  (deve ser ~0)\n', manipulabilidade(robot, qs));

%% 9. Exercicios Propostos
%
% -------------------------------------------------------------------------
% Exercicio 1 — Singularidade de Ombro
%   Encontre q tal que x=0 e y=0 no efetuador (sobre o eixo z0).
%   Verifique w(q) = 0 e identifique qual coluna de J e dependente.
%
% Exercicio 2 — Controle Proporcional de Posicao
%   xdot_d = Kp * (pd - pe). Simule convergencia para pd desejado.
%   Compare pinv e DLS. Qual converge mais suavemente?
%
% Exercicio 3 — Jacobiana de Posicao (Jv, 3x6)
%   Extraia so Jv e use pinv(Jv). O que acontece com a orientacao?
%
% Exercicio 4 — Otimizacao via Espaco Nulo
%   qdot = J^+ xdot_d + (I - J^+ J) * k * grad_q(w)
%   Maximize w(q) ao longo da trajetoria sem alterar a posicao do efetuador.
%
% Exercicio 5 — Cinematica Inversa por Newton-Raphson
%   q(k+1) = q(k) + J^+ * (pd - f(q(k)))
%   Compare convergencia com ikine() do Robotics System Toolbox.
% -------------------------------------------------------------------------

% Template: gradiente de manipulabilidade e projetor de espaco nulo
q_ex   = deg2rad([0, -30, 60, 0, 60, 0]).';
J_ex   = geometricJacobian(robot, q_ex, 'tcp');
gw     = gradManip(robot, q_ex);
N_ex   = eye(6) - pinv(J_ex) * J_ex;   % projetor no espaco nulo

fprintf('Gradiente de w(q):');
fprintf('  %.4f', gw.');
fprintf('\n');
fprintf('(Para robo nao redundante n=6 com posto 6: N = 0 => sem espaco nulo)\n');

%% Funcoes Auxiliares
% (Colocar em arquivos .m separados ou ao final do script)

function J = jacGeometrica(robot, q)
% JACGEOMETRICA  Jacobiana Geometrica 6xn por deducao direta.
% Ji = [zi-1 x (pe - pi-1); zi-1]  para junta de revolucao i.
    bodies = robot.BodyNames;
    n      = numel(bodies) - 1;    % exclui tcp (fixed)
    T_tcp  = getTransform(robot, q, 'tcp', robot.BaseName);
    pe     = T_tcp(1:3, 4);
    J      = zeros(6, n);
    for i = 1:n
        if i == 1
            T_prev = eye(4);
        else
            T_prev = getTransform(robot, q, bodies{i-1}, robot.BaseName);
        end
        z        = T_prev(1:3, 3);
        p        = T_prev(1:3, 4);
        J(1:3,i) = cross(z, pe - p);
        J(4:6,i) = z;
    end
end

function w = manipulabilidade(robot, q)
% MANIPULABILIDADE  Indice de Yoshikawa w = sqrt(det(J*J^T)).
    J = geometricJacobian(robot, q, 'tcp');
    w = sqrt(max(0, det(J * J.')));
end

function qdot = pinvControl(J, xdot, tol)
% PINVCONTROL  Controle por pseudo-inversa de Moore-Penrose.
    if nargin < 3, tol = 1e-6; end
    qdot = pinv(J, tol) * xdot;
end

function qdot = dlsControl(J, xdot, lam)
% DLSCONTROL  Controle Damped Least Squares.
    if nargin < 3, lam = 0.05; end
    m    = size(J, 1);
    qdot = J.' * ((J * J.' + lam^2 * eye(m)) \ xdot);
end

function Jn = jacNumerica(robot, q, ep)
% JACNUMERICA  Jacobiana numerica de posicao por diferencas finitas centrais.
    if nargin < 3, ep = 1e-7; end
    n  = numel(q);
    Jn = zeros(3, n);
    for i = 1:n
        dq       = zeros(n,1); dq(i) = ep;
        Tp       = getTransform(robot, q+dq, 'tcp', robot.BaseName);
        Tm       = getTransform(robot, q-dq, 'tcp', robot.BaseName);
        Jn(:,i)  = (Tp(1:3,4) - Tm(1:3,4)) / (2*ep);
    end
end

function g = gradManip(robot, q, ep)
% GRADMANIP  Gradiente numerico da manipulabilidade.
    if nargin < 3, ep = 1e-5; end
    n = numel(q);
    g = zeros(n, 1);
    for i = 1:n
        dq   = zeros(n,1); dq(i) = ep;
        g(i) = (manipulabilidade(robot, q+dq) - manipulabilidade(robot, q-dq)) / (2*ep);
    end
end
