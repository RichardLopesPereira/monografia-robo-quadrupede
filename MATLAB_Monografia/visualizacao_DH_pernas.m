%% Visualização individual das quatro pernas com eixos DH e Z do efetuador ajustado
clear; close all; clc;
% ---------- TABELA DH (comum às pernas) ----------
% Colunas: Theta, d, a, Alpha. Valores iniciais são para 'homing' (juntas em 0).
DH = [
    0,    0.035,  0.065,  -90;   % Define a transformação T(J1 -> J2)
    0,    0.000,  0.100,   180;  % Define a transformação T(J2 -> J3)
    90,   0.000,  0.165,   180;  % Define a transformação T(J3 -> J4)
];
% NOTA: A tabela DH define a transformação entre juntas. 
% [Código de Transformações base -> primeira junta de cada perna permanece a mesma]
T0_to_J1_perna1 =  [1  0  0  0.07; 0 -1  0 -0.083; 0  0 -1  0; 0  0  0  1];
T0_to_J1_perna2 = [-1  0  0 -0.07; 0  1  0  0.05; 0  0 -1  0; 0  0  0  1];
T0_to_J1_perna3 = [1  0  0  0.07; 0 -1  0  0.05; 0  0 -1  0; 0  0  0  1];
T0_to_J1_perna4 = [-1  0  0 -0.07; 0  1  0 -0.083; 0  0 -1  0; 0  0  0  1];
% ---------- Configurações visuais ----------
frame_scale = 0.02;
axis_lim = 0.35; 
% ---------- Monta e Plota Cada Perna Individualmente ----------
plotPernaIndividual(1, T0_to_J1_perna1, DH, frame_scale, axis_lim);
plotPernaIndividual(2, T0_to_J1_perna2, DH, frame_scale, axis_lim);
plotPernaIndividual(3, T0_to_J1_perna3, DH, frame_scale, axis_lim);
plotPernaIndividual(4, T0_to_J1_perna4, DH, frame_scale, axis_lim);
%% ----------------- Funções Locais -------------------------
function plotPernaIndividual(num_perna, T0_base, DH, s, axis_lim)
    
    junta_offset = (num_perna - 1) * 3; % 3 juntas por perna

    % Monta e ajusta o Z do efetuador
    [Ts, T_end, T0_adj] = montarPerna(T0_base, DH);
    
    % Cria nova figura para cada perna
    figure('Color','w');
    hold on; grid on; axis equal;
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    view(120,20);
    
    % Posicionamento e Escala dos Gráficos 
    xlim([-axis_lim axis_lim]); ylim([-axis_lim axis_lim]); zlim([-0.07 0.25]); 
    
    % Título da Perna
    title(sprintf('Perna %d - Eixos DH (J%d a J%d e Efetuador)', num_perna, junta_offset+1, junta_offset+3));
    
    primeira_junta_nome = sprintf('J%d', junta_offset + 1);
    plotFrameCustom(T0_adj, s, primeira_junta_nome); 
    
    % Desenha os frames DH e as linhas dos elos
    % i varia de 1 a 3 (para as 3 TFs DH)
    for i = 1:size(DH,1)
        
        % Nome da junta atual (J2/J5/J8/J11) e J3/J6/J9/J12
        num_junta = junta_offset + i + 1;
        nome_junta = sprintf('J%d', num_junta);
        
        % Desenha o frame DH (Ts(:,:,2) a Ts(:,:,4))
        plotFrameCustom(Ts(:,:,i+1), s, nome_junta);
        
        % Desenha o elo entre o frame anterior e o frame atual
        p0 = Ts(1:3,4,i);    % Posição do frame anterior (J1, J2, J3)
        p1 = Ts(1:3,4,i+1);  % Posição do frame atual (J2, J3, J4)
        plot3([p0(1) p1(1)], [p0(2) p1(2)], [p0(3) p1(3)], '-k', 'LineWidth', 1.5);
    end
    
    % Efetuador final (em Z=0) - Corresponde à posição da ÚLTIMA junta (J3/J6/J9/J12)
    plot3(T_end(1,4), T_end(2,4), 0, 'ro', 'MarkerFaceColor','r','MarkerSize',6);

    text(T_end(1,4), T_end(2,4), 0.01, sprintf('E - P%d', num_perna), 'Color','black','FontSize',9);
    
   % Plano Z=0 (chão)
    fill3([-axis_lim axis_lim axis_lim -axis_lim], [-axis_lim -axis_lim axis_lim axis_lim], [0 0 0 0], [0.9 0.9 0.9], 'FaceAlpha',0.2, 'EdgeColor','none');
    
    % Legenda de Cores
    % Armazena os handles das linhas de legenda
    h(1) = plot3(NaN, NaN, NaN, 'r', 'LineWidth', 2, 'DisplayName', 'Eixo X');
    h(2) = plot3(NaN, NaN, NaN, 'g', 'LineWidth', 2, 'DisplayName', 'Eixo Y');
    h(3) = plot3(NaN, NaN, NaN, 'b', 'LineWidth', 2, 'DisplayName', 'Eixo Z');
    
    % Passa explicitamente os handles para a função legend
    legend(h, 'Location', 'NorthEastOutside');
    
    % Salva imagem
    exportgraphics(gca, sprintf('Perna_%d_J%d_a_J%d.png', num_perna, junta_offset+1, junta_offset+3), 'Resolution', 300);
end
function [Ts,T_end,T0_adj] = montarPerna(T0_base,DH)
    nloc = size(DH,1);
    
    % 1. Cáculo das TFs sem ajuste (para encontrar o z_end)
    Ts_calc1 = zeros(4,4,nloc+1);
    Ts_calc1(:,:,1) = T0_base;
    Tprev = T0_base;
    for ii=1:nloc
        theta = DH(ii,1); d = DH(ii,2); a = DH(ii,3); alpha = DH(ii,4);
        A_ii = dhA(theta,d,a,alpha);
        Tprev = Tprev * A_ii;
        Ts_calc1(:,:,ii+1) = Tprev;
    end
    T_end_calc1 = Ts_calc1(:,:,end);
    
    % Ajuste z_end -> 0: Determina o deslocamento necessário na base
    z_end = T_end_calc1(3,4);
    T0_adj = T0_base;
    T0_adj(3,4) = T0_adj(3,4) - z_end;
    
    % 2. Cáculo das TFs com o T0_adj (Frame base ajustado J1)
    Ts = zeros(4,4,nloc+1);
    Tprev = T0_adj; 
    Ts(:,:,1) = T0_adj; 
    for ii=1:nloc
        theta = DH(ii,1); d = DH(ii,2); a = DH(ii,3); alpha = DH(ii,4);
        A_ii = dhA(theta,d,a,alpha);
        Tprev = Tprev * A_ii;
        Ts(:,:,ii+1) = Tprev; 
    end
    T_end = Ts(:,:,end); % Efetuador final (Frame da última junta)
end

function plotFrameCustom(T, s, name)
    o = T(1:3,4); x = o + s * T(1:3,1); y = o + s * T(1:3,2); z = o + s * T(1:3,3);
    plot3([o(1) x(1)], [o(2) x(2)], [o(3) x(3)], 'r', 'LineWidth', 2);
    plot3([o(1) y(1)], [o(2) y(2)], [o(3) y(3)], 'g', 'LineWidth', 2);
    plot3([o(1) z(1)], [o(2) z(2)], [o(3) z(3)], 'b', 'LineWidth', 2);
    plot3(o(1), o(2), o(3), 'ro', 'MarkerSize', 4, 'MarkerFaceColor','r');
    text(x(1), x(2), x(3), ['x ' name], 'FontSize', 9, 'Color', 'k');
    text(y(1), y(2), y(3), ['y ' name], 'FontSize', 9, 'Color', 'k');
    text(z(1), z(2), z(3), ['z ' name], 'FontSize', 9, 'Color', 'k');
end
function A = dhA(theta_deg, d, a, alpha_deg)
    th = deg2rad(theta_deg); al = deg2rad(alpha_deg);
    A = [ cos(th), -sin(th)*cos(al),  sin(th)*sin(al), a*cos(th);
          sin(th),  cos(th)*cos(al), -cos(th)*sin(al), a*sin(th);
          0,        sin(al),          cos(al),          d;
          0,        0,                0,                1];
end
function T_final = aplicarRotacaoZAntes(angulo_graus, T_original)
    theta = angulo_graus;
    Rz = [cosd(theta) -sind(theta) 0; sind(theta)  cosd(theta) 0; 0 0 1];
    T_rotZ = [Rz, zeros(3,1); 0 0 0 1];
    T_final = T_original * T_rotZ; 
end