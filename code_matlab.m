%% Simulation de la Trajectoire de Référence (Sujet Sain)
% Ce script résout l'équation dynamique : I*ddtheta + b*dtheta + k*theta = tau(t)
% Méthode numérique : Runge-Kutta 4 (RK4)

clear; clc; close all;

%% 1. Paramètres du Système (Modèle Biomécanique)
I = 0.25;      % Moment d'inertie (kg.m^2)
b = 0.15;      % Coefficient de frottement (N.m.s/rad)
k = 2.0;       % Raideur articulaire (N.m/rad)

% Paramètres temporels
dt = 0.01;              % Pas de temps (s)
t_final = 2;            % Durée du geste (s)
t = 0:dt:t_final;       % Vecteur temps
n = length(t);

%% 2. Définition du Couple Musculaire (Entrée tau)
% On simule un effort fluide et optimal (Sain)
tau = 0.8 * sin(pi * t / t_final); 

%% 3. Initialisation des vecteurs d'état
theta = zeros(1, n);        % Position angulaire (theta_clean)
theta_dot = zeros(1, n);    % Vitesse angulaire (dtheta/dt)

% Conditions Initiales (Le membre est au repos)
theta(1) = 0; 
theta_dot(1) = 0;

%% 4. Résolution par Runge-Kutta 4 (RK4)
% Fonction d'accélération : ddtheta = f(theta, theta_dot, tau)
f_accel = @(th, v, trq) (trq - b*v - k*th) / I;

for i = 1:(n-1)
    % Étape 1
    k1_v = f_accel(theta(i), theta_dot(i), tau(i));
    k1_th = theta_dot(i);
    
    % Étape 2
    k2_v = f_accel(theta(i) + 0.5*dt*k1_th, theta_dot(i) + 0.5*dt*k1_v, tau(i));
    k2_th = theta_dot(i) + 0.5*dt*k1_v;
    
    % Étape 3
    k3_v = f_accel(theta(i) + 0.5*dt*k2_th, theta_dot(i) + 0.5*dt*k2_v, tau(i));
    k3_th = theta_dot(i) + 0.5*dt*k2_v;
    
    % Étape 4
    k4_v = f_accel(theta(i) + dt*k3_th, theta_dot(i) + dt*k3_v, tau(i));
    k4_th = theta_dot(i) + dt*k3_v;
    
    % Mise à jour de la vitesse et de la position
    theta_dot(i+1) = theta_dot(i) + (dt/6) * (k1_v + 2*k2_v + 2*k3_v + k4_v);
    theta(i+1)     = theta(i)     + (dt/6) * (k1_th + 2*k2_th + 2*k3_th + k4_th);
end

%% 5. Calcul des dérivées supérieures (Accélération et Jerk)
% Accélération (theta_ddot)
theta_ddot = (tau - b*theta_dot - k*theta) / I;

% Jerk (Dérivée de l'accélération - indicateur de fluidité)
jerk = diff(theta_ddot)/dt;

%% 6. Affichage des Courbes de Référence
figure('Color', 'w', 'Position', [100, 100, 800, 600]);

subplot(3,1,1);
plot(t, theta, 'g', 'LineWidth', 2);
ylabel('\theta (rad)'); title('Trajectoire Propre (Position Angulaire)');
grid on;

subplot(3,1,2);
plot(t, theta_dot, 'b', 'LineWidth', 1.5);
ylabel('d\theta/dt (rad/s)'); title('Vitesse Angulaire');
grid on;

subplot(3,1,3);
plot(t, tau, 'k', 'LineWidth', 1.5);
ylabel('\tau (N.m)'); title('Couple Musculaire Optimal (Effort)');
xlabel('Temps (s)');
grid on;
%% 1. Chargement ou Récupération de theta_clean
% On suppose que 'theta' (du code précédent) est déjà dans le workspace
t = 0:0.01:2; % Vecteur temps (identique au précédent)

%% 2. Paramétrage du Bruit (Simulation du handicap)
% On définit l'intensité du bruit (ex: 5% de l'amplitude maximale)
amplitude_bruit = 0.05; 

% Génération d'un bruit blanc gaussien (Normal distribution)
% Ce bruit simule les imperfections motrices et les erreurs de mesure
bruit = amplitude_bruit * randn(size(theta)); 

% Création de la trajectoire bruitée (Données du Patient)
theta_noisy = theta + bruit; 

%% 3. Visualisation de la Trajectoire Altérée
figure('Color', 'w', 'Name', 'Données Patient vs Référence');
hold on;

% Affichage des données bruitées (en rouge)
plot(t, theta_noisy, 'r.', 'MarkerSize', 8, 'DisplayName', 'Mouvement Patient (Bruité)');

% Affichage de la référence saine (en vert)
plot(t, theta, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Référence Saine (Propre)');

xlabel('Temps (s)');
ylabel('\theta (rad)');
title('Simulation du mouvement altéré (Patient vs Sain)');
legend('Location', 'best');
grid on;

%% 4. Exportation des données pour le PINN (Python)
% On exporte uniquement le temps (t) et la position bruitée (theta_noisy)
% Car c'est la seule chose que le clinicien peut mesurer réellement
data_patient = [t', theta_noisy'];
writematrix(data_patient, 'data_patient_noisy.csv');

fprintf('Fichier "data_patient_noisy.csv" généré avec succès.\n');
fprintf('Prêt à être utilisé comme entrée pour le modèle PINN en Python.\n');
%% 5. Sauvegarde des données propres (Pour validation du PINN)
data_clean = [t', theta', theta_dot', theta_ddot', tau'];
save('mouvement_sain_ref.mat', 'data_clean');
writematrix(data_clean, 'theta_clean.csv');
%% Simulation de la Dynamique Musculaire et Exportation des Résultats
% Auteur : Karim Chibani
% Ce script génère les données de référence (Sain) et les données bruitées (Patient)
clear; clc; close all;

%% 1. Configuration du chemin d'exportation (Dossier Downloads)
if ispc % Pour Windows
    downloadPath = fullfile(char(java.lang.System.getProperty('user.home')), 'Downloads');
else % Pour Mac/Linux
    downloadPath = '~/Downloads';
end

%% 2. Paramètres du Système (Modèle de second ordre)
I = 0.25;      % Moment d'inertie (kg.m^2)
b = 0.15;      % Coefficient de frottement (N.m.s/rad)
k = 2.0;       % Raideur articulaire (N.m/rad)

dt = 0.01;              % Pas de temps (s)
t_final = 2;            % Durée totale de la simulation
t = 0:dt:t_final;       % Vecteur temps
n = length(t);

% Génération du Torque (Effort musculaire fluide)
tau = 0.8 * sin(pi * t / t_final); 

%% 3. Résolution Numérique via Runge-Kutta 4 (RK4)
theta = zeros(1, n); 
theta_dot = zeros(1, n);
f_accel = @(th, v, trq) (trq - b*v - k*th) / I;

for i = 1:(n-1)
    k1_v = f_accel(theta(i), theta_dot(i), tau(i));
    k1_th = theta_dot(i);
    
    k2_v = f_accel(theta(i) + 0.5*dt*k1_th, theta_dot(i) + 0.5*dt*k1_v, tau(i));
    k2_th = theta_dot(i) + 0.5*dt*k1_v;
    
    k3_v = f_accel(theta(i) + 0.5*dt*k2_th, theta_dot(i) + 0.5*dt*k2_v, tau(i));
    k3_th = theta_dot(i) + 0.5*dt*k2_v;
    
    k4_v = f_accel(theta(i) + dt*k3_th, theta_dot(i) + dt*k3_v, tau(i));
    k4_th = theta_dot(i) + dt*k3_v;
    
    theta_dot(i+1) = theta_dot(i) + (dt/6) * (k1_v + 2*k2_v + 2*k3_v + k4_v);
    theta(i+1)     = theta(i)     + (dt/6) * (k1_th + 2*k2_th + 2*k3_th + k4_th);
end

%% 🖼️ FIGURE 1 : GÉNÉRATION DES DONNÉES DE RÉFÉRENCE (SAIN)
fig1 = figure('Color', 'w', 'Position', [100, 100, 800, 800]);

subplot(3,1,1);
plot(t, theta, 'g', 'LineWidth', 2.5);
ylabel('\theta (rad)'); title('Trajectoire Propre (Position Angulaire)');
grid on;

subplot(3,1,2);
plot(t, theta_dot, 'b', 'LineWidth', 2);
ylabel('d\theta/dt (rad/s)'); title('Vitesse Angulaire');
grid on;

subplot(3,1,3);
plot(t, tau, 'k', 'LineWidth', 2.5);
ylabel('\tau (N.m)'); title('Couple Musculaire Optimal (Effort)');
xlabel('Temps (s)');
grid on;

% Sauvegarde automatique dans Downloads
saveas(fig1, fullfile(downloadPath, 'Reference_Propre_Matlab.png'));

%% 🖼️ FIGURE 2 : SIMULATION DU BRUIT GAUSSIEN (PATIENT)
amplitude_bruit = 0.05;
theta_noisy = theta + amplitude_bruit * randn(size(theta));

fig2 = figure('Color', 'w', 'Name', 'Analyse du Bruit');
plot(t, theta_noisy, 'r.', 'MarkerSize', 10); hold on;
plot(t, theta, 'g-', 'LineWidth', 3);
xlabel('Temps (s)'); ylabel('\theta (rad)');
title('Simulation du mouvement altéré (Patient vs Sain)');
legend('Données Bruitées (Patient)', 'Trajectoire Réelle (Saine)');
grid on;

% Sauvegarde automatique dans Downloads
saveas(fig2, fullfile(downloadPath, 'Simulation_Patient_Bruit.png'));

%% 📊 EXPORTATION DES DONNÉES CSV
writematrix([t', theta_noisy'], 'data_patient_noisy.csv');
writematrix([t', theta', theta_dot', tau'], 'theta_clean.csv');

fprintf('✅ Simulation terminée avec succès.\n');
fprintf('📁 Images exportées vers le dossier Downloads.\n');