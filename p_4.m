clear; close all; clc;

% === 0.1 Leer audio (ajusta la ruta a tus archivos) ===
% Ejemplos del repositorio de Michigan (cambia por tus rutas):
% [x, fs] = audioread('data/09_apex_holo_sys_mur_supine_bell.mp3');
[x, fs] = audioread('data/02_apex_split_s1_supine_bell.mp3');
% [x, fs] = audioread('data/07_apex_mid_sys_mur_supine_bell.mp3');
% [x, fs] = audioread('data/10_apex_sys_click__late_sys_mur_lld_bell.mp3');
% [x, fs] = audioread('data/23_pulm_eject_sys_mur__single_s2__eject_click_supine_diaph.mp3');
% === 0.2 Si es estéreo, convertir a mono ===
if size(x,2) > 1
    x = mean(x, 2);
end

% === 0.3 Normalizar a [-1, 1] ===
a = -1; b = 1;
x = (x - min(x)) / max(eps, (max(x) - min(x)));  % [0,1]
x = x * (b - a) + a;                              % [-1,1]

% === 0.4 Ventana de trabajo (2 s para empezar) ===
t = (0:length(x)-1)/fs;
max_time = 5;                         % segundos
sel = t <= max_time;
x = x(sel);  t = t(sel);

% === 1.1 Probabilidad normalizada p[n] ===
p = abs(x);
p = p ./ max(eps, max(p));

% === 1.2 Entropía/energía de Shannon puntual (log10) ===
E = -p .* log(p + eps);     % variante equivalente a -|x|log|x|

% === 1.3 Estandarización y normalización a [0,1] ===
E_z  = (E - mean(E)) / (std(E) + eps);
Env0 = (E_z - min(E_z)) / max(eps, (max(E_z) - min(E_z)));   % envolvente base [0,1]

% === 2.1 LPF sobre la envolvente (no sobre x) ===
fc = 12;                                         % corte ~10 Hz
[b,a] = butter(4, fc/(fs/2), 'low');
Env = filtfilt(b, a, Env0);

% === 2.2 Normalización final a [0,1] ===
Env = (Env - min(Env)) / max(eps, (max(Env) - min(Env)));

% Visual rápido
figure('Name','Señal & Envolvente (Shannon LPF)');
subplot(2,1,1); plot(t, x, 'k'); grid on;
title('Señal PCG (normalizada)'); xlabel('Tiempo (s)'); ylabel('Amplitud');

subplot(2,1,2); plot(t, Env, 'b'); grid on;
title('Envolvente de Shannon (LPF)'); xlabel('Tiempo (s)'); ylabel('Amplitud norm.');

% Visual rápido
figure('Name','Señal & Envolvente (Shannon LPF) superpuestas');
plot(t, x, 'k');  hold on;
plot(t, Env, 'b'); grid on;
legend('Señal PCG (normalizada)','Envolvente de Shannon (LPF)'); xlabel('Tiempo (s)'); ylabel('Amplitud norm.');
grid on;

% === 3.1 Derivada discreta y detección de cambios de signo ===
d = diff(Env);
idx_ext = []; tipo = [];  % tipo: 1 = min, 2 = max

for i = 1:length(d)-1
    if d(i) < 0 && d(i+1) > 0
        idx_ext(end+1) = i+1; tipo(end+1) = 1;  % mínimo
    elseif d(i) > 0 && d(i+1) < 0
        idx_ext(end+1) = i+1; tipo(end+1) = 2;  % máximo
    end
end

% Visual
figure('Name','Envolvente con min/máx');
plot(t, Env, 'b'); hold on; grid on;
plot(t(idx_ext(tipo==1)), Env(idx_ext(tipo==1)), 'go', 'DisplayName','Mínimos');
plot(t(idx_ext(tipo==2)), Env(idx_ext(tipo==2)), 'ro', 'DisplayName','Máximos');
legend; xlabel('Tiempo (s)'); ylabel('Amplitud norm.');
title('Extremos locales sobre la envolvente');

% === 4.1 Construcción de tripletes mín–máx–mín ===
tri_samp = []; tri_time = []; tri_amp = [];
for k = 1:length(tipo)-2
    if tipo(k)==1 && tipo(k+1)==2 && tipo(k+2)==1
        i1 = idx_ext(k); i2 = idx_ext(k+1); i3 = idx_ext(k+2);
        tri_samp(end+1,:) = [i1,i2,i3];                  %#ok<SAGROW>
        tri_time(end+1,:) = [t(i1), t(i2), t(i3)];       %#ok<SAGROW>
        tri_amp(end+1,:)  = [Env(i1), Env(i2), Env(i3)]; %#ok<SAGROW>
    end
end

% === 4.2 Área de cada triángulo ===
areas = zeros(size(tri_time,1),1);
for i = 1:size(tri_time,1)
    x1=tri_time(i,1); y1=tri_amp(i,1);
    x2=tri_time(i,2); y2=tri_amp(i,2);
    x3=tri_time(i,3); y3=tri_amp(i,3);
    areas(i) = 0.5 * abs( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) );
end

% Visual
figure('Name','Triángulos mín–máx–mín');
plot(t, Env, 'g', 'LineWidth', 1.2); hold on; grid on;
for i = 1:size(tri_time,1)
    tx = [tri_time(i,:), tri_time(i,1)];
    ty = [tri_amp(i,:),  tri_amp(i,1)];
    plot(tx, ty, 'r-', 'LineWidth', 1.0);
end
plot(t(idx_ext), Env(idx_ext), 'ko', 'MarkerSize', 4);
xlabel('Tiempo (s)'); ylabel('Amplitud norm.');
title('Envolvente y triángulos candidatos');

% === 5.1 Selección por área (umbral: promedio) ===
Amed = mean(areas);
Amed = Amed * 1.2;
mask_big = areas > Amed;


% === 5.2 Propuesta de ciclos: [min1, min2] de cada triángulo grande ===
ciclos_idx = tri_samp(mask_big, [1 3]);  % índices
ciclos_idx = sortrows(ciclos_idx, 1);

% === 5.3 Limpieza por duración fisiológica y no solape ===
minRR = round(0.20*fs);   % 0.30 s
maxRR = round(2.00*fs);   % 1.50 s
ciclos_ref = [];
for i = 1:size(ciclos_idx,1)
    L = ciclos_idx(i,2) - ciclos_idx(i,1);
    if L >= minRR && L <= maxRR
        if isempty(ciclos_ref) || ciclos_idx(i,1) > ciclos_ref(end,2)
            ciclos_ref = [ciclos_ref; ciclos_idx(i,:)]; %#ok<AGROW>
        end
    end
end

% Visual
figure('Name','Triángulos y original');
plot(t, x, 'k'); grid on; hold on;
plot(t, Env, 'g', 'LineWidth', 1.2); hold on; grid on;
for i = find(mask_big)'
    tx = [tri_time(i,:), tri_time(i,1)];
    ty = [tri_amp(i,:),  tri_amp(i,1)];
    plot(tx, ty, 'r-', 'LineWidth', 1.0);
end
plot(t(idx_ext), Env(idx_ext), 'ko', 'MarkerSize', 4);
xlabel('Tiempo (s)'); ylabel('Amplitud norm.');
title('Envolvente y triángulos candidatos');

% ---------- INPUTS (debe existir) ----------
% t, x, Env, tri_time, tri_amp, mask_big
% ---------- FIN INPUTS ----------

% --- 1) Detectar crestas (simplificado) ---
tri_idx = find(mask_big)';
if isempty(tri_idx)
    error('No se encontraron triángulos candidatos (mask_big está vacío).');
end

% Los picos SON el punto central (i2) de los triángulos "grandes"
peak_indices = tri_samp(mask_big, 2);

% Ordenar y eliminar duplicados (si los triángulos se solaparon)
peak_indices = unique(peak_indices);

peaks_time = t(peak_indices);
peaks_amp  = Env(peak_indices);

nPeaks = numel(peaks_time);

if nPeaks < 2
    error('No hay suficientes crestas detectadas para clasificar S1/S2.');
end

% --- 2) Calcular dt y separar en "corto" y "largo" (kmeans) ---
dt = diff(peaks_time);   % nPeaks-1

% Fallback por si kmeans no está disponible o falla
useKmeans = true;
try
    % requerirá Statistics Toolbox, replicates para estabilidad
    [clustIdx, C] = kmeans(dt(:), 2, 'Replicates', 10, 'MaxIter', 500);
catch
    useKmeans = false;
end

if useKmeans
    % cluster menor = intervalo corto (S1->S2)
    [~, shortCluster] = min(C);
    isShort = (clustIdx == shortCluster);
    short_mean = C(shortCluster);
    long_mean  = C(3-shortCluster);
else
    % fallback simple: umbral por mediana
    med = median(dt);
    isShort = dt < med;
    short_mean = mean(dt(isShort));
    long_mean  = mean(dt(~isShort));
end

fprintf('Intervalos: media corto=%.3f s, media largo=%.3f s, nPeaks=%d\n', short_mean, long_mean, nPeaks);

% --- 3) Etiquetado secuencial robusto: si dt(i) es corto -> i=S1, i+1=S2.
labels = strings(nPeaks,1);
i = 1;
while i <= nPeaks
    if i == nPeaks
        % último pico sin forward interval: decidir por vecino
        if i>1 && labels(i-1) == "S1"
            labels(i) = "S2";
        else
            labels(i) = "S1";
        end
        break;
    end

    if isShort(i)   % dt(i) = peaks_time(i+1)-peaks_time(i) corto -> forma par S1-S2
        labels(i)   = "S1";
        labels(i+1) = "S2";
        i = i + 2;
    else
        % dt(i) largo -> es muy probable que estemos en S2 -> S1 siguiente
        if labels(i) == ""
            labels(i) = "S2";
        end
        i = i + 1;
    end
end

% --- 4) Rellenar cualquier etiqueta vacía con heurística de vecino/amplitud ---
for i = 1:nPeaks
    if labels(i) == ""
        if i>1 && labels(i-1) == "S1"
            labels(i) = "S2";
        elseif i<nPeaks && labels(i+1) == "S2"
            labels(i) = "S1";
        else
            % fallback por amplitud relativa
            if peaks_amp(i) >= median(peaks_amp)
                labels(i) = "S1";
            else
                labels(i) = "S2";
            end
        end
    end
end

% --- 5) Calcular BPM desde S1 detectadas ---
S1_times = peaks_time(labels == "S1");
if numel(S1_times) >= 2
    rr = diff(S1_times);       % segundos entre S1
    bpm = 60 / median(rr);
else
    bpm = NaN;
end


% --- 6) Graficar resultados ---
figure('Name','Detección S1/S2 mejorada'); clf;
ax_s1s2 = gca; % <-- GUARDAR EL HANDLE DE LOS EJES

plot(ax_s1s2, t, x, 'k'); hold(ax_s1s2, 'on'); grid(ax_s1s2, 'on');
plot(ax_s1s2, t, Env, 'g', 'LineWidth', 1.2);
% triángulos
tri_idx = find(mask_big)'; % Asegúrate que tri_idx exista (de la corrección 1)
for iTri = tri_idx
    tx = [tri_time(iTri,:), tri_time(iTri,1)];
    ty = [tri_amp(iTri,:),  tri_amp(iTri,1)];
    plot(ax_s1s2, tx, ty, 'r-', 'LineWidth', 1.0);
end
% marcadores para S1 y S2
S1_idx = find(labels == "S1");
S2_idx = find(labels == "S2");
plot(ax_s1s2, peaks_time(S1_idx), peaks_amp(S1_idx), 'bv', 'MarkerFaceColor','b', 'MarkerSize',8); % S1 azul
plot(ax_s1s2, peaks_time(S2_idx), peaks_amp(S2_idx), 'rv', 'MarkerFaceColor','r', 'MarkerSize',8); % S2 rojo
% Etiquetas de texto
for k = 1:nPeaks
    txt = labels(k);
    text(ax_s1s2, peaks_time(k), peaks_amp(k)+0.02*range(Env), char(txt), 'HorizontalAlignment','center', 'FontWeight','bold');
end
xlabel(ax_s1s2, 'Tiempo (s)'); ylabel(ax_s1s2, 'Amplitud normalizada');
title(ax_s1s2, sprintf('Clasificación S1/S2 (BPM~%.1f)', bpm));
% LA LEYENDA SE MUEVE AL FINAL
% legend('Señal original','Envolvente','Triángulos','S1','S2','Location','best');

% --- 7) Información de depuración opcional ---
fprintf('Detectadas: %d S1, %d S2\n', sum(labels=="S1"), sum(labels=="S2"));

smooth_ms = 15;             % <-- ASEGÚRATE DE QUE ESTA LÍNEA EXISTA
search_pre  = 0.12;         % segundos antes del pico S1 para buscar (ajusta 0.05–0.15)
search_post = 0.03;         % pequeño margen después del pico
min_prominence = 0.01;      % prominencia mínima (ajustable)
min_depth_rel = 0.10;       % profundidad mínima relativa (10%)

% ... (definición de Fs, smooth_ms, etc.) ...
% --- Suavizar envolvente ---
smooth_n = max(3, round(fs * (smooth_ms/1000)));
% Env_s = smooth(Env, smooth_n); % REQUIERE CURVE FITTING TOOLBOX
Env_s = movmean(Env, smooth_n);   % Alternativa (MATLAB >= R2016a)
if length(Env_s) ~= length(Env)  % Asegurar compatibilidad
     Env_s = movmean(Env, smooth_n, 'Endpoints','fill');
end

% === Aquí comienza el código que falta ===

% --- Pre-alocar variables para los mínimos ---
S1_min_time = nan(size(S1_idx));
S1_min_amp  = nan(size(S1_idx));

% --- Bucle principal para encontrar mínimos izquierdos ---
for kk = 1:length(S1_idx)
    idxPeak = S1_idx(kk);
    tS1 = peaks_time(idxPeak);
    pAmp = peaks_amp(idxPeak); % Amplitud del pico (en Env original)

    % === Buscar ANTES del pico (mínimo izquierdo) ===
    i_start = find(t >= max(t(1), tS1 - search_pre), 1, 'first');
    i_end   = find(t <= min(t(end), tS1 + search_post), 1, 'last');
    
    if isempty(i_start) || isempty(i_end) || i_end <= i_start
        continue;
    end
    
    seg_t = t(i_start:i_end);
    seg_env = Env_s(i_start:i_end); % Usar la envolvente SUAVIZADA (Env_s)
    
    % Buscar picos en -Env (equivale a mínimos)
    % NOTA: Requiere Signal Processing Toolbox
    [negPeaks, negLocs] = findpeaks(-seg_env, 'MinPeakProminence', min_prominence);
    
    if ~isempty(negLocs)
        locs_time = seg_t(negLocs);
        % Elegir el último mínimo antes del pico S1 (izquierdo)
        before_mask = locs_time <= tS1;
        
        if any(before_mask)
            candIdx = find(before_mask);
            chosen = candIdx(end); % el más cercano antes del pico
            min_idx = i_start + negLocs(chosen) - 1;
        else
            % Si no hay ninguno antes (raro), usar el mínimo global de la ventana
            [~, im] = min(seg_env);
            min_idx = i_start + im - 1;
        end
        
        % Guardar tiempo (t) y amplitud de la envolvente ORIGINAL (Env)
        S1_min_time(kk) = t(min_idx);
        S1_min_amp(kk)  = Env(min_idx); 
    end
end
% === Aquí termina el código que faltaba ===

% --- Graficar ---
% Usa el handle 'ax_s1s2' guardado para dibujar en la figura correcta
plot(ax_s1s2, S1_min_time, S1_min_amp, 'mv', 'MarkerFaceColor','m', 'MarkerSize',8);
for jj = 1:length(S1_min_time)
    if ~isnan(S1_min_time(jj))
        text(ax_s1s2, S1_min_time(jj), S1_min_amp(jj) - 0.02*range(Env), ...
            'S1minL', 'Color','m', 'FontWeight','bold', ...
            'HorizontalAlignment','center','VerticalAlignment','top');
    end
end

% --- ACTUALIZAR LEYENDA (AHORA SÍ) ---
legend(ax_s1s2, 'Señal original','Envolvente','Triángulos','S1','S2','S1minL','Location','best');
hold(ax_s1s2, 'off'); % <- Buena práctica

fprintf('Mínimos izquierdos detectados: %d de %d S1\n', ...
        sum(~isnan(S1_min_time)), length(S1_idx));


% === 5.4 Visual: sombrear ciclos en la envolvente ===
figure('Name','Ciclos cardiacos detectados');
plot(t, Env, 'k'); grid on; hold on;
for i = 1:size(ciclos_ref,1)
    i1 = ciclos_ref(i,1); i2 = ciclos_ref(i,2);
    fill([t(i1) t(i2) t(i2) t(i1)], [0 0 1 1], ...
         'c', 'FaceAlpha', 0.18, 'EdgeColor','none');
end
legend('Envolvente','Ciclos'); xlabel('Tiempo (s)'); ylabel('Amplitud norm.');
title('Ciclos cardiacos (propuestos)');

% --- Comprobación final antes de guardar ---

% Comprobar que 'ciclos_ref' no esté vacío Y tenga 2 columnas
if ~isempty(ciclos_ref) && size(ciclos_ref, 2) == 2
    ciclos_t = [ t(ciclos_ref(:,1)) , t(ciclos_ref(:,2)) ];
    disp('>> Segmentación completada. Guardado: resultados_segmentacion.mat');
else
    % Si no se encontraron ciclos, 'ciclos_t' se crea como vacío
    ciclos_t = []; 
    disp('>> ADVERTENCIA: No se detectaron ciclos cardíacos válidos. Guardando resultados vacíos.');
end

% Guardar las variables (ciclos_t y ciclos_ref estarán vacíos si no se encontró nada)
save resultados_segmentacion.mat fs t x Env ciclos_ref ciclos_t
