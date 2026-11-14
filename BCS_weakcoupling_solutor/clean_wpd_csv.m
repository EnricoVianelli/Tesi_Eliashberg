function clean_wpd_csv(input_csv, output_csv, T_list)
% Converte il CSV multi-colonna di WebPlotDigitizer nel formato [T, nu, sigma]
% Esempio d'uso:
% clean_wpd_csv('raw_sigma1.csv', 'Real_cond_MgB2.csv', [3, 15, 20, 27]);

% 1. Legge i dati numerici saltando le prime 2 righe di intestazione
raw_data = readmatrix(input_csv, 'NumHeaderLines', 2);

T_col     = [];
nu_col    = [];
sigma_col = [];

% 2. Estrae ciascuna coppia di colonne (X, Y) per ogni temperatura
for i = 1:numel(T_list)
    col_x = 2*i - 1;
    col_y = 2*i;

    nu_i    = raw_data(:, col_x);
    sigma_i = raw_data(:, col_y);

    % Rimuove le righe vuote (NaN) dovute alle diverse lunghezze delle serie
    valid = ~isnan(nu_i) & ~isnan(sigma_i);
    nu_i    = nu_i(valid);
    sigma_i = sigma_i(valid);

    % Ordina i punti per frequenza (nu) crescente
    [nu_i, sort_idx] = sort(nu_i);
    sigma_i = sigma_i(sort_idx);

    % Accoda i dati puliti
    T_col     = [T_col; repmat(T_list(i), numel(nu_i), 1)];
    nu_col    = [nu_col; nu_i];
    sigma_col = [sigma_col; sigma_i];
end

% 3. Crea la tabella ordinata e la salva nel nuovo CSV
T_clean = table(T_col, nu_col, sigma_col, 'VariableNames', {'T', 'nu', 'sigma'});
writetable(T_clean, output_csv);

fprintf('File salvato con successo in "%s" (%d punti totali).\n', output_csv, height(T_clean));
end
