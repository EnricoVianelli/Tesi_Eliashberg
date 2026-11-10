img = imread('Figure4_conductivity_Nb.png');   % your file name
figure; imshow(img); hold on

% ---- calibration ----
disp('Click x = 0 then x = 30 cm^-1 on the bottom axis');   cx = ginput(2);
disp('PANEL a: click y = 0 then y = 2 on the left axis');   ca = ginput(2);
disp('PANEL b: click y = 0 then y = 8 on the left axis');   cb = ginput(2);

px2x  = @(px) 0 + (px - cx(1,1)) * 30 / (cx(2,1) - cx(1,1));
py2s1 = @(py) 0 + (py - ca(1,2)) * 2  / (ca(2,2) - ca(1,2));
py2s2 = @(py) 0 + (py - cb(1,2)) * 8  / (cb(2,2) - cb(1,2));

temps = [7.5 7 6 4.5 9];                      
rows = [];
for panel = 1:2
    for T = temps
        fprintf('Panel %d, T = %.1f K: click all markers, then press Enter\n', panel, T);
        p = ginput;                           % finish with Enter
        if isempty(p), continue, end
        nu = px2x(p(:,1));
        if panel == 1, s = py2s1(p(:,2)); else, s = py2s2(p(:,2)); end
        rows = [rows; repmat(panel,numel(nu),1), repmat(T,numel(nu),1), nu, s]; %#ok<AGROW>
    end
end
raw = array2table(rows, 'VariableNames', {'panel','T','nu','sigma'});
writetable(raw, 'nb_digitized_raw.csv');