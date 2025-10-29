%% Beta-wise Force–Displacement Curves
% File: beta_stiffness.xlsx
% Sheets: Sheet1 -> vertical, Sheet2 -> longitudinal, Sheet3 -> lateral
% Columns: for each beta header ('0 degree','10 degree','20 degree','30 degree')
%          the next two columns are [Displacement, Force] pairs.
% This script plots Displacement (X) vs Force (Y) curves for each beta.

clear; clc;

excelFile = 'beta_stiffness.xlsx';

% ---- Sheet mapping ----
sheetInfo = struct('name', {'vertical','longitudinal','lateral'}, 'idx', {1,2,3});

% ---- Beta labels to detect in the sheet headers (robust to spaces/case) ----
betaHeadersRaw = {'0 degree','10 degree','20 degree','30 degree'};
betaLegend     = {'\beta=0^\circ','\beta=10^\circ','\beta=20^\circ','\beta=30^\circ'};
betaKeys       = {'b0','b10','b20','b30'}; % struct field keys

% ---- Options ----
showMarkers = true;       % 선 위에 작은 점 표시
mark_U1 = true; U1 = 5;   % longitudinal plot에서 U1=5 mm 위치 마커
mark_U2 = true; U2 = 5;   % lateral plot에서  U2=5 mm 위치 마커

colors = lines(numel(betaKeys));

% ---------------- Parse all three sheets ----------------
Data = struct();
for i = 1:numel(sheetInfo)
    Data.(sheetInfo(i).name) = parse_beta_blocks(excelFile, sheetInfo(i).idx, betaHeadersRaw);
end

% ---------------- Plot curves ----------------
figure('Name','Beta: Force–Displacement Curves','Color','w');

dirs  = {'vertical','longitudinal','lateral'};
titles= {'Vertical','Longitudinal','Lateral'};

for i = 1:numel(dirs)
    subplot(1,3,i); hold on; grid on; box on;
    leg = {};
    for k = 1:numel(betaKeys)
        if isfield(Data.(dirs{i}), betaKeys{k})
            d = Data.(dirs{i}).(betaKeys{k});
            % 정렬 후 라인 플롯 (X=Displacement, Y=Force)
            [x, idx] = sort(d.disp(:));
            y = d.force(:); y = y(idx);
            if showMarkers
                plot(x, y, 'LineWidth', 2.0, 'Color', colors(k,:), ...
                     'Marker','.', 'MarkerSize',10);
            else
                plot(x, y, 'LineWidth', 2.0, 'Color', colors(k,:));
            end
            leg{end+1} = betaLegend{k}; %#ok<SAGROW>
        end
    end
    xlabel('Displacement (mm)'); ylabel('Force (N)');
    title([titles{i} ' F-\delta']);
    if ~isempty(leg), legend(leg, 'Location','best'); end

    % ---- lateral subplot만 y축 0부터 시작 ----
    if strcmp(dirs{i}, 'lateral')
        ylim([0 inf]);   % y축 하한 0, 상한 자동
    end

    % ---- longitudinal, lateral subplot에서 x축 [0,5 mm] 제한 ----
    if any(strcmp(dirs{i}, {'longitudinal','lateral'}))
        xlim([0 5]);
    end
end

% ---------------- Optional markers at U1/U2 = 5 mm ----------------
% Longitudinal: mark U1
if mark_U1 && isfield(Data, 'longitudinal')
    subplot(1,3,2); % longitudinal panel
    for k = 1:numel(betaKeys)
        if isfield(Data.longitudinal, betaKeys{k})
            d = Data.longitudinal.(betaKeys{k});
            [x, idx] = sort(d.disp(:));
            y = d.force(:); y = y(idx);
            if U1 >= min(x) && U1 <= max(x)
                F5 = interp1(x, y, U1, 'linear');
                plot(U1, F5, 'o', 'MarkerSize', 7, ...
                     'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k');
            end
        end
    end
end

% Lateral: mark U2
if mark_U2 && isfield(Data, 'lateral')
    subplot(1,3,3); % lateral panel
    for k = 1:numel(betaKeys)
        if isfield(Data.lateral, betaKeys{k})
            d = Data.lateral.(betaKeys{k});
            [x, idx] = sort(d.disp(:));
            y = d.force(:); y = y(idx);
            if U2 >= min(x) && U2 <= max(x)
                F5 = interp1(x, y, U2, 'linear');
                plot(U2, F5, 's', 'MarkerSize', 7, ...
                     'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k');
            end
        end
    end
end

%% ================== Helper functions ==================
function D = parse_beta_blocks(xfile, sheetIdx, betaHeadersRaw)
    % Read mixed-type cells
    C = readcell(xfile, 'Sheet', sheetIdx);
    [nRow, nCol] = size(C);

    % Find beta headers within top rows
    headers = struct();
    maxHeaderRow = min(8, nRow);
    for r = 1:maxHeaderRow
        for c = 1:nCol
            if ischar(C{r,c}) || isstring(C{r,c})
                lbl = strtrim(string(C{r,c}));
                if is_beta_header(lbl, betaHeadersRaw)
                    key = header_to_key_beta(lbl); % '0 degree' -> 'b0', etc.
                    headers.(key) = [r c];
                end
            end
        end
    end

    % Collect numeric pairs [Displacement, Force] for each beta block
    D = struct();
    hnames = fieldnames(headers);
    for i = 1:numel(hnames)
        hk  = hnames{i};         % e.g., 'b0','b10','b20','b30'
        pos = headers.(hk);      % [row col] of the header
        cD = pos(2);             % column of Displacement
        cF = pos(2) + 1;         % next column: Force

        % Find first numeric row (skip labels)
        firstDataRow = NaN;
        for r = pos(1)+1:nRow
            if is_num(C{r,cD}) && is_num(C{r,cF})
                firstDataRow = r; break;
            end
        end
        if isnan(firstDataRow), firstDataRow = min(pos(1)+2, nRow); end

        % Collect numeric pairs
        DD = []; FF = [];
        for r = firstDataRow:nRow
            if is_num(C{r,cD}) && is_num(C{r,cF})
                DD(end+1,1) = double(C{r,cD}); %#ok<AGROW>
                FF(end+1,1) = double(C{r,cF}); %#ok<AGROW>
            end
        end

        if ~isempty(DD)
            D.(hk).disp  = DD(:).';
            D.(hk).force = FF(:).';
        end
    end
end

function tf = is_beta_header(lbl, headerList)
    s = lower(strtrim(lbl));
    s = regexprep(s, '\s+', ' ');
    s = strrep(s, 'degrees', 'degree');
    if any(strcmp(s, lower(headerList)))
        tf = true; return;
    end
    % Also accept compact forms (e.g., '0deg','0°')
    tf = ~isempty(regexp(s, '^(0|10|20|30)\s*(degree|deg|°)$','once'));
end

function key = header_to_key_beta(lbl)
    s = lower(regexprep(char(lbl), '\s+', ''));
    s = strrep(s,'degrees','degree');
    tok = regexp(s, '^(\d+)(degree|deg|°)$', 'tokens','once');
    if isempty(tok)
        t2 = regexp(lower(char(lbl)), '^(\d+)\s*(degree|deg|°)$', 'tokens','once');
        if ~isempty(t2), tok = t2; end
    end
    if ~isempty(tok)
        deg = tok{1};
        key = ['b' deg];  % 'b0','b10','b20','b30'
    else
        key = matlab.lang.makeValidName(s);
    end
end

function tf = is_num(x)
    tf = isnumeric(x) && isfinite(x);
end
