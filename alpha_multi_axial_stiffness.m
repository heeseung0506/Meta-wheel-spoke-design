%% Curve & stiffness from Excel by sheet index (FIXED field names)
% Sheet 1 -> vertical, Sheet 2 -> longitudinal, Sheet 3 -> lateral

clear; clc;

excelFile = 'stiffness.xlsx';
sheetInfo = struct('name', {'vertical','longitudinal','lateral'}, 'idx', {1,2,3});

rlabels = {'r=1','r=3','r=4','r=6'};   % 표시용 라벨(legend 등)
rkShort  = {'r1','r3','r4','r6'};      % 구조체 필드명(합법)
colors   = lines(numel(rlabels));

% ------- 시트 파싱 -------
Data = struct();
for i = 1:numel(sheetInfo)                                             
    dirName = sheetInfo(i).name;
    Data.(dirName) = parse_r_blocks(excelFile, sheetInfo(i).idx, rlabels);
end

% ------- Force–Displacement 곡선 (축 반대로: Force = X, Disp = Y) -------
figure('Name','Force–Displacement Curves (Reversed Axes)','Color','w');
dirs = {'vertical','longitudinal','lateral'};
for i = 1:numel(dirs)
    subplot(1,3,i); hold on; grid on; box on;
    legends = {};
    for k = 1:numel(rkShort)
        if isfield(Data.(dirs{i}), rkShort{k})
            d = Data.(dirs{i}).(rkShort{k});
            % ⬇⬇⬇ force가 X축, disp가 Y축
            plot(d.force, d.disp, 'LineWidth', 1.8, 'Color', colors(k,:));
            legends{end+1} = rlabels{k}; %#ok<SAGROW>
        end
    end
    ylabel('Force (N)'); xlabel('Displacement (mm)');
    title([upper(dirs{i}(1)) dirs{i}(2:end) ' F-\delta']); % δ–F
    if ~isempty(legends), legend(legends, 'Location','best'); end

    % ---- longitudinal, lateral은 x축 0~5 mm 제한 ----
    if any(strcmp(dirs{i}, {'longitudinal','lateral'}))
        xlim([0 5]);
    end
end

% ------- 강성 계산 (변경 없음: 계산은 Disp→Force 기울기) -------
frac = 0.30;  % vertical/longitudinal 초기 구간 비율
S = struct();

% vertical & longitudinal: 초기 기울기
for dn = {'vertical','longitudinal'}
    dname = dn{1};
    for k = 1:numel(rkShort)
        if isfield(Data.(dname), rkShort{k})
            dd = Data.(dname).(rkShort{k});
            S.(dname).(rkShort{k}) = initial_slope(dd.disp, dd.force, frac);
        end
    end
end

% lateral: 등가 강성 ΔF/Δδ
for k = 1:numel(rkShort)
    if isfield(Data.lateral, rkShort{k})
        dd = Data.lateral.(rkShort{k});
        S.lateral.(rkShort{k}) = equiv_slope(dd.disp, dd.force);
    end
end

% ------- 강성 막대그래프 -------
figure('Name','Stiffness Summary','Color','w');

subplot(1,3,1); hold on; grid on; box on;
vals = ordered_vals(S, 'vertical', rkShort);
bar(categorical(rlabels(1:numel(vals))), vals, 'LineWidth', 1);
title('Vertical stiffness (N/mm)');

subplot(1,3,2); hold on; grid on; box on;
vals = ordered_vals(S, 'longitudinal', rkShort);
bar(categorical(rlabels(1:numel(vals))), vals, 'LineWidth', 1);
title('Longitudinal stiffness (N/mm)');

subplot(1,3,3); hold on; grid on; box on;
vals = ordered_vals(S, 'lateral', rkShort);
bar(categorical(rlabels(1:numel(vals))), vals, 'LineWidth', 1);
title('Lateral equivalent stiffness (N/mm)');

% ------- 콘솔 요약 -------
fprintf('\n===== Stiffness Summary (N/mm) =====\n');
print_line('Vertical',        S, 'vertical',     rkShort);
print_line('Longitudinal',    S, 'longitudinal', rkShort);
print_line('Lateral (ΔF/Δδ)', S, 'lateral',      rkShort);

%% ================== Helper functions ==================
function D = parse_r_blocks(xfile, sheetIdx, rlabels)
    C = readcell(xfile, 'Sheet', sheetIdx);
    [nRow, nCol] = size(C);

    % 상단부에서 r=1, r=3, r=4, r=6 헤더 탐색
    headers = struct();
    maxHeaderRow = min(8, nRow);
    for r = 1:maxHeaderRow
        for c = 1:nCol
            if ischar(C{r,c}) || isstring(C{r,c})
                label = strtrim(string(C{r,c}));
                if any(strcmp(label, rlabels))
                    key = header_to_key(label);   % 'r=1' -> 'r1'
                    headers.(key) = [r c];        % 유효 필드명으로 저장
                end
            end
        end
    end

    % 각 헤더 위치에서 Force/Displacement 열 수집
    D = struct();
    hnames = fieldnames(headers);
    for i = 1:numel(hnames)
        hk  = hnames{i};               % 예: 'r1'
        pos = headers.(hk);            % [row col]
        cF = pos(2); cD = pos(2) + 1;  % Force, Displacement 예상 열

        % 실제 데이터 시작 행(두 열 모두 숫자인 첫 행)
        firstDataRow = NaN;
        for r = pos(1)+1:nRow
            if is_num(C{r,cF}) && is_num(C{r,cD})
                firstDataRow = r; break;
            end
        end
        if isnan(firstDataRow)
            firstDataRow = min(pos(1)+2, nRow);
        end

        F = []; Dd = [];
        for r = firstDataRow:nRow
            if is_num(C{r,cF}) && is_num(C{r,cD})
                F(end+1,1)  = double(C{r,cF}); %#ok<AGROW>
                Dd(end+1,1) = double(C{r,cD}); %#ok<AGROW>
            end
        end

        if ~isempty(F)
            D.(hk).force = F(:).';
            D.(hk).disp  = Dd(:).';
        end
    end
end

function key = header_to_key(label)
    % 'r=1', 'r = 3', 'r=4mm' 등 변형을 안전하게 'r1','r3','r4'로 변환
    s = regexprep(char(label), '\s+', '');          % 공백 제거
    tok = regexp(s, '^r=([\d\.]+)', 'tokens','once');
    if ~isempty(tok)
        numpart = tok{1};
        numpart = regexprep(numpart,'[^0-9]','');   % 정수만 사용
        if isempty(numpart), numpart = '0'; end
        key = ['r' numpart];
    else
        % 최후 수단: 유효한 필드명으로 정규화
        key = matlab.lang.makeValidName(s);
    end
end

function tf = is_num(x)
    tf = isnumeric(x) && isfinite(x);
end

function s = initial_slope(disp, force, frac)
    % Disp→Force 기울기(초기 구간) 계산용
    disp = disp(:); force = force(:);
    if numel(disp) < 2 || numel(force) < 2
        s = NaN; return;
    end
    [disp, idx] = sort(disp); force = force(idx);
    dmin = disp(1); dmax = disp(end);
    cutoff = dmin + frac * (dmax - dmin);
    sel = disp <= cutoff;
    if nnz(sel) < 2, sel(1:min(2,numel(disp))) = true; end
    p = polyfit(disp(sel), force(sel), 1);
    s = p(1);
end

function s = equiv_slope(disp, force)
    disp = disp(:); force = force(:);
    dD = max(disp) - min(disp);
    dF = max(force) - min(force);
    if dD == 0, s = NaN; else, s = dF / dD; end
end

function vals = ordered_vals(S, fieldName, rkShort)
    vals = nan(1, numel(rkShort));
    if isfield(S, fieldName)
        for k = 1:numel(rkShort)
            if isfield(S.(fieldName), rkShort{k})
                vals(k) = S.(fieldName).(rkShort{k});
            end
        end
    end
end

function print_line(titleStr, S, fieldName, rkShort)
    labs = {'r=1','r=3','r=4','r=6'};
    vals = ordered_vals(S, fieldName, rkShort);
    fprintf('%-18s | ', titleStr);
    for k = 1:numel(vals)
        if isnan(vals(k)), fprintf('%s: %8s | ', labs{k}, 'NaN');
        else,              fprintf('%s: %8.4f | ', labs{k}, vals(k));
        end
    end
    fprintf('\n');
end
