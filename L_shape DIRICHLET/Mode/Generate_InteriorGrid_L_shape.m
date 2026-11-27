%% ==============================================================
% 由 F15.DAT L 型域邊界產生 1F80.DAT 場內點
% MATLAB R2024a
%% ==============================================================

clear; clc; close all;

%% ===== 使用者設定 =====
inputFile  = 'F15.DAT';
outputFile = '1F80.DAT';

nGrid = 100;    % x、y 方向網格點數

%% ===== 讀取 F15 邊界座標 =====
fid = fopen(inputFile, 'r');

if fid == -1
    error('無法開啟 %s。', inputFile);
end

% 跳過 F15 前兩行控制資料
fgetl(fid);
fgetl(fid);

xb = [];
yb = [];

while ~feof(fid)

    lineText = strtrim(fgetl(fid));

    if isempty(lineText)
        continue;
    end

    % 節點資料以 -1 結束
    if strcmp(lineText, '-1')
        break;
    end

    data = sscanf(lineText, '%f');

    % F15 第 5、6 欄為 x、y 座標
    if numel(data) >= 6
        xb(end+1,1) = data(5); %#ok<AGROW>
        yb(end+1,1) = data(6); %#ok<AGROW>
    end

end

fclose(fid);

if numel(xb) < 3
    error('邊界點不足，請檢查 F15.DAT。');
end

%% ===== 建立規則網格 =====
xGrid = linspace(min(xb), max(xb), nGrid);
yGrid = linspace(min(yb), max(yb), nGrid);

[X, Y] = meshgrid(xGrid, yGrid);

%% ===== 保留L 型域內點 =====
[inDomain, onBoundary] = inpolygon(X, Y, xb, yb);

% 排除邊界點，只保留場內點
inside = inDomain & ~onBoundary;

fieldPoints = sortrows([X(inside), Y(inside)], [2, 1]);

xin = fieldPoints(:,1);
yin = fieldPoints(:,2);

nField = size(fieldPoints,1);

%% ===== 輸出 1F80.DAT =====
fid = fopen(outputFile, 'w');

if fid == -1
    error('無法建立 %s。', outputFile);
end

for i = 1:nField

    % Fortran FORMAT(I6,2E13.5)
    fprintf(fid, '%6d%13.5E%13.5E\r\n', ...
        i, xin(i), yin(i));

end

fclose(fid);

%% ===== 繪製場內點分布 =====
figure;
hold on;

plot(xin, yin, '.', ...
    'MarkerSize', 7);

plot([xb; xb(1)], [yb; yb(1)], ...
    '-', 'LineWidth', 1.5);

plot(xb, yb, 'o', ...
    'MarkerSize', 4, ...
    'LineWidth', 1);

axis equal;
axis tight;

box on;
grid on;

xlabel('x');
ylabel('y');

title(sprintf( ...
    'L 型域 1F80 場內點分布（%d × %d 網格，共 %d 點）', ...
    nGrid, nGrid, nField));

legend('1F80 場內點', ...
       'L 型邊界', ...
       'F15 邊界點', ...
       'Location', 'best');

%% ===== 顯示結果 =====
fprintf('F15 邊界點數 = %d\n', numel(xb));
fprintf('1F80 場內點數 = %d\n', nField);
fprintf('%s 已成功建立。\n', outputFile);