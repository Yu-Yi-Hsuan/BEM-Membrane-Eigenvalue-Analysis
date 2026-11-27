%% ==============================================================
% 由 F15.DAT 的環形邊界點產生 1F80.DAT 場內點
%
% 流程：
% 1. 讀取 F15 邊界座標
% 2. 自動分離外圓與內圓
% 3. 建立規則網格並以遮罩保留環形區域
% 4. 輸出 1F80.DAT
% 5. 繪製場內點與邊界
%
% MATLAB R2024a
%% ==============================================================

clear; clc; close all;

%% ===== 使用者設定 =====
inputFile  = 'F15.DAT';
outputFile = '1F80.DAT';
nGrid      = 100;           % x、y 方向網格點數

%% ===== 讀取 F15 邊界點 =====
fid = fopen(inputFile, 'r');

if fid == -1
    error('無法開啟 %s。', inputFile);
end

% 前兩行為控制資料
fgetl(fid);
fgetl(fid);

xb = [];
yb = [];

while ~feof(fid)
    lineText = strtrim(fgetl(fid));

    if isempty(lineText)
        continue;
    end

    % 第一段節點資料以 -1 結束
    if strcmp(lineText, '-1')
        break;
    end

    data = sscanf(lineText, '%f').';

    % F15 第 5、6 欄為 x、y 座標
    if numel(data) >= 6
        xb(end+1,1) = data(5); %#ok<AGROW>
        yb(end+1,1) = data(6); %#ok<AGROW>
    end
end

fclose(fid);

if numel(xb) < 6
    error('邊界點不足，請檢查 F15.DAT。');
end

%% ===== 分離外圓與內圓 =====
% F15 節點依序為：外圓節點 → 內圓節點
% 兩組邊界交界處的相鄰點距離最大
jumpDistance = hypot(diff(xb), diff(yb));
[~, splitIndex] = max(jumpDistance);

xBoundary1 = xb(1:splitIndex);
yBoundary1 = yb(1:splitIndex);

xBoundary2 = xb(splitIndex+1:end);
yBoundary2 = yb(splitIndex+1:end);

if numel(xBoundary1) < 3 || numel(xBoundary2) < 3
    error('無法分離內、外圓邊界，請檢查 F15.DAT 節點排列。');
end

% 以平均半徑判斷內、外圓
rBoundary1 = mean(hypot( ...
    xBoundary1-mean(xBoundary1), ...
    yBoundary1-mean(yBoundary1)));

rBoundary2 = mean(hypot( ...
    xBoundary2-mean(xBoundary2), ...
    yBoundary2-mean(yBoundary2)));

if rBoundary1 > rBoundary2
    xOuter = xBoundary1;
    yOuter = yBoundary1;
    xInner = xBoundary2;
    yInner = yBoundary2;
else
    xOuter = xBoundary2;
    yOuter = yBoundary2;
    xInner = xBoundary1;
    yInner = yBoundary1;
end

%% ===== 建立規則網格 =====
xVector = linspace(min(xOuter), max(xOuter), nGrid);
yVector = linspace(min(yOuter), max(yOuter), nGrid);

[X, Y] = meshgrid(xVector, yVector);

%% ===== 建立環形區域遮罩 =====
[inOuter, onOuter] = inpolygon(X, Y, xOuter, yOuter);
[inInner, onInner] = inpolygon(X, Y, xInner, yInner);

% 保留外圓內部，排除內圓與兩條邊界
inside = inOuter & ~onOuter & ~inInner & ~onInner;

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

fprintf('已讀取 %d 個邊界點。\n', numel(xb));
fprintf('外圓邊界點數 = %d\n', numel(xOuter));
fprintf('內圓邊界點數 = %d\n', numel(xInner));
fprintf('共產生 %d 個環形場內點。\n', nField);
fprintf('已輸出 %s。\n', outputFile);

%% ===== 繪製場內點分布圖 =====
figure;
plot(xin, yin, '.', 'MarkerSize', 7);
hold on;

plot([xOuter; xOuter(1)], [yOuter; yOuter(1)], ...
    '-', 'LineWidth', 1.5);

plot([xInner; xInner(1)], [yInner; yInner(1)], ...
    '-', 'LineWidth', 1.5);

plot(xb, yb, 'o', ...
    'MarkerSize', 4, ...
    'LineWidth', 1);

axis equal;
axis tight;
grid on;
box on;

xlabel('x');
ylabel('y');

title(sprintf( ...
    '環形 1F80 場內點分布（%d × %d 網格，共 %d 點）', ...
    nGrid, nGrid, nField));

legend('1F80 場內點', ...
       '外圓邊界', ...
       '內圓邊界', ...
       'F15 邊界點', ...
       'Location', 'best');