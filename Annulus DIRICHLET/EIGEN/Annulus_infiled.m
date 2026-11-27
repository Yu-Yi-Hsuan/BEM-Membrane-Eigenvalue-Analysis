%% ==========================================================
% 環形邊界 F15.DAT 產生程式
% 外圓逆時針、內圓順時針，並繪製邊界佈點圖
% MATLAB R2024a
%% ==========================================================

clear; clc; close all;

%% ===== 使用者設定 =====
outputFile = 'F15.DAT';

R_outer = 1;      % 外圓半徑
R_inner = 0.43495;      % 內圓半徑

nOuter = 130;        % 外圓節點數
nInner = 57;        % 內圓節點數

centerX = 0.0;      % 圓心 x 座標
centerY = 0.0;      % 圓心 y 座標

%% ===== 參數檢查 =====
if R_inner <= 0 || R_outer <= R_inner
    error('需滿足 0 < R_inner < R_outer。');
end

if nOuter < 3 || nInner < 3
    error('內、外圓節點數皆須至少為 3。');
end

%% ===== 建立環形邊界節點 =====
% 外圓：逆時針
thetaOuter = 2*pi*(0:nOuter-1)/nOuter;
xOuter = centerX + R_outer*cos(thetaOuter);
yOuter = centerY + R_outer*sin(thetaOuter);

% 內圓：順時針
thetaInner = -2*pi*(0:nInner-1)/nInner;
xInner = centerX + R_inner*cos(thetaInner);
yInner = centerY + R_inner*sin(thetaInner);

x = [xOuter, xInner];
y = [yOuter, yInner];
nNode = nOuter + nInner;

%% ===== 輸出 F15.DAT =====
fid = fopen(outputFile, 'w');

if fid == -1
    error('無法建立 %s。', outputFile);
end

% 節點資料
fprintf(fid, '%6d\n', -1);
fprintf(fid, '%6d\n', 15);

for i = 1:nNode
    fprintf(fid, '%10d%10d%10d%10d%13.5E%13.5E%13.5E\n', ...
        i, 0, 0, 11, x(i), y(i), 0.0);
end

% 元素資料
fprintf(fid, '%6d\n', -1);
fprintf(fid, '%6d\n', -1);
fprintf(fid, '%6d\n', 71);

% 外圓元素
for i = 1:nOuter
    nextNode = mod(i, nOuter) + 1;

    fprintf(fid, '%10d%10d%10d%10d%10d%10d%10d\n', ...
        i, 1, 21, 1, 1, 7, 2);
    fprintf(fid, '%10d%10d\n', i, nextNode);
end

% 內圓元素
for j = 1:nInner
    currentNode = nOuter + j;
    nextNode = nOuter + mod(j, nInner) + 1;
    elementID = nOuter + j;

    fprintf(fid, '%10d%10d%10d%10d%10d%10d%10d\n', ...
        elementID, 1, 21, 1, 1, 7, 2);
    fprintf(fid, '%10d%10d\n', currentNode, nextNode);
end

fprintf(fid, '%6d\n', -1);
fclose(fid);

%% ===== 繪製邊界佈點圖 =====
figure;
hold on;
axis equal;
box on;
grid on;

% 合併成同一組 plot，使內、外圓使用相同的 MATLAB 預設顏色。
xPlot = [xOuter, xOuter(1), NaN, xInner, xInner(1)];
yPlot = [yOuter, yOuter(1), NaN, yInner, yInner(1)];

plot(xPlot, yPlot, '-o', ...
    'MarkerSize', 4, ...
    'LineWidth', 1);

xlabel('x');
ylabel('y');
title('同心圓邊界佈點圖');

xlim(centerX + 1.15*R_outer*[-1, 1]);
ylim(centerY + 1.15*R_outer*[-1, 1]);

%% ===== 顯示結果 =====
fprintf('%s 已成功建立。\n', outputFile);
fprintf('總節點數 = %d\n', nNode);
fprintf('總元素數 = %d\n', nNode);
