%% ==========================================================
% 圓形邊界 F15.DAT 產生程式
% 圓形邊界逆時針佈點，並繪製邊界佈點圖
% MATLAB R2024a
%% ==========================================================

clear; clc; close all;

%% ===== 使用者設定 =====
outputFile = 'F15.DAT';

R = 0.4938;            % 圓半徑
nNode = 40;         % 邊界節點數，也是邊界元素數

centerX = 0.0;      % 圓心 x 座標
centerY = 0.0;      % 圓心 y 座標

%% ===== 參數檢查 =====
if R <= 0
    error('圓半徑 R 必須大於 0。');
end

if nNode < 3 || fix(nNode) ~= nNode
    error('邊界節點數 nNode 必須是大於或等於 3 的整數。');
end

%% ===== 建立圓形邊界節點 =====
% 圓形外邊界：逆時針
theta = 2*pi*(0:nNode-1)/nNode;

x = centerX + R*cos(theta);
y = centerY + R*sin(theta);

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

for i = 1:nNode
    nextNode = mod(i, nNode) + 1;

    fprintf(fid, '%10d%10d%10d%10d%10d%10d%10d\n', ...
        i, 1, 21, 1, 1, 7, 2);

    fprintf(fid, '%10d%10d\n', i, nextNode);
end

fprintf(fid, '%6d\n', -1);
fclose(fid);

%% ===== 繪製邊界佈點圖 =====
figure;
hold on;
axis equal;
box on;
grid on;

% 將第一個節點接回最後一個節點，使圓形邊界閉合
xPlot = [x, x(1)];
yPlot = [y, y(1)];

plot(xPlot, yPlot, '-o', ...
    'MarkerSize', 4, ...
    'LineWidth', 1);

xlabel('x');
ylabel('y');
title('圓形邊界佈點圖');

xlim(centerX + 1.15*R*[-1, 1]);
ylim(centerY + 1.15*R*[-1, 1]);

%% ===== 顯示結果 =====
fprintf('%s 已成功建立。\n', outputFile);
fprintf('圓半徑 = %.6g\n', R);
fprintf('總節點數 = %d\n', nNode);
fprintf('總元素數 = %d\n', nNode);
fprintf('節點方向 = 逆時針\n');
