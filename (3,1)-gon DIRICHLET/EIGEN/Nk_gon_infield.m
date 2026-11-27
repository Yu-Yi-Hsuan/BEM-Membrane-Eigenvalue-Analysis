%% ==========================================================
% (N,k)-gon 星形多邊形邊界 F15.DAT 產生程式
% MATLAB R2024a
%% ==========================================================

clear; clc; close all;

%% ===== 使用者設定 =====
outputFile = 'F15.DAT';

N = 3;              % 外尖點數
k = 1;              % 星形跨步參數
R = 1.0;           % 外尖點至中心距離(外接圓半徑)
nNode = 60;        % 邊界總節點數

%% ===== 參數檢查 =====
if N < 3 || N ~= floor(N)
    error('N 必須為大於等於 3 的整數。');
end

if k < 1 || k ~= floor(k)
    error('k 必須為大於等於 1 的整數。');
end

if k > floor((N-1)/2)
    error('需滿足 k <= floor((N-1)/2)。');
end

if R <= 0
    error('R 必須大於 0。');
end

if mod(nNode, 2*N) ~= 0
    error('nNode 必須可以被 2*N 整除。');
end

nPerSeg = nNode / (2*N);

%% ===== 建立外尖點 P =====
theta = pi/2 + 2*pi*(0:N-1)/N;

Px = R*cos(theta);
Py = R*sin(theta);

%% ===== 建立凹點 Q =====
Qx = zeros(1, N);
Qy = zeros(1, N);

if k == 1
    for i = 1:N
        iNext = mod(i, N) + 1;

        Qx(i) = 0.5*(Px(i) + Px(iNext));
        Qy(i) = 0.5*(Py(i) + Py(iNext));
    end
else
    for i = 1:N
        i1 = i;
        i2 = mod(i-1+k, N) + 1;
        i3 = mod(i, N) + 1;
        i4 = mod(i-k+N, N) + 1;

        x1 = Px(i1);
        y1 = Py(i1);

        dx1 = Px(i2) - Px(i1);
        dy1 = Py(i2) - Py(i1);

        x3 = Px(i3);
        y3 = Py(i3);

        dx2 = Px(i4) - Px(i3);
        dy2 = Py(i4) - Py(i3);

        den = dx1*dy2 - dy1*dx2;

        if abs(den) < 1.0e-14
            error('建立 Q(%d) 時兩直線平行。', i);
        end

        s = ((x3-x1)*dy2 - (y3-y1)*dx2) / den;

        Qx(i) = x1 + s*dx1;
        Qy(i) = y1 + s*dy1;
    end
end

%% ===== 建立邊界控制點 =====
Bx = zeros(1, 2*N);
By = zeros(1, 2*N);

for i = 1:N
    Bx(2*i-1) = Px(i);
    By(2*i-1) = Py(i);

    Bx(2*i) = Qx(i);
    By(2*i) = Qy(i);
end

%% ===== 沿邊界均勻產生節點 =====
x = zeros(1, nNode);
y = zeros(1, nNode);

idx = 0;

for i = 1:2*N
    iNext = mod(i, 2*N) + 1;

    for j = 0:nPerSeg-1
        idx = idx + 1;

        t = j / nPerSeg;

        x(idx) = (1-t)*Bx(i) + t*Bx(iNext);
        y(idx) = (1-t)*By(i) + t*By(iNext);
    end
end

%% ===== 輸出 F15.DAT =====
fid = fopen(outputFile, 'w');

if fid == -1
    error('無法建立 %s。', outputFile);
end

fprintf(fid, '%6d\n', -1);
fprintf(fid, '%6d\n', 15);

for i = 1:nNode
    fprintf(fid, '%10d%10d%10d%10d%13.5E%13.5E%13.5E\n', ...
        i, 0, 0, 11, x(i), y(i), 0.0);
end

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

plot([x, x(1)], [y, y(1)], '-o', ...
    'MarkerSize', 4, ...
    'LineWidth', 1);

xlabel('x');
ylabel('y');
title(sprintf('(%d,%d)-gon 邊界佈點圖', N, k));

margin = 0.15*max([max(x)-min(x), max(y)-min(y)]);

xlim([min(x)-margin, max(x)+margin]);
ylim([min(y)-margin, max(y)+margin]);

%% ===== 顯示結果 =====
fprintf('%s 已成功建立。\n', outputFile);
fprintf('(N,k)-gon = (%d,%d)\n', N, k);
fprintf('總節點數 = %d\n', nNode);
