%% ==============================================================
% 由圓形 F15.DAT 的邊界點產生 1F80.DAT 場內點
%
% 重要：
% F15 可能只有少量離散邊界點，因此不能直接用 inpolygon 將其
% 視為多邊形來判斷內域。
%
% 本程式會：
% 1. 讀取 F15.DAT 第 5、6 欄的圓周邊界座標
% 2. 由 F15 邊界點自動擬合圓心與半徑
% 3. 在真圓外接正方形中建立 nGrid × nGrid 規則網格
% 4. 直接利用圓方程判斷真圓內部點
% 5. 輸出成 Fortran 可讀取的 1F80.DAT
% 6. 繪製真圓、F15 離散邊界點與場內點
%
% MATLAB R2024a
%% ==============================================================

clear; clc; close all;

%% ===== 可修改設定 =====
inputFile  = 'F15.DAT';    % 輸入的圓形邊界點檔案
outputFile = '1F80.DAT';   % 輸出的場內點檔案
nGrid      = 100;           % x、y 方向各建立多少個候選網格點

nCirclePlot = 500;         % 僅控制畫圖時真圓的平滑程度

%% ===== 讀取 F15 邊界點 =====
fid = fopen(inputFile, 'r');

if fid == -1
    error('無法開啟 %s，請確認檔案與程式位於同一資料夾。', inputFile);
end

% F15 前兩行為控制資料，不是邊界點，因此略過
fgetl(fid);
fgetl(fid);

xb = [];
yb = [];

while ~feof(fid)

    lineText = fgetl(fid);

    if ~ischar(lineText)
        break;
    end

    lineText = strtrim(lineText);

    if isempty(lineText)
        continue;
    end

    % F15 第一段邊界資料以 -1 作為結束標記
    if strcmp(lineText, '-1')
        break;
    end

    data = sscanf(lineText, '%f').';

    % F15 邊界座標位於第 5、6 欄
    if numel(data) >= 6
        xb(end+1,1) = data(5); %#ok<AGROW>
        yb(end+1,1) = data(6); %#ok<AGROW>
    end
end

fclose(fid);

if numel(xb) < 3
    error('讀取到的邊界點不足，請檢查 F15.DAT 格式。');
end

fprintf('已讀取 %d 個 F15 邊界點。\n', numel(xb));

%% ===== 由 F15 邊界點自動擬合真圓 =====
%
% 圓方程：
%   (x-xc)^2 + (y-yc)^2 = R^2
%
% 展開後：
%   x^2 + y^2 + A*x + B*y + C = 0
%
% 利用最小平方法求 A、B、C，再得到圓心與半徑。

M = [xb, yb, ones(size(xb))];
rhs = -(xb.^2 + yb.^2);

coef = M \ rhs;

A = coef(1);
B = coef(2);
C = coef(3);

centerX = -A/2;
centerY = -B/2;

R2 = centerX^2 + centerY^2 - C;

if R2 <= 0
    error('由 F15 邊界點無法得到有效圓半徑，請檢查邊界資料。');
end

R = sqrt(R2);

fprintf('自動擬合圓心 = (%.10f, %.10f)\n', centerX, centerY);
fprintf('自動擬合半徑 R = %.10f\n', R);

%% ===== 建立真圓外接正方形規則網格 =====
xVector = linspace(centerX-R, centerX+R, nGrid);
yVector = linspace(centerY-R, centerY+R, nGrid);

[X, Y] = meshgrid(xVector, yVector);

%% ===== 利用真圓方程判斷場內點 =====
%
% 不再使用：
%   inpolygon(X,Y,xb,yb)
%
% 而是直接使用：
%   (x-xc)^2 + (y-yc)^2 < R^2
%
% 因此即使 F15 只有 10 個邊界點，
% 1F80 的內域仍然是圓，不會變成正十邊形。

radiusSquared = (X-centerX).^2 + (Y-centerY).^2;

inside = radiusSquared < R^2;

xin = X(inside);
yin = Y(inside);

%% ===== 場內點排序 =====
% 先依 y 由小到大，再依 x 由小到大
fieldPoints = sortrows([xin, yin], [2, 1]);

xin = fieldPoints(:,1);
yin = fieldPoints(:,2);

nField = length(xin);

fprintf('共產生 %d 個真圓場內點。\n', nField);

%% ===== 輸出 1F80.DAT =====
fid = fopen(outputFile, 'w');

if fid == -1
    error('無法建立輸出檔案 %s。', outputFile);
end

for i = 1:nField

    % 對應 Fortran：
    % FORMAT(I6,2E13.5)
    fprintf(fid, '%6d%13.5E%13.5E\r\n', ...
        i, xin(i), yin(i));
end

fclose(fid);

fprintf('已輸出 %s。\n', outputFile);

%% ===== 繪製場內點分布圖 =====
figure;
hold on;

% 1F80 真圓場內點
hField = plot(xin, yin, '.', 'MarkerSize', 7);

% 畫真正的圓，不使用 F15 節點直線連接來代表圓
thetaFine = linspace(0, 2*pi, nCirclePlot);

xCircle = centerX + R*cos(thetaFine);
yCircle = centerY + R*sin(thetaFine);

hCircle = plot(xCircle, yCircle, ...
    'k-', 'LineWidth', 1.5);

% F15 離散邊界點
hNode = plot(xb, yb, ...
    'ro', 'MarkerSize', 5, 'LineWidth', 1);

axis equal;

xlim(centerX + 1.10*R*[-1, 1]);
ylim(centerY + 1.10*R*[-1, 1]);

grid on;
box on;

xlabel('x');
ylabel('y');

title(sprintf( ...
    '圓形 1F80 場內點分布（%d × %d 網格，共 %d 點）', ...
    nGrid, nGrid, nField));

legend([hField,hCircle,hNode], ...
    '1F80 場內點', ...
    '真圓邊界', ...
    'F15 邊界點', ...
    'Location', 'best');

%% ===== 最後資訊 =====
fprintf('============================================\n');
fprintf('F15 邊界點數 = %d\n', numel(xb));
fprintf('圓心          = (%.10f, %.10f)\n', centerX, centerY);
fprintf('圓半徑 R      = %.10f\n', R);
fprintf('nGrid         = %d\n', nGrid);
fprintf('1F80 場內點數 = %d\n', nField);
fprintf('內域判斷方式  = 真圓方程\n');
fprintf('============================================\n');
