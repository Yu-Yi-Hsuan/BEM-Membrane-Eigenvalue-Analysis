%% ==============================================================
% 由雙聯通 (N,k)-gon 的 F15.DAT 產生 1F80.DAT 場內點
%
% 適用：
%   外邊界 = 一組 (N,k)-gon
%   內邊界 = 另一組 (N,k)-gon
%
% 流程：
% 1. 讀取 F15.DAT 節點座標
% 2. 讀取元素連接關係，自動辨識外、內兩個閉合邊界
% 3. 在外邊界外接矩形建立 nGrid × nGrid 規則網格
% 4. 保留「外邊界內、內邊界外」的雙聯通區域
% 5. 輸出 Fortran 可讀取的 1F80.DAT
% 6. 繪製場內點與內外 (N,k)-gon 邊界
%
% MATLAB R2024a
%% ==============================================================

clear; clc; close all;

%% ===== 使用者設定 =====
inputFile  = 'F15.DAT';
outputFile = '1F80.DAT';

nGrid = 100;    % x、y 方向各建立多少個候選網格點

%% ==============================================================
% 讀取 F15.DAT
%% ==============================================================

fid = fopen(inputFile,'r');

if fid == -1
    error('無法開啟 %s。',inputFile);
end

%% ===== 讀取節點資料 =====
% F15 前兩行：
%   -1
%   15
fgetl(fid);
fgetl(fid);

nodeID = [];
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

    % 節點資料結束
    if strcmp(lineText,'-1')
        break;
    end

    data = sscanf(lineText,'%f').';

    if numel(data) >= 6
        nodeID(end+1,1) = data(1); %#ok<AGROW>
        xb(end+1,1)     = data(5); %#ok<AGROW>
        yb(end+1,1)     = data(6); %#ok<AGROW>
    end
end

nNode = numel(nodeID);

if nNode < 6
    fclose(fid);
    error('讀取到的邊界節點不足，請檢查 F15.DAT。');
end

%% ===== 找到元素資料開頭 =====
% 節點後面通常為：
%   -1
%   71
% 或
%   -1
%   -1
%   71
%
% 所以繼續讀到 71 為止。

found71 = false;

while ~feof(fid)

    lineText = strtrim(fgetl(fid));

    if strcmp(lineText,'71')
        found71 = true;
        break;
    end
end

if ~found71
    fclose(fid);
    error('找不到 F15 元素資料控制碼 71。');
end

%% ===== 讀取元素連接關係 =====
elemID = [];
node1List = [];
node2List = [];

while ~feof(fid)

    headerLine = fgetl(fid);

    if ~ischar(headerLine)
        break;
    end

    headerLine = strtrim(headerLine);

    if isempty(headerLine)
        continue;
    end

    if strcmp(headerLine,'-1')
        break;
    end

    headerData = sscanf(headerLine,'%f').';

    if isempty(headerData)
        continue;
    end

    % 下一行為該元素兩個節點編號
    if feof(fid)
        break;
    end

    connLine = strtrim(fgetl(fid));
    connData = sscanf(connLine,'%f').';

    if numel(connData) >= 2
        elemID(end+1,1)    = headerData(1); %#ok<AGROW>
        node1List(end+1,1) = connData(1);   %#ok<AGROW>
        node2List(end+1,1) = connData(2);   %#ok<AGROW>
    end
end

fclose(fid);

nElem = numel(elemID);

if nElem < 6
    error('讀取到的元素數不足，請檢查 F15.DAT。');
end

%% ==============================================================
% 自動辨識外、內兩個閉合邊界
%% ==============================================================

% 對這種 F15 產生方式：
%
% 外邊界最後一個元素：
%   nOuter -> 1
%
% 內邊界最後一個元素：
%   nTotal -> nOuter+1
%
% 因此找第一個「node2 < node1」的閉合元素，
% 即可得到外邊界節點數。

closureIndex = find(node2List < node1List);

if numel(closureIndex) < 2
    error(['無法辨識兩個閉合邊界。' ...
           '請確認 F15.DAT 為外邊界 + 內邊界的雙聯通格式。']);
end

% 第一個閉合環 = 外邊界
outerEndElem = closureIndex(1);
nOuter = node1List(outerEndElem);

% 內邊界從 nOuter+1 開始，到最後一個節點
outerIDs = (1:nOuter).';
innerIDs = (nOuter+1:nNode).';

if numel(innerIDs) < 3
    error('內邊界節點不足。');
end

xOuter = xb(outerIDs);
yOuter = yb(outerIDs);

xInner = xb(innerIDs);
yInner = yb(innerIDs);

fprintf('已讀取總節點數 = %d\n',nNode);
fprintf('已讀取總元素數 = %d\n',nElem);
fprintf('外邊界節點數   = %d\n',numel(xOuter));
fprintf('內邊界節點數   = %d\n',numel(xInner));

%% ==============================================================
% 建立規則網格
%% ==============================================================

xVector = linspace(min(xOuter),max(xOuter),nGrid);
yVector = linspace(min(yOuter),max(yOuter),nGrid);

[X,Y] = meshgrid(xVector,yVector);

%% ==============================================================
% 建立雙聯通區域遮罩
%% ==============================================================

% 外邊界
[inOuter,onOuter] = inpolygon(X,Y,xOuter,yOuter);

% 內邊界
[inInner,onInner] = inpolygon(X,Y,xInner,yInner);

% 物理區域：
%   外邊界內
%   且不在內邊界內
%   且排除兩條邊界本身
inside = inOuter & ~onOuter & ~inInner & ~onInner;

%% ===== 整理場內點 =====
fieldPoints = sortrows([X(inside),Y(inside)],[2,1]);

xin = fieldPoints(:,1);
yin = fieldPoints(:,2);

nField = size(fieldPoints,1);

fprintf('共產生 %d 個雙聯通場內點。\n',nField);

%% ==============================================================
% 輸出 1F80.DAT
%% ==============================================================

fid = fopen(outputFile,'w');

if fid == -1
    error('無法建立 %s。',outputFile);
end

for i = 1:nField

    % Fortran FORMAT(I6,2E13.5)
    fprintf(fid,'%6d%13.5E%13.5E\r\n', ...
        i,xin(i),yin(i));
end

fclose(fid);

fprintf('已輸出 %s。\n',outputFile);

%% ==============================================================
% 繪製場內點分布圖
%% ==============================================================

figure;
hold on;

% 場內點
hField = plot(xin,yin,'.','MarkerSize',7);

% 外邊界
hOuter = plot([xOuter;xOuter(1)], ...
              [yOuter;yOuter(1)], ...
              '-','LineWidth',1.5);

% 內邊界
hInner = plot([xInner;xInner(1)], ...
              [yInner;yInner(1)], ...
              '-','LineWidth',1.5);

% F15 離散邊界點
hNode = plot(xb,yb,'o', ...
    'MarkerSize',4, ...
    'LineWidth',1);

axis equal;
axis tight;

grid on;
box on;

xlabel('x');
ylabel('y');

title(sprintf( ...
    'Double-connected (N,k)-gon interior points (%d x %d grid, %d points)', ...
    nGrid,nGrid,nField));

legend([hField,hOuter,hInner,hNode], ...
    '1F80 場內點', ...
    '外邊界', ...
    '內邊界', ...
    'F15 邊界點', ...
    'Location','best');
