%% ================================================================
% (7,3)-gon 重根模態的二維特徵空間：線性組合與繪圖
%
% 讀取 Fortran 已輸出的兩組獨立「域內」模態 xspmoU2*.dat
%% ================================================================
clear; clc; close all;

%% ===== 使用者設定 =====
anglesDeg    = [45];                   % 可自由增減角度（度）
nGrid        = 500;                    % 內插繪圖網格密度
scaleXY      = 5;                      % 三維 XY 視覺放大倍率
scaleZ       = 270;                    % 三維 Z 視覺放大倍率
saveFigures  = true;                   % 是否輸出 TIFF 圖片
outputDir    = 'Mode_Superposition';   % 新檔儲存資料夾
boundaryFile = 'F15.DAT';              % 七角星 (7,3)-gon 邊界

% ===== 沿用原始 (N,k)-gon 的零節線與圖形設定 =====
contourNumber     = 10;     % 一般等高線數量
contourWidth      = 0.4;    % 一般等高線線寬
nodalWidth        = 2.5;    % 黑色內部零節線線寬
boundaryWidth     = 2.0;    % 幾何邊界線寬
boundaryCells     = 1;      % 排除緊鄰邊界的網格，避免尖角旁的假節線
windowCells       = 11;     % 局部振幅視窗大小（奇數）
rangeRatio        = 0.010;  % 低振幅平坦區過濾比例
centerProtectRatio = 0.05;  % 保留中心附近可能存在的真正節線交會
pauseBeforeSave   = false;  % 若要手動調整三維視角再存圖，改為 true
% 疊加模態不人工延伸垂直節線，避免畫出原始資料不存在的線。

% ===== 讀取兩組內部模態 =====
mode1File = 'xspmoU2.dat';
mode2File = 'xspmoU2_2.dat';

A = load(mode1File);
B = load(mode2File);

if size(A,2)<3 || size(B,2)<3 || isempty(A) || isempty(B)
    error('兩個模態檔都必須含至少三欄：[x, y, u]。');
end
A = A(:,1:3); B = B(:,1:3);
if size(A,1) ~= size(B,1) || any(~isfinite(A(:))) || any(~isfinite(B(:)))
    error('兩個模態的場點數必須相同，且不可含 NaN / Inf。');
end
if max(abs(A(:,1:2)-B(:,1:2)),[],'all') > 1e-8
    error('兩組模態的場點座標或排列不同，不能直接逐點線性組合。');
end

xy = A(:,1:2);
u1 = A(:,3);
u2 = B(:,3);

%% ===== 檢查模態是否為非零向量 =====
normU1 = norm(u1,2);
normU2 = norm(u2,2);
if normU1 == 0 || normU2 == 0
    error('其中一組內部模態全為零，無法進行有效的模態疊加。');
end



%% ===== 模態相關性與 MAC =====
% 僅為檢驗重新排列場點，不改變原始模態及後續疊加結果。
[~, p] = sortrows(xy, [1 2]);

u1check = u1(p);
u2check = u2(p);

% 與獨立性檢查程式採用相同的正規化及內積計算。
u1n = u1check / norm(u1check);
u2n = u2check / norm(u2check);

eta = abs(u1n' * u2n);
MAC12 = eta^2;

fprintf('\n===== 模態相關性分析 =====\n');
fprintf('正規化離散內積 eta = %.10e\n', eta);
fprintf('模態保證準則 MAC12 = %.10e\n', MAC12);

% 僅檢查線性獨立性，不對模態正交化或正規化。
s = svd([u1,u2],'econ');
if s(2) <= 1e-10*s(1)
    error('原始兩組模態接近線性相關，請先檢查 Fortran 的重根輸出。');
end

if ~exist(outputDir,'dir'), mkdir(outputDir); end
fprintf('第一組模態：%s\n',mode1File);
fprintf('第二組模態：%s\n',mode2File);
fprintf('共同場點：%d 個\n',size(xy,1));
fprintf('原始模態內積：%.5e\n',u1'*u2);
fprintf('奇異值：%.9g、%.9g\n',s(1),s(2));

% 只把兩份原始模態的數值副本存入輸出資料夾，不修改任何原檔。
write_mode(fullfile(outputDir,'BEM_original_u1.dat'),xy,u1);
write_mode(fullfile(outputDir,'BEM_original_u2.dat'),xy,u2);

%% ===== 建立對應的 (7,3)-gon 繪圖邊界 =====
% 無正確七角星 F15 時仍可輸出組合數值，但不會錯畫成其他幾何。
canPlot = false;
if isfile(boundaryFile)
    [xb,yb] = read_F15(boundaryFile);
    if numel(xb) >= 3
        % 與原始場點一起檢查，避免誤用 L 型或五角星的 F15.DAT
        fractionInside = mean(inpolygon(xy(:,1),xy(:,2),xb,yb));
        if fractionInside >= 0.995
            canPlot = true;
        else
            warning('F15.DAT 與模態場點不符（內含率 %.1f%%），跳過圖片，仍輸出數值檔。',100*fractionInside);
        end
    end
else
    warning('未找到 F15.DAT：只存數值檔。要畫七角星請提供對應 F15.DAT。');
end

if canPlot
    xv = linspace(min(xb),max(xb),nGrid);
    yv = linspace(min(yb),max(yb),nGrid);
    [X,Y] = meshgrid(xv,yv);
    inside = inpolygon(X,Y,xb,yb);
end

%% ===== 直接使用原始 u1、u2 作線性疊加 =====
% u_theta = cos(theta)*u1 + sin(theta)*u2
% 不正交化、不正規化、不另外產生 u+、u-。
% 0° 與 90° 分別對應原始 u1 與 u2；其餘為混合模態。
for theta = anglesDeg
    uTheta = cosd(theta)*u1 + sind(theta)*u2;
    tag = sprintf('BEM_theta_%03ddeg',round(theta));
    write_mode(fullfile(outputDir,[tag '.dat']),xy,uTheta);
    if canPlot
        plot_mode_pair(xy,uTheta,xb,yb,X,Y,inside,scaleXY,scaleZ, ...
            sprintf('theta = %.0f deg',theta),fullfile(outputDir,tag),saveFigures, ...
            contourNumber,contourWidth,nodalWidth,boundaryWidth, ...
            boundaryCells,windowCells,rangeRatio,centerProtectRatio,pauseBeforeSave);
    end
    fprintf('theta = %6.1f°：已輸出 %s.dat\n',theta,tag);
end

fprintf('完成：所有新結果存於 %s\n',outputDir);
if ~canPlot
    fprintf('提醒：提供七角星 F15.DAT 後才會額外輸出模態圖。\n');
end

%% ================================================================
% 本檔案使用的區域函數
%% ================================================================
function write_mode(filename,xy,u)
    fid = fopen(filename,'w');
    if fid<0, error('無法輸出 %s',filename); end
    fprintf(fid,'%.10e %.10e %.10e\n',[xy,u].');
    fclose(fid);
end

function [xb,yb] = read_F15(filename)
    fid=fopen(filename,'r');
    if fid<0, error('無法讀取 %s',filename); end
    cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fgetl(fid); fgetl(fid);  % 跳過兩行控制資料
    xb=[]; yb=[];
    while ~feof(fid)
        str=fgetl(fid);
        if ~ischar(str), break; end
        str=strtrim(str);
        if isempty(str), continue; end
        if strcmp(str,'-1'), break; end
        vals=sscanf(str,'%f');
        if numel(vals)>=6
            xb(end+1,1)=vals(5); %#ok<AGROW>
            yb(end+1,1)=vals(6); %#ok<AGROW>
        end
    end
    if numel(xb)>1 && hypot(xb(end)-xb(1),yb(end)-yb(1))<1e-12
        xb(end)=[]; yb(end)=[];
    end
end

function plot_mode_pair(xy,u,xb,yb,X,Y,inside,scaleXY,scaleZ, ...
    plotName,outPrefix,saveFigures, ...
    contourNumber,contourWidth,nodalWidth,boundaryWidth, ...
    boundaryCells,windowCells,rangeRatio,centerProtectRatio,pauseBeforeSave)

    %% ===== 建立內插函數（與原始模態繪圖程式一致） =====
    % 彩色模態於固定邊界補 u=0，僅用於繪圖；原始 u 不變。
    Fsurface = scatteredInterpolant( ...
        [xy(:,1);xb],[xy(:,2);yb],[u;zeros(size(xb))], ...
        'natural','none');
    % 內部節線只依原始場內位移建立，避免固定邊界產生整圈假節線。
    Fnodal = scatteredInterpolant(xy(:,1),xy(:,2),u,'natural','none');
    uq = Fsurface(X,Y);
    uqNodal = Fnodal(X,Y);
    uq(~inside) = NaN;
    uqNodal(~inside) = NaN;

    %% ===== 二維模態：使用原始零交越篩選，去除邊界假節線 =====
    fig2=figure('Name',[plotName ' contour'],'Color','w');
    contourf(X,Y,uq,60,'LineStyle','none');
    hold on;
    contour(X,Y,uq,contourNumber, ...
        'Color',[0.35 0.35 0.35],'LineWidth',contourWidth);

    draw_nodal_line_2d_zero(X,Y,uqNodal,xb,yb,nodalWidth, ...
        boundaryCells,windowCells,rangeRatio,centerProtectRatio);

    plot([xb;xb(1)],[yb;yb(1)],'k-','LineWidth',boundaryWidth);
    axis equal; axis tight; grid on; box on;
    xlabel('x'); ylabel('y');
    colorbar; colormap(jet);
    set(gca,'FontSize',12);
    if saveFigures
        exportgraphics(fig2,[outPrefix '_Contour.tif'],'Resolution',600);
    end

    %% ===== 三維模態：恢復原本三維節線、光照、視角與比例設定 =====
    fig3=figure('Name',[plotName ' surface'],'Color','w');
    surf(scaleXY*X,scaleXY*Y,scaleZ*uq,'EdgeColor','none');
    hold on;

    draw_nodal_line_3d_zero(X,Y,uqNodal,Fsurface,xb,yb, ...
        scaleXY,scaleZ,nodalWidth, ...
        boundaryCells,windowCells,rangeRatio,centerProtectRatio);

    plot3(scaleXY*[xb;xb(1)],scaleXY*[yb;yb(1)], ...
        zeros(numel(xb)+1,1),'k-','LineWidth',boundaryWidth);
    axis tight;
    daspect([1 1 1]);
    grid on; box on;
    xlabel('x'); ylabel('y'); zlabel('Scaled u(x,y)');
    view(45,35);
    camproj('perspective');
    colorbar; colormap(jet);
    lighting gouraud;
    camlight headlight;
    material dull;
    set(gca,'FontSize',12);

    if pauseBeforeSave && saveFigures
        fprintf('請調整 %s 的三維視角，回 Command Window 按 Enter 儲存...\n',plotName);
        pause;
    end
    if saveFigures
        exportgraphics(fig3,[outPrefix '_Surface.tif'],'Resolution',600);
    end
end

%% ===== 以下零交越函數沿用原始七角星程式 =====
function segments = draw_nodal_line_2d_zero( ...
    x,y,u,xb,yb,lineWidth,boundaryCells,windowCells,rangeRatio,centerProtectRatio)

    segments = zeros(0,4);

    values = u(isfinite(u));
    if isempty(values) || min(values) > 0 || max(values) < 0
        return;
    end

    dx = mean(diff(x(1,:)));
    dy = mean(diff(y(:,1)));
    boundaryTol = boundaryCells*max(abs(dx),abs(dy));

    % 排除靠近真實幾何邊界的網格，避免固定邊界附近的數值毛邊。
    distanceGrid = point_to_polygon_distance(x,y,xb,yb);
    u(distanceGrid <= boundaryTol) = NaN;

    % 排除近乎全零的低振幅平坦區，只作用於節線擷取副本。
    u = suppress_flat_zero_regions( ...
        x,y,u,xb,yb,windowCells,rangeRatio,centerProtectRatio);

    % 中心保護完成後再次套用邊界遮罩，避免邊緣雜線被保留。
    u(distanceGrid <= boundaryTol) = NaN;

    segments = extract_zero_crossing_segments(x,y,u);
    plot_zero_segments_2d(segments,lineWidth);

end

%% ==========================================================
% 三維內域零節線
%% ==========================================================
function draw_nodal_line_3d_zero( ...
    x,y,u,F,xb,yb,scaleXY,scaleZ,lineWidth, ...
    boundaryCells,windowCells,rangeRatio,centerProtectRatio)

    values = u(isfinite(u));
    if isempty(values) || min(values) > 0 || max(values) < 0
        return;
    end

    dx = mean(diff(x(1,:)));
    dy = mean(diff(y(:,1)));
    boundaryTol = boundaryCells*max(abs(dx),abs(dy));

    distanceGrid = point_to_polygon_distance(x,y,xb,yb);
    u(distanceGrid <= boundaryTol) = NaN;
    u = suppress_flat_zero_regions( ...
        x,y,u,xb,yb,windowCells,rangeRatio,centerProtectRatio);

    % 中心保護完成後再次套用邊界遮罩，避免邊緣雜線被保留。
    u(distanceGrid <= boundaryTol) = NaN;

    segments = extract_zero_crossing_segments(x,y,u);
    plot_zero_segments_3d(segments,F,scaleXY,scaleZ,lineWidth);

end

%% ==========================================================
% 以 Marching Squares 概念擷取每個網格中的零交越線段
%
% segments 每一列：
% [x1, y1, x2, y2]
%% ==========================================================
function segments = extract_zero_crossing_segments(x,y,u)

    [nRow,nCol] = size(u);
    segments = zeros(0,4);

    % 最多每個 cell 產生兩條線段，先配置較大空間以降低重複擴充成本。
    maxSegments = 2*(nRow-1)*(nCol-1);
    buffer = nan(maxSegments,4);
    nSegment = 0;

    for i = 1:nRow-1
        for j = 1:nCol-1

            % 四個角點依逆時針排列：
            % 4 ----- 3
            % |       |
            % 1 ----- 2
            px = [x(i,j),x(i,j+1),x(i+1,j+1),x(i+1,j)];
            py = [y(i,j),y(i,j+1),y(i+1,j+1),y(i+1,j)];
            v  = [u(i,j),u(i,j+1),u(i+1,j+1),u(i+1,j)];

            if any(~isfinite(v))
                continue;
            end

            % 全正或全負時沒有零交越。
            if all(v >= 0) || all(v < 0)
                continue;
            end

            % 邊順序：下、右、上、左。
            edgePair = [1 2;2 3;3 4;4 1];
            edgePoint = nan(4,2);
            edgeHasCrossing = false(1,4);

            for e = 1:4
                a = edgePair(e,1);
                b = edgePair(e,2);

                % 將精確零值歸入非負側，可穩定處理恰好落在網格線上的節線。
                sideA = v(a) >= 0;
                sideB = v(b) >= 0;

                if sideA ~= sideB
                    edgeHasCrossing(e) = true;
                    edgePoint(e,:) = interpolate_zero_point( ...
                        px(a),py(a),v(a),px(b),py(b),v(b));
                end
            end

            crossedEdges = find(edgeHasCrossing);
            nCrossing = numel(crossedEdges);

            if nCrossing == 2
                p1 = edgePoint(crossedEdges(1),:);
                p2 = edgePoint(crossedEdges(2),:);
                [buffer,nSegment] = append_segment(buffer,nSegment,p1,p2);

            elseif nCrossing == 4
                % 對角正負分布會產生四個交點。
                % 使用 cell 中心的雙線性近似值決定連接方式。
                centerValue = mean(v);
                code = (v(1)>=0) + 2*(v(2)>=0) + ...
                       4*(v(3)>=0) + 8*(v(4)>=0);

                if code == 5
                    if centerValue >= 0
                        pairs = [1 2;3 4];
                    else
                        pairs = [4 1;2 3];
                    end
                elseif code == 10
                    if centerValue >= 0
                        pairs = [4 1;2 3];
                    else
                        pairs = [1 2;3 4];
                    end
                else
                    % 理論上四交點只會出現在 code 5 或 10。
                    pairs = [1 2;3 4];
                end

                for k = 1:2
                    p1 = edgePoint(pairs(k,1),:);
                    p2 = edgePoint(pairs(k,2),:);
                    [buffer,nSegment] = append_segment(buffer,nSegment,p1,p2);
                end
            end
        end
    end

    segments = buffer(1:nSegment,:);

    if isempty(segments)
        return;
    end

    % 移除長度接近零的退化線段。
    lengthSegment = hypot(segments(:,3)-segments(:,1), ...
                          segments(:,4)-segments(:,2));
    tolerance = 1e-12*max(1,max(abs([x(:);y(:)])));
    segments = segments(lengthSegment > tolerance,:);
end

%% ==========================================================
% 線性內插求邊上的零交越位置
%% ==========================================================
function point = interpolate_zero_point(x1,y1,u1,x2,y2,u2)

    denominator = u1-u2;

    if abs(denominator) <= eps(max(abs([u1,u2,1])))
        t = 0.5;
    else
        t = u1/denominator;
    end

    t = max(0,min(1,t));
    point = [x1+t*(x2-x1),y1+t*(y2-y1)];
end

%% ==========================================================
% 將線段加入預配置陣列
%% ==========================================================
function [buffer,nSegment] = append_segment(buffer,nSegment,p1,p2)

    if any(~isfinite([p1,p2]))
        return;
    end

    nSegment = nSegment+1;
    buffer(nSegment,:) = [p1,p2];
end

%% ==========================================================
% 一次畫出全部二維零交越線段
%% ==========================================================
function plot_zero_segments_2d(segments,lineWidth)

    if isempty(segments)
        return;
    end

    % 將每個網格內的小線段依共同端點串成完整折線。
    % 避免大量獨立短線段在螢幕反鋸齒下看起來像虛線。
    polylines = connect_zero_segments(segments);

    for k = 1:numel(polylines)
        xy = polylines{k};
        if size(xy,1) >= 2
            plot(xy(:,1),xy(:,2),'k-', ...
                'LineWidth',lineWidth);
        end
    end
end

%% ==========================================================
% 將 Marching Squares 的短線段串成連續折線
%% ==========================================================
function polylines = connect_zero_segments(segments)

    if isempty(segments)
        polylines = {};
        return;
    end

    points = [segments(:,1:2);segments(:,3:4)];
    scale = max(1,max(abs(points),[],'all'));
    tol = 1e-8*scale;

    % 將數值上幾乎相同的端點量化成同一個節點。
    key = round(points/tol);
    [~,~,nodeId] = unique(key,'rows','stable');

    nSeg = size(segments,1);
    edgeNodes = [nodeId(1:nSeg),nodeId(nSeg+1:end)];
    nNode = max(nodeId);

    nodeXY = zeros(nNode,2);
    for i = 1:nNode
        nodeXY(i,:) = mean(points(nodeId==i,:),1);
    end

    adjacency = cell(nNode,1);
    for e = 1:nSeg
        a = edgeNodes(e,1);
        b = edgeNodes(e,2);
        adjacency{a}(end+1,:) = [e,b]; %#ok<AGROW>
        adjacency{b}(end+1,:) = [e,a]; %#ok<AGROW>
    end

    used = false(nSeg,1);
    polylines = {};

    % 先從度數不等於 2 的端點開始，可正確處理開放節線。
    degree = cellfun(@(a)size(a,1),adjacency);
    startNodes = find(degree~=2 & degree>0);

    for s = startNodes(:).'
        while any(~used(adjacency{s}(:,1)))
            [lineNodes,used] = trace_polyline(s,adjacency,used);
            if numel(lineNodes)>=2
                polylines{end+1} = nodeXY(lineNodes,:); %#ok<AGROW>
            end
        end
    end

    % 剩下的是封閉迴圈。
    while any(~used)
        e = find(~used,1);
        s = edgeNodes(e,1);
        [lineNodes,used] = trace_polyline(s,adjacency,used);
        if numel(lineNodes)>=2
            polylines{end+1} = nodeXY(lineNodes,:); %#ok<AGROW>
        end
    end
end

function [lineNodes,used] = trace_polyline(startNode,adjacency,used)

    lineNodes = startNode;
    current = startNode;
    previousEdge = 0;

    while true
        choices = adjacency{current};
        if isempty(choices)
            break;
        end

        available = choices(~used(choices(:,1)),:);
        if previousEdge~=0
            available = available(available(:,1)~=previousEdge,:);
        end

        if isempty(available)
            break;
        end

        edge = available(1,1);
        nextNode = available(1,2);
        used(edge) = true;
        lineNodes(end+1) = nextNode; %#ok<AGROW>
        previousEdge = edge;
        current = nextNode;

        if current==startNode
            break;
        end
    end
end

%% ==========================================================
% 一次畫出全部三維零交越線段
%% ==========================================================
function plot_zero_segments_3d(segments,F,scaleXY,scaleZ,lineWidth)

    if isempty(segments)
        return;
    end

    x1 = segments(:,1);
    y1 = segments(:,2);
    x2 = segments(:,3);
    y2 = segments(:,4);

    z1 = scaleZ*F(x1,y1);
    z2 = scaleZ*F(x2,y2);

    keep = isfinite(z1) & isfinite(z2);
    if ~any(keep)
        return;
    end

    x1 = x1(keep);
    y1 = y1(keep);
    x2 = x2(keep);
    y2 = y2(keep);
    z1 = z1(keep);
    z2 = z2(keep);

    % 稍微抬高節線，避免與曲面完全重合造成閃爍。
    zRange = max([z1;z2])-min([z1;z2]);
    zOffset = 1e-5*max(1,zRange);

    n = numel(x1);
    xPlot = [scaleXY*x1,scaleXY*x2,nan(n,1)].';
    yPlot = [scaleXY*y1,scaleXY*y2,nan(n,1)].';
    zPlot = [z1+zOffset,z2+zOffset,nan(n,1)].';

    plot3(xPlot(:),yPlot(:),zPlot(:), ...
        'k-','LineWidth',lineWidth);
end

%% ==========================================================
% 判斷 Zero-crossing 線段中是否存在指定位置附近的垂直節線
%% ==========================================================
function uOut = suppress_flat_zero_regions( ...
    x,y,uIn,xb,yb,windowCells,rangeRatio,centerProtectRatio)

    uOut = uIn;
    finiteValues = uIn(isfinite(uIn));

    if isempty(finiteValues)
        return;
    end

    globalRange = max(finiteValues)-min(finiteValues);
    if ~isfinite(globalRange) || globalRange <= eps
        return;
    end

    % 視窗大小固定為至少 3 的奇數。
    windowCells = max(3,round(windowCells));
    if mod(windowCells,2) == 0
        windowCells = windowCells+1;
    end

    % 局部最大、最小與振幅範圍。
    localMax = movmax(uIn,windowCells,1,'omitnan');
    localMax = movmax(localMax,windowCells,2,'omitnan');

    localMin = movmin(uIn,windowCells,1,'omitnan');
    localMin = movmin(localMin,windowCells,2,'omitnan');

    localRange = localMax-localMin;

    % 原本要刪除的低振幅平坦區域。
    flatMask = localRange < rangeRatio*globalRange;

    % 局部同時存在正值與負值，表示可能有真實零交越。
    hasSignChange = localMin < 0 & localMax > 0;

    % 只保護幾何中心附近，避免把尖角和邊界內插雜訊也保留下來。
    centerX = mean(xb);
    centerY = mean(yb);

    geometryWidth  = max(xb)-min(xb);
    geometryHeight = max(yb)-min(yb);
    protectRadius  = centerProtectRatio*min(geometryWidth,geometryHeight);

    centerMask  = hypot(x-centerX,y-centerY) <= protectRadius;
    protectMask = centerMask & hasSignChange;

    % 中心交會節線保留；其他低振幅區域仍刪除。
    flatMask(protectMask) = false;
    uOut(flatMask) = NaN;
end

%% ==========================================================
% 計算各網格點到多邊形邊界的最短距離
%% ==========================================================
function dMin = point_to_polygon_distance(xp,yp,xb,yb)

    xb = xb(:);
    yb = yb(:);

    x1 = xb;
    y1 = yb;
    x2 = [xb(2:end);xb(1)];
    y2 = [yb(2:end);yb(1)];

    dMin = inf(size(xp));

    for i = 1:numel(x1)
        vx = x2(i)-x1(i);
        vy = y2(i)-y1(i);
        denominator = vx^2+vy^2;

        if denominator <= eps
            distance = hypot(xp-x1(i),yp-y1(i));
        else
            t = ((xp-x1(i))*vx+(yp-y1(i))*vy)/denominator;
            t = max(0,min(1,t));

            xProjection = x1(i)+t*vx;
            yProjection = y1(i)+t*vy;
            distance = hypot(xp-xProjection,yp-yProjection);
        end

        dMin = min(dMin,distance);
    end
end

