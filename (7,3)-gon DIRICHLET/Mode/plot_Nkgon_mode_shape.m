%% ==========================================================
% (N,k)-gon 模態繪圖－Zero-crossing 版本
% MATLAB R2024a
%
% 說明：
% 1. contourf、contour 僅用於彩色模態與一般等高線。
% 2. 黑色內域節線不使用 contourc，而是以網格零交越擷取。
% 3. 不平滑、不強制拉直，也不修改原始模態資料。
% 4. 只保護幾何中心附近的真實節線交會，邊界區仍會清除雜線。
% 5. 只有原始 Zero-crossing 確實存在垂直節線時，才延伸到邊界。
%% ==========================================================

clc;
clear;
close all;

%% ===== 可修改設定 =====
nGrid          = 500;    % 二維、三維內插網格密度
scaleXY        = 5;      % 三維 x、y 放大倍率
scaleZ         = 270;    % 三維 z 放大倍率

contourNumber  = 10;     % 二維一般等高線數量
contourWidth   = 0.4;    % 二維一般等高線線寬
nodalWidth     = 2.5;    % 內域零節線線寬
boundaryWidth  = 2.0;    % 幾何邊界線寬

boundaryCells  = 1;      % 邊界向內排除幾個網格，避免邊緣假節線
windowCells    = 11;     % 局部振幅視窗大小，建議使用奇數
rangeRatio        = 0.015; % 局部振幅低於全域振幅此比例時，不擷取節線
centerProtectRatio = 0.05;  % 只保護幾何中心附近的真實節線交會區

% 指定垂直節線延伸到幾何邊界
extendVerticalToBoundary = true; % true：將指定垂直零節線補到上下邊界
targetX                  = 0.0;  % 要延伸的垂直節線位置
targetTol                = 0.10; % 判斷節線是否位於 targetX 附近
verticalMinSpanRatio     = 0.15; % 原節線至少占幾何高度此比例才進行延伸
verticalBoundaryGap      = 0.20; % 垂直節線與上下邊界保留的實際距離

%% ===== 讀取 (N,k)-gon 邊界 =====
[xb,yb] = read_F15('F15.DAT');

%% ===== 第一模態 =====
plot_mode( ...
    'xspedU2.DAT', ...
    'xspmoU2.dat', ...
    'Mode 1', ...
    xb,yb,nGrid,scaleXY,scaleZ, ...
    contourNumber,contourWidth,nodalWidth,boundaryWidth, ...
    boundaryCells,windowCells,rangeRatio,centerProtectRatio, ...
    extendVerticalToBoundary,targetX,targetTol,verticalMinSpanRatio,verticalBoundaryGap);

%% ===== 第二模態：檔案存在且內容不全為零才繪製 =====
boundaryFile2 = 'xspedU2_2.DAT';
modeFile2     = 'xspmoU2_2.dat';

if isfile(boundaryFile2) && isfile(modeFile2)
    M2 = load(modeFile2);
    check_columns(M2,3,modeFile2);

    u2 = M2(:,3);
    u2 = u2(isfinite(u2));

    if ~isempty(u2) && any(abs(u2) > 1e-12)
        plot_mode( ...
            boundaryFile2,modeFile2,'Mode 2', ...
            xb,yb,nGrid,scaleXY,scaleZ, ...
            contourNumber,contourWidth,nodalWidth,boundaryWidth, ...
            boundaryCells,windowCells,rangeRatio,centerProtectRatio, ...
            extendVerticalToBoundary,targetX,targetTol,verticalMinSpanRatio,verticalBoundaryGap);
    else
        disp('第二模態內容全為 0，不繪製。');
    end
else
    disp('找不到完整的第二模態資料，只繪製第一模態。');
end

%% ==========================================================
% 繪製單一模態
%% ==========================================================
function plot_mode(boundaryFile,modeFile,modeName, ...
    xb,yb,nGrid,scaleXY,scaleZ, ...
    contourNumber,contourWidth,nodalWidth,boundaryWidth, ...
    boundaryCells,windowCells,rangeRatio,centerProtectRatio, ...
    extendVerticalToBoundary,targetX,targetTol,verticalMinSpanRatio,verticalBoundaryGap)

    %% 讀取資料
    if ~isfile(boundaryFile) || ~isfile(modeFile)
        error('找不到 %s 或 %s。',boundaryFile,modeFile);
    end

    B = load(boundaryFile);   % 第 1 欄：邊界編號；第 2 欄：邊界密度
    M = load(modeFile);       % 第 1、2、3 欄：x、y、u

    check_columns(B,2,boundaryFile);
    check_columns(M,3,modeFile);

    valid = all(isfinite(M(:,1:3)),2);
    M = M(valid,:);

    if isempty(M)
        error('%s 沒有可用的有限數值資料。',modeFile);
    end

    modeTag = regexprep(modeName,'\s+','');

    %% 邊界密度圖
    figure;
    plot(B(:,1),B(:,2),'-o', ...
        'LineWidth',1.5,'MarkerSize',5);

    bmax = max(abs(B(:,2)),[],'omitnan');
    if isfinite(bmax) && bmax > 0
        ylim(1.1*[-bmax bmax]);
    end

    grid on;
    box on;
    xlabel('Element index');
    ylabel('Boundary traction');
    set(gca,'FontSize',12);

    exportgraphics(gcf,[modeTag '_Boundary.tif'], ...
    'Resolution',1200);

    %% 建立內插函數
    % 彩色模態曲面加入固定邊界 u = 0，使曲面在邊界回到零。
    Fsurface = scatteredInterpolant( ...
        [M(:,1);xb], ...
        [M(:,2);yb], ...
        [M(:,3);zeros(size(xb))], ...
        'natural','none');

    % 節線只使用場內資料，避免固定邊界的零值形成整圈假節線。
    Fnodal = scatteredInterpolant( ...
        M(:,1),M(:,2),M(:,3), ...
        'natural','none');

    %% 建立 Cartesian 規則網格
    xVector = linspace(min(xb),max(xb),nGrid);
    yVector = linspace(min(yb),max(yb),nGrid);
    [xq,yq] = meshgrid(xVector,yVector);

    uq      = Fsurface(xq,yq);
    uqNodal = Fnodal(xq,yq);

    % 只保留多邊形內域。
    inside = inpolygon(xq,yq,xb,yb);
    uq(~inside)      = NaN;
    uqNodal(~inside) = NaN;

    %% 二維模態圖
    figure;
    contourf(xq,yq,uq,60,'LineStyle','none');
    hold on;

    contour(xq,yq,uq,contourNumber, ...
        'Color',[0.35 0.35 0.35], ...
        'LineWidth',contourWidth);

    segments2D = draw_nodal_line_2d_zero( ...
        xq,yq,uqNodal,xb,yb,nodalWidth, ...
        boundaryCells,windowCells,rangeRatio,centerProtectRatio);

    % 邊界直接由 F15.DAT 節點連線，不經過內插。
    plot([xb;xb(1)],[yb;yb(1)], ...
        'k-','LineWidth',boundaryWidth);

    % 只有 Zero-crossing 本身確實存在 targetX 附近的垂直節線時，
    % 才將該節線補到幾何邊界；避免其他模態被人工加入垂直線。
    hasVerticalNodal = has_vertical_nodal_segment( ...
        segments2D,targetX,targetTol,verticalMinSpanRatio,xb,yb);

    if extendVerticalToBoundary && hasVerticalNodal
        hVertical2D = draw_vertical_boundary_extension_2d( ...
            targetX,xb,yb,nodalWidth,verticalBoundaryGap);
        if isgraphics(hVertical2D)
            uistack(hVertical2D,'top');
        end
    end

    axis equal;
    axis tight;
    grid on;
    box on;
    xlabel('x');
    ylabel('y');
    colorbar;
    colormap(jet);
    set(gca,'FontSize',12);

    exportgraphics(gcf,[modeTag '_Contour.tif'], ...
    'Resolution',600);

    %% 三維模態圖
    figure;
    surf(scaleXY*xq,scaleXY*yq,scaleZ*uq, ...
        'EdgeColor','none');
    hold on;

    draw_nodal_line_3d_zero( ...
        xq,yq,uqNodal,Fsurface,xb,yb, ...
        scaleXY,scaleZ,nodalWidth, ...
        boundaryCells,windowCells,rangeRatio,centerProtectRatio);

    % 固定邊界 u = 0，直接畫在 z = 0。
    plot3(scaleXY*[xb;xb(1)], ...
          scaleXY*[yb;yb(1)], ...
          zeros(numel(xb)+1,1), ...
          'k-','LineWidth',boundaryWidth);

    % 三維圖沿用二維 Zero-crossing 的判斷結果。
    if extendVerticalToBoundary && hasVerticalNodal
        draw_vertical_boundary_extension_3d( ...
            targetX,xb,yb,scaleXY,nodalWidth,scaleZ,verticalBoundaryGap);
    end

    axis tight;
    daspect([1 1 1]);
    grid on;
    box on;
    xlabel('x');
    ylabel('y');
    zlabel('Scaled u(x,y)');
    view(45,35);
    camproj('perspective');
    colorbar;
    colormap(jet);
    lighting gouraud;
    camlight headlight;
    material dull;
    set(gca,'FontSize',12);

    disp(['請調整 ' modeName ' 的三維視角，完成後回 Command Window 按 Enter 儲存...']);
    pause;

    exportgraphics(gcf,[modeTag '_Surface.tif'], ...
    'Resolution',600);
end

%% ==========================================================
% 二維內域零節線
%% ==========================================================
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
function hasVertical = has_vertical_nodal_segment( ...
    segments,targetX,targetTol,minSpanRatio,xb,yb)

    hasVertical = false;

    if isempty(segments)
        return;
    end

    xMid = 0.5*(segments(:,1)+segments(:,3));
    dx   = abs(segments(:,3)-segments(:,1));
    dy   = abs(segments(:,4)-segments(:,2));

    % 只挑選 targetX 附近且明顯偏垂直的 Zero-crossing 小線段。
    % 使用 dy >= 2*dx，避免中心附近的斜線或水平節線被誤判。
    candidate = abs(xMid-targetX) <= targetTol & dy >= 2*dx;

    if ~any(candidate)
        return;
    end

    % 使用候選線段端點的完整 y 範圍，而不是只使用中點。
    yCandidate = [segments(candidate,2);segments(candidate,4)];
    geometryHeight = max(yb)-min(yb);
    verticalSpan = max(yCandidate)-min(yCandidate);

    % 還要求候選垂直線段至少分布在多個不同高度，避免單一短線誤判。
    yTolerance = 2*eps(max(1,max(abs(yCandidate))));
    nDistinctLevels = numel(unique(round(yCandidate/max(yTolerance,1e-12))));

    hasVertical = verticalSpan >= minSpanRatio*geometryHeight && ...
                  nDistinctLevels >= 4;
end

%% ==========================================================
% 將指定垂直節線補畫到多邊形上下邊界（二維）
%% ==========================================================
function hLine = draw_vertical_boundary_extension_2d( ...
    targetX,xb,yb,lineWidth,boundaryGap)

    hLine = gobjects(0);
    yHit = vertical_line_polygon_intersections(targetX,xb,yb);

    if numel(yHit) < 2
        warning('x = %.6g 與多邊形邊界的交點不足，無法延伸垂直節線。',targetX);
        return;
    end

    yStart = min(yHit)+boundaryGap;
    yEnd   = max(yHit)-boundaryGap;

    if yStart >= yEnd
        warning('verticalBoundaryGap 太大，垂直節線沒有可繪製的範圍。');
        return;
    end

    hLine = plot([targetX,targetX],[yStart,yEnd], ...
        'k-','LineWidth',lineWidth,'Clipping','on');
end

%% ==========================================================
% 將指定垂直節線補畫到多邊形上下邊界（三維，零節線位於 z=0）
%% ==========================================================
function draw_vertical_boundary_extension_3d( ...
    targetX,xb,yb,scaleXY,lineWidth,scaleZ,boundaryGap)

    yHit = vertical_line_polygon_intersections(targetX,xb,yb);

    if numel(yHit) < 2
        warning('x = %.6g 與多邊形邊界的交點不足，無法延伸三維垂直節線。',targetX);
        return;
    end

    yStart = min(yHit)+boundaryGap;
    yEnd   = max(yHit)-boundaryGap;

    if yStart >= yEnd
        warning('verticalBoundaryGap 太大，三維垂直節線沒有可繪製的範圍。');
        return;
    end

    zOffset = max(eps,1e-5*abs(scaleZ));
    plot3(scaleXY*[targetX,targetX], ...
          scaleXY*[yStart,yEnd], ...
          zOffset*[1,1], ...
          'k-','LineWidth',lineWidth,'Clipping','on');
end

%% ==========================================================
% 求垂直線 x = targetX 與多邊形邊界的所有交點
%% ==========================================================
function yHit = vertical_line_polygon_intersections(targetX,xb,yb)

    xb = xb(:);
    yb = yb(:);

    x1 = xb;
    y1 = yb;
    x2 = [xb(2:end);xb(1)];
    y2 = [yb(2:end);yb(1)];

    yHit = [];

    for i = 1:numel(x1)
        edgeDx = x2(i)-x1(i);

        if abs(edgeDx) <= eps
            if abs(targetX-x1(i)) <= 1e-12
                yHit = [yHit;y1(i);y2(i)]; %#ok<AGROW>
            end
            continue;
        end

        t = (targetX-x1(i))/edgeDx;

        if t >= -1e-12 && t <= 1+1e-12
            yHit(end+1,1) = y1(i)+t*(y2(i)-y1(i)); %#ok<AGROW>
        end
    end

    if isempty(yHit)
        return;
    end

    yHit = sort(yHit);
    tolerance = 1e-10*max(1,max(abs(yHit)));
    yHit = yHit([true;diff(yHit)>tolerance]);
end

%% ==========================================================
% 排除低振幅、近乎全零的平坦區域
% 只作用於節線擷取副本，不修改彩色模態或原始資料
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

%% ==========================================================
% 讀取 F15.DAT 邊界節點
%% ==========================================================
function [xb,yb] = read_F15(inputFile)

    fid = fopen(inputFile,'r');
    if fid == -1
        error('無法開啟 %s。',inputFile);
    end

    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    % 跳過前兩行控制資料。
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

        % 節點資料以 -1 結束。
        if strcmp(lineText,'-1')
            break;
        end

        data = sscanf(lineText,'%f');

        % F15 第 5、6 欄為 x、y。
        if numel(data) >= 6
            xb(end+1,1) = data(5); %#ok<AGROW>
            yb(end+1,1) = data(6); %#ok<AGROW>
        end
    end

    if numel(xb) < 3
        error('F15.DAT 的邊界點不足。');
    end

    % 若最後一點與第一點重複，只保留一份。
    if hypot(xb(end)-xb(1),yb(end)-yb(1)) < 1e-12
        xb(end) = [];
        yb(end) = [];
    end
end

%% ==========================================================
% 檢查資料欄數
%% ==========================================================
function check_columns(A,nColumn,fileName)

    if size(A,2) < nColumn
        error('%s 的資料欄數不足，至少需要 %d 欄。', ...
            fileName,nColumn);
    end
end
