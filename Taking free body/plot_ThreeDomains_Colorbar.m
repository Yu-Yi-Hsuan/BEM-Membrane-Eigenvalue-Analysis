%% 三個區域模態比較：單根 1 張、重根 2 張；三圖皆畫等高線與零節線
% 每張圖包含原始完整域 / 內圓 / 穿孔域，不重新計算 Fortran 模態
% MATLAB R2024a
clc; clear; close all;

%% ===== 使用者設定 =====
nGrid      = 500;       % 二維內插網格：nGrid × nGrid
nLevels    = 100;       % 彩色等高線層數
boundaryLW = 2.0;       % 幾何邊界線寬
nodalLW    = 2.5;       % 三個區域的 u=0 零節線線寬
contourLW  = 0.4;       % 一般等高線線寬
contourNum = 10;        % 每張圖的一般等高線數量
saveTIFF   = true;      % 是否輸出 600 dpi TIFF

% 星形邊界零節線修正：true 啟用；false 不啟用（圓形案例使用 false）
starBoundaryFix = true;
% 以下設定僅在 starBoundaryFix = true 時使用，參考原 (N,k)-gon 繪圖程式
boundaryCells = 2;       % 排除緊貼星形外邊界的網格層數
windowCells = 11;        % 局部振幅檢查視窗（奇數）
rangeRatio = 0.010;      % 排除局部振幅變化過小的區域
centerProtectRatio = 0.05; % 保護幾何中心附近的真實零節線交會

% raw：保留原始數值，以三組資料的共同最大絕對值設定色階
% normalized：各模態分別正規化至最大絕對值 1，再共用 [-1,1]
colorMode = 'raw';

% 獨立模態可整體反號；此設定僅對齊正負方向，不調整振幅大小
phaseSign = [1, -1, 1; ...  % Mode 1：[完整域、內圓、穿孔域]
             1, -1, 1];     % Mode 2：可依模態對應方向調整

% 三組模態與其對應的原始 F15 邊界資料
% 星形情況只對第 1、3 個區域處理星形外邊界；第 2 個內圓照原本方式繪圖。
modeFiles = { ...
    'xspmoU2(Full (5,2)-gon).dat',   'xspmoU2(Inner Circle).dat',   'xspmoU2(Perforated (5,2)-gon).dat'; ...
    'xspmoU2_2(Full (5,2)-gon).dat', 'xspmoU2_2(Inner Circle).dat', 'xspmoU2_2(Perforated (5,2)-gon).dat'};
f15Files = { ...
    'F15(Full (5,2)-gon).DAT', ...
    'F15(Inner Circle).DAT', ...
    'F15(Perforated (5,2)-gon).DAT'};

folder = fileparts(mfilename('fullpath'));
if isempty(folder)
    folder = pwd;
end

%% ===== 讀取三組 F15 幾何（Mode 1、Mode 2 共用） =====
G = cell(1,3);
for j = 1:3
    f15Path = fullfile(folder,f15Files{j});
    if ~isfile(f15Path)
        error('缺少邊界檔：%s',f15Path);
    end
    [bx,by] = read_f15_points(f15Path);
    if j == 3
        [xo,yo,xi,yi] = split_annulus(bx,by);
        G{j} = struct('outerX',xo,'outerY',yo,'innerX',xi,'innerY',yi);
    else
        G{j} = struct('outerX',bx,'outerY',by,'innerX',[],'innerY',[]);
    end
end

%% ===== 判斷是否有第二模態（重根） =====
% 第一模態必須齊全；第二模態三份都存在且非零才另外繪製。
secondExists = false(1,3);
for j = 1:3
    secondExists(j) = isfile(fullfile(folder,modeFiles{2,j}));
end
if any(secondExists) && ~all(secondExists)
    error('第二模態資料不完整：必須同時有 Full、Inner、Annulus 的 xspmoU2_2 檔案。');
end
nModes = 1 + all(secondExists);

%% ===== 共用座標與網格，保留內圓的實際相對尺寸 =====
fullX = G{1}.outerX;
fullY = G{1}.outerY;
pad = 0.04*max(max(fullX)-min(fullX),max(fullY)-min(fullY));
plotXLim = [min(fullX)-pad, max(fullX)+pad];
plotYLim = [min(fullY)-pad, max(fullY)+pad];
[X,Y] = meshgrid(linspace(plotXLim(1),plotXLim(2),nGrid), ...
                 linspace(plotYLim(1),plotYLim(2),nGrid));

%% ===== 逐一繪製 Mode 1、Mode 2（各為三圖共用色階） =====
for modeIdx = 1:nModes
    D = cell(1,3);
    zeroMode = false(1,3);

    for j = 1:3
        modePath = fullfile(folder,modeFiles{modeIdx,j});
        if ~isfile(modePath)
            error('缺少模態檔：%s',modePath);
        end
        A = readmatrix(modePath);
        if size(A,2) < 3
            error('%s 至少需要三欄 x、y、u。',modeFiles{modeIdx,j});
        end
        A = A(:,1:3);
        A = A(all(isfinite(A),2),:);
        if isempty(A)
            error('%s 沒有有效模態資料。',modeFiles{modeIdx,j});
        end
        zeroMode(j) = ~any(abs(A(:,3)) > 1e-12);

        % 僅調整正負方向；重根的兩個獨立模態各有自己的 phaseSign。
        A(:,3) = phaseSign(modeIdx,j)*A(:,3);
        if strcmpi(colorMode,'normalized') && ~zeroMode(j)
            A(:,3) = A(:,3)/max(abs(A(:,3)));
        elseif ~strcmpi(colorMode,'raw') && ~strcmpi(colorMode,'normalized')
            error('colorMode 必須為 raw 或 normalized。');
        end
        D{j} = A;
    end

    % 單根時，Fortran 可能輸出全零的 _2 檔案，避免畫空白第二模態。
    if modeIdx == 2 && all(zeroMode)
        fprintf('第二模態三份資料全為零，判定為單根：不繪製 Mode 2。\n');
        break;
    end
    if any(zeroMode)
        error('Mode %d 的三份資料中有非零與全零混合，無法形成完整三域比較。',modeIdx);
    end

    % 每一個 Mode 都取三組模態的共同振幅界限，不各自縮放。
    if strcmpi(colorMode,'normalized')
        colorLimit = [-1,1];
    else
        maxAbs = max(cellfun(@(a)max(abs(a(:,3))),D));
        colorLimit = [-maxAbs,maxAbs];
    end
    levels = linspace(colorLimit(1),colorLimit(2),nLevels);

    fig = figure('Color','w','Position',[80,80,1450,520]);
    tile = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

    for j = 1:3
        ax = nexttile(tile,j);
        A = D{j};
        g = G{j};

        % 使用 F15 的真實幾何邊界（固定邊界 u=0）輔助內插。
        boundaryX = [g.outerX;g.innerX];
        boundaryY = [g.outerY;g.innerY];
        F = scatteredInterpolant( ...
            [A(:,1);boundaryX],[A(:,2);boundaryY], ...
            [A(:,3);zeros(numel(boundaryX),1)],'natural','none');
        U = F(X,Y);
        inside = inpolygon(X,Y,g.outerX,g.outerY);
        if ~isempty(g.innerX)
            inside = inside & ~inpolygon(X,Y,g.innerX,g.innerY);
        end
        U(~inside) = NaN;

        contourf(ax,X,Y,U,levels,'LineStyle','none');
        hold(ax,'on');

        % 三個區域都畫一般等高線；排除 0，避免與黑色零節線重疊。
        finiteU = U(isfinite(U));
        if ~isempty(finiteU) && max(finiteU) > min(finiteU)
            lineLevels = linspace(min(finiteU),max(finiteU),contourNum+2);
            lineLevels = lineLevels(2:end-1);
            lineLevels(abs(lineLevels) < 1e-10*max(abs(finiteU))) = [];
            if ~isempty(lineLevels)
                contour(ax,X,Y,U,lineLevels, ...
                    'LineColor',[0.35 0.35 0.35],'LineWidth',contourLW);
            end
        end

        % 三個區域皆畫 u=0 節線；星形外邊界可選擇特殊清理。
        if ~isempty(finiteU) && min(finiteU)<0 && max(finiteU)>0
            if starBoundaryFix && (j==1 || j==3)
                % 只用內域數值擷取節線，避免 F15 固定邊界 u=0 造成假零線。
                Fnodal = scatteredInterpolant( ...
                    A(:,1),A(:,2),A(:,3),'natural','none');
                nodalU = Fnodal(X,Y);
                nodalU(~inside) = NaN;

                % 按真實「星形外邊界」距離排除尖角與邊緣假節線，
                % 不使用圓半徑近似；內圓孔洞保留原有的 inside 遮罩。
                dx = mean(diff(X(1,:)));
                dy = mean(diff(Y(:,1)));
                boundaryTol = boundaryCells*max(abs(dx),abs(dy));
                distOuter = distance_to_polygon(X,Y,g.outerX,g.outerY);
                nodalU(distOuter <= boundaryTol) = NaN;

                % 抑制近乎全零的局部雜訊，但保留中心附近的零交越。
                nodalU = suppress_flat_nodal( ...
                    X,Y,nodalU,g.outerX,g.outerY, ...
                    windowCells,rangeRatio,centerProtectRatio);
                nodalU(distOuter <= boundaryTol) = NaN;
            else
                % 圓形：維持原來的處理，不進行星形邊界修正。
                cx = 0.5*(min(g.outerX)+max(g.outerX));
                cy = 0.5*(min(g.outerY)+max(g.outerY));
                radial = hypot(X-cx,Y-cy);
                outerR = max(hypot(g.outerX-cx,g.outerY-cy));
                innerR = 0;
                if ~isempty(g.innerX)
                    innerR = mean(hypot(g.innerX-cx,g.innerY-cy));
                end
                margin = 0.02*(outerR-innerR);
                nodalU = U;
                nodalU(radial >= outerR-margin) = NaN;
            end

            v = nodalU(isfinite(nodalU));
            if ~isempty(v) && min(v)<0 && max(v)>0
                contour(ax,X,Y,nodalU,[0,0],'k','LineWidth',nodalLW);
            end
        end

        % 幾何邊界沿用各域自己的 F15 節點。
        plot(ax,[g.outerX;g.outerX(1)],[g.outerY;g.outerY(1)], ...
             'k-','LineWidth',boundaryLW);
        if ~isempty(g.innerX)
            plot(ax,[g.innerX;g.innerX(1)],[g.innerY;g.innerY(1)], ...
                 'k-','LineWidth',boundaryLW);
        end
        axis(ax,'equal');
        xlim(ax,plotXLim);
        ylim(ax,plotYLim);
        clim(ax,colorLimit);
        colormap(ax,jet(256));
        xlabel(ax,'x');
        ylabel(ax,'y');
        set(ax,'FontSize',11,'Box','on','Layer','top');
        hold(ax,'off');
    end

    % 每組三圖只放一個 colorbar；不顯示額外圖標題。
    cb = colorbar(ax,'eastoutside');
    cb.Layout.Tile = 'east';
    cb.FontSize = 11;
    if strcmpi(colorMode,'normalized')
        cb.Ticks = [-1,-0.5,0,0.5,1];
    end

    if saveTIFF
        outName = sprintf('ThreeDomains_Mode%d_SharedColorbar.tif',modeIdx);
        exportgraphics(fig,fullfile(folder,outName),'Resolution',600);
    end
    fprintf('Mode %d：三域共用 %s 色階 [%.6g, %.6g]。\n', ...
        modeIdx,colorMode,colorLimit(1),colorLimit(2));
end

%% ===== F15：跳過前兩行，讀取第 5、6 欄座標，-1 表示結束 =====
function [xb,yb] = read_f15_points(filePath)
    fid = fopen(filePath,'r');
    if fid < 0
        error('無法讀取 %s',filePath);
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fgetl(fid);
    fgetl(fid);

    xb = [];
    yb = [];
    while ~feof(fid)
        lineText = fgetl(fid);
        if ~ischar(lineText), break; end
        lineText = strtrim(lineText);
        if isempty(lineText), continue; end
        if strcmp(lineText,'-1'), break; end

        data = sscanf(lineText,'%f');
        if numel(data) >= 6
            xb(end+1,1) = data(5); %#ok<AGROW>
            yb(end+1,1) = data(6); %#ok<AGROW>
        end
    end

    if numel(xb) < 3
        error('%s 邊界節點不足。',filePath);
    end
    if hypot(xb(end)-xb(1),yb(end)-yb(1)) < 1e-12
        xb(end) = [];
        yb(end) = [];
    end
end

%% ===== 分離環形域 F15 的外圈與內圈 =====
function [xo,yo,xi,yi] = split_annulus(xb,yb)
    % F15 通常先列完一圈再接另一圈；跨圈連接距離最大
    jumps = hypot(diff(xb),diff(yb));
    [~,idx] = max(jumps);
    x1 = xb(1:idx);       y1 = yb(1:idx);
    x2 = xb(idx+1:end);   y2 = yb(idx+1:end);

    if numel(x1) < 3 || numel(x2) < 3
        error('環形域 F15 無法辨識內外邊界，請確認原始格式。');
    end
    % 以多邊形面積判別內外圈
    if polyarea(x1,y1) >= polyarea(x2,y2)
        xo = x1; yo = y1; xi = x2; yi = y2;
    else
        xo = x2; yo = y2; xi = x1; yi = y1;
    end
end

%% ===== 星形外邊界至網格點的最短距離（依真實 F15 多邊形） =====
function dMin = distance_to_polygon(X,Y,xb,yb)
    xb = xb(:); yb = yb(:);
    dMin = inf(size(X));
    for i = 1:numel(xb)
        k = mod(i,numel(xb))+1;
        vx = xb(k)-xb(i); vy = yb(k)-yb(i);
        den = vx^2+vy^2;
        if den <= eps
            d = hypot(X-xb(i),Y-yb(i));
        else
            t = max(0,min(1,((X-xb(i))*vx+(Y-yb(i))*vy)/den));
            d = hypot(X-(xb(i)+t*vx),Y-(yb(i)+t*vy));
        end
        dMin = min(dMin,d);
    end
end

%% ===== 星形零節線雜訊抑制（只更動擷取節線的副本） =====
function U = suppress_flat_nodal(X,Y,U,xb,yb,windowCells,rangeRatio,centerProtectRatio)
    vals = U(isfinite(U));
    if isempty(vals), return; end
    ampRange = max(vals)-min(vals);
    if ampRange <= eps(max(abs(vals))), return; end
    windowCells = max(3,round(windowCells));
    if mod(windowCells,2)==0, windowCells=windowCells+1; end

    localMax = movmax(movmax(U,windowCells,1,'omitnan'),windowCells,2,'omitnan');
    localMin = movmin(movmin(U,windowCells,1,'omitnan'),windowCells,2,'omitnan');
    flat = (localMax-localMin) < rangeRatio*ampRange;

    centerX = 0.5*(min(xb)+max(xb));
    centerY = 0.5*(min(yb)+max(yb));
    radius = centerProtectRatio*min(max(xb)-min(xb),max(yb)-min(yb));
    protect = hypot(X-centerX,Y-centerY)<=radius & localMin<0 & localMax>0;
    flat(protect) = false;
    U(flat) = NaN;
end
