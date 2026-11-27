clc; clear; close all;

%% ===== 可修改設定 =====
nGrid         = 500;   % 二維內插網格密度
scaleXY       = 5;     % 三維 x、y 放大倍率
scaleZ        = 120;    % 三維 z 放大倍率
nodalWidth    = 2.5;   % 內域節線線寬
boundaryWidth = 2.0;   % 內、外邊界節線線寬
contourNumber = 10;    % 二維等高線數量
contourWidth  = 0.4;   % 二維等高線線寬

%% ===== 讀取環形邊界 =====
[xOuter,yOuter,xInner,yInner] = read_F15('F15.DAT');

%% ===== 第一模態 =====
plot_mode('xspedU2.DAT','xspmoU2.dat','Mode 1', ...
    xOuter,yOuter,xInner,yInner, ...
    nGrid,scaleXY,scaleZ,nodalWidth,boundaryWidth, ...
    contourNumber,contourWidth);

%% ===== 第二模態：檔案存在且資料不全為零才繪製 =====
boundaryFile2 = 'xspedU2_2.DAT';
modeFile2     = 'xspmoU2_2.dat';

if isfile(boundaryFile2) && isfile(modeFile2)
    M2 = load(modeFile2);
    check_columns(M2,3,modeFile2);

    u2 = M2(:,3);
    if any(abs(u2(isfinite(u2))) > 1e-12)
        plot_mode(boundaryFile2,modeFile2,'Mode 2', ...
            xOuter,yOuter,xInner,yInner, ...
            nGrid,scaleXY,scaleZ,nodalWidth,boundaryWidth, ...
    contourNumber,contourWidth);
    else
        disp('第二模態為 0，不繪製。');
    end
else
    disp('找不到完整的第二模態資料，只繪製第一模態。');
end

%% ==========================================================
% 繪製單一環形模態
%% ==========================================================
function plot_mode(boundaryFile,modeFile,modeName, ...
    xOuter,yOuter,xInner,yInner, ...
    nGrid,scaleXY,scaleZ,nodalWidth,boundaryWidth, ...
    contourNumber,contourWidth)

    %% 固定設定
    nRadius = 250;   % 三維徑向網格數
    nTheta  = 500;   % 三維圓周網格數

    %% 讀取資料
    if ~isfile(boundaryFile) || ~isfile(modeFile)
        error('找不到 %s 或 %s。',boundaryFile,modeFile);
    end

    B = load(boundaryFile);   % 編號、邊界密度
    M = load(modeFile);       % x、y、u

    check_columns(B,2,boundaryFile);
    check_columns(M,3,modeFile);

    modeTag = regexprep(modeName,'\s+','');

    %% 邊界密度
    figure;
    plot(B(:,1),B(:,2),'-o','LineWidth',1.5,'MarkerSize',5);

    bmax = max(abs(B(:,2)));
    if bmax > 0
        ylim(1.1*[-bmax bmax]);
    end

    grid on; box on;
    xlabel('Element index');
    ylabel('Boundary traction');
    set(gca,'FontSize',12);

    exportgraphics(gcf,[modeTag '_Boundary.tif'], ...
    'Resolution',1200);

    %% 建立內插函數：內、外固定邊界皆補入 u = 0
    F = scatteredInterpolant( ...
        [M(:,1); xOuter; xInner], ...
        [M(:,2); yOuter; yInner], ...
        [M(:,3); zeros(size(xOuter)); zeros(size(xInner))], ...
        'natural','none');

    %% 二維模態
    xVector = linspace(min(xOuter),max(xOuter),nGrid);
    yVector = linspace(min(yOuter),max(yOuter),nGrid);
    [xq,yq] = meshgrid(xVector,yVector);

    uq = F(xq,yq);

    inOuter = inpolygon(xq,yq,xOuter,yOuter);
    inInner = inpolygon(xq,yq,xInner,yInner);
    inside  = inOuter & ~inInner;
    uq(~inside) = NaN;

    figure;
    contourf(xq,yq,uq,60,'LineStyle','none');
    hold on;

    % 灰色等高線（顯示內域模態變化）
    contour(xq,yq,uq,contourNumber, ...
        'Color',[0.35 0.35 0.35], ...
        'LineWidth',contourWidth);

    % 黑色內域節線 u = 0
    draw_nodal_line_2d( ...
        xq,yq,uq,xOuter,yOuter,xInner,yInner,nodalWidth);

    % 直接利用 F15 邊界點畫內、外固定邊界節線
    plot([xOuter; xOuter(1)],[yOuter; yOuter(1)], ...
        'k-','LineWidth',boundaryWidth);
    plot([xInner; xInner(1)],[yInner; yInner(1)], ...
        'k-','LineWidth',boundaryWidth);

    axis equal;
    xlim([min(xOuter) max(xOuter)]);
    ylim([min(yOuter) max(yOuter)]);
    grid on; box on;
    xlabel('x'); ylabel('y');
    colorbar;
    colormap(jet);
    set(gca,'FontSize',12);

    exportgraphics(gcf,[modeTag '_Contour.tif'], ...
    'Resolution',600);

    %% 三維模態
    centerX = mean(xOuter);
    centerY = mean(yOuter);
    rOuter  = mean(hypot(xOuter-centerX,yOuter-centerY));
    rInner  = mean(hypot(xInner-centerX,yInner-centerY));

    radius = linspace(rInner,rOuter,nRadius);
    angle  = linspace(0,2*pi,nTheta);
    [A,R]  = meshgrid(angle,radius);

    xp = centerX + R.*cos(A);
    yp = centerY + R.*sin(A);
    up = F(xp,yp);

    inOuter3 = inpolygon(xp,yp,xOuter,yOuter);
    inInner3 = inpolygon(xp,yp,xInner,yInner);
    up(~(inOuter3 & ~inInner3)) = NaN;

    figure;
    surf(scaleXY*xp,scaleXY*yp,scaleZ*up,'EdgeColor','none');
    hold on;

    draw_nodal_line_3d( ...
        xq,yq,uq,F,xOuter,yOuter,xInner,yInner, ...
        scaleXY,scaleZ,nodalWidth);

    % 直接利用 F15 邊界點在 z = 0 畫內、外固定邊界節線
    plot3(scaleXY*[xOuter; xOuter(1)], ...
          scaleXY*[yOuter; yOuter(1)], ...
          zeros(numel(xOuter)+1,1), ...
          'k-','LineWidth',boundaryWidth);

    plot3(scaleXY*[xInner; xInner(1)], ...
          scaleXY*[yInner; yInner(1)], ...
          zeros(numel(xInner)+1,1), ...
          'k-','LineWidth',boundaryWidth);

    axis tight;
    daspect([1 1 1]);
    grid on; box on;
    xlabel('x'); ylabel('y'); zlabel('Scaled u(x,y)');
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
% 二維節線
%% ==========================================================
function draw_nodal_line_2d( ...
    x,y,u,xOuter,yOuter,xInner,yInner,lineWidth)

    values = u(isfinite(u));
    if isempty(values) || min(values) > 0 || max(values) < 0
        return;
    end

    u = remove_boundary_zero( ...
        x,y,u,xOuter,yOuter,xInner,yInner);

    contour(x,y,u,[0 0],'k','LineWidth',lineWidth);
end

%% ==========================================================
% 三維節線
%% ==========================================================
function draw_nodal_line_3d( ...
    x,y,u,F,xOuter,yOuter,xInner,yInner, ...
    scaleXY,scaleZ,lineWidth)

    values = u(isfinite(u));
    if isempty(values) || min(values) > 0 || max(values) < 0
        return;
    end

    u = remove_boundary_zero( ...
        x,y,u,xOuter,yOuter,xInner,yInner);

    C = contourc(x(1,:),y(:,1),u,[0 0]);

    zOffset = 1e-5*scaleZ*max(abs(values));
    if ~isfinite(zOffset) || zOffset == 0
        zOffset = eps;
    end

    k = 1;
    while k < size(C,2)
        n = C(2,k);
        xLine = C(1,k+1:k+n);
        yLine = C(2,k+1:k+n);
        zLine = scaleZ*F(xLine,yLine);

        valid = isfinite(zLine);
        if any(valid)
            plot3(scaleXY*xLine(valid), ...
                  scaleXY*yLine(valid), ...
                  zLine(valid)+zOffset, ...
                  'k-','LineWidth',lineWidth);
        end

        k = k+n+1;
    end
end

%% ==========================================================
% 排除內、外固定邊界附近的零值線，只保留內域節線
%% ==========================================================
function u = remove_boundary_zero( ...
    x,y,u,xOuter,yOuter,xInner,yInner)

    centerX = mean(xOuter);
    centerY = mean(yOuter);

    r = hypot(x-centerX,y-centerY);
    rOuter = mean(hypot(xOuter-centerX,yOuter-centerY));
    rInner = mean(hypot(xInner-centerX,yInner-centerY));

    margin = 0.02*(rOuter-rInner);

    u(r >= rOuter-margin | r <= rInner+margin) = NaN;
end


%% ==========================================================
% 讀取 F15 並分離內、外邊界
%% ==========================================================
function [xOuter,yOuter,xInner,yInner] = read_F15(inputFile)

    fid = fopen(inputFile,'r');
    if fid == -1
        error('無法開啟 %s。',inputFile);
    end

    cleaner = onCleanup(@() fclose(fid));

    fgetl(fid);
    fgetl(fid);

    xb = [];
    yb = [];

    while ~feof(fid)
        lineText = strtrim(fgetl(fid));

        if isempty(lineText)
            continue;
        end

        if strcmp(lineText,'-1')
            break;
        end

        data = sscanf(lineText,'%f').';
        if numel(data) >= 6
            xb(end+1,1) = data(5); %#ok<AGROW>
            yb(end+1,1) = data(6); %#ok<AGROW>
        end
    end

    if numel(xb) < 6
        error('F15.DAT 的邊界點不足。');
    end

    jumpDistance = hypot(diff(xb),diff(yb));
    [~,splitIndex] = max(jumpDistance);

    x1 = xb(1:splitIndex);
    y1 = yb(1:splitIndex);
    x2 = xb(splitIndex+1:end);
    y2 = yb(splitIndex+1:end);

    if numel(x1) < 3 || numel(x2) < 3
        error('無法由 F15.DAT 分離內、外邊界。');
    end

    if abs(polyarea(x1,y1)) >= abs(polyarea(x2,y2))
        xOuter = x1; yOuter = y1;
        xInner = x2; yInner = y2;
    else
        xOuter = x2; yOuter = y2;
        xInner = x1; yInner = y1;
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
