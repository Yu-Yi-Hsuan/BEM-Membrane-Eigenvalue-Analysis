%% =========================================================
% 星形域薄膜模態動畫（清晰版）
%
% 改進：
% 1. 正確讀取 F15.DAT
% 2. 自動判斷第二模態，第二模態全零時略過
% 3. 在固定邊界加入 u = 0，避免三維邊界出現尖刺
% 4. 二維動畫使用對比增強，讓低振幅區域更容易辨識
% 5. 不繪製節線
% 6. 僅輸出 GIF，不輸出 MP4
%
% 模態振動：
% u(x,y,t) = phi(x,y) cos(omega*t),  omega = c*k
%
% MATLAB R2024a
%% =========================================================

clear;
clc;
close all;

%% ================= 使用者設定 =================
boundaryFile = 'F15.DAT';

modeFiles = {
    'xspmoU2.dat'
    'xspmoU2_2.dat'
};

modeNames = {
    '第一模態'
    '第二模態'
};

nGrid  = 250;        % 插值網格密度

% ===== 真實/相對振動速度設定 =====
kValue        = 4.8790;  % 此次要製作動畫的特徵波數 k
cWave         = 1.0;     % 正規化波速；若有實際材料參數，可改成 sqrt(T/rho_s)
fps           = 30;      % 所有模態固定相同播放影格率
videoDuration = 5.0;     % 每支動畫固定相同影片時間（秒）
nFrame        = round(fps*videoDuration);  % 總影格數
omega         = cWave*kValue;              % 自然圓頻率 omega = c*k

% 若只比較不同 k 的相對快慢，保持 cWave = 1 即可。
% 此時 k 越大 -> omega 越大 -> 相同影片時間內振動次數越多。

amplitudeScale = 1.7;
viewAngle      = [45,30];

showFigure = false;

nContour = 80;       % 二維色階層數
contrastGamma = 0.50;
% contrastGamma < 1：加強低振幅區域顏色
% contrastGamma = 1：完全依原始振幅顯示

enforceZeroBoundary = true;
% 固定邊界問題建議設為 true，使邊界位移為 0

%% ================= 讀取 F15 邊界 =================
[xb,yb] = readF15Boundary(boundaryFile);

fprintf('已讀取邊界：%s\n',boundaryFile);
fprintf('邊界點數：%d\n',numel(xb)-1);
fprintf('x範圍：%.6g ～ %.6g\n',min(xb),max(xb));
fprintf('y範圍：%.6g ～ %.6g\n',min(yb),max(yb));

if ~(isfinite(kValue) && kValue > 0)
    error('kValue 必須為正的有限數值。');
end
if ~(isfinite(cWave) && cWave > 0)
    error('cWave 必須為正的有限數值。');
end
if ~(isfinite(videoDuration) && videoDuration > 0)
    error('videoDuration 必須大於 0。');
end

fprintf('動畫設定：k = %.6f, c = %.6f, omega = %.6f rad/s\n', ...
    kValue,cWave,omega);
fprintf('影片時間：%.2f s，fps = %d，總影格數 = %d\n', ...
    videoDuration,fps,nFrame);

%% ================= 處理模態 =================
foundMode = false;
validMode = false;

for iMode = 1:numel(modeFiles)

    modeFile = modeFiles{iMode};

    if ~isfile(modeFile)
        fprintf('\n未找到 %s，略過%s。\n',modeFile,modeNames{iMode});
        continue;
    end

    foundMode = true;

    fprintf('\n========================================\n');
    fprintf('正在處理%s：%s\n',modeNames{iMode},modeFile);
    fprintf('========================================\n');

    modeData = readmatrix(modeFile);
    modeData = modeData(:,all(~isnan(modeData),1));
    modeData = modeData(all(isfinite(modeData),2),:);

    if isempty(modeData)
        fprintf('%s 沒有有效資料，略過。\n',modeFile);
        continue;
    end

    if size(modeData,2) >= 4
        x   = modeData(:,2);
        y   = modeData(:,3);
        phi = modeData(:,4);
        fprintf('資料格式：四欄 [編號, x, y, u]\n');
    elseif size(modeData,2) == 3
        x   = modeData(:,1);
        y   = modeData(:,2);
        phi = modeData(:,3);
        fprintf('資料格式：三欄 [x, y, u]\n');
    else
        fprintf('%s 欄數不足，略過。\n',modeFile);
        continue;
    end

    fprintf('模態點數：%d\n',numel(phi));
    fprintf('x範圍：%.6g ～ %.6g\n',min(x),max(x));
    fprintf('y範圍：%.6g ～ %.6g\n',min(y),max(y));
    fprintf('u範圍：%.6g ～ %.6g\n',min(phi),max(phi));

    phiMaxOriginal = max(abs(phi));

    if phiMaxOriginal <= 1.0e-14
        fprintf('%s 的模態值全部為零，略過此模態。\n',modeFile);
        continue;
    end

    validMode = true;
    phi = phi ./ phiMaxOriginal;

    %% 固定邊界加入 u = 0
    if enforceZeroBoundary
        xbOpen = xb(1:end-1);
        ybOpen = yb(1:end-1);

        xInterp   = [x; xbOpen];
        yInterp   = [y; ybOpen];
        phiInterp = [phi; zeros(size(xbOpen))];
    else
        xInterp   = x;
        yInterp   = y;
        phiInterp = phi;
    end

    %% 建立規則網格
    xq = linspace(min(xb),max(xb),nGrid);
    yq = linspace(min(yb),max(yb),nGrid);
    [Xq,Yq] = meshgrid(xq,yq);

    inside = inpolygon(Xq,Yq,xb,yb);

    %% 插值
    F = scatteredInterpolant( ...
        xInterp,yInterp,phiInterp, ...
        'natural','none');

    PhiGrid = F(Xq,Yq);
    PhiGrid(~inside) = NaN;

    % 只補非常少量的域內插值空白
    missingInside = inside & isnan(PhiGrid);

    if any(missingInside,'all')
        Fnearest = scatteredInterpolant( ...
            xInterp,yInterp,phiInterp, ...
            'nearest','nearest');

        PhiNearest = Fnearest(Xq,Yq);
        PhiGrid(missingInside) = PhiNearest(missingInside);
    end

    validInside = inside & isfinite(PhiGrid);

    if ~any(validInside,'all')
        error('%s 插值後沒有有效域內資料。',modeFile);
    end

    maxPhi = max(abs(PhiGrid(validInside)));

    colorLimit = amplitudeScale*[-maxPhi,maxPhi];
    zLimit     = 1.10*colorLimit;

    %% 輸出名稱
    modeTag = sprintf('Mode%d',iMode);

    gif2D = [modeTag,'_2D_Clear.gif'];
    gif3D = [modeTag,'_3D_Clear.gif'];

    if isfile(gif2D), delete(gif2D); end
    if isfile(gif3D), delete(gif3D); end

    %% 二維動畫
    create2DAnimation( ...
        Xq,Yq,PhiGrid,xb,yb,modeNames{iMode}, ...
        gif2D,nFrame,fps,omega, ...
        amplitudeScale,colorLimit,nContour, ...
        contrastGamma,showFigure);

    %% 三維動畫
    create3DAnimation( ...
        Xq,Yq,PhiGrid,xb,yb,modeNames{iMode}, ...
        gif3D,nFrame,fps,omega, ...
        amplitudeScale,colorLimit,zLimit, ...
        viewAngle,showFigure);

    fprintf('%s輸出完成（k=%.6f, omega=%.6f rad/s）：\n', ...
        modeNames{iMode},kValue,omega);
    fprintf('  %s\n',gif2D);
    fprintf('  %s\n',gif3D);
end

if ~foundMode
    error('找不到任何模態檔案。');
end

if validMode
    fprintf('\n全部有效模態動畫製作完成。\n');
else
    fprintf('\n所有模態皆為零或格式無效。\n');
end

%% =========================================================
% 讀取 F15 邊界
%% =========================================================
function [xb,yb] = readF15Boundary(filename)

    if ~isfile(filename)
        error('找不到邊界檔案：%s',filename);
    end

    fid = fopen(filename,'r');

    if fid == -1
        error('無法開啟：%s',filename);
    end

    cleaner = onCleanup(@() fclose(fid));

    fgetl(fid); % -1
    fgetl(fid); % 15

    xb = [];
    yb = [];

    while ~feof(fid)

        line = strtrim(fgetl(fid));

        if isempty(line)
            continue;
        end

        values = sscanf(line,'%f').';

        if numel(values) == 1 && values(1) == -1
            break;
        end

        if numel(values) >= 7
            xb(end+1,1) = values(5); %#ok<AGROW>
            yb(end+1,1) = values(6); %#ok<AGROW>
        end
    end

    if numel(xb) < 3
        error('從 %s 讀到的邊界點不足。',filename);
    end

    keep = [true; hypot(diff(xb),diff(yb)) > 1.0e-14];
    xb = xb(keep);
    yb = yb(keep);

    if hypot(xb(1)-xb(end),yb(1)-yb(end)) > 1.0e-14
        xb(end+1) = xb(1);
        yb(end+1) = yb(1);
    end
end

%% =========================================================
% 二維動畫
%% =========================================================
function create2DAnimation( ...
    X,Y,Phi,xb,yb,modeName, ...
    gifFile,nFrame,fps,omega, ...
    amplitudeScale,colorLimit,nContour, ...
    contrastGamma,showFigure)

    visibility = getVisibility(showFigure);

    fig = figure( ...
        'Color','w', ...
        'Position',[100,100,900,760], ...
        'Visible',visibility);

    ax = axes(fig);

    gifIndex = 0;
    maxAbsPhi = max(abs(Phi(:)),[],'omitnan');

    for iFrame = 1:nFrame

        t = (iFrame-1)/fps;
        phase = omega*t;
        U = amplitudeScale*Phi*cos(phase);

            % 僅用於顏色顯示的對比增強
            Udisplay = sign(U).* ...
                (abs(U)./(amplitudeScale*maxAbsPhi)).^contrastGamma;

            cla(ax);

            contourf(ax,X,Y,Udisplay,nContour,'LineColor','none');
            hold(ax,'on');

            % 星形邊界
            plot(ax,xb,yb,'k-','LineWidth',2.0);

            axis(ax,'equal');
            axis(ax,'tight');
            box(ax,'on');

            xlabel(ax,'x');
            ylabel(ax,'y');

            title(ax,{ ...
                [modeName,'：二維視角動畫'], ...
                sprintf('t = %.2f s，phase = %.0f°',t,mod(rad2deg(phase),360))});

            colormap(ax,turbo);
            colorbar(ax);
            clim(ax,[-1,1]);

            drawnow;

            frame = getframe(fig);

        gifIndex = gifIndex+1;
        writeGifFrame(frame,gifFile,gifIndex,fps);
    end

    close(fig);
end

%% =========================================================
% 三維動畫
%% =========================================================
function create3DAnimation( ...
    X,Y,Phi,xb,yb,modeName, ...
    gifFile,nFrame,fps,omega, ...
    amplitudeScale,colorLimit,zLimit, ...
    viewAngle,showFigure)

    visibility = getVisibility(showFigure);

    fig = figure( ...
        'Color','w', ...
        'Position',[100,100,920,760], ...
        'Visible',visibility);

    ax = axes(fig);

    gifIndex = 0;

    for iFrame = 1:nFrame

        t = (iFrame-1)/fps;
        phase = omega*t;
        U = amplitudeScale*Phi*cos(phase);

            cla(ax);

            surf(ax,X,Y,U,U, ...
                'EdgeColor','none', ...
                'FaceColor','interp');

            hold(ax,'on');

            plot3(ax,xb,yb,zeros(size(xb)), ...
                'k-', ...
                'LineWidth',2.0);

            view(ax,viewAngle);

            axis(ax,'equal');
            box(ax,'on');

            xlim(ax,[min(xb),max(xb)]);
            ylim(ax,[min(yb),max(yb)]);
            zlim(ax,zLimit);
            clim(ax,colorLimit);

            xlabel(ax,'x');
            ylabel(ax,'y');
            zlabel(ax,'u(x,y,t)');

            title(ax,{ ...
                [modeName,'：三維視角動畫'], ...
                sprintf('t = %.2f s，phase = %.0f°',t,mod(rad2deg(phase),360))});

            colormap(ax,turbo);
            colorbar(ax);

            camlight(ax,'headlight');
            lighting(ax,'gouraud');

            drawnow;

            frame = getframe(fig);

        gifIndex = gifIndex+1;
        writeGifFrame(frame,gifFile,gifIndex,fps);
    end

    close(fig);
end

%% =========================================================
% GIF 寫入
%% =========================================================
function writeGifFrame(frame,gifFile,frameIndex,fps)

    [indexedImage,colorMap] = rgb2ind(frame.cdata,256);

    if frameIndex == 1
        imwrite(indexedImage,colorMap,gifFile, ...
            'gif','LoopCount',inf,'DelayTime',1/fps);
    else
        imwrite(indexedImage,colorMap,gifFile, ...
            'gif','WriteMode','append','DelayTime',1/fps);
    end
end

%% =========================================================
% 顯示控制
%% =========================================================
function visibility = getVisibility(showFigure)

    if showFigure
        visibility = 'on';
    else
        visibility = 'off';
    end
end
