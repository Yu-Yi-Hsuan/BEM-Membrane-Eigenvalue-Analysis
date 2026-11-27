%% 正三角形 Dual BEM：特徵波數收斂分析（整理版）
% 用途：
% 1. stage = 'prepare'：建立不同邊界元素數的 F15.DAT
% 2. 在各資料夾執行 Fortran，產生 Complex_U.DAT
% 3. stage = 'analyze'：讀取掃根結果，繪製收斂圖與相對誤差圖（TIF）

clear; clc; close all;

%% ===== 0. 使用者設定 =====
stage = 'analyze';                    % 'prepare' 或 'analyze'
R = 1.0;                              % 正三角形外接圓半徑
nNodeList = [15 30 45 60 75 90];      % 總邊界元素數，必須是 3 的倍數
rootFolder = fullfile(pwd,'Triangle_Convergence');
sigmaFile = 'Complex_U.DAT';          % Fortran 掃描結果：[k, sigma_min]

% 根的辨識
manualExactK = 6.3985;                    % []：使用解析公式；或自行指定，例如 4.1888
selectedRoot = [];                    % []：自動辨識；或指定 1、2、3...
maxRootIndex = 30;                    % 最多列出前幾組不同的解析根（不計重複）

% 谷值搜尋控制
searchHalfWidth = 0.20;               % 只在目標解析根 ± 此範圍內找谷值
maxSearchError = 0.10;                % 找到的根若偏離解析值太多，標示可疑
useQuadraticFit = false;              % false：直接取掃描谷值；true：三點二次擬合

% 圖形設定
numXTicks = 6;                        % X 軸主要刻度數量
numYTicks = 6;                        % Y 軸主要刻度數量
extraYMargin = 0.12;                  % Y 軸上下留白比例
plotError = true;                     % 是否另外畫相對誤差圖
exportTif = true;                     % 是否輸出 tif
exportResolution = 2400;              % tif 解析度（參考你原本的作圖程式）

%% ===== 1. 基本檢查 =====
assert(R>0 && isfinite(R),'R 必須大於 0。');
assert(all(nNodeList>0 & mod(nNodeList,3)==0), ...
    'nNodeList 每一項都必須是 3 的正整數倍。');
assert(numel(unique(nNodeList))==numel(nNodeList), ...
    'nNodeList 不能有重複值。');
assert(isempty(selectedRoot) || (isscalar(selectedRoot) && selectedRoot>=1 && selectedRoot==floor(selectedRoot)), ...
    'selectedRoot 必須是正整數或 []。');
assert(isempty(manualExactK) || (isscalar(manualExactK) && isfinite(manualExactK) && manualExactK>0), ...
    'manualExactK 必須是正數或 []。');

switch lower(stage)
    case 'prepare'
        runPrepare(rootFolder,nNodeList,R);

    case 'analyze'
        runAnalyze(rootFolder,nNodeList,R,sigmaFile,manualExactK,selectedRoot,maxRootIndex, ...
            searchHalfWidth,maxSearchError,useQuadraticFit, ...
            numXTicks,numYTicks,extraYMargin,plotError,exportTif,exportResolution);

    otherwise
        error('stage 只能填 prepare 或 analyze。');
end

%% ========================================================================
function runPrepare(rootFolder,nNodeList,R)
    if ~exist(rootFolder,'dir')
        mkdir(rootFolder);
    end

    for t = 1:numel(nNodeList)
        nb = nNodeList(t);
        folder = fullfile(rootFolder,sprintf('N%04d',nb));
        if ~exist(folder,'dir')
            mkdir(folder);
        end
        filePath = fullfile(folder,'F15.DAT');
        writeTriangleF15(filePath,nb,R);
        fprintf('已建立 %s（邊界元素 %d）\n',filePath,nb);
    end

    fprintf('\n下一步：在每個 N 資料夾中執行 Fortran，產生 %s。\n', 'Complex_U.DAT');
    fprintf('之後把 stage 改成 analyze，再重新執行此程式。\n');
end

%% ========================================================================
function runAnalyze(rootFolder,nNodeList,R,sigmaFile,manualExactK,selectedRoot,maxRootIndex, ...
    searchHalfWidth,maxSearchError,useQuadraticFit, ...
    numXTicks,numYTicks,extraYMargin,plotError,exportTif,exportResolution)

    [kAnalytical,multiplicity] = triangleRoots(R,maxRootIndex);

    nCases = numel(nNodeList);
    scans = cell(nCases,1);
    lowerBound = -inf;
    upperBound = inf;

    for t = 1:nCases
        file = fullfile(rootFolder,sprintf('N%04d',nNodeList(t)),sigmaFile);
        if ~isfile(file)
            warning('找不到 %s，這組資料將略過。',file);
            continue;
        end

        A = readmatrix(file);
        if size(A,2) < 2
            warning('資料欄不足：%s',file);
            continue;
        end

        valid = isfinite(A(:,1)) & isfinite(A(:,2)) & A(:,1)>0 & A(:,2)>=0;
        A = A(valid,1:2);
        A = sortrows(A,1);
        [~,ia] = unique(A(:,1),'stable');
        A = A(ia,:);

        if size(A,1) < 3
            warning('有效資料少於 3 列：%s',file);
            continue;
        end

        scans{t} = A;
        lowerBound = max(lowerBound,min(A(:,1)));
        upperBound = min(upperBound,max(A(:,1)));
    end

    validCases = find(~cellfun(@isempty,scans));
    assert(~isempty(validCases),'沒有可用的掃描資料，請先確認各資料夾中的 Complex_U.DAT。');
    assert(lowerBound < upperBound,'各組掃描範圍沒有交集，請檢查 EIGMIN / EIGMAX。');

    scanStep = median(diff(scans{validCases(end)}(:,1)));
    guard = max(2*scanStep,1e-10);
    candidates = find(kAnalytical > lowerBound + guard & kAnalytical < upperBound - guard);

    if ~isempty(selectedRoot)
        assert(selectedRoot <= numel(kAnalytical),'selectedRoot 超過解析根列表長度。');
        rootIndex = selectedRoot;
        assert(kAnalytical(rootIndex) > lowerBound && kAnalytical(rootIndex) < upperBound, ...
            '指定的解析根 k_%d 不在共同掃描範圍 [%.6f, %.6f]。',rootIndex,lowerBound,upperBound);
    else
        if numel(candidates) == 1
            rootIndex = candidates;
        elseif isempty(candidates)
            error('共同掃描範圍 [%.6f, %.6f] 內沒有完整涵蓋任何解析根。',lowerBound,upperBound);
        else
            fprintf('目前共同掃描範圍包含多個解析根：\n');
            disp(table(candidates,kAnalytical(candidates),multiplicity(candidates), ...
                'VariableNames',{'RootIndex','kExact','Multiplicity'}));
            error('請指定 selectedRoot，例如 1、2、3。');
        end
    end

    formulaExact = kAnalytical(rootIndex);
    if isempty(manualExactK)
        referenceK = formulaExact;
    else
        referenceK = manualExactK;
        assert(abs(referenceK-formulaExact) < max(0.02,0.005*formulaExact), ...
            '手動解析值與公式解析值差異過大，請確認輸入。');
    end

    fprintf('自動/指定辨識：第 %d 組不同的特徵波數 k_%d，重數 %d\n',rootIndex,rootIndex,multiplicity(rootIndex));
    fprintf('公式解析值：%.12f；作圖參考值：%.12f\n',formulaExact,referenceK);
    fprintf('共同掃描範圍：[%.6f, %.6f]\n',lowerBound,upperBound);

    kNumerical = nan(nCases,1);
    minSigma = nan(nCases,1);
    status = repmat({''},nCases,1);

    for t = validCases(:).'
        A = scans{t};
        [kNumerical(t),minSigma(t),status{t}] = trackMinimum(A(:,1),A(:,2), ...
            formulaExact,searchHalfWidth,maxSearchError,useQuadraticFit);
        if ~strcmp(status{t},'OK')
            warning('N = %d：%s',nNodeList(t),status{t});
        end
    end

    relativeError = abs(kNumerical-referenceK)/referenceK*100;
    Result = table(nNodeList(:),kNumerical,repmat(referenceK,nCases,1),relativeError,minSigma,status, ...
        'VariableNames',{'BoundaryElements','k_BEM','k_Analytical','RelativeError_percent','MinimumSigma','Status'});
    disp(Result);

    % ===== 圖 1：特徵波數收斂圖 =====
    fig1 = figure('Color','w','Name',sprintf('k%d convergence',rootIndex));
    hold on;
    yline(referenceK,'-','Color',[0.18 0.42 0.77],'LineWidth',1.3, ...
        'DisplayName','Analytical solution');
    plot(nNodeList,kNumerical,'--o', ...
        'Color',[0.5 0.5 0.5], ...
        'MarkerFaceColor',[0.90 0.10 0.10], ...
        'MarkerEdgeColor',[0.90 0.10 0.10], ...
        'MarkerSize',6, ...
        'LineWidth',1.0, ...
        'DisplayName','Dual BEM');
    xlabel('Number of boundary elements');
    ylabel(sprintf('k_%d',rootIndex));
    formatAxis(nNodeList,numXTicks,numYTicks,extraYMargin,[referenceK; kNumerical(isfinite(kNumerical))]);
    legend('Location','best');
    box on;
    grid off;

    if exportTif
        figure(fig1);
        drawnow;
        disp('請先調整第一張收斂圖的位置、視窗大小與圖例，完成後回 Command Window 按 Enter 儲存...');
        pause;
        exportgraphics(fig1,fullfile(rootFolder,sprintf('Triangle_k%d_Convergence.tif',rootIndex)), ...
            'Resolution',exportResolution);
        fprintf('第一張收斂圖已儲存，接著產生相對誤差圖。\n');
    end

    % ===== 圖 2：相對誤差圖 =====
    if plotError
        fig2 = figure('Color','w','Name',sprintf('k%d relative error',rootIndex));
        plot(nNodeList,relativeError,'-o', ...
            'Color',[0.18 0.42 0.77], ...
            'MarkerFaceColor',[0.18 0.42 0.77], ...
            'LineWidth',1.4);
        xlabel('Number of boundary elements');
        ylabel('Relative error (%)');
        xticks(nNodeList);
        xlim([min(nNodeList)-3, max(nNodeList)+3]);
        yData = relativeError(isfinite(relativeError));
        if isempty(yData)
            yData = [0 1];
        end
        yMin = min(yData);
        yMax = max(yData);
        yPad = max((yMax-yMin)*extraYMargin,1e-6);
        if yMin == yMax
            yMin = yMin - yPad;
            yMax = yMax + yPad;
        else
            yMin = max(0,yMin-yPad);
            yMax = yMax + yPad;
        end
        ylim([yMin yMax]);
        yticks(linspace(yMin,yMax,max(2,numYTicks)));
        grid on;
        box on;

        if exportTif
            figure(fig2);
            drawnow;
            disp('請先調整第二張相對誤差圖的位置、視窗大小與圖例，完成後回 Command Window 按 Enter 儲存...');
            pause;
            exportgraphics(fig2,fullfile(rootFolder,sprintf('Triangle_k%d_RelativeError.tif',rootIndex)), ...
                'Resolution',exportResolution);
            fprintf('第二張相對誤差圖已儲存。\n');
        end
    end
end

%% ========================================================================
function formatAxis(nNodeList,numXTicks,numYTicks,extraYMargin,yValues)
    xticks(nNodeList);
    if numel(nNodeList) > numXTicks && numXTicks >= 2
        xTickValues = round(linspace(min(nNodeList),max(nNodeList),numXTicks));
        xticks(unique(xTickValues));
    end
    xlim([min(nNodeList)-3, max(nNodeList)+3]);

    yValues = yValues(isfinite(yValues));
    yMin = min(yValues);
    yMax = max(yValues);
    yPad = max((yMax-yMin)*extraYMargin,1e-6);
    if yMin == yMax
        yMin = yMin - yPad;
        yMax = yMax + yPad;
    else
        yMin = yMin - yPad;
        yMax = yMax + yPad;
    end
    ylim([yMin yMax]);
    yticks(linspace(yMin,yMax,max(2,numYTicks)));
end

%% ========================================================================
function writeTriangleF15(filename,nNode,R)
    theta = pi/2 + (0:2)'*2*pi/3;
    V = R*[cos(theta), sin(theta)];
    perEdge = nNode/3;
    points = zeros(nNode,2);
    count = 0;

    for s = 1:3
        a = V(s,:);
        b = V(mod(s,3)+1,:);
        for j = 0:perEdge-1
            count = count + 1;
            t = j/perEdge;
            points(count,:) = (1-t)*a + t*b;
        end
    end

    fid = fopen(filename,'wt');
    assert(fid>0,'無法寫入 %s',filename);
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid,'%6d\n%6d\n',-1,15);
    for i = 1:nNode
        fprintf(fid,'%10d%10d%10d%10d%13.5E%13.5E%13.5E\n', ...
            i,0,0,11,points(i,1),points(i,2),0);
    end
    fprintf(fid,'%6d\n%6d\n%6d\n',-1,-1,71);
    for i = 1:nNode
        fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d\n',i,1,21,1,1,7,2);
        fprintf(fid,'%10d%10d\n',i,mod(i,nNode)+1);
    end
    fprintf(fid,'%6d\n',-1);
end

%% ========================================================================
function [rootK,rootMultiplicity] = triangleRoots(R,numRoots)
    side = sqrt(3)*R;
    lim = max(20,ceil(sqrt(numRoots)*5));
    vals = [];
    mult = [];

    for m = 1:lim
        for n = m:lim
            kv = 4*pi/(3*side)*sqrt(m*m + m*n + n*n);
            vals(end+1,1) = kv; %#ok<AGROW>
            mult(end+1,1) = 1 + (m<n); %#ok<AGROW>
        end
    end

    [vals,order] = sort(vals);
    mult = mult(order);

    rootK = [];
    rootMultiplicity = [];
    for j = 1:numel(vals)
        if isempty(rootK) || abs(vals(j)-rootK(end)) > 1e-9*max(1,vals(j))
            rootK(end+1,1) = vals(j); %#ok<AGROW>
            rootMultiplicity(end+1,1) = mult(j); %#ok<AGROW>
        else
            rootMultiplicity(end) = rootMultiplicity(end) + mult(j);
        end

        if numel(rootK) >= numRoots && j < numel(vals) && ...
                vals(j+1) > rootK(end) + 1e-9*max(1,rootK(end))
            break;
        end
    end

    rootK = rootK(1:numRoots);
    rootMultiplicity = rootMultiplicity(1:numRoots);
end

%% ========================================================================
function [kBest,sBest,status] = trackMinimum(k,s,kRef,halfWidth,maxError,fitEnabled)
    kBest = nan;
    sBest = nan;
    status = '搜尋範圍沒有足夠資料';

    idx = find(abs(k-kRef) <= halfWidth);
    if numel(idx) < 3
        return;
    end

    ii = idx(2:end-1);
    ii = ii(s(ii) <= s(ii-1) & s(ii) <= s(ii+1));
    if isempty(ii)
        status = '沒有內部局部谷值（可能根不在掃描範圍）';
        return;
    end

    [~,best] = min(abs(k(ii)-kRef));
    i = ii(best);
    kBest = k(i);
    sBest = s(i);

    if fitEnabled && all(s(i-1:i+1) > 0)
        x = k(i-1:i+1) - k(i);
        p = polyfit(x,log(max(s(i-1:i+1),realmin)),2);
        if p(1) > 0
            dx = -p(2)/(2*p(1));
            if dx > x(1) && dx < x(3)
                kBest = k(i) + dx;
            end
        end
    end

    if abs(kBest-kRef) > maxError
        status = '候選谷值與解析根偏離較大';
    else
        status = 'OK';
    end
end
