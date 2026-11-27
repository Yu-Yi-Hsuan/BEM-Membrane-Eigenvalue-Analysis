%% 環形 Dual BEM：特徵波數收斂分析（內外圓解析根）
% 用途：
% 1. stage = 'prepare'：建立不同邊界元素數的 F15.DAT
% 2. 在各資料夾執行 Fortran，產生 Complex_U.DAT
% 3. stage = 'analyze'：讀取掃根結果，繪製收斂圖與相對誤差圖（TIF）

clear; clc; close all;

%% ===== 0. 使用者設定 =====
stage = 'analyze';                    % 'prepare'：產生 F15.DAT；'analyze'：分析與繪圖
R_outer = 1.0;                        % 環形外半徑
R_inner = 0.43565; % 解析節線半徑（亦可自行輸入 0.43565）
nOuterList = [70 90 110 130 150 170 190 210]; % 各組外圓元素數（自行設定）
nInnerList = [30 39 48 57 65 74 83 91];  % 各組內圓元素數（自行設定）
nNodeList = nOuterList + nInnerList;   % 自動計算總元素數（收斂圖橫軸）
rootFolder = fullfile(pwd,'Annulus_Convergence'); % 避免與圓形資料混用
sigmaFile = 'Complex_U_L.DAT';          % Fortran 掃根結果：[k, sigma_min]

% ===== 根的辨識：環形固定邊界的 Bessel 特徵方程 =====
manualExactK = [];         % []：用解析根；或直接輸入目標解析波數
selectedRoot = 1;                     % 1：環形第一組不同解析根；[]：依共同掃描範圍自動辨識
maxRootIndex = 30;                    % 列出前幾組不同解析根，包含各階角向模態
% 說明：manualExactK 有值時，會優先用它尋找對應的解析根。
%       manualExactK=[] 且 selectedRoot=[]，若掃描範圍有多個根，則列出清單供指定。
% 注意：更改 R_outer/R_inner 後，必須重新 prepare、重跑 Fortran 再 analyze。

% 谷值搜尋控制
searchHalfWidth = 0.20;               % 只在目標解析根 ± 此範圍內找谷值
maxSearchError = 0.10;                % 找到的根若偏離解析值太多，標示可疑
useQuadraticFit = false;              % false：直接取掃描谷值；true：三點二次擬合

% 圖形設定
numXTicks = 8;                        % X 軸主要刻度數量
numYTicks = 5;                        % Y 軸主要刻度數量
extraYMargin = 0.12;                  % Y 軸上下留白比例
plotError = false;                     % 是否另外畫相對誤差圖
exportTif = true;                     % 是否輸出 tif
exportResolution = 2400;              % tif 解析度（參考你原本的作圖程式）

%% ===== 1. 基本檢查 =====
assert(R_inner>0 && R_outer>R_inner && all(isfinite([R_inner,R_outer])), ...
    '需滿足 0 < R_inner < R_outer，且半徑必須為有限值。');
assert(numel(nOuterList)==numel(nInnerList) && ~isempty(nOuterList), ...
    'nOuterList 與 nInnerList 必須有相同且非零的組數。');
assert(all(nOuterList>=3 & nOuterList==fix(nOuterList)) && ...
       all(nInnerList>=3 & nInnerList==fix(nInnerList)), ...
    '內外圓各組元素數都必須是至少為 3 的整數。');
assert(numel(unique(nNodeList))==numel(nNodeList), ...
    'nNodeList 不能有重複值。');

assert(isempty(selectedRoot) || (isscalar(selectedRoot) && ...
    selectedRoot>=1 && selectedRoot==floor(selectedRoot)), ...
    'selectedRoot 必須是正整數或 []。');
assert(isempty(manualExactK) || (isscalar(manualExactK) && ...
    isfinite(manualExactK) && manualExactK>0), ...
    'manualExactK 必須是正數或 []。');
assert(isscalar(maxRootIndex) && maxRootIndex>=1 && ...
    maxRootIndex==floor(maxRootIndex),'maxRootIndex 必須是正整數。');

switch lower(stage)
    case 'prepare'
        runPrepare(rootFolder,nOuterList,nInnerList,R_outer,R_inner);
    case 'analyze'
        runAnalyze(rootFolder,nNodeList,nOuterList,nInnerList,R_outer,R_inner,sigmaFile,manualExactK,selectedRoot,maxRootIndex, ...
            searchHalfWidth,maxSearchError,useQuadraticFit, ...
            numXTicks,numYTicks,extraYMargin,plotError,exportTif,exportResolution);
    otherwise
        error('stage 只能填 prepare 或 analyze。');
end

%% ========================================================================
function runPrepare(rootFolder,nOuterList,nInnerList,R_outer,R_inner)
    nNodeList = nOuterList + nInnerList;
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
        nOuter = nOuterList(t);
        nInner = nInnerList(t);
        writeAnnulusF15(filePath,nOuter,nInner,R_outer,R_inner);
        fprintf('已建立 %s（總元素 %d：外圓 %d + 內圓 %d）\n', ...
            filePath,nb,nOuter,nInner);
    end

    fprintf('\n下一步：在每個 N 資料夾中執行 Fortran，產生 %s。\n', 'Complex_U.DAT');
    fprintf('之後把 stage 改成 analyze，再重新執行此程式。\n');
end

%% ========================================================================
function runAnalyze(rootFolder,nNodeList,nOuterList,nInnerList,R_outer,R_inner,sigmaFile,manualExactK,selectedRoot,maxRootIndex, ...
    searchHalfWidth,maxSearchError,useQuadraticFit, ...
    numXTicks,numYTicks,extraYMargin,plotError,exportTif,exportResolution)

    [kAnalytical,besselOrder,zeroIndex,multiplicity] = annulusRoots(R_outer,R_inner,maxRootIndex);

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
    candidates = find(kAnalytical>lowerBound+guard & ...
                      kAnalytical<upperBound-guard);

    if ~isempty(selectedRoot)
        assert(selectedRoot<=numel(kAnalytical), ...
            'selectedRoot 超過已列出的解析根數量。');
        rootIndex = selectedRoot;
    elseif ~isempty(manualExactK)
        % 手動給定解析值時，直接尋找最接近的 Bessel 解析根
        [delta,rootIndex] = min(abs(kAnalytical-manualExactK));
        assert(delta < max(0.02,0.005*manualExactK), ...
            'manualExactK 與前 %d 組環形解析根不符，請增加 maxRootIndex 或檢查數值。',maxRootIndex);
    elseif numel(candidates)==1
        rootIndex = candidates;
    elseif isempty(candidates)
        error('共同掃描範圍 [%.6f, %.6f] 沒有完整包含前 %d 組解析根。', ...
            lowerBound,upperBound,maxRootIndex);
    else
        fprintf('共同掃描範圍中包含多組解析根，請設定 selectedRoot 或 manualExactK：\n');
        disp(table(candidates,kAnalytical(candidates),besselOrder(candidates), ...
            zeroIndex(candidates),multiplicity(candidates), ...
            'VariableNames',{'RootIndex','kExact','BesselOrder','ZeroIndex','Multiplicity'}));
        error('請指定 selectedRoot，例如 4；或填入 manualExactK。');
    end

    formulaExact = kAnalytical(rootIndex);
    if ~isempty(manualExactK)
        assert(abs(manualExactK-formulaExact)<max(0.02,0.005*formulaExact), ...
            'selectedRoot 與 manualExactK 對應不同解析根，請確認設定。');
    end
    if isempty(manualExactK)
        referenceK = formulaExact;
    else
        referenceK = manualExactK;
    end
    assert(referenceK>lowerBound && referenceK<upperBound, ...
        '解析參考值 %.8f 不在共同掃描範圍 [%.6f, %.6f]。', ...
        referenceK,lowerBound,upperBound);

    fprintf('\n===== 環形解析根辨識 =====\n');
    fprintf('R_outer = %.12f；R_inner = %.12f\n',R_outer,R_inner);
    fprintf('第 %d 組不同解析根：角向階數 n=%d，徑向序號 m=%d；物理重數 %d\n', ...
        rootIndex,besselOrder(rootIndex),zeroIndex(rootIndex),multiplicity(rootIndex));
    fprintf('公式解析值 = %.12f；作圖參考值 = %.12f\n',formulaExact,referenceK);
    fprintf('共同掃描範圍：[%.6f, %.6f]\n',lowerBound,upperBound);

    kNumerical = nan(nCases,1);
    minSigma = nan(nCases,1);
    status = repmat({''},nCases,1);

    for t = validCases(:).'
        A = scans{t};
        [kNumerical(t),minSigma(t),status{t}] = trackMinimum(A(:,1),A(:,2), ...
            referenceK,searchHalfWidth,maxSearchError,useQuadraticFit);
        if ~strcmp(status{t},'OK')
            warning('N = %d：%s',nNodeList(t),status{t});
        end
    end

    relativeError = abs(kNumerical-referenceK)/referenceK*100;
    outerCounts=nOuterList(:);
    innerCounts=nInnerList(:);
    Result = table(nNodeList(:),outerCounts,innerCounts,kNumerical, ...
        repmat(referenceK,nCases,1),relativeError,minSigma,status, ...
        'VariableNames',{'BoundaryElements','OuterElements','InnerElements', ...
        'k_BEM','k_Analytical','RelativeError_percent','MinimumSigma','Status'});
    disp(Result);

    % ===== 圖 1：特徵波數收斂圖 =====
    fig1 = figure('Color','w','Name','Annulus convergence');
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
    ylabel('k');
    formatAxis(nNodeList,numXTicks,numYTicks,extraYMargin,[referenceK; kNumerical(isfinite(kNumerical))]);
    legend('Location','best');
    box on;
    grid off;

    if exportTif
        figure(fig1);
        drawnow;
        disp('請先調整第一張收斂圖的位置、視窗大小與圖例，完成後回 Command Window 按 Enter 儲存...');
        pause;
        exportgraphics(fig1,fullfile(rootFolder,sprintf('Annulus_k%d_Convergence.tif',rootIndex)), ...
            'Resolution',exportResolution);
        fprintf('第一張收斂圖已儲存，接著產生相對誤差圖。\n');
    end

    % ===== 圖 2：相對誤差圖 =====
    if plotError
        fig2 = figure('Color','w','Name','Annulus relative error');
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
            exportgraphics(fig2,fullfile(rootFolder,sprintf('Annulus_k%d_RelativeError.tif',rootIndex)), ...
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
function writeAnnulusF15(filename,nOuter,nInner,R_outer,R_inner)
    % 依 Annulus_infiled(1).m：外圓逆時針、內圓順時針。
    thetaOuter = 2*pi*(0:nOuter-1)'/nOuter;
    thetaInner = -2*pi*(0:nInner-1)'/nInner;
    x = [R_outer*cos(thetaOuter); R_inner*cos(thetaInner)];
    y = [R_outer*sin(thetaOuter); R_inner*sin(thetaInner)];
    nTotal = nOuter+nInner;

    fid = fopen(filename,'wt');
    assert(fid>0,'無法寫入 %s',filename);
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid,'%6d\n%6d\n',-1,15);
    for i=1:nTotal
        fprintf(fid,'%10d%10d%10d%10d%13.5E%13.5E%13.5E\n', ...
            i,0,0,11,x(i),y(i),0);
    end

    fprintf(fid,'%6d\n%6d\n%6d\n',-1,-1,71);
    for i=1:nOuter
        fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d\n',i,1,21,1,1,7,2);
        fprintf(fid,'%10d%10d\n',i,mod(i,nOuter)+1);
    end
    for j=1:nInner
        currentNode = nOuter+j;
        nextNode = nOuter+mod(j,nInner)+1;
        fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d\n',currentNode,1,21,1,1,7,2);
        fprintf(fid,'%10d%10d\n',currentNode,nextNode);
    end
    fprintf(fid,'%6d\n',-1);
end

%% ========================================================================
function [rootK,besselOrder,zeroIndex,multiplicity] = annulusRoots(R_outer,R_inner,numRoots)
    % 同心環形 Dirichlet 特徵方程：
    % J_n(k*R_inner)*Y_n(k*R_outer) - J_n(k*R_outer)*Y_n(k*R_inner)=0。
    % 使用正波數 k 掃描括根後以 fzero 求根；不依賴外部 besselzero。
    % n=0 為單重模態；n>=1 為角向雙重模態。
    kMax = max(15/R_outer, 2*pi/(R_outer-R_inner));
    while true
        kGrid = linspace(0.1/R_outer,kMax,max(1500,ceil(100*kMax*R_outer)))';
        values = zeros(0,2);
        % 對 Dirichlet 環形，角向階數 n 高於 k*R_outer 時不會有特徵根。
        nMax = ceil(kMax*R_outer)+1;
        for n=0:nMax
            y = annulusCharacteristic(kGrid,n,R_outer,R_inner);
            crossings = find(isfinite(y(1:end-1)) & isfinite(y(2:end)) & ...
                y(1:end-1).*y(2:end)<0);
            for q=1:numel(crossings)
                i=crossings(q);
                z=fzero(@(k)annulusCharacteristic(k,n,R_outer,R_inner), ...
                    [kGrid(i),kGrid(i+1)]);
                if isfinite(z) && z>0
                    values(end+1,:)=[z,n]; %#ok<AGROW>
                end
            end
        end
        values=sortrows(values,1);
        uniqueRows=zeros(0,2);
        for p=1:size(values,1)
            if isempty(uniqueRows) || ...
                    abs(values(p,1)-uniqueRows(end,1))>1e-8*max(1,values(p,1))
                uniqueRows(end+1,:)=values(p,:); %#ok<AGROW>
            end
        end
        if size(uniqueRows,1)>=numRoots
            break;
        end
        kMax=1.5*kMax;
        assert(kMax<500/R_outer,'環形解析根搜尋範圍過大，請檢查半徑及 maxRootIndex。');
    end

    uniqueRows=uniqueRows(1:numRoots,:);
    rootK=uniqueRows(:,1);
    besselOrder=uniqueRows(:,2);
    zeroIndex=zeros(numRoots,1);
    for n=unique(besselOrder)'
        ii=find(besselOrder==n);
        zeroIndex(ii)=(1:numel(ii))';
    end
    multiplicity=ones(numRoots,1);
    multiplicity(besselOrder>=1)=2;
end

%% ========================================================================
function f=annulusCharacteristic(k,n,R_outer,R_inner)
    f=besselj(n,k.*R_inner).*bessely(n,k.*R_outer) - ...
      besselj(n,k.*R_outer).*bessely(n,k.*R_inner);
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
