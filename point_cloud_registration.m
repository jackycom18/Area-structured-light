% =========================================================================
% 点云匹配（拼接）程序
% 步骤：
%   1. 基于标志点的粗匹配（SVD 刚体变换求解）
%   2. 基于点到面 ICP 的精匹配
% =========================================================================

clear; close all; clc;

%% ======================== 参数设置 ========================
% 数据文件路径（请根据实际情况修改）
sourcePath = 'source_cloud.ply';   % 源点云（待变换的点云）
targetPath = 'target_cloud.ply';   % 目标点云（固定参考点云）

% ICP 参数
icpMaxIter   = 50;       % 最大迭代次数
icpTolerance = 1e-6;     % 收敛阈值（变换矩阵变化量）
icpMaxDist   = inf;      % 对应点最大距离阈值（inf 表示不限制）

% 标志点数量（需要在两个点云中分别选取相同数量的对应点）
numMarkers = 4;          % 至少 3 个标志点

%% ======================== 加载点云 ========================
fprintf('加载点云数据...\n');

% ---- 方式1：从文件加载 ----
% sourcePts = pcread(sourcePath);
% targetPts = pcread(targetPath);
% sourceCloud = sourcePts.Location;
% targetCloud = targetPts.Location;

% ---- 方式2：生成示例数据（供测试用）----
fprintf('使用示例数据（请替换为实际数据加载代码）...\n');
% 生成一个兔子形状的模拟源点云
rng(42);
theta = linspace(0, 2*pi, 500)';
sourceCloud = [5*cos(theta).*sin(theta*2), 5*sin(theta).*sin(theta), ...
               3*cos(theta*3)*0.5 + 2];

% 对源点云做一个已知的刚体变换，生成目标点云
R_true = eul2rotm([0.3, 0.5, 0.2]);  % 真实旋转（用于测试）
t_true = [2; -1; 3]';                % 真实平移（用于测试）
targetCloud = (R_true * sourceCloud' + t_true')';

fprintf('源点云点数: %d, 目标点云点数: %d\n', size(sourceCloud,1), size(targetCloud,1));

%% ======================== 可视化原始点云 ========================
figure('Name', '原始点云', 'NumberTitle', 'off');
pcshow(sourceCloud, [1 0 0]);   % 红色 - 源点云
hold on;
pcshow(targetCloud, [0 0 1]);   % 蓝色 - 目标点云
title('原始点云（红：源点云  蓝：目标点云）');
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; grid on; view(3);
legend('Source', 'Target');

%% ======================== 步骤1：基于标志点的粗匹配 ========================
fprintf('\n========== 步骤1：基于标志点的粗匹配 ==========\n');
fprintf('请在两个点云中分别选取 %d 个对应的标志点\n', numMarkers);
fprintf('操作方法：在图上点击选取，先选源点云上的标志点，再选目标点云上对应的标志点\n');

% ---- 方式1：手动交互选取 ----
% figure('Name', '选取标志点', 'NumberTitle', 'off');
% pcshow(sourceCloud, [1 0 0]);
% hold on; pcshow(targetCloud, [0 0 1]);
% title('请先在源点云（红色）选取标志点，按回车确认后，再在目标点云（蓝色）选取对应点');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% axis equal; grid on; view(3);
% 
% % 选取源点云标志点
% fprintf('请选取源点云（红色）上的 %d 个标志点...\n', numMarkers);
% [srcMarkers, ~] = selectMarkers(sourceCloud, numMarkers, [1 0 0]);
% 
% % 选取目标点云标志点
% fprintf('请选取目标点云（蓝色）上对应的 %d 个标志点...\n', numMarkers);
% [tgtMarkers, ~] = selectMarkers(targetCloud, numMarkers, [0 0 1]);

% ---- 方式2：自动使用最近点模拟标志点（仅用于测试）----
fprintf('【测试模式】自动选取标志点...\n');
% 在目标点云中选几个点作为标志点
markerIdx = round(linspace(1, size(targetCloud,1), numMarkers+2));
markerIdx = markerIdx(2:end-1);  % 去掉首尾
tgtMarkers = targetCloud(markerIdx, :);
% 通过真实变换反算源点云中的对应标志点
srcMarkers = ((tgtMarkers - t_true) / R_true);

fprintf('标志点选取完成。\n');
fprintf('源点云标志点:\n'); disp(srcMarkers);
fprintf('目标点云标志点:\n'); disp(tgtMarkers);

% 计算粗匹配变换矩阵（SVD 方法）
[R_coarse, t_coarse, rmsError_coarse] = computeRigidTransform(srcMarkers, tgtMarkers);

fprintf('\n粗匹配结果:\n');
fprintf('旋转矩阵 R:\n'); disp(R_coarse);
fprintf('平移向量 t: '); disp(t_coarse');
fprintf('标志点 RMS 误差: %.6f\n', rmsError_coarse);

% 应用粗匹配变换到源点云
sourceCoarseAligned = (R_coarse * sourceCloud' + t_coarse')';

% 可视化粗匹配结果
figure('Name', '粗匹配结果', 'NumberTitle', 'off');
pcshow(sourceCoarseAligned, [1 0.6 0]);  % 橙色 - 粗匹配后的源点云
hold on;
pcshow(targetCloud, [0 0 1]);            % 蓝色 - 目标点云
title('粗匹配结果（橙：变换后源点云  蓝：目标点云）');
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; grid on; view(3);
legend('Source (Coarse)', 'Target');

%% ======================== 步骤2：点到面 ICP 精匹配 ========================
fprintf('\n========== 步骤2：点到面 ICP 精匹配 ==========\n');

[R_icp, t_icp, icpHistory] = pointToPlaneICP(sourceCoarseAligned, targetCloud, ...
    'MaxIterations', icpMaxIter, ...
    'Tolerance', icpTolerance, ...
    'MaxDistance', icpMaxDist, ...
    'Verbose', true);

% 组合粗匹配和精匹配的变换
R_final = R_icp * R_coarse;
t_final = R_icp * t_coarse' + t_icp';
t_final = t_final';

% 应用最终变换
sourceFinalAligned = (R_final * sourceCloud' + t_final')';

fprintf('\n最终变换结果:\n');
fprintf('总旋转矩阵 R:\n'); disp(R_final);
fprintf('总平移向量 t: '); disp(t_final);

%% ======================== 结果可视化与评估 ========================
% 最终对齐结果
figure('Name', '最终匹配结果', 'NumberTitle', 'off');
pcshow(sourceFinalAligned, [0 1 0]);  % 绿色 - 最终对齐的源点云
hold on;
pcshow(targetCloud, [0 0 1]);         % 蓝色 - 目标点云
title('最终匹配结果（绿：最终对齐源点云  蓝：目标点云）');
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; grid on; view(3);
legend('Source (Final)', 'Target');

% 对比：原始 vs 粗匹配 vs 精匹配
figure('Name', '匹配过程对比', 'NumberTitle', 'off');
subplot(1,3,1);
pcshow(sourceCloud, [1 0 0]); hold on; pcshow(targetCloud, [0 0 1]);
title('原始位置'); axis equal; view(3);
subplot(1,3,2);
pcshow(sourceCoarseAligned, [1 0.6 0]); hold on; pcshow(targetCloud, [0 0 1]);
title('粗匹配后'); axis equal; view(3);
subplot(1,3,3);
pcshow(sourceFinalAligned, [0 1 0]); hold on; pcshow(targetCloud, [0 0 1]);
title('ICP精匹配后'); axis equal; view(3);

% ICP 收敛曲线
if ~isempty(icpHistory.RMSE)
    figure('Name', 'ICP 收敛曲线', 'NumberTitle', 'off');
    plot(0:length(icpHistory.RMSE), [icpHistory.RMSE(1); icpHistory.RMSE(:)], ...
         'b-o', 'LineWidth', 1.5);
    xlabel('迭代次数'); ylabel('RMSE');
    title('ICP 迭代收敛曲线');
    grid on;
end

% 计算并显示误差统计
distances = sqrt(sum((sourceFinalAligned - targetCloud(...
    knnsearch(targetCloud, sourceFinalAligned), :)).^2, 2));
fprintf('\n========== 最终误差统计 ==========\n');
fprintf('平均距离误差: %.6f\n', mean(distances));
fprintf('最大距离误差: %.6f\n', max(distances));
fprintf('RMS 误差:     %.6f\n', rms(distances));
fprintf('标准差:       %.6f\n', std(distances));

fprintf('\n点云匹配完成！\n');

%% ======================== 保存结果 ========================
% 保存对齐后的点云
% ptCloudAligned = pointCloud(sourceFinalAligned);
% pcwrite(ptCloudAligned, 'aligned_cloud.ply');

% 保存变换矩阵
% save('transform_matrix.mat', 'R_final', 't_final', 'R_coarse', 't_coarse', 'R_icp', 't_icp');


% =========================================================================
%                         辅助函数
% =========================================================================

%% 函数：手动选取标志点
function [markers, indices] = selectMarkers(pointCloud, numPoints, color)
    % 在点云可视化窗口中手动点击选取标志点
    % 输入：
    %   pointCloud - Nx3 点云坐标
    %   numPoints  - 需要选取的点数
    %   color      - 显示颜色（用于区分）
    % 输出：
    %   markers    - 选取的标志点坐标
    %   indices    - 标志点在原始点云中的索引
    
    markers = zeros(numPoints, 3);
    indices = zeros(numPoints, 1);
    
    % 使用 datacursormode 或者 ginput3d 方式选取
    % 这里使用旋转视图后按点的索引（简化方式）
    % 实际使用中建议使用 ginput 或专门的选点工具
    view(3);
    rotate3d on;
    
    for i = 1:numPoints
        fprintf('  选取第 %d/%d 个点...\n', i, numPoints);
        rotate3d off;
        [x, y] = ginput(1);  % 获取2D屏幕坐标
        % 注意：ginput 获取的是2D坐标，需要通过视图变换反算3D坐标
        % 这里简化为选择最近投影点
        % 获取当前视图
        ax = gca;
        pos = ax.CurrentPoint;
        
        % 使用 data tips 或基于投影选择最近点
        % 简化版本：将2D坐标近似匹配
        camPos = campos;
        camTarget = camtarget;
        camVec = camTarget - camPos;
        camVec = camVec / norm(camVec);
        
        % 计算各点到射线的距离，选最近的点
        ptProj = pointCloud - camPos;
        projLen = ptProj * camVec';
        ptOnRay = camPos + projLen * camVec;
        distToRay = sqrt(sum((pointCloud - ptOnRay).^2, 2));
        [~, idx] = min(distToRay);
        
        markers(i, :) = pointCloud(idx, :);
        indices(i) = idx;
        
        % 在图上标记选中的点
        hold on;
        plot3(markers(i,1), markers(i,2), markers(i,3), ...
              'o', 'MarkerSize', 12, 'MarkerFaceColor', 'g', ...
              'MarkerEdgeColor', 'k', 'LineWidth', 2);
        text(markers(i,1), markers(i,2), markers(i,3), ...
             sprintf(' %d', i), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'g');
        drawnow;
    end
    rotate3d on;
end

%% 函数：SVD 求解刚体变换（Kabsch / Umeyama 算法）
function [R, t, rmsError] = computeRigidTransform(P, Q)
    % 使用 SVD 方法计算两点集之间的最优刚体变换（旋转 + 平移）
    % 输入：
    %   P - 源点集 (Nx3)
    %   Q - 目标点集 (Nx3)，与 P 逐点对应
    % 输出：
    %   R        - 3x3 旋转矩阵
    %   t        - 3x1 平移向量
    %   rmsError - 变换后对应点的 RMS 误差
    
    n = size(P, 1);
    
    % 1. 计算质心
    centroidP = mean(P, 1);
    centroidQ = mean(Q, 1);
    
    % 2. 中心化
    Pc = P - centroidP;
    Qc = Q - centroidQ;
    
    % 3. 计算协方差矩阵 H = Pc' * Qc
    H = Pc' * Qc;
    
    % 4. SVD 分解
    [U, ~, V] = svd(H);
    
    % 5. 计算旋转矩阵（处理反射情况）
    R = V * U';
    
    % 确保 det(R) = 1（不是反射），如果为 -1 则修正
    if det(R) < 0
        V(:, 3) = -V(:, 3);
        R = V * U';
    end
    
    % 6. 计算平移向量
    t = centroidQ' - R * centroidP';
    
    % 7. 计算 RMS 误差
    P_transformed = (R * P' + t)';
    errors = sqrt(sum((P_transformed - Q).^2, 2));
    rmsError = rms(errors);
end

%% 函数：点到面 ICP（Point-to-Plane ICP）
function [R, t, history] = pointToPlaneICP(source, target, varargin)
    % 点到面 ICP 精配准算法
    % 基于 Low (2004) 的线性近似方法
    %
    % 输入：
    %   source  - 源点云 (Nx3)，已经过粗匹配
    %   target  - 目标点云 (Mx3)，固定参考
    % 参数（名值对）：
    %   'MaxIterations' - 最大迭代次数（默认 50）
    %   'Tolerance'     - 收敛容差（默认 1e-6）
    %   'MaxDistance'   - 对应点最大距离（默认 inf）
    %   'Verbose'       - 是否打印迭代信息（默认 false）
    % 输出：
    %   R       - 3x3 旋转矩阵
    %   t       - 1x3 平移向量
    %   history - 迭代历史记录
    
    % 解析参数
    p = inputParser;
    addParameter(p, 'MaxIterations', 50, @(x) isscalar(x) && x > 0);
    addParameter(p, 'Tolerance', 1e-6, @(x) isscalar(x) && x > 0);
    addParameter(p, 'MaxDistance', inf, @(x) isscalar(x) && x > 0);
    addParameter(p, 'Verbose', false, @islogical);
    addParameter(p, 'NormalK', 10, @(x) isscalar(x) && x > 0);              % 法向量估计和边缘检测共用邻域大小
    addParameter(p, 'EnableEdgeFilter', true, @islogical);                   % 是否启用边缘点过滤
    addParameter(p, 'EdgeAngleThreshold', 90, @(x) isscalar(x) && x > 0);    % 边缘判定角度阈值（度）
    addParameter(p, 'EnableNormalFilter', true, @islogical);                  % 是否启用法向量一致性过滤（部分重叠场景推荐开启）
    addParameter(p, 'NormalAngleThreshold', 45, @(x) isscalar(x) && x > 0);   % 法向量夹角阈值（度，论文推荐30~45）
    addParameter(p, 'EnableTrimmed', true, @islogical);                       % 是否启用Trimmed ICP（部分重叠场景推荐开启）
    addParameter(p, 'OverlapRatio', 0, @(x) isscalar(x) && x >= 0 && x <= 1); % 重叠率（0=自动估计，论文Chetverikov et al.）
    parse(p, varargin{:});
    
    maxIter      = p.Results.MaxIterations;
    tol          = p.Results.Tolerance;
    maxDist      = p.Results.MaxDistance;
    verbose      = p.Results.Verbose;
    normalK      = p.Results.NormalK;
    enableEdge   = p.Results.EnableEdgeFilter;
    edgeAngleTh  = p.Results.EdgeAngleThreshold;
    enableNormal = p.Results.EnableNormalFilter;         % ← 新增
    normalAngleTh = p.Results.NormalAngleThreshold;      % ← 新增
    enableTrimmed = p.Results.EnableTrimmed;             % ← 新增
    overlapRatio  = p.Results.OverlapRatio;              % ← 新增
    
    % 初始化
    R = eye(3);
    t = [0, 0, 0];
    
    % 历史记录
    history.RMSE = [];
    history.TransformChange = [];
    
    % 将源点云复制，后续迭代中不断更新
    sourceTransformed = source;
    
    % 预计算目标点云的法向量（使用 PCA 方法）
    if verbose
        fprintf('  计算目标点云法向量...\n');
    end
    targetNormals = estimateNormals(target, normalK);
    
    % 预计算源点云法向量（用于法向量一致性过滤）
    sourceNormals = [];            % ← 新增
    if enableNormal
        if verbose
            fprintf('  计算源点云法向量...\n');
        end
        sourceNormalsOrig = estimateNormals(source, normalK);
    end
    
    % 构建目标点云的 KD-tree（使用 MATLAB 内置函数或暴力搜索）
    % 对于大规模点云，建议使用 createns 或者 KDTreeSearcher
    if size(target, 1) > 10000
        kdtree = createns(target, 'NSMethod', 'kdtree');
    else
        kdtree = createns(target, 'NSMethod', 'exhaustive');
    end
    
    if verbose
        fprintf('  开始 ICP 迭代...\n');
        fprintf('  %-6s %-14s %-14s\n', 'Iter', 'RMSE', 'Delta');
        fprintf('  %s\n', repmat('-', 1, 40));
    end
        
    % 预计算目标点云边缘点（迭代前一次性算好，使用与法向量估计相同的 normalK）
    if enableEdge
        if verbose
            fprintf('  检测目标点云边缘点 (k=%d, 角度阈值=%.1f°)...\n', normalK, edgeAngleTh);
        end
        targetIsEdge = detectEdgePoints(target, targetNormals, normalK, edgeAngleTh);
        if verbose
            fprintf('  边缘点: %d / %d (%.1f%%)\n', ...
                sum(targetIsEdge), size(target,1), 100*sum(targetIsEdge)/size(target,1));
        end
    end
    
    % 自动估计重叠率（Chetverikov TrICP）
    if enableTrimmed && overlapRatio == 0
        if verbose
            fprintf('  自动估计点云重叠率...\n');
        end
        overlapRatio = estimateOverlap(source, target, kdtree, normalK);
        if verbose
            fprintf('  估计重叠率: %.1f%%\n', overlapRatio * 100);
        end
    end
    
    for iter = 1:maxIter
        % 1. 查找最近对应点
        [indices, distances] = knnsearch(kdtree, sourceTransformed);
        
        % 2. 根据距离阈值过滤
        validMask = distances <= maxDist;
        
        % 2.5 过滤掉匹配到边缘目标点的对应关系
        if enableEdge
            targetIdx = indices(validMask);
            edgeMask = targetIsEdge(targetIdx);
            validMask(validMask) = ~edgeMask;
        end
        
        % 2.6 法向量一致性过滤（论文：PCL CorrespondenceRejectorSurfaceNormal）
        if enableNormal
            % 将源法向量旋转到当前姿态
            srcNrm = (R * sourceNormalsOrig(validMask, :)')';
            tgtNrm = targetNormals(indices(validMask), :);
            % 计算法向量夹角（度）
            cosAngles = abs(sum(srcNrm .* tgtNrm, 2));
            angles = acosd(min(cosAngles, 1));  % acosd: 反余弦返回度数
            % 只保留夹角小于阈值的配对
            validMask(validMask) = angles <= normalAngleTh;
        end
        
        validSrc  = sourceTransformed(validMask, :);
        validTgt  = target(indices(validMask), :);
        validNrm  = targetNormals(indices(validMask), :);
        validDist = distances(validMask);
        
        if sum(validMask) < 6
            warning('有效对应点太少 (%d)，ICP 提前终止。', sum(validMask));
            break;
        end
        
        % 2.7 Trimmed ICP：只保留距离最小的前 overlapRatio*N 对（Chetverikov et al.）
        if enableTrimmed && overlapRatio < 1
            nKeep = max(6, round(sum(validMask) * overlapRatio));
            [sortedDist, sortIdx] = sort(validDist);
            keepIdx = sortIdx(1:nKeep);
            validSrc  = validSrc(keepIdx, :);
            validTgt  = validTgt(keepIdx, :);
            validNrm  = validNrm(keepIdx, :);
            validDist = validDist(keepIdx);
        end
        
        % 3. 计算当前 RMSE（仅记录，不用于收敛判断）
        currentRMSE = sqrt(mean(validDist.^2));
        history.RMSE = [history.RMSE; currentRMSE];
        
        if verbose
            fprintf('  %-6d %-14.6f\n', iter, currentRMSE);
        
        end
        % 4. 构建点到面的线性系统并求解
        %    对于每个对应点对 (s_i, d_i, n_i)：
        %    最小化 Σ||(R·s_i + t - d_i)·n_i||²
        %    
        %    使用小角度近似：R ≈ I + α^
        %    其中 α = [α_x, α_y, α_z] 是小旋转角
        %    α^ 是 α 的反对称矩阵：
        %        [ 0   -α_z  α_y ]
        %        [ α_z  0   -α_x]
        %        [-α_y  α_x  0   ]
        %
        %    线性化后：(s_i × n_i)·α + n_i·t = n_i·(d_i - s_i)
        %    即：A * x = b, 其中 x = [α; t]
        
        nPts = size(validSrc, 1);
        A = zeros(nPts, 6);  % [α_x, α_y, α_z, t_x, t_y, t_z]

        % 鲁棒核尺度：取对应点距离的中位数（Cauchy 核）
        sigma = max(median(validDist), 1e-6);
        b = zeros(nPts, 1);
        
        for i = 1:nPts
            si = validSrc(i, :)';
            di = validTgt(i, :)';
            ni = validNrm(i, :)';
            
            % 权重：Cauchy 鲁棒核，抑制大距离/噪声点对解的支配
            wi = 1 / (1 + (validDist(i) / sigma)^2);
            
            % s_i × n_i  (叉积)
            crossProd = cross(si, ni);
            
            A(i, 1:3) = wi * crossProd';     % (s_i × n_i) 对于 α
            A(i, 4:6) = wi * ni';            % n_i 对于 t
            
            b(i) = wi * (ni' * (di - si));     % n_i · (d_i - s_i)
        end
        
        % 求解线性最小二乘问题 A * x = b
        x = A \ b;
        
        alpha = x(1:3);  % 旋转增量（小角度）
        dt    = x(4:6)'; % 平移增量

        % 大旋转警告：小角度线性化假设失效
        if norm(alpha) > deg2rad(5)
            warning('点面ICP：单步旋转增量 %.1f°，超出小角度线性化假设，结果可能不可靠。', rad2deg(norm(alpha)));
        end
        
        % 5. 将小角度旋转向量转换为旋转矩阵（Rodrigues 公式）
        dR = rotationVectorToMatrix(alpha);
        
        % 6. 更新累积变换
        R = dR * R;
        t = t * dR' + dt;
        
        % 7. 更新源点云位置
        sourceTransformed = (dR * sourceTransformed' + dt')';
        
        % 5. 收敛判断：变换变化量（比 RMSE 变化更可靠）
        deltaTransform = norm(alpha) + norm(dt);
        history.TransformChange = [history.TransformChange; deltaTransform];
        if deltaTransform < tol
            if verbose
                fprintf('  ICP 收敛于第 %d 次迭代（变换变化量 %.3e）。\n', iter, deltaTransform);
            end
            break;
        end
    end
    
    if verbose && iter == maxIter
        fprintf('  ICP 达到最大迭代次数 %d。\n', maxIter);
    end
end

%% 函数：估计点云法向量（PCA 局部平面拟合）
function normals = estimateNormals(points, k)
    % 使用 PCA 方法估计每个点的法向量
    % 输入：
    %   points - Nx3 点云
    %   k      - 邻域点数
    % 输出：
    %   normals - Nx3 法向量
    
    n = size(points, 1);
    normals = zeros(n, 3);
    
    % 构建 KD-tree 用于邻域搜索
    if n > 10000
        kdtree = createns(points, 'NSMethod', 'kdtree');
    else
        kdtree = createns(points, 'NSMethod', 'exhaustive');
    end
    
    % 对每个点计算法向量
    for i = 1:n
        % 查找 k 个最近邻
        [indices, ~] = knnsearch(kdtree, points(i,:), 'K', min(k+1, n));
        % 去掉自身
        neighbors = points(indices(2:end), :);
        
        % 去中心化
        centroid = mean(neighbors, 1);
        centered = neighbors - centroid;
        
        % PCA：协方差矩阵的最小特征值对应的特征向量即为法向量
        covMat = centered' * centered;
        [V, D] = eig(covMat);
        [~, minIdx] = min(diag(D));
        normal = V(:, minIdx);
        
        % 归一化
        normals(i, :) = normal' / norm(normal);
    end
    
    % 统一法向量方向（使它们指向同一侧）
    % 简单做法：所有法向量指向质心的外侧
    centroid = mean(points, 1);
    for i = 1:n
        dirToCentroid = centroid - points(i, :);
        if dot(normals(i, :), dirToCentroid) < 0
            normals(i, :) = -normals(i, :);
        end
    end
end

%% 函数：旋转向量 → 旋转矩阵（Rodrigues 公式）
function R = rotationVectorToMatrix(omega)
    % 将旋转向量（轴角表示）转换为 3x3 旋转矩阵
    % omega = θ * k，其中 k 是单位旋转轴，θ 是旋转角
    %
    % R = I + sin(θ)/θ * ω^ + (1-cos(θ))/θ² * ω^²
    
    theta = norm(omega);
    
    if theta < 1e-12
        R = eye(3);
        return;
    end
    
    k = omega / theta;  % 单位旋转轴
    
    % 反对称矩阵 K = k^
    K = [0,    -k(3),  k(2);
         k(3),  0,    -k(1);
        -k(2),  k(1),  0];
    
    R = eye(3) + sin(theta) * K + (1 - cos(theta)) * (K * K);
end

%% 函数：边缘点检测（基于切平面投影的角度间隙）
function isEdge = detectEdgePoints(points, normals, k, angleThreshold)
    % 检测点云中的边缘点
    % 原理：将 k 个邻域点投影到切平面，计算投影方向的角度间隙，
    %       若最大间隙超过阈值则为边缘点
    % 输入：
    %   points         - Nx3 点云
    %   normals        - Nx3 预计算的法向量（切平面的法向）
    %   k              - 邻域点数（与法向量估计使用相同的 k）
    %   angleThreshold - 角度阈值（度），默认 90°
    % 输出：
    %   isEdge - Nx1 逻辑数组，true 表示边缘点
    
    n = size(points, 1);
    isEdge = false(n, 1);
    
    % 构建 KD-tree
    if n > 10000
        kdtree = createns(points, 'NSMethod', 'kdtree');
    else
        kdtree = createns(points, 'NSMethod', 'exhaustive');
    end
    
    % 对每个点进行边缘检测
    for i = 1:n
        % 1. 查找 k 个最近邻（含自身）
        [idx, ~] = knnsearch(kdtree, points(i,:), 'K', min(k+1, n));
        neighbors = points(idx(2:end), :);  % 去掉自身
        
        % 2. 获取该点的法向量（切平面法向）
        ni = normals(i, :)';
        
        % 3. 将邻域点投影到以该点为中心的切平面上
        %    构建切平面的两个正交基向量
        if abs(ni(3)) < 0.9
            u = cross(ni, [0; 0; 1]);
        else
            u = cross(ni, [1; 0; 0]);
        end
        u = u / norm(u);
        v = cross(ni, u);  % u, v, ni 构成正交坐标系
        
        % 邻域点相对于中心点的向量，投影到切平面 (u,v)
        projU = (neighbors - points(i,:)) * u;
        projV = (neighbors - points(i,:)) * v;
        
        % 4. 计算每个投影向量的方向角 [0, 2π)
        angles = atan2(projV, projU);
        angles = mod(angles, 2*pi);
        
        % 5. 排序角度
        angles = sort(angles);
        
        % 6. 计算相邻角度差（含首尾闭环）
        gaps = diff(angles);
        gaps = [gaps; angles(1) + 2*pi - angles(end)];
        
        % 7. 最大间隙超过阈值则为边缘点
        maxGap = max(gaps) * 180 / pi;  % 转为度数
        
        if maxGap > angleThreshold
            isEdge(i) = true;
        end
    end
end

%% 函数：自动估计点云重叠率（Chetverikov TrICP 论文方法）
function xi = estimateOverlap(source, target, kdtree, k)
    % 基于 Chetverikov et al. "Trimmed ICP" 自动估计两片点云的重叠率
    % 原理：尝试多个候选重叠率 ξ，计算 ψ(ξ) = e(ξ) / ξ^(1+λ)
    %       选使 ψ(ξ) 最小的 ξ（λ=2 为论文推荐值）
    
    lambda = 2;  % Chetverikov 推荐
    
    % 对源点云随机采样（加速）
    nSample = min(500, size(source, 1));
    rng(0);  % 固定随机种子，保证可重复
    sampleIdx = randperm(size(source, 1), nSample);
    sourceSample = source(sampleIdx, :);
    
    % 查找每个采样点的最近目标点距离
    [~, distances] = knnsearch(kdtree, sourceSample);
    distances = sort(distances);
    
    % 候选重叠率：0.4 ~ 1.0，步长 0.05
    xiCandidates = 0.4:0.05:1.0;
    bestXi = 0.5;
    bestPsi = inf;
    
    for xi = xiCandidates
        nTrim = max(6, round(nSample * xi));
        trimmedDist = distances(1:nTrim);
        e_xi = mean(trimmedDist.^2);  % trimmed MSE
        
        % ψ(ξ) = e(ξ) / ξ^(1+λ)
        psi = e_xi / (xi^(1 + lambda));
        
        if psi < bestPsi
            bestPsi = psi;
            bestXi = xi;
        end
    end
    
    xi = bestXi;
end
