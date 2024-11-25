# ALF模型训练流程文档

## 1. 程序概述
此程序用于训练ALF（Adjustable Load Fleet）模型，主要通过逆向优化方法来学习可转移负荷的参数。

## 2. 主要步骤

### 2.1 数据准备
- 从预设的数据集中加载数据（支持cement、steelpowder、steelmaking等数据集）
- 加载历史电价和用电量数据
- 设置训练集大小（NOFTRAIN = 21天）

### 2.2 参数设置

```matlab
TimeLimit = 120;  % 每次迭代的最大求解时间（秒）
max_itr = 100;    % 最大迭代次数
delta_t = 1;      % 时间间隔长度
```

### 2.3 模型训练循环
程序对1-3个可转移负荷模型分别进行训练：

1. **变量初始化**
   - 调用`add_varDef_and_initVal`初始化变量和初始值

2. **迭代训练过程**
   - 随机选择训练集中的一天数据
   - 获取该天的电价和实际用电量数据
   - 执行反向优化求解

3. **参数更新**
   - 使用自适应学习率(α)更新模型参数：
     - 当迭代次数≤40时：α = 0.9
     - 当迭代次数>40时：α = 1/√(迭代次数-40)
   - 更新以下参数：
     - P_max_i：最大功率限制
     - P_min_i：最小功率限制
     - E_max_i：最大能量限制
     - E_min_i：最小能量限制

4. **结果记录**
   - 记录每次迭代的损失函数值
   - 记录更新后的参数值

### 2.4 收敛判断
满足以下条件时停止迭代：
- 迭代次数超过训练集天数
- 最近NOFTRAIN天的损失函数平均绝对值小于误差阈值（err = 1e-3）

### 2.5 结果保存
- 将训练结果保存到`results/data_[数据集名称][可转移负荷数量]ALs.mat`文件中

## 3. 输出结果
训练过程会记录以下数据：
- J_theta：损失函数值
- P_max_i：最大功率限制
- P_min_i：最小功率限制
- E_max_i：最大能量限制
- E_min_i：最小能量限制

## 4. 注意事项
- 程序使用YALMIP工具箱进行优化求解
- 每个模型训练完成后会清理YALMIP环境
- 学习率会随着迭代次数的增加而自适应调整，以确保收敛性

## 5. 代码示例

### 5.1 数据加载
```matlab
load("data_set/dataset_" + data_set_name + ".mat");
data_set.NOFTRAIN = 21;
```

### 5.2 参数更新示例
```matlab
% 学习率计算
alpha = 0.9;
if idx_itr > 40
    alpha = 1 / (idx_itr - 40)^0.5;
end

% 参数更新
P_max_i_val = (1 - alpha) * P_max_i_val + alpha * value(P_max_i);
assign(P_max_i, P_max_i_val);
```

### 5.3 收敛判断示例
```matlab
err = 1e-3;
if idx_itr > data_set.NOFTRAIN && ...
        mean(abs(result.J_theta(end - data_set.NOFTRAIN + 1:end))) < err
    break;
end
```
