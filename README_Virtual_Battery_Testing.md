# 电解铝负荷降维可行域（虚拟储能）模型测试说明

## 概述

本目录包含了电解铝负荷降维可行域（虚拟储能）模型的测试实现，参考了 `test_reduced_constraints_ALF.m` 的结构，但针对虚拟储能模型进行了适配。

## 文件说明

### 1. 主要测试文件

#### `test_reduced_constraints_virtual_battery.m`
- **功能**: 电解铝虚拟储能模型的主要测试文件
- **参考**: `test_reduced_constraints_ALF.m`
- **特点**: 
  - 测试4种不同批次大小的虚拟储能模型（1-4天）
  - 在10天测试集上评估模型性能
  - 计算多种性能指标（MAE、MAPE、RMSE）
  - 记录求解时间和状态

#### `primal_problem_virtual_battery.m`
- **功能**: 虚拟储能模型的原始问题求解
- **参考**: `primal_problem.m`
- **特点**:
  - 实现虚拟储能模型的状态转移约束
  - 包含功率和能量状态约束
  - 使用Gurobi求解器优化

#### `test_virtual_battery_simple.m`
- **功能**: 简单测试脚本，验证文件完整性和基本功能
- **用途**: 在运行完整测试前检查环境

### 2. 依赖文件

- `add_varDef_and_initVal_virtual_battery.m`: 虚拟储能模型变量定义
- `inverse_optimization_virtual_battery.m`: 虚拟储能模型逆优化
- `model_training_eal.m`: 虚拟储能模型训练脚本
- `data_set/EAL_Implementation/dataset_eal.mat`: 电解铝数据集

## 使用方法

### 1. 环境准备

确保已安装以下工具：
- MATLAB
- YALMIP工具箱
- Gurobi求解器

### 2. 训练模型（如果尚未训练）

```matlab
% 运行虚拟储能模型训练
model_training_eal
```

这将生成4个训练结果文件：
- `results/data_eal_virtual_battery_1batch.mat`
- `results/data_eal_virtual_battery_2batch.mat`
- `results/data_eal_virtual_battery_3batch.mat`
- `results/data_eal_virtual_battery_4batch.mat`

### 3. 运行测试

#### 方法1：简单测试（推荐先运行）
```matlab
test_virtual_battery_simple
```

#### 方法2：完整测试
```matlab
test_reduced_constraints_virtual_battery
```

### 4. 查看结果

测试结果保存在 `results/` 目录下：
- `test_virtual_battery_eal_1batch.mat`
- `test_virtual_battery_eal_2batch.mat`
- `test_virtual_battery_eal_3batch.mat`
- `test_virtual_battery_eal_4batch.mat`

每个结果文件包含：
- `E_reduced_constraints`: 模型预测的功率消耗
- `E_actual`: 实际功率消耗
- `optimization_times`: 求解时间
- `optimization_status`: 求解状态
- `mae`, `mape`, `rmse`: 性能指标
- `model_params`: 模型参数

## 与ALF模型的对比

| 特性 | ALF模型 | 虚拟储能模型 |
|------|---------|-------------|
| 模型数量 | 1-3个 | 1个（固定） |
| 约束类型 | 功率+能量总量 | 功率+能量状态+状态转移 |
| 状态耦合 | 无 | 有（通过状态转移方程） |
| 物理意义 | 可调负荷 | 虚拟储能系统 |
| 参数数量 | 4个（P_max, P_min, E_max, E_min） | 6个（+theta, w） |

## 虚拟储能模型特点

### 1. 状态转移方程
```
e_{t+1} = θ * e_t + p_t * Δt + w * Δt
```
其中：
- `θ`: 自放电系数（0.8-1.0）
- `w`: 基础功率消耗
- `e_t`: 能量状态
- `p_t`: 功率消耗

### 2. 约束条件
- 功率约束：`p_t ∈ [p_min, p_max]`
- 能量约束：`e_t ∈ [e_min, e_max]`
- 状态转移约束：上述状态转移方程

### 3. 物理意义
- 将280个电解槽视为一个虚拟储能系统
- 考虑电解槽的自放电特性
- 包含维持运行的基础功率消耗

## 性能指标说明

- **MAE (Mean Absolute Error)**: 平均绝对误差
- **MAPE (Mean Absolute Percentage Error)**: 平均绝对百分比误差
- **RMSE (Root Mean Square Error)**: 均方根误差
- **求解时间**: 每次优化的计算时间
- **求解状态**: 0表示成功，其他值表示失败

## 注意事项

1. 确保数据集文件路径正确
2. 训练过程可能需要较长时间（建议先用简单测试验证）
3. 如果求解失败，检查Gurobi许可证和参数设置
4. 虚拟储能模型比ALF模型计算复杂度更高

## 故障排除

### 常见问题

1. **文件缺失错误**
   - 检查文件路径是否正确
   - 确保已运行模型训练

2. **求解器错误**
   - 检查Gurobi许可证
   - 验证YALMIP安装

3. **内存不足**
   - 减少批次大小
   - 增加系统内存

4. **求解失败**
   - 检查模型参数合理性
   - 调整求解器设置
