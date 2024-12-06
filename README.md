# Data-Driven-Dimension-Reduction

## Project Structure

### Data Sets (`data_set/`)
Contains simulation data for industrial processes modeled using:
- Cement plant (STN-based modeling)
- Steel powder plant (STN-based modeling)
- Steel plant (RTN-based modeling)

Each dataset includes detailed documentation and is used for ALF model training.

### Results (`results/`)
Contains analysis outputs and model performance data generated using the datasets.

### Visualization (`visualization/`)
Tools and scripts for analyzing and visualizing results.

### OVB (`OVB/`)
Contains the reproduced OVB methods. Please refer to the documentation inside for detailed explanations.

## Core Components

### Training Pipeline
1. `main_training.m`
   - Main entry point for ALF model training
   - Handles training across all three datasets

2. `model_training.m`
   - Implements iterative training algorithm for ALF models
   - Core model optimization logic

### Model Configuration
3. `add_varDef_and_initVal.m`
   - Defines and initializes variables for inverse optimization
   - Separated from constraint construction for flexibility

### Optimization Modules
4. `inverse_optimization.m`
   - Fits ALF model parameters using electricity meter data
   - Implements inverse optimization approach

5. `primal_problem.m`
   - Solves the original optimization problem
   - Determines optimal energy consumption based on ALF

### Testing and Validation
6. `test_reduced_constraints_ALF.m`
   - Evaluates optimal energy consumption on test set
   - Uses identified ALF model parameters

7. `test_reduced_constraints_SLF.m`
   - Computes optimal energy consumption using simplified model (SAL)
   - Provides baseline comparison

## Getting Started
1. Ensure all dependencies are installed:
   - Gurobi
   - YALMIP
2. Place your data in the `data_set` folder
3. Run `main_training.m` to begin the training process
4. View results in the `results` folder
5. Use visualization tools to analyze outputs

---

# Data-Driven-Dimension-Reduction (中文)

## 项目结构

### 数据集 (`data_set/`)
包含使用以下方式建模的工业过程仿真数据：
- 水泥厂（基于STN建模）
- 钢粉厂（基于STN建模）
- 钢铁厂（基于RTN建模）

每个数据集都包含详细文档，用于ALF模型训练。

### 结果 (`results/`)
包含使用数据集生成的分析输出和模型性能数据。

### 可视化 (`visualization/`)
用于分析和可视化结果的工具和脚本。

### OVB (`OVB/`)
包含复现的OVB方法。详细说明请参见内部文档。

## 核心组件

### 训练流程
1. `main_training.m`
   - ALF模型训练的主入口点
   - 处理所有三个数据集的训练

2. `model_training.m`
   - 实现ALF模型的迭代训练算法
   - 核心模型优化逻辑

### 模型配置
3. `add_varDef_and_initVal.m`
   - 定义和初始化逆优化变量
   - 与约束构建分离以提高灵活性

### 优化模块
4. `inverse_optimization.m`
   - 使用电表数据拟合ALF模型参数
   - 实现逆优化方法

5. `primal_problem.m`
   - 求解原始优化问题
   - 基于ALF确定最优能耗

### 测试和验证
6. `test_reduced_constraints_ALF.m`
   - 评估测试集上的最优能耗
   - 使用已识别的ALF模型参数

7. `test_reduced_constraints_SLF.m`
   - 使用简化模型(SAL)计算最优能耗
   - 提供基准比较

## 开始使用
1. 确保安装所有依赖项：
   - Gurobi
   - YALMIP
2. 将数据放入`data_set`文件夹
3. 运行`main_training.m`开始训练过程
4. 在`results`文件夹中查看结果
5. 使用可视化工具分析输出

---

