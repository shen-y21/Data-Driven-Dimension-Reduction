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

# Frequently Asked Questions (FAQ)

## 1. Why choose the Adjustable Load Fleet (ALF) model over the Virtual Battery (VB) model?

The ALF model offers several advantages over the VB model:

- **Modeling Philosophy**: While the VB model focuses on load deviations from a baseline, the ALF model directly represents the actual composition of industrial loads, making it more intuitive for industrial applications
- **Parameter Flexibility**: The VB model uses only 3 parameters, whereas the ALF model offers 4I parameters (where I is the number of ALs), allowing for more accurate representation of complex industrial load characteristics
- **Physical Interpretation**: The ALF model better reflects the physical reality of industrial loads, where multiple pieces of equipment operate simultaneously with their own power and energy constraints
- **Baseline Independence**: Unlike the VB model, the ALF model operates independently of any baseline, making it more robust to changes in operating conditions

## 2. How does the ALF model handle coupling relationships between different devices?

In practical production, device coupling typically falls into three scenarios:

1. **Large Buffer Scenario**: When intermediate product buffer areas are sufficiently large, production stages are effectively decoupled and can be modeled as parallel ALs
2. **Small Buffer Scenario**: When buffer areas are very small, adjacent stages must open/close simultaneously and can be modeled as a single AL
3. **Medium Buffer Scenario**: For cases between these extremes, additional constraints can be introduced to the ALF model to describe device coupling relationships

## 3. How is feasibility guaranteed after model simplification?

Our approach doesn't pursue strict inner approximation (ensuring RC solutions are always feasible for OC) for several reasons:

1. In actual grid interaction scenarios, demand-side resources typically allow for a 10-20% error margin, and our method's error (under 10%) falls within this acceptable range
2. Inner approximation methods often exhibit strong conservatism that can significantly limit the benefits of industrial user participation in grid interactions
3. In practical applications, industrial users often have other flexible resources (e.g., energy storage, electric vehicles) that can help achieve desired outcomes

## 4. How is the number of Adjustable Loads (ALs) determined?

The selection of AL numbers involves balancing several factors:

- More ALs allow RCs to better approximate OCs but increase computational costs
- Too many ALs may lead to overfitting, especially with limited training data
- The choice should be based on specific application scenarios
- Optimal AL numbers can be selected from training data through methods like cross-validation

## 5. What are the main innovations of the D3R framework?

The D3R framework's key innovations include:

1. Bypassing the convexity requirement of traditional analytical methods by working with optimal solutions rather than constraints themselves
2. Handling both continuous and integer variables simultaneously
3. Learning reduced constraints directly from data, making it more adaptable to different types of industrial loads
4. Providing a flexible framework where different forms of reduced constraints can be chosen based on needs

## 6. What are the limitations of the current ALF model?

The current ALF model has several limitations:

1. It may not fully capture temporal coupling between production processes
2. The model relaxes integer variables, which might affect accuracy in some scenarios
3. It uses a simplified representation of complex industrial processes
4. The decomposition of total energy consumption into individual devices requires additional consideration

## 7. How is the training convergence ensured?

Training convergence is achieved through:

1. Careful selection of learning rates (e.g., 0.05 instead of 0.9)
2. Implementation of different iteration strategies (single-day or multi-day batch training)
3. Monitoring of key parameters during training until stabilization (relative change < 0.02%)
4. Use of appropriate stopping criteria based on parameter convergence
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

# Frequently Asked Questions (FAQ)

## 1. 为什么选择使用可调负荷群(ALF)模型而不是虚拟电池(VB)模型?

ALF模型相比VB模型有以下优势:

- **建模理念不同**: VB模型关注负荷与基准负荷的偏差,而ALF模型直接表示工业负荷的实际组成,更适合工业应用场景
- **参数灵活性**: VB模型仅使用3个参数,而ALF模型有4I个参数(I为可调负荷数量),能更准确地表示复杂的工业负荷特性
- **物理解释性**: ALF模型更好地反映了工业负荷的物理现实,即多个设备同时运行且各自具有功率和能量约束
- **基准独立性**: 与VB模型不同,ALF模型不需要定义基准负荷,使其对运行条件变化更具鲁棒性

## 2. ALF模型如何处理不同设备之间的耦合关系?

实际生产中设备耦合主要有三种情况:

1. **大缓冲区场景**: 当中间产品缓冲区足够大时,各生产阶段基本解耦,可以用并行的AL建模
2. **小缓冲区场景**: 当缓冲区很小时,相邻阶段必须同时开/关,可以作为单个AL建模
3. **中等缓冲区场景**: 对于介于两者之间的情况,可以在ALF模型中引入额外约束来描述设备间耦合关系

## 3. 模型简化后如何保证可行性?

本方法不追求严格的内部近似(确保RC解对OC始终可行)的原因是:

1. 实际电网互动场景中,需求侧资源通常允许10-20%的误差范围,而本方法误差(低于10%)在可接受范围内
2. 内部近似方法往往过于保守,会显著限制工业用户参与电网互动的收益
3. 实际应用中,工业用户通常有其他灵活资源(如储能、电动汽车等)可以帮助实现期望的效果

## 4. 如何选择可调负荷(AL)的数量?

AL数量的选择需要权衡以下因素:

- 更多的AL允许RC更准确地近似OC,但会增加计算成本
- AL数量过多可能导致过拟合,特别是在训练数据有限的情况下
- 需要根据具体应用场景选择合适的AL数量
- 可以通过交叉验证等方法从训练数据中选择最优AL数量

## 5. D3R框架的主要创新点是什么?

D3R框架的主要创新在于:

1. 绕过了传统分析方法对凸性的要求,通过处理最优解而不是约束本身
2. 能同时处理连续变量和整数变量
3. 直接从数据中学习简化约束,使其更适应不同类型的工业负荷
4. 提供了一个灵活的框架,可以根据需要选择不同的简化约束形式
---

