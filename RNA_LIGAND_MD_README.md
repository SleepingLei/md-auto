# RNA-Ligand MD - RNA-配体复合物分子动力学模拟

## 功能概述

`rna_ligand_md.py` 专门用于 **RNA-配体复合物** 的MD模拟，完全参考 `md.py` 的设计风格和参数设置。

### 主要特点

- ✅ 支持所有 RNA 力场（OL3, Shaw, YIL, ROC, LJbb）
- ✅ 使用与 `md.py` 相同的离子浓度（Mg²⁺ + K⁺）
- ✅ 高温采样策略（默认 350K）用于 RNA 构象探索
- ✅ 相同的盒子大小（15Å）和时间步长（2fs）
- ✅ 自动配体参数化（AM1-BCC 电荷）
- ✅ 批量处理多个对接 pose

## 与 md.py 的一致性

| 参数 | md.py (RNA) | rna_ligand_md.py (RNA-配体) |
|------|-------------|---------------------------|
| RNA力场 | OL3/Shaw/YIL/ROC/LJbb | ✅ 相同 |
| 离子 | Mg²⁺ + K⁺ | ✅ 相同 (0.15M + 0.1M) |
| 盒子边界 | 15 Å | ✅ 相同 |
| 时间步长 | 2 fs | ✅ 相同 |
| 温度 | 350 K | ✅ 相同 |
| 限制原子 | @P | ✅ 相同 |
| MD协议 | EM→Heat→RS→HighT | ✅ 相同 |

## 输入文件要求

```
input_dir/
├── ligand.smi          # 配体SMILES（单行）
├── receptor.mol2       # RNA受体结构（MOL2格式）
├── pose_0.pdb          # 对接pose 1
├── pose_1.pdb          # 对接pose 2
└── pose_N.pdb          # 对接pose N
```

## 依赖环境

### Python依赖
```bash
pip install rdkit
```

### AMBER工具
- `antechamber` - 配体参数化
- `parmchk2` - 力场参数检查
- `tleap` - 系统构建
- `parmed` - 拓扑修改
- `pmemd.cuda` - MD引擎

## 使用方法

### 基本用法

```bash
python rna_ligand_md.py auto-complex-example
```

### 高级用法

```bash
# 使用Shaw力场，更高温度，更长时间
python rna_ligand_md.py auto-complex-example \
    --rna_ff shaw \
    --temp 370 \
    --ns 50 \
    --box_d 18

# 调整离子浓度
python rna_ligand_md.py auto-complex-example \
    --mg_m 0.20 \
    --k_m 0.15
```

### 参数说明

| 参数 | 默认值 | 说明 | 来自md.py |
|------|--------|------|-----------|
| `--mg_m` | 0.15 | Mg²⁺离子浓度 (M) | ✅ |
| `--k_m` | 0.1 | K⁺离子浓度 (M) | ✅ |
| `--box_d` | 15 | 溶剂盒子边界距离 (Å) | ✅ |
| `--temp` | 350 | 模拟温度 (K) | ✅ |
| `--ns` | 20 | 生产MD时长 (ns) | ✅ |
| `--rna_ff` | ol3 | RNA力场 | ✅ |
| `--ligand_ff` | gaff2 | 配体力场 | - |
| `--wat_ff` | tip3p | 水模型 | ✅ |

## 支持的力场

### RNA 力场（与 md.py 完全一致）

| 力场 | 参数值 | LEaP命令 |
|------|--------|---------|
| AMBER OL3 | `ol3` | `source leaprc.RNA.OL3` |
| DESRES Shaw | `shaw` | `source leaprc.RNA.Shaw` |
| YIL | `yil` | `source leaprc.RNA.YIL` |
| ROC | `roc` | `source leaprc.RNA.ROC` |
| LJbb | `ljbb` | `source leaprc.RNA.LJbb` |

### 配体力场

| 力场 | 参数值 |
|------|--------|
| GAFF2 (推荐) | `gaff2` |
| GAFF | `gaff` |

### 水模型

| 模型 | 参数值 |
|------|--------|
| TIP3P | `tip3p` |
| OPC | `opc` |
| OPC3 | `opc3` |

## 工作流程

### 1. 配体参数化

```bash
# 自动从SMILES计算净电荷
net_charge = get_charge_from_smiles(smiles)

# Antechamber参数化
antechamber -i pose_0.pdb -fi pdb -o pose_0.mol2 -fo mol2 \
    -c bcc -nc {net_charge} -rn LIG -at gaff2

# 生成缺失力场参数
parmchk2 -i pose_0.mol2 -f mol2 -o pose_0.frcmod
```

### 2. 系统构建（LEaP）

```tcl
# 加载力场
source leaprc.RNA.OL3
source leaprc.gaff2
source leaprc.water.tip3p

# 加载配体参数
loadamberparams pose_0.frcmod
LIG = loadmol2 pose_0.mol2

# 加载RNA受体
RNA = loadmol2 receptor.mol2

# 组装复合物
complex = combine {RNA LIG}

# 添加溶剂和离子（使用八面体盒子）
solvateoct complex TIP3PBOX 15

# 添加Mg2+和K+离子（根据水分子数计算）
addionsrand complex MG <mg_count>
addionsrand complex K+ <k_count>
addionsrand complex K+ 0
addionsrand complex Cl- 0
```

### 3. MD模拟协议（与 md.py 完全一致）

| 步骤 | 输入文件 | 时长 | 温度 | 系综 | 限制 |
|------|----------|------|------|------|------|
| 1. 能量最小化 | em.in | 200步 | - | - | 0.1 kcal/mol/Å² (@P) |
| 2. 加热 | heat.in | 400 ps | 0→300K | NVT | 0.1 kcal/mol/Å² (@P) |
| 3. 限制性平衡 | rs.in | 2 ns | 300K | NPT | 0.1 kcal/mol/Å² (@P) |
| 4. 高温限制平衡 | 350K_rs.in | 1 ns | 300→350K | NVT | 0.1 kcal/mol/Å² (@P) |
| 5. 高温生产MD | 350K.in | 20 ns | 350K | NVT | 无 |

**关键设置（参考 md.py）**:
- 限制原子：磷原子 (`restraintmask='@P'`)
- Langevin动力学：`ntt=3, gamma_ln=1.0`
- 截断距离：10 Å
- SHAKE约束：`ntf=2, ntc=2`
- 时间步长：2 fs

### 4. 输出文件结构

```
md_pose_0/
├── receptor.mol2              # RNA受体副本
├── ligand.pdb                 # 配体副本
├── ligand.smi                 # SMILES副本
├── pose_0.mol2                # 参数化的配体
├── pose_0.frcmod              # 配体力场参数
├── complex.tleap              # LEaP输入
├── complex.top                # AMBER拓扑（完整）
├── complex.hmr.top            # HMR拓扑
├── complex.nosol.pdb          # 无溶剂复合物
├── complex.pdb                # 溶剂化复合物
├── em.in                      # 能量最小化输入
├── heat.in                    # 加热输入
├── rs.in                      # 限制性平衡输入
├── 350K_rs.in                 # 高温限制平衡输入
├── 350K.in                    # 生产MD输入
├── md_0.slurm                 # SLURM脚本
└── [输出轨迹和日志文件]
    ├── system.em.0.rst
    ├── system.heat.0.nc
    ├── system.rs.0.nc
    ├── system.350K_rs.0.nc
    └── system.350K.0.nc       # 生产轨迹
```

## Python API使用

### 单个pose处理

```python
from rna_ligand_md import construct_rna_ligand_top, run_rna_ligand_md

# 构建拓扑
top = construct_rna_ligand_top(
    receptor_mol2="receptor.mol2",
    ligand_pdb="pose_0.pdb",
    ligand_smi="ligand.smi",
    work_pwd="md_pose_0",
    mg_m=0.15,
    k_m=0.1,
    box_d=15,
    rna_ff="ol3",
    ligand_ff="gaff2",
    wat_ff="tip3p"
)

# 运行MD
run_rna_ligand_md(
    work_pwd="md_pose_0",
    top=top,
    replica=0,
    temp=350,
    ns=20
)
```

### 批量处理

```python
from rna_ligand_md import batch_run_rna_ligand_poses

batch_run_rna_ligand_poses(
    input_dir="auto-complex-example",
    mg_m=0.15,
    k_m=0.1,
    box_d=15,
    temp=350,
    ns=20,
    rna_ff="ol3"
)
```

## 与 md.py 和 complex_md.py 的对比

| 特性 | md.py | rna_ligand_md.py | complex_md.py |
|------|-------|------------------|---------------|
| 受体类型 | RNA | RNA | 蛋白质 |
| 力场 | RNA专用 | RNA + GAFF | Protein + GAFF |
| 离子 | Mg²⁺ + K⁺ | Mg²⁺ + K⁺ | Na⁺ |
| 盒子类型 | 八面体 | 八面体 | 八面体 |
| 盒子大小 | 15 Å | 15 Å | 10 Å |
| 温度 | 350 K | 350 K | 300 K |
| 时间步长 | 2 fs | 2 fs | 4 fs (HMR) |
| 距离限制 | 二级结构 | 无 | 无 |
| 用途 | RNA ensemble | RNA-配体对接 | 蛋白-配体对接 |

## 高温采样的原理

### 为什么使用 350K？（参考 md.py）

1. **增强采样**: 高温增加构象转换频率，快速探索RNA折叠空间
2. **克服势垒**: 帮助系统跨越能量势垒，发现多个亚稳态
3. **配体解离**: 允许配体从非最优结合模式中逃逸
4. **符合RNA研究习惯**: 许多RNA MD研究使用高温采样

### 降温策略（如需要）

如果需要在常规温度下精修，可以使用两步策略：

```bash
# 第一步：高温采样（探索构象空间）
python rna_ligand_md.py data --temp 350 --ns 20

# 第二步：从高温轨迹提取结构，降温精修（需手动实现）
# 使用 cpptraj 聚类 → 提取代表性结构 → 300K再模拟
```

## 实际应用案例

### 案例1：适配体-小分子复合物

```bash
# 研究RNA适配体与小分子配体的结合
python rna_ligand_md.py aptamer_screening \
    --rna_ff ol3 \
    --temp 350 \
    --ns 50 \
    --mg_m 0.1 \
    --k_m 0.15
```

### 案例2：核糖开关构象变化

```bash
# 研究配体诱导的核糖开关折叠
python rna_ligand_md.py riboswitch \
    --rna_ff shaw \
    --temp 370 \
    --ns 100 \
    --box_d 18
```

### 案例3：虚拟筛选验证

```bash
# 对接产生100个pose，验证结合稳定性
python rna_ligand_md.py virtual_screening \
    --temp 350 \
    --ns 10
# 输出：100个轨迹，分析RMSD和结合自由能
```

## 轨迹分析

### 使用 cpptraj

```bash
# 计算RMSD
cpptraj -p complex.hmr.top <<EOF
trajin system.350K.0.nc
rms RNA :1-20 # RNA残基
rms LIG :LIG # 配体
run
EOF

# 氢键分析
cpptraj -p complex.hmr.top <<EOF
trajin system.350K.0.nc
hbond RNA-LIG :1-20 :LIG avgout hbond.dat
run
EOF
```

### 使用 MMPBSA.py

```bash
# 计算结合自由能
MMPBSA.py -O -i mmpbsa.in \
    -o FINAL_RESULTS.dat \
    -sp complex.hmr.top \
    -cp complex.hmr.top \
    -rp receptor.top \
    -lp ligand.top \
    -y system.350K.0.nc
```

## 注意事项

1. **RNA受体格式**: 必须是MOL2格式，包含正确的原子类型
2. **配体位置**: pose PDB中的配体应已对接到RNA结合口袋
3. **盒子大小**: RNA通常需要更大的盒子（15Å），避免周期性影响
4. **离子类型**: RNA偏好Mg²⁺和K⁺，不同于蛋白质的Na⁺
5. **高温注意**: 350K会增加构象灵活性，不适合精确结合能计算
6. **限制策略**: 仅限制磷原子，允许碱基和配体自由运动

## 故障排查

### LEaP报错: "Unknown residue"
- 检查RNA受体MOL2中的残基名称（应为A/U/G/C）
- 确认使用了正确的RNA力场

### 高温模拟不稳定
- 检查初始结构质量
- 适当增加平衡时间
- 考虑降低温度到340K或330K

### 离子参数化失败
- 确认安装了正确的AMBER力场文件
- Shaw力场与TIP4P-D配合使用

## SLURM配置

默认SLURM参数（在脚本中可修改）：

```bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=480:00:00    # 20天
#SBATCH --partition=4090
```

## 示例数据准备

```bash
# 创建测试目录
mkdir rna_aptamer_example
cd rna_aptamer_example

# 准备配体SMILES（例如：茶碱）
echo "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" > ligand.smi

# 准备RNA受体（使用Chimera/PyMOL转换为MOL2）
# receptor.mol2 - 茶碱适配体结构

# 对接产生pose（使用AutoDock/rDock等）
# pose_0.pdb, pose_1.pdb, ...
```

## 参考文献

如果使用此脚本，请引用相关力场和方法：
- **OL3力场**: Zgarbová et al., J. Chem. Theory Comput. 2011
- **Shaw力场**: Tan et al., J. Phys. Chem. B 2018
- **GAFF**: Wang et al., J. Comput. Chem. 2004
- **AM1-BCC电荷**: Jakalian et al., J. Comput. Chem. 2002

## 开发路线图

未来计划添加：
- [ ] 自动聚类分析
- [ ] 结合自由能计算流程
- [ ] 支持二级结构距离限制（参考md.py）
- [ ] 多温度副本交换（REMD）
- [ ] 配体RMSD和结合模式分析脚本

## 联系与支持

基于 `md.py` 开发，保持最大兼容性。如有问题请参考AMBER手册或联系开发者。
