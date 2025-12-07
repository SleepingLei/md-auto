# Complex MD - 蛋白质-配体复合物分子动力学模拟

## 功能概述

`complex_md.py` 用于自动化处理蛋白质-配体复合物的MD模拟，参考 `md.py` 的设计风格。

### 主要特点

- 自动从SMILES计算配体电荷
- 使用Antechamber进行配体参数化（AM1-BCC电荷）
- 支持批量处理多个对接pose
- 本地准备系统，SLURM提交集群任务
- 氢质量重分配（HMR）支持4fs时间步长

## 输入文件要求

在输入文件夹中需要包含：

```
input_dir/
├── ligand.smi          # 配体SMILES（单行）
├── receptor.mol2       # 受体结构（MOL2格式）
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
python complex_md.py auto-complex-example
```

### 高级用法

```bash
python complex_md.py auto-complex-example \
    --temp 310 \
    --ns 20 \
    --na_m 0.15 \
    --box_d 12 \
    --protein_ff ff19SB \
    --ligand_ff gaff2 \
    --wat_ff opc
```

### 参数说明

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--na_m` | 0.15 | Na+离子浓度 (M) |
| `--box_d` | 10 | 溶剂盒子边界距离 (Å) |
| `--temp` | 300 | 模拟温度 (K) |
| `--ns` | 10 | 生产MD时长 (ns) |
| `--protein_ff` | ff14SB | 蛋白质力场 (ff14SB/ff19SB) |
| `--ligand_ff` | gaff2 | 配体力场 (gaff/gaff2) |
| `--wat_ff` | tip3p | 水模型 (tip3p/opc/tip4pew) |

## 工作流程

### 1. 配体参数化

```python
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
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p

# 加载配体参数
loadamberparams pose_0.frcmod
LIG = loadmol2 pose_0.mol2

# 加载受体
REC = loadmol2 receptor.mol2

# 组装复合物
complex = combine {REC LIG}

# 添加溶剂和离子
solvatebox complex TIP3PBOX 10
addionsrand complex Na+ <calculated>
addionsrand complex Na+ 0
addionsrand complex Cl- 0
```

### 3. MD模拟协议

| 步骤 | 输入文件 | 时长 | 系综 | 限制 | 说明 |
|------|----------|------|------|------|------|
| 1. 能量最小化 | 01_em.in | 5000步 | - | 10 kcal/mol/Å² | 限制蛋白质和配体 |
| 2. 加热 | 02_heat.in | 50 ps | NVT | 5 kcal/mol/Å² | 0K → 300K |
| 3. 密度平衡 | 03_equil_npt.in | 500 ps | NPT | 2 kcal/mol/Å² | Berendsen控压 |
| 4. 无限制平衡 | 04_equil_free.in | 2 ns | NPT | 无 | 完全弛豫 |
| 5. 生产MD | 05_prod.in | 10 ns | NVT | 无 | 数据采集 |

### 4. 输出文件结构

```
md_pose_0/
├── receptor.mol2              # 受体副本
├── ligand.pdb                 # 配体副本
├── ligand.smi                 # SMILES副本
├── pose_0.mol2                # 参数化的配体
├── pose_0.frcmod              # 配体力场参数
├── complex.tleap              # LEaP输入
├── complex.top                # AMBER拓扑（完整）
├── complex.hmr.top            # HMR拓扑
├── complex.inpcrd             # 初始坐标
├── complex.pdb                # 溶剂化系统PDB
├── 01_em.in ~ 05_prod.in      # MD输入文件
├── run_md_0.slurm             # SLURM脚本
└── [输出轨迹和日志文件]
```

## Python API使用

### 单个pose处理

```python
from complex_md import construct_complex_top, run_complex_md

# 构建拓扑
top = construct_complex_top(
    receptor_mol2="receptor.mol2",
    ligand_pdb="pose_0.pdb",
    ligand_smi="ligand.smi",
    work_pwd="md_pose_0",
    na_m=0.15,
    box_d=10,
    protein_ff="ff14SB",
    ligand_ff="gaff2",
    wat_ff="tip3p"
)

# 运行MD
run_complex_md(
    work_pwd="md_pose_0",
    top=top,
    replica=0,
    temp=300,
    ns=10
)
```

### 批量处理

```python
from complex_md import batch_run_poses

batch_run_poses(
    input_dir="auto-complex-example",
    na_m=0.15,
    box_d=10,
    temp=300,
    ns=10
)
```

## 与 md.py 的主要区别

| 特性 | md.py (RNA) | complex_md.py (蛋白质-配体) |
|------|-------------|----------------------------|
| 系统类型 | RNA单链 | 蛋白质-配体复合物 |
| 输入格式 | PDB + 序列 + 二级结构 | MOL2 + PDB + SMILES |
| 力场 | RNA力场(OL3/Shaw等) | 蛋白质力场 + GAFF |
| 配体处理 | 无 | Antechamber参数化 |
| 距离限制 | 基于二级结构 | 无 |
| 时间步长 | 2 fs | 4 fs (HMR) |
| 温度 | 高温采样(350K) | 常规MD(300K) |

## SLURM配置

默认SLURM参数（在脚本中可修改）：

```bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=480:00:00
#SBATCH --partition=4090
```

## 注意事项

1. **电荷计算**: 确保SMILES格式正确，脚本会自动计算形式电荷
2. **受体格式**: 必须是MOL2格式（包含原子类型）
3. **配体构象**: pose PDB中的配体应已对接到受体结合位点
4. **盒子大小**: box_d=10Å 适用于大多数情况，大蛋白质可增加
5. **离子浓度**: 默认0.15M NaCl，根据实验条件调整

## 故障排查

### Antechamber失败
- 检查配体PDB是否包含氢原子
- 检查SMILES是否正确
- 查看 `ANTECHAMBER*` 临时文件排查错误

### LEaP报错
- 检查受体MOL2原子类型是否正确
- 确认力场文件已安装（`$AMBERHOME/dat/leap/cmd/`）
- 查看 `tleap.log` 详细信息

### 任务未提交
- 检查SLURM是否可用：`squeue`
- 验证模块加载路径：`module avail amber`

## 示例数据准备

```bash
# 创建测试目录
mkdir auto-complex-example
cd auto-complex-example

# 准备配体SMILES
echo "CC(C)Cc1ccc(cc1)C(C)C(=O)O" > ligand.smi

# 从对接软件导出pose（例如AutoDock Vina）
# pose_0.pdb, pose_1.pdb, ...

# 准备受体（使用Chimera/PyMOL转换为MOL2）
# receptor.mol2
```

## 联系与支持

如有问题请参考AMBER手册或联系开发者。
