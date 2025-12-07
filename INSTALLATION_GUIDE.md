# 安装和使用指南 - 完整版

## 问题修复说明

### 问题：LEaP无法识别MOL2原子类型

**错误信息：**
```
Error! For atom (.R<XXX>.A<YYY>) could not find vdW (or other) parameters for type (O.3)
Error! For atom (.R<XXX>.A<YYY>) could not find vdW (or other) parameters for type (C.3)
Error! For atom (.R<XXX>.A<YYY>) could not find vdW (or other) parameters for type (N.ar)
```

**原因：**
- 对接软件输出的MOL2文件使用Sybyl原子类型（`O.3`, `C.3`, `N.ar`等）
- LEaP无法直接识别这些原子类型，需要AMBER原子类型
- 可能缺少氢原子或存在结构问题

**完整解决方案（已实现）：**
脚本现已实现标准的AMBER系统准备流程，包含三个关键步骤：

1. **OpenBabel转换** - MOL2 → PDB格式转换
2. **reduce加氢** - 智能添加氢原子，优化氢键网络
3. **pdb4amber清理** - AMBER兼容性处理

## 依赖安装

### 1. Python环境

```bash
conda create -n amber_md python=3.10
conda activate amber_md
```

### 2. 必需的Python包

```bash
# RDKit（用于SMILES解析和格式转换）
conda install -c conda-forge rdkit

# OpenBabel（用于MOL2到PDB转换，强烈推荐）
conda install -c conda-forge openbabel
```

### 3. AMBER工具

确保已安装AMBER24并加载：
```bash
module load amber24/24
```

需要的AMBER工具：
- `antechamber` - 配体参数化
- `parmchk2` - 力场参数检查
- `tleap` - 系统构建
- `parmed` - 拓扑修改
- `pmemd.cuda` - MD引擎
- **`pdb4amber`** - PDB清理工具（AMBER自带）

### 4. 额外推荐工具

```bash
# reduce - 智能加氢工具（强烈推荐）
conda install -c conda-forge reduce

# 或从源码编译
git clone https://github.com/rlabduke/reduce.git
cd reduce && make && make install
```

## 完整的受体准备流程

脚本中的 `prepare_receptor()` 函数实现了标准的AMBER准备流程：

```python
def prepare_receptor(receptor_mol2, output_pdb):
    """
    步骤1: MOL2 → PDB转换
      - 使用OpenBabel或RDKit
      - 移除Sybyl原子类型

    步骤2: reduce添加氢原子
      - 智能添加氢原子
      - 优化氢键网络
      - 处理His/Asn/Gln侧链方向

    步骤3: pdb4amber清理
      - 删除替代构象
      - 重命名残基为AMBER格式
      - 检测二硫键
      - 添加缺失的OXT末端原子
    """
```

### 详细步骤说明

#### 步骤1: OpenBabel转换

```bash
obabel receptor.mol2 -O receptor_temp.pdb
```

**作用：**
- 转换Sybyl原子类型为PDB标准原子名
- 保留坐标信息

**备选方案：** 如果OpenBabel不可用，使用RDKit

#### 步骤2: reduce加氢

```bash
reduce -BUILD receptor_temp.pdb > receptor_reduced.pdb
```

**作用：**
- **智能加氢**：考虑氢键网络、pH值、静电相互作用
- **侧链优化**：自动翻转His、Asn、Gln侧链以优化氢键
- **质子化状态**：根据环境自动确定组氨酸质子化

**为什么比OpenBabel -h更好？**
| 特性 | OpenBabel -h | reduce -BUILD |
|------|-------------|---------------|
| 氢键优化 | ❌ | ✅ |
| 侧链翻转 | ❌ | ✅ |
| pH考虑 | ❌ | ✅ |
| 生化准确性 | 一般 | 高 |

#### 步骤3: pdb4amber清理

```bash
pdb4amber -i receptor_reduced.pdb -o receptor_prepared.pdb -y --nohyd --noter
```

**参数说明：**
- `-y`: 自动确认所有提示
- `--nohyd`: 保留reduce添加的氢原子
- `--noter`: 不显示TER记录（简化输出）

**作用：**
- **残基重命名**：HIE/HID/HIP → 正确的His形式
- **删除非标准**：移除HETATM、HOH（后续LEaP会重新添加）
- **检测二硫键**：自动识别CYS-CYS二硫键
- **末端处理**：添加缺失的OXT原子
- **替代构象**：选择第一个构象，删除其他

## 脚本更新内容

### v2.0 (完整预处理版本)

两个脚本（`rna_ligand_md.py` 和 `complex_md.py`）均已更新：

#### 新增功能

1. **三步受体准备流程** `prepare_receptor()`:
   ```
   MOL2 → [OpenBabel] → PDB → [reduce] → H-PDB → [pdb4amber] → Clean-PDB
   ```

2. **容错机制**:
   - OpenBabel失败 → 尝试RDKit
   - reduce不可用 → 跳过但警告
   - pdb4amber失败 → 使用未清理版本

3. **详细日志输出**:
   ```
   ============================================================
   受体准备流程: receptor.mol2
   ============================================================

   [1/3] MOL2转PDB...
     ✓ OpenBabel转换成功

   [2/3] 使用reduce添加氢原子...
     ✓ reduce完成

   [3/3] 使用pdb4amber清理PDB...
     ✓ pdb4amber完成

   ✓ 受体准备完成: receptor_prepared.pdb
   ============================================================
   ```

### LEaP输入优化

```tcl
# 旧版（简单转换，可能出错）
RNA = loadmol2 receptor.mol2

# 新版（完整准备，高质量）
RNA = loadpdb receptor_prepared.pdb
```

## 使用方法（无变化）

```bash
# RNA-配体复合物
python rna_ligand_md.py auto-complex-example

# 蛋白质-配体复合物
python complex_md.py auto-complex-example
```

## 工作流程对比

### 旧版流程（v1.0）

```
receptor.mol2
    ↓
[OpenBabel -h]
    ↓
receptor.pdb → LEaP (可能报错)
```

### 新版流程（v2.0）

```
receptor.mol2
    ↓
[OpenBabel] → receptor_temp.pdb
    ↓
[reduce -BUILD] → receptor_reduced.pdb
    ↓
[pdb4amber -y --nohyd] → receptor_prepared.pdb
    ↓
LEaP (高质量输入) ✓
```

## 输出文件结构

```
工作目录/
├── receptor.mol2               # 原始MOL2（复制）
├── receptor_prepared.pdb       # 最终准备好的PDB ← 用于LEaP
├── ligand.pdb
├── ligand.smi
├── pose_0.mol2                 # 参数化的配体
├── pose_0.frcmod
├── complex.tleap               # LEaP输入脚本
├── complex.nosol.pdb           # 无溶剂复合物
├── complex.pdb                 # 溶剂化复合物
├── complex.top                 # AMBER拓扑
├── complex.hmr.top             # HMR拓扑
└── tleap.log                   # LEaP日志
```

## 常见问题 FAQ

### Q1: reduce未安装，脚本会失败吗？

**答：** 不会失败，但**强烈建议安装**。

```bash
# 如果reduce不可用，脚本会：
[2/3] 使用reduce添加氢原子...
  ⚠ reduce不可用或失败，跳过此步骤
  建议安装: conda install -c conda-forge reduce

# 然后继续使用OpenBabel加的氢（质量较低）
```

**安装reduce：**
```bash
conda install -c conda-forge reduce
# 或
conda install -c bioconda reduce
```

### Q2: pdb4amber报错怎么办？

**常见错误：**
```
ERROR: Unknown residue XXX
```

**原因：** 受体包含非标准残基（如修饰核苷酸、非天然氨基酸）

**解决方案：**
```bash
# 方法1: 手动编辑PDB，删除或重命名非标准残基

# 方法2: 使用pdb4amber的容忍模式
# （脚本已实现，会自动跳过清理步骤）
```

### Q3: 如何验证受体准备是否成功？

**检查清单：**

1. **查看日志**：所有步骤都显示 ✓
2. **检查氢原子数**：
   ```bash
   grep "^ATOM.*H" receptor_prepared.pdb | wc -l
   # 应有合理数量的氢原子
   ```

3. **验证残基名**：
   ```bash
   grep "^ATOM" receptor_prepared.pdb | awk '{print $4}' | sort -u
   # RNA: A, U, G, C, RA, RU, RG, RC
   # 蛋白质: ALA, VAL, LEU等标准名称
   ```

4. **LEaP加载测试**：
   ```bash
   tleap -f - <<EOF
   source leaprc.RNA.OL3
   mol = loadpdb receptor_prepared.pdb
   check mol
   quit
   EOF
   ```

### Q4: RNA残基名称不标准

**问题：** 对接软件可能输出 `A` 而不是 `RA`

**pdb4amber自动处理：**
```bash
# 输入: A, U, G, C
# 输出: RA, RU, RG, RC (AMBER标准)
```

**手动处理（如需要）：**
```bash
sed -i 's/ A  / RA /g' receptor.pdb
sed -i 's/ U  / RU /g' receptor.pdb
sed -i 's/ G  / RG /g' receptor.pdb
sed -i 's/ C  / RC /g' receptor.pdb
```

### Q5: 组氨酸质子化状态如何确定？

**reduce自动处理：**
```
HIS → HIE (ε质子化)
HIS → HID (δ质子化)
HIS → HIP (双质子化)
```

**手动指定（如需要）：**
```bash
# 编辑pdb4amber输出，手动重命名
vim receptor_prepared.pdb
# HIS 123 → HIE 123
```

### Q6: 配体的氢原子如何处理？

**配体氢原子由Antechamber处理：**
```bash
antechamber -i ligand.pdb -fi pdb -o ligand.mol2 \
    -c bcc -nc 0 -at gaff2
# 自动添加氢原子并计算AM1-BCC电荷
```

**注意：** 配体不需要reduce/pdb4amber处理

## 性能优化建议

### 1. 批量处理加速

```python
# 并行处理多个pose（需要修改脚本）
from multiprocessing import Pool

def process_pose(pose_pdb):
    # 准备系统并提交任务
    ...

with Pool(4) as p:
    p.map(process_pose, pose_files)
```

### 2. 跳过已准备的受体

```python
# 如果多个pose共用同一受体
if os.path.exists("receptor_prepared.pdb"):
    print("跳过受体准备，使用已有文件")
    receptor_pdb = "receptor_prepared.pdb"
else:
    receptor_pdb = prepare_receptor("receptor.mol2")
```

### 3. reduce参数优化

```bash
# 默认: pH 7.0
reduce -BUILD receptor.pdb

# 自定义pH
reduce -BUILD -pH7.4 receptor.pdb

# 跳过氢键翻转（更快但质量低）
reduce -BUILD -NOFLIP receptor.pdb
```

## 验证安装

### 检查所有工具

```bash
# 1. OpenBabel
obabel --version
# 期望: Open Babel 3.x.x

# 2. reduce
reduce -version
# 期望: reduce version 3.x

# 3. pdb4amber
pdb4amber --version
# 期望: (part of AmberTools)

# 4. RDKit
python -c "from rdkit import Chem; print(Chem.__version__)"
# 期望: 2023.x.x

# 5. AMBER
module load amber24/24
which pmemd.cuda
# 期望: /path/to/amber24/bin/pmemd.cuda
```

### 完整测试脚本

```bash
#!/bin/bash
# test_receptor_prep.sh

# 创建测试MOL2（仅示例）
cat > test_receptor.mol2 << 'EOF'
@<TRIPOS>MOLECULE
RNA_TEST
 5 4 1 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 P 0.0 0.0 0.0 P.3 1 A 0.0
2 O5' 1.5 0.0 0.0 O.3 1 A 0.0
3 C5' 2.0 1.0 0.0 C.3 1 A 0.0
4 H5'1 2.5 1.0 0.5 H 1 A 0.0
5 H5'2 2.5 1.0 -0.5 H 1 A 0.0
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
3 3 4 1
4 3 5 1
EOF

# 测试准备流程
python -c "
from rna_ligand_md import prepare_receptor
prepare_receptor('test_receptor.mol2', 'test_prepared.pdb')
print('✓ 测试成功！')
"
```

## 工具对比

| 工具 | 功能 | 必需性 | 替代方案 |
|------|------|--------|----------|
| **OpenBabel** | MOL2→PDB转换 | 必需 | RDKit（备选） |
| **reduce** | 智能加氢 | 强烈推荐 | OpenBabel -h（质量低） |
| **pdb4amber** | AMBER兼容性 | 推荐 | 手动编辑PDB |
| **antechamber** | 配体参数化 | 必需 | 无 |
| **parmchk2** | 力场参数 | 必需 | 无 |

## 更新日志

**v2.0 (2025-12-07) - 完整预处理版本**
- ✅ 添加reduce智能加氢
- ✅ 添加pdb4amber清理步骤
- ✅ 实现三步受体准备流程
- ✅ 容错机制（工具缺失时降级）
- ✅ 详细日志输出
- ✅ 更新两个脚本（RNA和蛋白质版本）

**v1.1 (2025-12-07)**
- 添加基本MOL2到PDB转换
- 支持OpenBabel和RDKit

**v1.0 (2025-12-07)**
- 初始版本

## 推荐工作流程

### 标准流程（全自动）

```bash
# 1. 准备输入文件
mkdir my_complex
cd my_complex
cp /path/to/receptor.mol2 .
cp /path/to/ligand.smi .
cp /path/to/pose_*.pdb .

# 2. 运行脚本（自动完成所有步骤）
python /path/to/rna_ligand_md.py .

# 3. 检查日志
tail -f md_pose_0/tleap.log

# 4. 监控任务
squeue -u $USER
```

### 高级流程（手动控制）

```bash
# 1. 单独准备受体（可复用）
python -c "
from rna_ligand_md import prepare_receptor
prepare_receptor('receptor.mol2', 'receptor_prepared.pdb')
"

# 2. 检查质量
pymol receptor_prepared.pdb
# 检查氢原子、结构完整性

# 3. 如有问题，手动调整
reduce -BUILD -pH7.4 -NOFLIP receptor_temp.pdb > receptor_custom.pdb
pdb4amber -i receptor_custom.pdb -o receptor_final.pdb

# 4. 使用自定义受体（修改脚本或手动拷贝）
cp receptor_final.pdb md_pose_0/receptor_prepared.pdb
```

## 技术支持

遇到问题时的检查顺序：

1. **查看脚本输出** - 每个步骤都有明确提示
2. **检查tleap.log** - LEaP的详细错误信息
3. **验证工具安装** - `which reduce pdb4amber obabel`
4. **检查输入文件** - PDB/MOL2格式是否正确
5. **参考AMBER手册** - AmberTools用户指南
6. **联系开发者** - 提供完整错误日志

---

## 快速参考卡片

```
受体准备三步曲：
1️⃣  MOL2 → PDB  [obabel]
2️⃣  加氢优化      [reduce]
3️⃣  AMBER清理    [pdb4amber]

必装工具：
✅ OpenBabel
✅ reduce
✅ AmberTools (含pdb4amber)
✅ RDKit

一键运行：
python rna_ligand_md.py <dir>
```
