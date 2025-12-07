# RNA残基名规范化 - 修复说明

## 问题描述

对接软件输出的MOL2文件中，RNA残基名可能包含非标准命名：

```
错误的残基名：C25, G25, A25, U25
LEaP期望的残基名：C, G, A, U (或 RC, RG, RA, RU)
```

### LEaP错误信息

```
FATAL: Atom .R<C25 1>.A<P 1> does not have a type.
FATAL: Atom .R<G25 3>.A<C4 15> does not have a type.
```

## 解决方案（已实现）

### 新增功能：四步受体准备流程

```
[1/4] MOL2 → PDB 转换
[2/4] reduce 智能加氢
[3/4] pdb4amber 清理
[4/4] RNA残基名规范化 ← 新增
```

### 第4步：残基名规范化

**功能：** `normalize_rna_residue_names()`

**处理的命名模式：**

| 输入残基名 | 输出残基名 | 说明 |
|-----------|-----------|------|
| `C25` | `C` | 提取首字母 |
| `G25` | `G` | 提取首字母 |
| `A_25` | `A` | 下划线分隔也支持 |
| `CYT` | `C` | 全名映射 |
| `GUA` | `G` | 全名映射 |
| `ADE` | `A` | 全名映射 |
| `URA` | `U` | 全名映射 |
| `THY` | `T` | DNA胸腺嘧啶 |

## 实现细节

### 检测逻辑

```python
# 1. 读取PDB，收集所有残基名
residues_found = {'C25', 'G25', 'A', 'U', 'LIG'}

# 2. 识别非标准RNA残基
standard_rna = {'A', 'U', 'G', 'C', 'RA', 'RU', 'RG', 'RC', ...}
non_standard = [r for r in residues_found
                if r not in standard_rna and r[0] in 'AUGCT']
# 结果: ['C25', 'G25']
```

### 修复策略

**策略1：直接映射**
```python
RNA_RESIDUE_MAP = {
    'CYT': 'C', 'GUA': 'G', 'ADE': 'A',
    'URA': 'U', 'THY': 'T'
}
```

**策略2：正则提取首字母**
```python
# 匹配: C25, C_25, C-25 等
pattern = r'^([AUGCT])\d+$|^([AUGCT])_'
# C25 → 提取 'C'
```

### PDB格式处理

```python
# PDB ATOM行格式（列位置）
# ATOM      1  P    C25 A   1      ...
#                  ^^^
#                17-20列：残基名（右对齐）

# 修复前
line = "ATOM      1  P    C25 A   1  ..."

# 修复后
line = "ATOM      1  P      C A   1  ..."
#                      ^^^
#                  右对齐3个字符
```

## 输出示例

### 正常情况（无需修复）

```
[4/4] 规范化RNA残基名...
  检测到的残基: ['A', 'C', 'G', 'U']
  ✓ 所有残基名已是标准格式
```

### 检测到非标准残基

```
[4/4] 规范化RNA残基名...
  检测到的残基: ['A25', 'C25', 'G25', 'U']
  ⚠ 检测到非标准RNA残基: ['A25', 'C25', 'G25']
  → 将规范化为标准AMBER格式
  ✓ 修复了 180 个原子的残基名
```

## 测试验证

### 1. 快速测试

```python
# 测试脚本
from rna_ligand_md import normalize_rna_residue_names

# 创建测试PDB
cat > test.pdb << 'EOF'
ATOM      1  P    C25 A   1       1.000   2.000   3.000  1.00  0.00           P
ATOM      2  O5'  C25 A   1       2.000   3.000   4.000  1.00  0.00           O
ATOM      3  P    G25 A   2       5.000   6.000   7.000  1.00  0.00           P
EOF

# 运行规范化
normalize_rna_residue_names('test.pdb')

# 检查结果
grep "^ATOM" test.pdb
# 期望输出：
# ATOM      1  P      C A   1  ...
# ATOM      2  O5'    C A   1  ...
# ATOM      3  P      G A   2  ...
```

### 2. 完整流程测试

```bash
# 准备测试数据
mkdir test_residue_fix
cd test_residue_fix
cp /path/to/receptor.mol2 .  # 含C25, G25的MOL2
cp /path/to/ligand.smi .
cp /path/to/pose_0.pdb .

# 运行脚本（会自动修复）
python /path/to/rna_ligand_md.py .

# 查看日志
tail -100 md_pose_0/tleap.log
# 应该不再有 "does not have a type" 错误
```

### 3. 手动检查修复后的PDB

```bash
# 查看受体准备后的残基名
grep "^ATOM" md_pose_0/receptor_prepared.pdb | awk '{print $4}' | sort -u
# 期望输出: A, C, G, U (而不是 A25, C25, G25, U25)
```

## 兼容性

### 支持的残基命名变体

✅ 数字后缀：`C25`, `G99`, `A1`, `U42`
✅ 下划线分隔：`C_25`, `G_99`
✅ 全名：`CYT`, `GUA`, `ADE`, `URA`, `THY`
✅ AMBER标准：`C`, `G`, `A`, `U`, `T`
✅ 带R前缀：`RC`, `RG`, `RA`, `RU` (已是标准，不修改)
✅ DNA：`DA`, `DT`, `DG`, `DC` (已是标准，不修改)

### 不会误修复的情况

❌ 蛋白质残基：`ALA`, `VAL`, `LEU` 等（首字母不是AUGCT）
❌ 配体残基：`LIG`, `XXX` 等
❌ 水分子：`WAT`, `HOH`
❌ 离子：`Na+`, `Cl-`, `MG`

## 工作流程集成

### 脚本自动调用

```python
# rna_ligand_md.py 和 complex_md.py 中
def prepare_receptor(receptor_mol2, output_pdb):
    # 步骤1: MOL2 → PDB
    # 步骤2: reduce 加氢
    # 步骤3: pdb4amber 清理
    # 步骤4: 残基名规范化 ← 自动调用
    normalize_rna_residue_names(pdb4amber_pdb)
```

### 对用户透明

用户**无需任何额外操作**，脚本会自动：
1. 检测非标准残基名
2. 规范化为AMBER格式
3. 显示修复统计

## 调试信息

### 查看详细输出

脚本会自动打印诊断信息：

```
============================================================
受体准备流程: receptor.mol2
============================================================

[1/4] MOL2转PDB...
  ✓ OpenBabel转换成功

[2/4] 使用reduce添加氢原子...
  ✓ reduce完成

[3/4] 使用pdb4amber清理PDB...
  ✓ pdb4amber完成

[4/4] 规范化RNA残基名...
  检测到的残基: ['A25', 'C25', 'G25', 'U']
  ⚠ 检测到非标准RNA残基: ['A25', 'C25', 'G25']
  → 将规范化为标准AMBER格式
  ✓ 修复了 180 个原子的残基名

✓ 受体准备完成: receptor_prepared.pdb
============================================================
```

## 常见问题

### Q: 如果是蛋白质-配体复合物会误修复吗？

**A:** 不会。函数只修复首字母为 `A/U/G/C/T` 的残基。

```python
# 蛋白质残基（不会被修改）
'ALA', 'VAL', 'LEU', 'ILE', ...  # 首字母A/V/L/I，但不匹配模式

# RNA残基（会被修复）
'A25', 'C25', 'G25', 'U25'  # 匹配 ^[AUGCT]\d+$
```

### Q: DNA会被正确处理吗？

**A:** 会。DNA使用 `DA`, `DT`, `DG`, `DC`，这些是AMBER标准格式，不会被修改。

### Q: 如果有DNA+RNA混合呢？

**A:** 也能正确处理：
- `DA`, `DT`, `DG`, `DC` → 保持不变（标准格式）
- `A25`, `C25` 等 → 修复为 `A`, `C`

### Q: 修复后残基编号会变吗？

**A:** **不会**。只修改残基类型名（17-20列），残基编号（23-26列）保持不变。

```
修复前：ATOM   1  P    C25 A   5  ...
修复后：ATOM   1  P      C A   5  ...
                    ^^^       ^
                 残基类型    残基编号(不变)
```

## 版本更新

**v2.1 (2025-12-07) - 残基名规范化**
- ✅ 添加 `normalize_rna_residue_names()` 函数
- ✅ 支持多种非标准命名模式
- ✅ 智能检测与诊断
- ✅ 集成到4步受体准备流程
- ✅ 更新 rna_ligand_md.py 和 complex_md.py

## 技术细节

### PDB残基名格式规范

根据PDB格式规范（ATOM记录）：

```
列号    内容          格式      说明
17-20   残基名        str(3)    右对齐，最多3个字符
23-26   残基编号      int(4)    残基序号
```

### 正则表达式模式

```python
# 完整模式
pattern = r'^([AUGCT])\d+$|^([AUGCT])_'

# 分解：
^([AUGCT])\d+$     # 匹配 C25, G99
    |              # 或
^([AUGCT])_        # 匹配 C_25, G_99
```

### 测试用例

```python
test_cases = [
    ('C25', 'C'),    # 数字后缀
    ('G_25', 'G'),   # 下划线
    ('CYT', 'C'),    # 全名
    ('RA', 'RA'),    # 标准AMBER（不修改）
    ('ALA', 'ALA'),  # 蛋白质（不修改）
]
```

## 参考文献

- PDB格式规范: https://www.wwpdb.org/documentation/file-format
- AMBER命名约定: AmberTools用户手册
- reduce文档: http://kinemage.biochem.duke.edu/software/reduce/
