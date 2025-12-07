import os
import subprocess
import shutil
import glob
import re
from rdkit import Chem


def normalize_rna_residue_names(pdb_file):
    """
    规范化RNA残基名称

    处理对接软件常见的非标准残基命名：
    - C25, C_25, CYT → C
    - G25, G_25, GUA → G
    - A25, A_25, ADE → A
    - U25, U_25, URA → U
    - T25, T_25, THY → T

    保留残基编号，仅修改残基类型
    """
    # RNA残基映射表
    RNA_RESIDUE_MAP = {
        'CYT': 'C', 'GUA': 'G', 'ADE': 'A', 'URA': 'U', 'THY': 'T',
    }

    # 读取PDB
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    # 诊断：检测所有残基名
    residues_found = set()
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            resname = line[17:20].strip()
            if resname:
                residues_found.add(resname)

    print(f"  检测到的残基: {sorted(residues_found)}")

    # 检测非标准RNA残基
    standard_rna = {'A', 'U', 'G', 'C', 'T', 'RA', 'RU', 'RG', 'RC', 'RT', 'DA', 'DT', 'DG', 'DC'}
    non_standard = [r for r in residues_found if r not in standard_rna and r[0] in 'AUGCT']

    if non_standard:
        print(f"  ⚠ 检测到非标准RNA残基: {non_standard}")
        print(f"  → 将规范化为标准AMBER格式")
    else:
        print(f"  ✓ 所有残基名已是标准格式")
        return pdb_file

    # 修复残基名
    fixed_lines = []
    fixes_count = 0

    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            resname = line[17:20].strip()
            new_resname = resname

            # 策略1: 直接映射（CYT → C）
            if resname in RNA_RESIDUE_MAP:
                new_resname = RNA_RESIDUE_MAP[resname]
                fixes_count += 1

            # 策略2: 提取首字母（C25 → C, G_25 → G）
            elif resname and resname[0] in 'AUGCT':
                # 匹配模式：字母 + 数字/下划线
                match = re.match(r'^([AUGCT])\d+$|^([AUGCT])_', resname)
                if match:
                    new_resname = match.group(1) or match.group(2)
                    fixes_count += 1

            # 重写PDB行（残基名在17-20列，右对齐）
            if new_resname != resname:
                line = line[:17] + f'{new_resname:>3}' + line[20:]

        fixed_lines.append(line)

    # 写回PDB
    with open(pdb_file, 'w') as f:
        f.writelines(fixed_lines)

    if fixes_count > 0:
        print(f"  ✓ 修复了 {fixes_count} 个原子的残基名")

    return pdb_file


def prepare_receptor(receptor_mol2, output_pdb="receptor_prepared.pdb"):
    """
    完整的受体准备流程（适用于RNA和蛋白质）

    步骤：
    1. MOL2 → PDB转换（OpenBabel）
    2. reduce添加氢原子（优化氢键网络）
    3. pdb4amber清理PDB（AMBER兼容性）
    4. RNA残基名规范化（修复C25→C等问题）

    参数：
        receptor_mol2: 输入MOL2文件
        output_pdb: 输出PDB文件
    """
    print(f"\n{'='*60}")
    print(f"受体准备流程: {receptor_mol2}")
    print(f"{'='*60}")

    # 步骤1: MOL2 → PDB（不加氢，下一步用reduce加）
    temp_pdb = "receptor_temp.pdb"
    print(f"\n[1/4] MOL2转PDB...")
    cmd = f"obabel {receptor_mol2} -O {temp_pdb}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0 or not os.path.exists(temp_pdb):
        # 备选: RDKit
        print("  OpenBabel失败，尝试RDKit...")
        try:
            from rdkit import Chem
            mol = Chem.MolFromMol2File(receptor_mol2, removeHs=False)
            if mol is not None:
                Chem.MolToPDBFile(mol, temp_pdb)
                print("  ✓ RDKit转换成功")
            else:
                raise RuntimeError("RDKit无法读取MOL2文件")
        except Exception as e:
            raise RuntimeError(f"MOL2转换失败: {e}\n请安装: conda install -c conda-forge openbabel")
    else:
        print(f"  ✓ OpenBabel转换成功")

    # 步骤2: reduce添加氢原子
    reduced_pdb = "receptor_reduced.pdb"
    print(f"\n[2/4] 使用reduce添加氢原子...")
    cmd = f"reduce -BUILD {temp_pdb} > {reduced_pdb} 2>/dev/null"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0 or not os.path.exists(reduced_pdb) or os.path.getsize(reduced_pdb) == 0:
        print("  ⚠ reduce不可用或失败，跳过此步骤")
        print("  建议安装: conda install -c conda-forge reduce")
        reduced_pdb = temp_pdb
    else:
        print(f"  ✓ reduce完成")

    # 步骤3: pdb4amber清理
    pdb4amber_pdb = "receptor_pdb4amber.pdb"
    print(f"\n[3/4] 使用pdb4amber清理PDB...")
    cmd = f"pdb4amber -i {reduced_pdb} -o {pdb4amber_pdb} -y --nohyd --noter"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0 or not os.path.exists(pdb4amber_pdb):
        print("  ⚠ pdb4amber失败，使用未清理的PDB")
        print(f"  警告信息: {result.stderr}")
        shutil.copy(reduced_pdb, pdb4amber_pdb)
    else:
        print(f"  ✓ pdb4amber完成")
        # 打印pdb4amber的信息（可能有警告）
        if result.stdout:
            print(f"\n  pdb4amber输出:\n{result.stdout}")

    # 步骤4: RNA残基名规范化
    print(f"\n[4/4] 规范化RNA残基名...")
    normalize_rna_residue_names(pdb4amber_pdb)
    shutil.copy(pdb4amber_pdb, output_pdb)

    # 清理临时文件
    for tmp in [temp_pdb, reduced_pdb, pdb4amber_pdb, "reduce_info.txt"]:
        if os.path.exists(tmp) and tmp != output_pdb:
            os.remove(tmp)

    print(f"\n✓ 受体准备完成: {output_pdb}")
    print(f"{'='*60}\n")
    return output_pdb


def get_charge_from_smiles(smiles):
    """
    从SMILES计算分子净电荷
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"无法解析SMILES: {smiles}")

    charge = Chem.GetFormalCharge(mol)
    return charge


def prepare_ligand(ligand_pdb, ligand_smi, work_pwd):
    """
    配体参数化流程
    1. 从SMILES获取电荷
    2. 使用Antechamber进行AM1-BCC电荷计算
    3. 使用Parmchk2生成力场参数
    """
    start_pwd = os.getcwd()
    os.chdir(work_pwd)

    # 读取SMILES并计算电荷
    with open(ligand_smi, "r") as f:
        smiles = f.read().strip()

    net_charge = get_charge_from_smiles(smiles)
    print(f"配体净电荷: {net_charge}")

    # Antechamber参数化
    ligand_name = os.path.splitext(os.path.basename(ligand_pdb))[0]
    cmd = f"antechamber -i {ligand_pdb} -fi pdb -o {ligand_name}.mol2 -fo mol2 -c bcc -nc {net_charge} -rn LIG -at gaff2 -pf y"
    print(f"运行: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Antechamber警告/错误:\n{result.stderr}")

    # Parmchk2生成缺失参数
    cmd = f"parmchk2 -i {ligand_name}.mol2 -f mol2 -o {ligand_name}.frcmod -a Y"
    print(f"运行: {cmd}")
    subprocess.run(cmd, shell=True)

    os.chdir(start_pwd)
    return f"{ligand_name}.mol2", f"{ligand_name}.frcmod"


def construct_rna_ligand_top(receptor_mol2, ligand_pdb, ligand_smi, work_pwd,
                             mg_m=0.15, k_m=0.1, box_d=15,
                             rna_ff="ol3", ligand_ff="gaff2", wat_ff="tip3p"):
    """
    构建RNA-配体复合物拓扑
    参考 md.py 的 RNA 力场设置
    """
    start_pwd = os.getcwd()
    os.makedirs(work_pwd, exist_ok=True)
    print(f"工作目录: {work_pwd}")

    # 复制文件到工作目录
    shutil.copy(receptor_mol2, os.path.join(work_pwd, "receptor.mol2"))
    shutil.copy(ligand_pdb, os.path.join(work_pwd, "ligand.pdb"))
    shutil.copy(ligand_smi, os.path.join(work_pwd, "ligand.smi"))

    os.chdir(work_pwd)

    # 受体准备流程（MOL2 → PDB + reduce + pdb4amber）
    receptor_pdb = prepare_receptor("receptor.mol2", "receptor_prepared.pdb")

    # 配体参数化
    ligand_mol2, ligand_frcmod = prepare_ligand("ligand.pdb", "ligand.smi", work_pwd)

    # RNA力场设置 (参考 md.py)
    if rna_ff.lower() == "shaw":
        load_rna_ff = "source leaprc.RNA.Shaw"
    elif rna_ff.lower() == "ol3":
        load_rna_ff = "source leaprc.RNA.OL3"
    elif rna_ff.lower() == "yil":
        load_rna_ff = "source leaprc.RNA.YIL"
    elif rna_ff.lower() == "roc":
        load_rna_ff = "source leaprc.RNA.ROC"
    elif rna_ff.lower() == "ljbb":
        load_rna_ff = "source leaprc.RNA.LJbb"
    else:
        raise ValueError(f"不支持的RNA力场: {rna_ff}")

    # 配体力场
    if ligand_ff.lower() == "gaff2":
        load_ligand_ff = "source leaprc.gaff2"
    elif ligand_ff.lower() == "gaff":
        load_ligand_ff = "source leaprc.gaff"
    else:
        raise ValueError(f"不支持的配体力场: {ligand_ff}")

    # 水模型设置 (参考 md.py)
    if wat_ff.lower() == "tip3p":
        load_wat_ff = "source leaprc.water.tip3p"
        wat_box = "TIP3PBOX"
    elif wat_ff.lower() == "opc":
        load_wat_ff = "source leaprc.water.opc"
        wat_box = "OPCBOX"
    elif wat_ff.lower() == "opc3":
        load_wat_ff = "source leaprc.water.opc3"
        wat_box = "OPC3BOX"
    elif wat_ff.lower() == "" and rna_ff.lower() == "shaw":
        load_wat_ff = ""
        wat_box = "TIP4PDBOX"
    else:
        raise ValueError(f"不支持的水模型: {wat_ff}")

    def generate_tleap(mg_ions=0, k_ions=0):
        """生成tleap输入文件"""
        # 离子添加逻辑 (参考 md.py)
        if mg_ions == 0 and k_ions > 0:
            add_ions = f"""
addionsrand complex K+ {k_ions}
addionsrand complex K+ 0
addionsrand complex Cl- 0"""
        elif mg_ions > 0 and k_ions == 0:
            add_ions = f"""
addionsrand complex MG {mg_ions}
addionsrand complex MG 0
addionsrand complex Cl- 0"""
        else:
            add_ions = f"""
addionsrand complex MG {mg_ions}
addionsrand complex K+ {k_ions}
addionsrand complex K+ 0
addionsrand complex Cl- 0"""

        tleap_input = f"""{load_rna_ff}
{load_ligand_ff}
{load_wat_ff}

# 加载配体参数
loadamberparams {ligand_frcmod}
LIG = loadmol2 {ligand_mol2}

# 加载RNA受体（使用PDB格式避免MOL2原子类型问题）
RNA = loadpdb {receptor_pdb}

# 合并形成复合物
complex = combine {{RNA LIG}}

# 保存无溶剂系统
saveamberparm complex complex.nosol.top complex.nosol.inpcrd
savepdb complex complex.nosol.pdb

# 添加溶剂和离子
solvateoct complex {wat_box} {box_d}
{add_ions}

# 检查电荷
charge complex

# 保存溶剂化系统
saveamberparm complex complex.top complex.inpcrd
savepdb complex complex.pdb

quit
"""
        with open("complex.tleap", "w") as f:
            f.write(tleap_input)

    # 第一次运行tleap获取水分子数
    print("第一次运行tleap...")
    generate_tleap()
    subprocess.run("tleap -f complex.tleap > tleap.log", shell=True)

    # 计算需要的离子数 (参考 md.py)
    with open("complex.pdb", "r") as f:
        wat_ids = set()
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if line[17:20] == "WAT":
                    wat_ids.add(line[20:27])

    total_water = len(wat_ids)
    print(f"水分子数: {total_water}")
    mg_ions = int(total_water / 55.56 * mg_m)
    k_ions = int(total_water / 55.56 * k_m)
    print(f"添加Mg2+离子数: {mg_ions}")
    print(f"添加K+离子数: {k_ions}")

    # 第二次运行tleap添加离子
    print("第二次运行tleap...")
    generate_tleap(mg_ions=mg_ions, k_ions=k_ions)
    subprocess.run("tleap -f complex.tleap", shell=True)

    # HMR (氢质量重分配) - 参考 md.py
    parmed_input = """HMassRepartition
parmout complex.hmr.top
go"""
    with open("complex.parmed", "w") as f:
        f.write(parmed_input)

    subprocess.run("parmed -p complex.top -i complex.parmed -O", shell=True)

    os.chdir(start_pwd)
    return "complex.hmr.top"


def run_rna_ligand_md(work_pwd, top, replica=0, temp=350, ns=20, equilibration_ns=2):
    """
    运行RNA-配体复合物MD模拟
    参考 md.py 的高温采样策略

    参数:
        temp: 模拟温度，默认350K用于RNA构象采样（参考md.py）
        ns: 生产MD时长
    """
    TIME_STEP = 0.002  # 2fs (参考 md.py)
    start_pwd = os.getcwd()
    system_name = os.path.basename(work_pwd.rstrip('/'))
    os.chdir(work_pwd)

    # 计算帧数
    frames = int(ns / TIME_STEP * 1000)

    # 获取原子数（用于轨迹输出限制）
    with open("complex.nosol.pdb", "r") as f:
        for line in f:
            words = line.split()
            if words[0] == "ATOM":
                atom_number = int(words[1])

    # 能量最小化
    with open("em.in", "w") as f:
        f.write("""EM
 &cntrl
  imin=1,
  ntx=1, irest=0,
  maxcyc=200, ncyc=100,
  ntpr=10, ntwr=10, ntwx=0,
  cut=10.0,
  ntr=1,
  restraint_wt=0.1,
  restraintmask='@P',
 /
""")

    # 加热 (参考 md.py: 400ps, 0K -> 300K)
    with open("heat.in", "w") as f:
        f.write(f"""Heat
 &cntrl
  imin=0,
  ntx=1, irest=0,
  nstlim=200000, dt=0.002,
  ntf=2, ntc=2,
  tempi=0.0,
  temp0=300.0,
  ntpr=20000, ntwx=20000, ntwr=20000,
  ioutfm=1,
  cut=10.0,
  ntb=1, ntp=0,
  ntt=3, gamma_ln=1.0,
  nmropt=1,
  ntr=1,
  restraint_wt=0.1,
  restraintmask='@P',
 /
&wt type='TEMP0', istep1=5000, istep2=195000, value1=0.0, value2=300.0 /
&wt type='TEMP0', istep1=195001, istep2=200000, value1=300.0, value2=300.0 /
&wt type='END' /
""")

    # 限制性平衡 (参考 md.py: 2ns @ 300K)
    with open("rs.in", "w") as f:
        f.write(f"""Restraint
 &cntrl
  imin=0,
  ntx=1, irest=0,
  nstlim=1000000, dt=0.002,
  ntf=2, ntc=2,
  tempi=300.0,
  temp0=300.0,
  ntpr=20000, ntwx=20000, ntwr=20000,
  ioutfm=1,
  cut=10.0,
  ntb=2, ntp=1,
  ntt=3, gamma_ln=1.0,
  ntr=1,
  restraint_wt=0.1,
  restraintmask='@P',
 /
""")

    # 高温限制性平衡 (参考 md.py: 1ns @ 350K)
    with open(f"{temp}K_rs.in", "w") as f:
        f.write(f"""Restraint
 &cntrl
  imin=0,
  ntx=1, irest=0,
  nstlim=500000, dt={TIME_STEP},
  ntf=2, ntc=2,
  temp0={temp}.0,
  tempi=300.0,
  ntpr=20000, ntwx=20000, ntwr=20000,
  ioutfm=1,
  cut=10.0,
  ntb=1, ntp=0,
  ntt=3, gamma_ln=1.0,
  ntr=1,
  restraint_wt=0.1,
  restraintmask='@P',
 /
""")

    # 生产MD (参考 md.py: 高温NVT采样)
    with open(f"{temp}K.in", "w") as f:
        f.write(f"""High Temperature Production
 &cntrl
  imin=0,
  ntx=1, irest=0,
  nstlim={frames}, dt={TIME_STEP},
  ntf=2, ntc=2,
  temp0={temp}.0,
  tempi=300.0,
  ntpr=50000, ntwx=50000, ntwr=50000,
  ioutfm=1,
  cut=10.0,
  ntwprt={atom_number},
  ntb=1, ntp=0,
  ntt=3, gamma_ln=1.0,
  nmropt=1,
 /
&wt type='END' /
""")

    # SLURM提交脚本 (参考 md.py)
    with open(f"md_{replica}.slurm", "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH --output={system_name}_{replica}.out
#SBATCH --error={system_name}_{replica}.err
#SBATCH --job-name={system_name}_{replica}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=480:00:00
#SBATCH --partition=4090

source /mnt/beegfs/opt/module/tools/modules/init/bash
module load amber24/24

pmemd -O -i em.in -p {top} -c complex.inpcrd -r system.em.{replica}.rst \\
-o system.em.{replica}.log -ref complex.inpcrd

pmemd.cuda -O -i heat.in -p {top} -c system.em.{replica}.rst -r system.heat.{replica}.rst \\
-o system.heat.{replica}.log -x system.heat.{replica}.nc -ref complex.inpcrd

pmemd.cuda -O -i rs.in -p {top} -c system.heat.{replica}.rst -r system.rs.{replica}.rst \\
-o system.rs.{replica}.log -x system.rs.{replica}.nc -ref complex.inpcrd

pmemd.cuda -O -i {temp}K_rs.in -p {top} -c system.rs.{replica}.rst -r system.{temp}K_rs.{replica}.rst \\
-o system.{temp}K_rs.{replica}.log -x system.{temp}K_rs.{replica}.nc -ref complex.inpcrd

pmemd.cuda -O -i {temp}K.in -p {top} -c system.{temp}K_rs.{replica}.rst -r system.{temp}K.{replica}.rst \\
-o system.{temp}K.{replica}.log -x system.{temp}K.{replica}.nc -ref complex.inpcrd

echo "RNA-ligand MD simulation completed!"
""")

    # 提交任务
    subprocess.run(f"sbatch md_{replica}.slurm", shell=True)
    print(f"任务已提交: {system_name}_{replica}")

    os.chdir(start_pwd)


def batch_run_rna_ligand_poses(input_dir, mg_m=0.15, k_m=0.1, box_d=15,
                               temp=350, ns=20, rna_ff="ol3",
                               ligand_ff="gaff2", wat_ff="tip3p"):
    """
    批量处理多个RNA-配体复合物pose的MD模拟

    参数:
        input_dir: 包含ligand.smi, receptor.mol2, pose_*.pdb的文件夹
        mg_m: Mg2+离子浓度 (M)
        k_m: K+离子浓度 (M)
        box_d: 溶剂盒子边界距离 (Å)，RNA默认15
        temp: 模拟温度 (K)，RNA采样默认350K
        ns: 生产MD时长 (ns)
        rna_ff: RNA力场 (ol3/shaw/yil/roc/ljbb)
        ligand_ff: 配体力场 (gaff/gaff2)
        wat_ff: 水模型 (tip3p/opc/opc3)
    """
    input_dir = os.path.abspath(input_dir)
    receptor_mol2 = os.path.join(input_dir, "receptor.mol2")
    ligand_smi = os.path.join(input_dir, "ligand.smi")

    # 检查必需文件
    if not os.path.exists(receptor_mol2):
        raise FileNotFoundError(f"未找到受体文件: {receptor_mol2}")
    if not os.path.exists(ligand_smi):
        raise FileNotFoundError(f"未找到SMILES文件: {ligand_smi}")

    # 找到所有pose
    pose_files = sorted(glob.glob(os.path.join(input_dir, "pose_*.pdb")))
    if not pose_files:
        raise FileNotFoundError(f"未找到pose文件: {input_dir}/pose_*.pdb")

    print(f"找到 {len(pose_files)} 个pose")

    for i, pose_pdb in enumerate(pose_files):
        pose_name = os.path.splitext(os.path.basename(pose_pdb))[0]
        work_pwd = os.path.join(input_dir, f"md_{pose_name}")

        print(f"\n{'='*60}")
        print(f"处理 {pose_name} ({i+1}/{len(pose_files)})")
        print(f"{'='*60}")

        # 构建拓扑
        top = construct_rna_ligand_top(
            receptor_mol2=receptor_mol2,
            ligand_pdb=pose_pdb,
            ligand_smi=ligand_smi,
            work_pwd=work_pwd,
            mg_m=mg_m,
            k_m=k_m,
            box_d=box_d,
            rna_ff=rna_ff,
            ligand_ff=ligand_ff,
            wat_ff=wat_ff
        )

        # 运行MD
        run_rna_ligand_md(
            work_pwd=work_pwd,
            top=top,
            replica=i,
            temp=temp,
            ns=ns
        )

    print(f"\n所有任务已提交!")


if __name__ == "__main__":
    # 示例用法
    import sys

    if len(sys.argv) < 2:
        print("用法: python rna_ligand_md.py <input_dir> [options]")
        print("\n必需文件:")
        print("  - ligand.smi: 配体SMILES")
        print("  - receptor.mol2: RNA受体结构")
        print("  - pose_*.pdb: 对接好的配体构象")
        print("\n可选参数:")
        print("  --mg_m FLOAT: Mg2+浓度 (默认: 0.15 M)")
        print("  --k_m FLOAT: K+浓度 (默认: 0.1 M)")
        print("  --box_d FLOAT: 溶剂盒子边界 (默认: 15 Å)")
        print("  --temp INT: 温度 (默认: 350 K，用于RNA构象采样)")
        print("  --ns INT: 生产MD时长 (默认: 20 ns)")
        print("  --rna_ff STR: RNA力场 (默认: ol3, 可选: shaw/yil/roc/ljbb)")
        print("  --ligand_ff STR: 配体力场 (默认: gaff2)")
        print("  --wat_ff STR: 水模型 (默认: tip3p)")
        print("\n示例:")
        print("  python rna_ligand_md.py auto-complex-example")
        print("  python rna_ligand_md.py auto-complex-example --temp 350 --ns 20 --rna_ff ol3")
        sys.exit(1)

    input_dir = sys.argv[1]

    # 解析参数
    kwargs = {
        'mg_m': 0.15,
        'k_m': 0.1,
        'box_d': 15,
        'temp': 350,
        'ns': 20,
        'rna_ff': 'ol3',
        'ligand_ff': 'gaff2',
        'wat_ff': 'tip3p'
    }

    i = 2
    while i < len(sys.argv):
        if sys.argv[i] == '--mg_m':
            kwargs['mg_m'] = float(sys.argv[i+1])
            i += 2
        elif sys.argv[i] == '--k_m':
            kwargs['k_m'] = float(sys.argv[i+1])
            i += 2
        elif sys.argv[i] == '--box_d':
            kwargs['box_d'] = float(sys.argv[i+1])
            i += 2
        elif sys.argv[i] == '--temp':
            kwargs['temp'] = int(sys.argv[i+1])
            i += 2
        elif sys.argv[i] == '--ns':
            kwargs['ns'] = int(sys.argv[i+1])
            i += 2
        elif sys.argv[i] == '--rna_ff':
            kwargs['rna_ff'] = sys.argv[i+1]
            i += 2
        elif sys.argv[i] == '--ligand_ff':
            kwargs['ligand_ff'] = sys.argv[i+1]
            i += 2
        elif sys.argv[i] == '--wat_ff':
            kwargs['wat_ff'] = sys.argv[i+1]
            i += 2
        else:
            print(f"未知参数: {sys.argv[i]}")
            sys.exit(1)

    batch_run_rna_ligand_poses(input_dir, **kwargs)
