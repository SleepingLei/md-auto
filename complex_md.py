import os
import subprocess
import shutil
import glob
from rdkit import Chem


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


def construct_complex_top(receptor_mol2, ligand_pdb, ligand_smi, work_pwd,
                          na_m=0.15, box_d=10, protein_ff="ff14SB",
                          ligand_ff="gaff2", wat_ff="tip3p"):
    """
    构建蛋白质-配体复合物拓扑
    """
    start_pwd = os.getcwd()
    os.makedirs(work_pwd, exist_ok=True)
    print(f"工作目录: {work_pwd}")

    # 复制文件到工作目录
    shutil.copy(receptor_mol2, os.path.join(work_pwd, "receptor.mol2"))
    shutil.copy(ligand_pdb, os.path.join(work_pwd, "ligand.pdb"))
    shutil.copy(ligand_smi, os.path.join(work_pwd, "ligand.smi"))

    os.chdir(work_pwd)

    # 配体参数化
    ligand_mol2, ligand_frcmod = prepare_ligand("ligand.pdb", "ligand.smi", work_pwd)

    # 设置力场
    if protein_ff.lower() == "ff14sb":
        load_protein_ff = "source leaprc.protein.ff14SB"
    elif protein_ff.lower() == "ff19sb":
        load_protein_ff = "source leaprc.protein.ff19SB"
    else:
        load_protein_ff = f"source leaprc.protein.{protein_ff}"

    if ligand_ff.lower() == "gaff2":
        load_ligand_ff = "source leaprc.gaff2"
    elif ligand_ff.lower() == "gaff":
        load_ligand_ff = "source leaprc.gaff"
    else:
        raise ValueError(f"不支持的配体力场: {ligand_ff}")

    if wat_ff.lower() == "tip3p":
        load_wat_ff = "source leaprc.water.tip3p"
        wat_box = "TIP3PBOX"
    elif wat_ff.lower() == "opc":
        load_wat_ff = "source leaprc.water.opc"
        wat_box = "OPCBOX"
    elif wat_ff.lower() == "tip4pew":
        load_wat_ff = "source leaprc.water.tip4pew"
        wat_box = "TIP4PEWBOX"
    else:
        raise ValueError(f"不支持的水模型: {wat_ff}")

    def generate_tleap(na_ions=0):
        """生成tleap输入文件"""
        if na_ions == 0:
            add_ions = """
addionsrand complex Na+ 0
addionsrand complex Cl- 0"""
        else:
            add_ions = f"""
addionsrand complex Na+ {na_ions}
addionsrand complex Na+ 0
addionsrand complex Cl- 0"""

        tleap_input = f"""{load_protein_ff}
{load_ligand_ff}
{load_wat_ff}

# 加载配体参数
loadamberparams {ligand_frcmod}
LIG = loadmol2 {ligand_mol2}

# 加载受体
REC = loadmol2 receptor.mol2

# 合并形成复合物
complex = combine {{REC LIG}}

# 保存无溶剂系统
saveamberparm complex complex.nosol.top complex.nosol.inpcrd
savepdb complex complex.nosol.pdb

# 添加溶剂和离子
solvatebox complex {wat_box} {box_d}
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

    # 计算需要的离子数
    with open("complex.pdb", "r") as f:
        wat_ids = set()
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if "WAT" in line or "HOH" in line:
                    wat_ids.add(line[22:26].strip())

    total_water = len(wat_ids)
    print(f"水分子数: {total_water}")
    na_ions = int(total_water / 55.56 * na_m)
    print(f"添加Na+离子数: {na_ions}")

    # 第二次运行tleap添加离子
    print("第二次运行tleap...")
    generate_tleap(na_ions=na_ions)
    subprocess.run("tleap -f complex.tleap", shell=True)

    # HMR (氢质量重分配)
    parmed_input = """HMassRepartition
parmout complex.hmr.top
go"""
    with open("complex.parmed", "w") as f:
        f.write(parmed_input)

    subprocess.run("parmed -p complex.top -i complex.parmed -O", shell=True)

    os.chdir(start_pwd)
    return "complex.hmr.top"


def run_complex_md(work_pwd, top, replica=0, temp=300, ns=10, equilibration_ns=2):
    """
    运行蛋白质-配体复合物MD模拟
    包括: EM -> Heat -> Equilibration -> Production
    """
    TIME_STEP = 0.004  # 4fs with HMR
    start_pwd = os.getcwd()
    system_name = os.path.basename(work_pwd.rstrip('/'))
    os.chdir(work_pwd)

    # 计算帧数
    eq_frames = int(equilibration_ns / TIME_STEP * 1000)
    prod_frames = int(ns / TIME_STEP * 1000)

    # 能量最小化
    with open("01_em.in", "w") as f:
        f.write("""Energy Minimization
 &cntrl
  imin=1,
  ntx=1, irest=0,
  maxcyc=5000, ncyc=2500,
  ntpr=100, ntwr=500, ntwx=0,
  cut=10.0,
  ntb=1,
  ntr=1,
  restraint_wt=10.0,
  restraintmask='!@H= & !:WAT,Na+,Cl-',
 /
""")

    # 加热 (50ps, 0K -> temp K)
    with open("02_heat.in", "w") as f:
        f.write(f"""Heating
 &cntrl
  imin=0,
  ntx=1, irest=0,
  nstlim=12500, dt=0.004,
  ntf=2, ntc=2,
  tempi=0.0, temp0={temp}.0,
  ntpr=500, ntwx=500, ntwr=500,
  ioutfm=1, ntxo=2,
  cut=10.0,
  ntb=1, ntp=0,
  ntt=3, gamma_ln=2.0,
  nmropt=1,
  ntr=1,
  restraint_wt=5.0,
  restraintmask='!@H= & !:WAT,Na+,Cl-',
 /
 &wt type='TEMP0', istep1=0, istep2=12500, value1=0.0, value2={temp}.0 /
 &wt type='END' /
""")

    # 密度平衡 (NTP, 500ps)
    with open("03_equil_npt.in", "w") as f:
        f.write(f"""Density Equilibration
 &cntrl
  imin=0,
  ntx=1, irest=0,
  nstlim=125000, dt=0.004,
  ntf=2, ntc=2,
  temp0={temp}.0,
  ntpr=5000, ntwx=5000, ntwr=5000,
  ioutfm=1, ntxo=2,
  cut=10.0,
  ntb=2, ntp=1, barostat=2,
  ntt=3, gamma_ln=2.0,
  ntr=1,
  restraint_wt=2.0,
  restraintmask='!@H= & !:WAT,Na+,Cl-',
 /
""")

    # 无限制平衡 (NTP)
    with open("04_equil_free.in", "w") as f:
        f.write(f"""Free Equilibration
 &cntrl
  imin=0,
  ntx=5, irest=1,
  nstlim={eq_frames}, dt={TIME_STEP},
  ntf=2, ntc=2,
  temp0={temp}.0,
  ntpr=10000, ntwx=10000, ntwr=10000,
  ioutfm=1, ntxo=2,
  cut=10.0,
  ntb=2, ntp=1, barostat=2,
  ntt=3, gamma_ln=2.0,
 /
""")

    # 生产MD (NVT)
    with open("05_prod.in", "w") as f:
        f.write(f"""Production MD
 &cntrl
  imin=0,
  ntx=5, irest=1,
  nstlim={prod_frames}, dt={TIME_STEP},
  ntf=2, ntc=2,
  temp0={temp}.0,
  ntpr=25000, ntwx=25000, ntwr=25000,
  ioutfm=1, ntxo=2,
  cut=10.0,
  ntb=1, ntp=0,
  ntt=3, gamma_ln=2.0,
 /
""")

    # SLURM提交脚本
    with open(f"run_md_{replica}.slurm", "w") as f:
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

# 能量最小化
pmemd.cuda -O -i 01_em.in -p {top} -c complex.inpcrd -r 01_em.rst \\
-o 01_em.{replica}.log -ref complex.inpcrd

# 加热
pmemd.cuda -O -i 02_heat.in -p {top} -c 01_em.rst -r 02_heat.rst \\
-o 02_heat.{replica}.log -x 02_heat.{replica}.nc -ref 01_em.rst

# NPT平衡
pmemd.cuda -O -i 03_equil_npt.in -p {top} -c 02_heat.rst -r 03_equil_npt.rst \\
-o 03_equil_npt.{replica}.log -x 03_equil_npt.{replica}.nc -ref 02_heat.rst

# 无限制平衡
pmemd.cuda -O -i 04_equil_free.in -p {top} -c 03_equil_npt.rst -r 04_equil_free.rst \\
-o 04_equil_free.{replica}.log -x 04_equil_free.{replica}.nc

# 生产MD
pmemd.cuda -O -i 05_prod.in -p {top} -c 04_equil_free.rst -r 05_prod.{replica}.rst \\
-o 05_prod.{replica}.log -x 05_prod.{replica}.nc

echo "MD simulation completed!"
""")

    # 提交任务
    subprocess.run(f"sbatch run_md_{replica}.slurm", shell=True)
    print(f"任务已提交: {system_name}_{replica}")

    os.chdir(start_pwd)


def batch_run_poses(input_dir, na_m=0.15, box_d=10, temp=300, ns=10,
                    protein_ff="ff14SB", ligand_ff="gaff2", wat_ff="tip3p"):
    """
    批量处理多个pose的MD模拟

    参数:
        input_dir: 包含ligand.smi, receptor.mol2, pose_*.pdb的文件夹
        na_m: Na+离子浓度 (M)
        box_d: 溶剂盒子边界距离 (Å)
        temp: 模拟温度 (K)
        ns: 生产MD时长 (ns)
        protein_ff: 蛋白质力场
        ligand_ff: 配体力场
        wat_ff: 水模型
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
        top = construct_complex_top(
            receptor_mol2=receptor_mol2,
            ligand_pdb=pose_pdb,
            ligand_smi=ligand_smi,
            work_pwd=work_pwd,
            na_m=na_m,
            box_d=box_d,
            protein_ff=protein_ff,
            ligand_ff=ligand_ff,
            wat_ff=wat_ff
        )

        # 运行MD
        run_complex_md(
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
        print("用法: python complex_md.py <input_dir> [options]")
        print("\n必需文件:")
        print("  - ligand.smi: 配体SMILES")
        print("  - receptor.mol2: 受体结构")
        print("  - pose_*.pdb: 对接好的配体构象")
        print("\n可选参数:")
        print("  --na_m FLOAT: Na+浓度 (默认: 0.15 M)")
        print("  --box_d FLOAT: 溶剂盒子边界 (默认: 10 Å)")
        print("  --temp INT: 温度 (默认: 300 K)")
        print("  --ns INT: 生产MD时长 (默认: 10 ns)")
        print("  --protein_ff STR: 蛋白质力场 (默认: ff14SB)")
        print("  --ligand_ff STR: 配体力场 (默认: gaff2)")
        print("  --wat_ff STR: 水模型 (默认: tip3p)")
        print("\n示例:")
        print("  python complex_md.py auto-complex-example")
        print("  python complex_md.py auto-complex-example --temp 310 --ns 20")
        sys.exit(1)

    input_dir = sys.argv[1]

    # 解析参数
    kwargs = {
        'na_m': 0.15,
        'box_d': 10,
        'temp': 300,
        'ns': 10,
        'protein_ff': 'ff14SB',
        'ligand_ff': 'gaff2',
        'wat_ff': 'tip3p'
    }

    i = 2
    while i < len(sys.argv):
        if sys.argv[i] == '--na_m':
            kwargs['na_m'] = float(sys.argv[i+1])
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
        elif sys.argv[i] == '--protein_ff':
            kwargs['protein_ff'] = sys.argv[i+1]
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

    batch_run_poses(input_dir, **kwargs)
