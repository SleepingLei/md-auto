import os
import subprocess
from collections import deque
import string
import shutil
import math


def get_atom_selection(res):
    if res == "A":
        return "N1"
    elif res == "U":
        return "N3"
    elif res == "G":
        return "N1"
    else:
        return "N3"
    
    
def generate_rst(pdb_pwd, seq, rna_ss, cpptraj_pwd, rst_pwd):
    # 基于二级结构生成距离限制文件
    symbols = ["()", "[]", "{}"]
    # 获取所有的小写字母和大写字母
    lowercase_letters = string.ascii_lowercase
    uppercase_letters = string.ascii_uppercase
    for i in range(len(lowercase_letters)):
        symbols.append(uppercase_letters[i] + lowercase_letters[i])
    all_symbols = "".join(symbols) + "."
    # 将rna单链转换为图
    assert "&" not in rna_ss
    hbonds = []
    stacks = [deque() for i in range(len(symbols))]
    for i in range(len(rna_ss)):
        rna_s = rna_ss[i]
        assert rna_s in all_symbols
        if rna_s == ".":
            continue
        for j in range(len(stacks)):
            if symbols[j][0] == rna_s:
                stacks[j].append(i)
                break
            elif symbols[j][1] == rna_s:
                hbonds.append((stacks[j].pop(), i))
                break
    distance_rst = ""
    for hbond in hbonds:
        res1 = seq[hbond[0]]
        res2 = seq[hbond[1]]
        atom1 = get_atom_selection(res1)
        atom2 = get_atom_selection(res2)
        if res1 in "AU" and res2 in "AU" and res1 != res2:
            distance_rst += f"rst :{hbond[0]+1}@{atom1} :{hbond[1]+1}@{atom2} r1 2 r2 2.5 r3 3.1 r4 5 rk2 10.0 rk3 10.0 out {rst_pwd}\n"
        elif res1 in "GC" and res2 in "GC" and res1 != res2:
            distance_rst += f"rst :{hbond[0]+1}@{atom1} :{hbond[1]+1}@{atom2} r1 2 r2 2.5 r3 3.1 r4 5 rk2 10.0 rk3 10.0 out {rst_pwd}\n"
        else:
            distance_rst += f"rst :{hbond[0]+1}@{atom1} :{hbond[1]+1}@{atom2} r1 2 r2 2.9 r3 4.1 r4 6 rk2 10.0 rk3 10.0 out {rst_pwd}\n"
    with open(cpptraj_pwd, "w") as ifile:
        ifile.write(f"""parm {pdb_pwd}
{distance_rst}                    
run""")
    subprocess.run(f"cpptraj -i {cpptraj_pwd}", shell=True)
    
    
    
def construct_top(pdb_pwd, seq, rna_ss, work_pwd, mg_m=0.15, k_m=0.1, box_d=15, ff="ol3", wat_ff="tip3p"):
    """ M
    """
    start_pwd = os.getcwd()
    os.makedirs(work_pwd, exist_ok=True)
    print(work_pwd)
    shutil.copy(pdb_pwd, os.path.join(work_pwd, "system.pdb"))
    os.chdir(work_pwd)
    # tleap
    if ff.lower() == "shaw":
        load_ff = "source leaprc.RNA.Shaw"
    elif ff.lower() == "ol3":
        load_ff = "source leaprc.RNA.OL3"
    elif ff.lower() == "yil":
        load_ff = "source leaprc.RNA.YIL"
    elif ff.lower() == "roc":
        load_ff = "source leaprc.RNA.ROC"
    elif ff.lower() == "ljbb":
        load_ff = "source leaprc.RNA.LJbb"
    else:
        raise
    if wat_ff.lower() == "tip3p":
        load_wat_ff = "source leaprc.water.tip3p"
        wat_box = "TIP3PBOX"
    elif wat_ff.lower() == "opc":
        load_wat_ff = "source leaprc.water.opc"
        wat_box = "OPCBOX"
    elif wat_ff.lower() == "opc3":
        load_wat_ff = "source leaprc.water.opc3"
        wat_box = " OPC3BOX"
    elif wat_ff.lower() == "" and ff.lower() == "shaw":
        load_wat_ff = ""
        wat_box = " TIP4PDBOX"
    else:
        raise
    def generate_tleap(pdb_pwd, d=box_d, mg_ions=0, k_ions=0):
        if mg_ions == 0 and k_ions > 0:
            add_ions = f"""
addionsrand complex K+ {k_ions}
addionsrand complex K+ 0
addionsrand complex Cl- 0"""
        elif mg_ions >= 0 and k_ions == 0:
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
        pbs = \
    f'''{load_ff}
{load_wat_ff}
complex = loadpdb system.pdb
saveamberparm complex system.nosol.top system.nosol.inpcrd
savepdb complex system.check.nosol.pdb
solvateoct complex {wat_box} {d}
{add_ions}
charge complex
saveamberparm complex system.top system.inpcrd
savepdb complex system.check.pdb
quit
'''
        with open(f"system.tleap", "w") as ifile:
            ifile.write(pbs)

    generate_tleap("system_amber.pdb")
    subprocess.run(f"tleap -f system.tleap > tleap.log", shell=True)
    with open("system.check.pdb", "r") as ifile:
        wat_ids = set()
        for line in ifile:
            if line.startswith("ATOM") and line[17:20] == "WAT":
                wat_ids.add(line[20:27])
    total_water = len(wat_ids)
    print("Water molecules:", total_water)
    mg_ions = int(total_water / 55.56 * mg_m)
    k_ions = int(total_water / 55.56 * k_m)
    generate_tleap("system_amber.pdb", mg_ions=mg_ions, k_ions=k_ions)
    subprocess.run(f"tleap -f system.tleap", shell=True)
    # parmed
    parmed_in = \
        f"""HMassRepartition
parmout system.hmr.top
go"""
    with open(f"system.parmed", "w") as ifile:
        ifile.write(parmed_in)
    subprocess.run(f"parmed -p system.top -i system.parmed -O", shell=True)
    generate_rst("system.check.nosol.pdb", seq, rna_ss, "rst.cpptraj", "distance.rst")
    os.chdir(start_pwd)
    top = "system.hmr.top"
    return top


def run_md(work_pwd, top, replica, temp=350, ns=20):
    TIME_STEP = 0.002
    start_pwd = os.getcwd()
    system_name = work_pwd.split("/")[-1]
    os.chdir(work_pwd)
    frames = int(ns / TIME_STEP * 1000)
    with open(f"system.check.nosol.pdb", "r") as ifile:
        for line in ifile:
            words = line.split()
            if words[0] == "ATOM":
                atom_number = int(words[1])   
                    
    with open(f"em.in", "w") as ifile:
        ifile.write(f"""EM
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
        
    with open(f"heat.in", "w") as ifile:
        ifile.write(f"""Heat
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
    
    with open(f"rs.in", "w") as ifile:
        ifile.write(f"""Restraint
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
    with open(f"{temp}K_rs.in", "w") as ifile:
        ifile.write(f"""Restraint
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
    with open(f"{temp}K.in", "w") as ifile:
        ifile.write(f"""High
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
DISANG=distance.rst
""")
    
    with open(f"md_{replica}.slurm", "w") as ifile:
        ifile.write(f"""#!/bin/bash
#SBATCH --output={system_name}_{replica}.out
#SBATCH --error={system_name}_{replica}.err
#SBATCH --job-name={system_name}_{replica}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8           # 每个任务分配的CPU核心数
#SBATCH --gres=gpu:1
#SBATCH --time=480:00:00             # 运行时间限制 (HH:MM:SS)
#SBATCH --partition=4090        # 分区名

source /mnt/beegfs/opt/module/tools/modules/init/bash
module load amber24/24

pmemd -O -i em.in -p {top} -c system.inpcrd -r system.em.{replica}.rst \
-o system.em.{replica}.log -ref system.inpcrd

pmemd.cuda -O -i heat.in -p {top} -c system.em.{replica}.rst -r system.heat.{replica}.rst \
-o system.heat.{replica}.log -x system.heat.{replica}.nc -ref system.inpcrd

pmemd.cuda -O -i rs.in -p {top} -c system.heat.{replica}.rst -r system.rs.{replica}.rst \
-o system.rs.{replica}.log -x system.rs.{replica}.nc -ref system.inpcrd

pmemd.cuda -O -i {temp}K_rs.in -p {top} -c system.rs.{replica}.rst -r system.{temp}K_rs.{replica}.rst \
-o system.{temp}K_rs.{replica}.log -x system.{temp}K_rs.{replica}.nc -ref system.inpcrd

pmemd.cuda -O -i {temp}K.in -p {top} -c system.{temp}K_rs.{replica}.rst -r system.{temp}K.{replica}.rst \
-o system.{temp}K.{replica}.log -x system.{temp}K.{replica}.nc  -ref system.inpcrd
"""
)
    subprocess.run(f"sbatch md_{replica}.slurm", shell=True)
    os.chdir(start_pwd)
    
def extend_md(work_pwd, top, replica, start_rst, temp=350, ns=20):
    TIME_STEP = 0.002
    start_pwd = os.getcwd()
    system_name = work_pwd.split("/")[-1]
    os.chdir(work_pwd)
    frames = int(ns / TIME_STEP * 1000)
    with open(f"system.check.nosol.pdb", "r") as ifile:
        for line in ifile:
            words = line.split()
            if words[0] == "ATOM":
                atom_number = int(words[1])   
                    
    with open(f"heat.in", "w") as ifile:
        ifile.write(f"""Heat
 &cntrl
  imin=0,
  ntx=1, irest=0,
  nstlim=200000, dt=0.002,
  ntf=2, ntc=2,
  tempi=0.0,
  temp0=300.0,
  ntpr=20000, ntwx=20000, ntwr=20000,
  ioutfm=1,
  cut=9.0,
  ntb=1, ntp=0,
  ntt=3, gamma_ln=2.0,
  nmropt=1,
  ntr=1,
  restraint_wt=0.1,
  restraintmask='@P',
 /
&wt type='TEMP0', istep1=5000, istep2=195000, value1=0.0, value2=300.0 /
&wt type='TEMP0', istep1=195001, istep2=200000, value1=300.0, value2=300.0 /
&wt type='END' /
""")

    with open(f"{temp}K.in", "w") as ifile:
        ifile.write(f"""High
 &cntrl
  imin=0,
  ntx=1, irest=0,
  nstlim={frames}, dt={TIME_STEP},
  ntf=2, ntc=2,
  temp0={temp}.0,
  tempi=300.0,
  ntpr=20000, ntwx=20000, ntwr=20000,
  ioutfm=1,
  cut=9.0,
  ntwprt={atom_number},
  ntb=1, ntp=0,
  ntt=3, gamma_ln=2.0,
  nmropt=1,             
 /
&wt type='END' /
DISANG=distance.rst
""")
    
    with open(f"md_{replica}.slurm", "w") as ifile:
        ifile.write(f"""#!/bin/bash
#SBATCH --output={system_name}_{replica}.out
#SBATCH --error={system_name}_{replica}.err
#SBATCH --job-name={system_name}_{replica}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8           # 每个任务分配的CPU核心数
#SBATCH --gres=gpu:1
#SBATCH --time=480:00:00             # 运行时间限制 (HH:MM:SS)
#SBATCH --partition=4090        # 分区名

source /mnt/beegfs/opt/module/tools/modules/init/bash
module load amber24/24

pmemd.cuda -O -i heat.in -p {top} -c {start_rst} -r system.heat.{replica}.rst \
-o system.heat.{replica}.log -x system.heat.{replica}.nc -ref {start_rst}

pmemd.cuda -O -i {temp}K.in -p {top} -c system.heat.{replica}.rst -r system.{temp}K.{replica}.rst \
-o system.{temp}K.{replica}.log -x system.{temp}K.{replica}.nc -ref {start_rst}
"""
)
    subprocess.run(f"sbatch md_{replica}.slurm", shell=True)
    os.chdir(start_pwd)
    

def dpeaks_top0(work_pwd, parm_pwd, traj_pwds, out_pdb_pwd):
    # cluster setting
    density_cut = 3.0
    distance_cut = 3.0
    align_selection = "@C3'"
    # file name
    start_pwd = os.getcwd()
    os.chdir(work_pwd)
    dvd_file = f"{work_pwd}/dpeaks.dvd"
    trajins = ""
    for traj_pwd in traj_pwds:
        trajins += f"trajin {traj_pwd} 1 last 10\n"
    with open("dpeaks.cpptraj", "w") as ifile:
        ifile.write(f"""
parm {parm_pwd}
{trajins}
autoimage
rmsd {align_selection}
cluster dpeaks epsilon {density_cut} dvdfile {dvd_file} rms {align_selection} nofit
run
""")
    subprocess.run(f"cpptraj -i dpeaks.cpptraj", shell=True)
    with open(dvd_file, "r") as ifile:
        lines = ifile.readlines()[::-1]
    rank = 0
    for line in lines[:-1]:
        words = line.split()
        density = int(words[0])
        distance = float(words[1])
        frame = int(words[2][1:-1])
        break
                
    with open(f"dpeaks_top0.cpptraj", "w") as ifile:
        ifile.write(f"""
parm {parm_pwd}
{trajins}
autoimage
rmsd {align_selection}
trajout {out_pdb_pwd} onlyframes {frame}
""")
    subprocess.run(f"cpptraj -i dpeaks_top0.cpptraj", shell=True)
    os.chdir(start_pwd)
    return out_pdb_pwd


def dpeaks_topN(work_pwd, parm_pwd, traj_pwds, out_pdb_prefix, N=1):
    # cluster setting
    density_cut = 3.0
    distance_cut = 3.0
    align_selection = "@C3'"
    # file name
    start_pwd = os.getcwd()
    os.chdir(work_pwd)
    dvd_file = f"{work_pwd}/dpeaks.dvd"
    trajins = ""
    for traj_pwd in traj_pwds:
        trajins += f"trajin {traj_pwd} 1 last 10\n"
    with open("dpeaks.cpptraj", "w") as ifile:
        ifile.write(f"""
parm {parm_pwd}
{trajins}
autoimage
rmsd {align_selection}
cluster dpeaks epsilon {density_cut} dvdfile {dvd_file} rms {align_selection} nofit
run
""")
    subprocess.run(f"cpptraj -i dpeaks.cpptraj", shell=True)
    with open(dvd_file, "r") as ifile:
        lines = ifile.readlines()[::-1]
    rank = 0
    out_frames = []
    for line in lines[:-1]:
        words = line.split()
        density = int(words[0])
        distance = float(words[1])
        frame = int(words[2][1:-1])
        if len(out_frames) < N and distance > distance_cut:
            out_frames.append(frame)
        elif len(out_frames) >= N:
            break
    for i, frame in enumerate(out_frames):            
        with open(f"dpeaks_top{i}.cpptraj", "w") as ifile:
            ifile.write(f"""
parm {parm_pwd}
{trajins}
autoimage
strip :WAT,Na+,Cl-,MG
rmsd {align_selection}
trajout {out_pdb_prefix}_top{i}.pdb onlyframes {frame}
    """)
        subprocess.run(f"cpptraj -i dpeaks_top{i}.cpptraj", shell=True)
    os.chdir(start_pwd)
    return


def kmeans(work_pwd, parm_pwd, traj_pwds, out_pdb_prefix, clusters, interval=10):
    # cluster setting
    align_selection = "@C3'"
    # file name
    start_pwd = os.getcwd()
    os.chdir(work_pwd)
    trajins = ""
    for traj_pwd in traj_pwds:
        trajins += f"trajin {traj_pwd} 1 last {interval}\n"
    with open("kmeans.cpptraj", "w") as ifile:
        ifile.write(f"""
parm {parm_pwd}
{trajins}
autoimage
rmsd {align_selection}
cluster kmeans clusters {clusters} repout {out_pdb_prefix} repfmt pdb rms {align_selection} nofit
run
""")
    subprocess.run(f"cpptraj -i kmeans.cpptraj", shell=True)
    os.chdir(start_pwd)
    
