# 用rosetta对初始结构建模
import os
import subprocess

def pdb_modify_name(pdb_pwd, out_pdb_pwd):
    with open(pdb_pwd, "r") as ifile:
        lines = ifile.readlines()
        out_lines = []
        for line in lines:
            if not line.startswith("ATOM"):
                continue
            resid = line[22:26]
            res_atom = line[12:16].strip()
            atom_type = line[77]
            if int(resid) == 1 and res_atom in ["P", "OP1", "OP2"]:
                continue
            elif atom_type == "H":
                continue
            out_lines.append(line)
    with open(out_pdb_pwd, "w") as ifile:
        for line in out_lines:
            ifile.write(line)

def sbatch_farfar2(work_dir, seq, rna_ss, out_pdb):
    start_pwd = os.getcwd()
    os.makedirs(work_dir, exist_ok=True)
    os.chdir(work_dir)
    with open(f"farfar2.slurm", "w") as ifile:
        ifile.write(f"""#!/bin/bash
#SBATCH --job-name=farfar2
#SBATCH --output=farfar2.out
#SBATCH --error=farfar2.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=240:00:00
#SBATCH --partition=4090

source /mnt/beegfs/opt/module/tools/modules/init/bash
module load rosetta/3.14-mpi

rna_denovo.cxx11threadmpiserialization.linuxgccrelease -sequence "{seq.lower()}" -secstruct "{rna_ss}" -nstruct 1 -out:file:silent test.out -minimize_rna true -overwrite true
rna_minimize.cxx11threadmpiserialization.linuxgccrelease -in:file:silent test.out -out:pdb true
cp S_000001_minimize.pdb {out_pdb}
"""
)
    subprocess.run(f"sbatch farfar2.slurm", shell=True)
    os.chdir(start_pwd)


input_txt = "/gene/home/ldw/Desktop/RNA_ensemble/input/mir_list_20251202.txt"
save_dir = "/gene/home/ldw/Desktop/RNA_ensemble/input/mir_list_20251202"
rna_infos = {}
with open(input_txt, "r") as ifile:
    for line in ifile.readlines():
        words = line.split()
        rna_infos[words[0]] = {"seq": words[1], "rna_ss": words[2]}

for name in rna_infos:
    seq = rna_infos[name]['seq']
    rna_ss = rna_infos[name]['rna_ss']
    if len(seq) != len(rna_ss):
        print(name, len(seq), len(rna_ss), seq, rna_ss)
    work_dir = os.path.join(save_dir, name)
    out_pdb = os.path.join(save_dir, f"{name}.pdb")
    sbatch_farfar2(work_dir, seq, rna_ss, out_pdb)
