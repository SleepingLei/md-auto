import os
import subprocess
import numpy as np
import math
import sys


def hie_pocket_topN(work_pwd, parm_pwd, traj_pwds, out_pwd, out_prefix, frame_per_ns, frame_interval, start_ns, end_ns, pocket_resids, N=1, nofit=False):
    # cluster setting
    if len(pocket_resids) > 0:
        align_resids = ":"
        for start_pos, end_pos in pocket_resids:
            align_resids += f"{start_pos}-{end_pos},"
        align_selection = align_resids + "@C3'"
    else:
        align_selection = "@C3'"
    if nofit:
        fit_parm = "nofit"
    else:
        fit_parm = ""
    # file name
    start_pwd = os.getcwd()
    os.chdir(work_pwd)
    start_frame = int(frame_per_ns * start_ns)
    end_frame = int(frame_per_ns * end_ns)
    trajins = ""
    for traj_pwd in traj_pwds:
        if os.path.exists(traj_pwd):
            trajins += f"trajin {traj_pwd} {start_frame} {end_frame} {frame_interval}\n"
    with open(f"hie_pocket_{start_ns}-{end_ns}ns_N{N}_i{frame_interval}.cpptraj", "w") as ifile:
        ifile.write(f"""
parm {parm_pwd}
{trajins}
autoimage
rmsd {align_selection}
cluster clusters {N} rms {align_selection} {fit_parm} repout {out_pwd}/{out_prefix}_pocket_{start_ns}-{end_ns}ns_hie repfmt pdb summary {out_pwd}/{out_prefix}_pocket_{start_ns}-{end_ns}ns_hie.dat
run
""")
    subprocess.run(f"cpptraj -i hie_pocket_{start_ns}-{end_ns}ns_N{N}_i{frame_interval}.cpptraj > cpptraj.log", shell=True)
    os.chdir(start_pwd)
    return

def find_all_positions(text, pattern):
    positions = []
    start = 0
    while start < len(text):
        idx = text.find(pattern, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1  # 从下一个位置继续搜索
    return positions

def get_pocket_resid(seq, pocket_seq):
    if pocket_seq == None:
        return []
    pocket_seqs = pocket_seq.strip("...").split("...")
    pocket_resids = []
    for sub_seq in pocket_seqs:
        positions = find_all_positions(seq, sub_seq)
        if len(positions) != 1:
            raise
        pocket_resids.append((positions[0], positions[0]+len(sub_seq)))
    return pocket_resids


traj_ns = 500  # 500 ns
num_of_replica = 3  # 平行轨迹数
box_d = 12  # 水盒子大小
k_m = 0.12  # 0.12 M K+
mg_m = 0.02  # 0.02 M Mg2+
temp = 400  # 400 K
ff = "ljbb"  # RNA force field
wat_ff = "opc3"  # water force field
input_txt = "/gene/home/ldw/Desktop/RNA_ensemble/input/mir_list_20251202.txt"
save_dir = "/gene/home/ldw/Desktop/RNA_ensemble/input/mir_list_20251202"  # 初始pdb结构路径
md_pwd = "/gene/home/ldw/Desktop/RNA_ensemble/md/mir_list_20251202"  # md路径
out_pwd = "/gene/home/ldw/Desktop/RNA_ensemble/analysis/mir_list_20251202"  # 聚类输出路径
frame_per_ns = 10
start_ns = 0
end_ns = traj_ns
N = 10  # 聚类输出结构数目
frame_interval = 2

os.makedirs(out_pwd, exist_ok=True)
rna_infos = {}
with open(input_txt, "r") as ifile:
    for line in ifile:
        words = line.split()
        rna_infos[words[0]] = {"seq": words[1], "rna_ss": words[2]}
        pocket_resids = []
        for pocket_resid in words[3:]:
            pocket_resid = pocket_resid.split("-")
            pocket_resids.append((int(pocket_resid[0]), int(pocket_resid[1])))
        rna_infos[words[0]]["pocket_resids"] = pocket_resids

for rna_name in list(rna_infos.keys())[:]:
    md_name = f"{rna_name}_{ff}_{wat_ff}_{temp}K_{traj_ns}ns_{k_m}k_{mg_m}mg"
    work_pwd = os.path.join(md_pwd, md_name)
    parm_pwd = "system.nosol.top"
    traj_pwds = [f"system.{temp}K.{i}.nc" for i in range(num_of_replica)]
    pocket_resids = rna_infos[rna_name]["pocket_resids"]
    print(rna_name, pocket_resids)
    hie_pocket_topN(work_pwd, parm_pwd, traj_pwds, out_pwd, md_name, frame_per_ns, frame_interval, start_ns, end_ns, pocket_resids, N=N, nofit=False)
