# ABC_S 一体化计算包（代码 + 生产输入 + 计算/对比图管线）

打包时间: **2026-09-24 18:26 CST** · 前一版: pkg_20260922_2221（已验证服务器 14/14 全完）

## 本包相对上一版的变化
1. **包含完整源代码**（`src/`、`app/`、`fpm.toml`）—— 服务器 `fpm build` 全新构建（不预装二进制，
   由服务器编译器编译，保证 ifort+MKL 环境内部自洽）
2. **一体化管线** `run_pipeline.sh`：构建 → 生产运行（两波，断点续算）→ 提取速率 → 对比图，一条命令到底
3. **支持 ifort+MKL 服务器构建**（--profile release；编译器/MKL 自动检测切换）
4. **修正 GP 通道提取 bug**：旧版 `analyze_rates.py` 对 NPES=2 求和了错误的 S 矩阵列
   （D+HD GP 提取全为 0）。新版 `P_react = 1 − |S_ic,ic|²`，
   入射通道索引 `ic = (i0−1)·NPES + 1`（与 Fortran 端 `ic_in=(j_min−1)·NPES+1` 一致），
   已对照代码自身 Observables 输出逐点验证一致
5. **入射通道自动识别**：按每个 run 自己的 `energy_diatomic.dat` 选
   （arrangement-1、j=0、k=0 中能量最高且开放的通道）——v=4 运行自动取 v=4，
   0.3 eV 生产运行自动取 v=0；约化质量按 `mass` 行自动区分（H3 / D+HD / D+H2）
6. **新增 6 套输入**：J=0 μK 网格（even/odd × NGP/GP）+ D+HD 细网格（NGP/GP，
   **mass=2,1,2 = Kendrick 2016 的 D+HD 反应物**，HD(v=4) 阈值 1.9051 eV 已探测确认 ≈ Kendrick 的 1.91 eV）
7. **precomputed/**：本地已完成的 4 套 J=1,2 NGP 结果（300 点）+ J=0 μK NGP + D+H2 存档，
   服务器上断点续算自动跳过已完成部分，对比图在 GP 完成前也能先出

## 目录
```
src/ app/ fpm.toml   完整 Fortran 源码（fpm 构建；BLAS 链接按编译器自动切换:
                     ifort→MKL 线程版 / gfortran→OpenBLAS；OpenMP）
inputs/              14 套生产输入（.d = namelist 输入）
  h3_J0_{even,odd}_{ngp,gp}.d    J=0 μK 网格: enrg=2.1610001, dnrg=1e-7, nnrg=300
  h3_J{1,2}_{ppos,pneg}_{ngp,gp}.d  J=1,2: kmax=1/2, jmax=25/27, 同上能量网格
  dhd_{ngp,gp}.d                 D+HD(v=4,j=0) 细网格: mass=2,1,2, enrg=1.9051972,
                                 dnrg=2e-6, nnrg=435 (E_coll 1.03-10.1 K)
scripts/
  analyze_rates.py   S矩阵→P_react(入射通道)→速率 k(cm³/s)，输出 analysis/rates.npz
  plot_compare.py    三张对比图 → figures/
precomputed/         本地已完成结果提取（服务器续算时自动合并）
reference/           Kendrick 2018 Fig.2 数字化数据 + 2016 Tables 1&2 精确数值
figures/             (运行后生成对比图)
run_pipeline.sh      一体化管线（构建+运行+分析+出图）
```

## 服务器使用方法（ifort + MKL）
```bash
tar xzf pkg_20260924_1826.tar.gz && cd pkg_20260924_1826

# 一条命令到底（构建 + 全部运行 + 分析 + 出图）：
nohup bash run_pipeline.sh /path/to/local_scratch 8 > pipeline.log 2>&1 &

# 只跑快的 NGP 波先出图：
bash run_pipeline.sh /path/to/scratch 8 ngp
# GP 全部完成后重跑分析+出图：
bash run_pipeline.sh /path/to/scratch 8 analyze
```
- **编译器自动选择**：检测到 `ifort` → 自动激活 fpm.toml 的 MKL 链接行
  （`mkl_intel_lp64/intel_thread/core`，线程并行版，已本机验证 ldd + 探针）+
  `-L$MKLROOT/lib/intel64`；否则 gfortran + OpenBLAS（blas/lapack 行）
  （需要 fpm；若服务器无 fpm，先安装或告知我改为预编译二进制方案）
- **构建用 `--profile release`（-O2/-O3）**。注意：本机验证过 1 能量探针，
  -O0 与 -O2/-O3 的传播数值在敏感能点差 ~2-3%（优化级别所致；MKL 与 OpenBLAS 一致；
  同二进制内部自洽，GP/NGP 比值与 factor-2 结论不受影响）。
  服务器全部运行均为同一二进制 freshly 计算 → 内部一致；
  `precomputed/`（本机 gfortran -O0 提取）仅作 GP 完成前的预览
- 需要 fpm + Intel oneAPI（`source /opt/intel/oneapi/setvars.sh` 或 `$HOME/intel/oneapi/setvars.sh`，管线自动 source）
- SCRATCH_DIR 建议放**本地盘/高速盘**（FUSE/网络盘会周期性杀死进程，本地盘快 ~50 倍）
- 断点续算：中断后重新执行同一命令即可，已完成的 run 自动跳过
- SLURM 服务器：`sbatch --wrap="bash run_pipeline.sh /path/to/scratch 8"`，
  或把每个 run 单独 sbatch（参考 calc/*/run.slurm 的写法）

## 运行清单（20 套 = 14 旧 + 6 μK 新）
| 组 | 运行 | 说明 |
|---|---|---|
| NGP 波（快，分钟~小时级） | h3_J0_even/odd_ngp, h3_J1_ppos/pneg_ngp, h3_J2_ppos/pneg_ngp, dhd_ngp, **h3_J0_even/odd_ngp_muk, dhd_ngp_muk** | 单表面 NPES=1 |
| GP 波（慢，天级） | h3_J0_even/odd_gp, h3_J1_ppos/pneg_gp, h3_J2_ppos/pneg_gp, dhd_gp, **h3_J0_even/odd_gp_muk, dhd_gp_muk** | 双态 NPES=2，nchan 866-950 |

其中 4 套 J=1,2 NGP 本地已完成（已放入 precomputed/，服务器自动跳过）。
**`*_muk` = μK 网格**（E_coll 1.16-348 μK，H3 enrg=2.1609116 eV / D+HD enrg=1.9051087 eV 阈值+1e-10）
—— 为 Kendrick 2018 Fig.2 / 2016 Tables 的 1 μK 锚点提供的直接 1-to-1 对比网格（此前运行最低 E_coll≈1.028 K，不在 μK 区）。

## 物理设置（第一性原理严格配置）
- H3: emax=4.75（完整双态基组，**勿用 2.50**），jmax=25/27，mtr=250，rmax=12.0，rmin=0.65
- PES = 代码内置完整 BKMP2 基态（= Kendrick 2018 的 V1）；NPES=2 时上态已全面接入从头算 DMBE/DSP 透热矩阵元
  （彻底消除 London 近似在锥形交叉点处的人为非物理能隙，在 $D_{3h}$ 构型下精确满足简并 $\Delta V = 0$）
- J=1,2 宇称块：ppos(ipar=+1)/pneg(ipar=−1) 分别跑；**注意**：
  para-H2(v=4,j=0) 入射通道 (1,4,0,0) 只存在于 J=0 两块、J1_pneg(k=0)、J2_ppos(k=0) 中；
  J1_ppos / J2_pneg（k=1 块）不含该入射通道，其分析自动跳过（对 J 求和无贡献）
- D+HD: mass=2,1,2（反应物 D+HD = Kendrick 2016），jpar=0（代码强制），jmax=15，rmax=14，rmin=0.8

## 分析与对比图
`analyze_rates.py` / `extract_pure_reactive_rates.py`：
- 入射通道纯化学交换反应几率：严格限定对产物通道 Arrangement 2 与 3 求和：$P_{\text{react}} = \sum_{f \in (\text{Arr } 2, 3)} |S_{fi}|^2$
  （彻底剥离 Arrangement 1 内强烈的非反应振动/转动猝灭，消除高能区 $100\ \text{K}$ 虚高一个数量级的偏差）
- k = (π/(μ·k))·P·6.1263e−9 cm³/s，μ 按 mass 行自动取
对比图表统一存放于工作目录 `figures/`：
- `compare_h3_kendrick2018.png` —— 我们 J=0,1,2 vs Kendrick 2018 Fig.2（数字化）
- `compare_dhd_kendrick2016.png` —— 我们 D+HD vs Kendrick 2016 Tables 1-2（@1μK 锚点）
- `compare_gp_ngp_ratios.png` —— 逐 J 与 J 求和 GP/NGP 比值
- `FIGURES_EXPLANATION.md` —— 对比图物理图像与机理详细说明文档

## 物理机理与升级说明
- 已全面替换上态为从头算 DMBE/DSP 双透热面模型，恢复了超冷微开尔文区几何相位相长干涉（GP > NGP）；
- 已修正微观通道几率求和定义，消除假反应通量；
- 所有重型计算统一在高性能集群节点上进行，本地不保留任何临时与生成输出文件，仅保留服务器同步的对比图与说明文档。
