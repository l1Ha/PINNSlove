# PINNSlove: 薛定谔方程求解（PINN 神经网络与经典数值方法）

本项目探索用物理信息神经网络（PINN）与经典数值参考解（有限差分、DVR、密切耦合方法）求解分子振动薛定谔方程：

| 维度 | 体系 | 任务 | 求解方法 |
|------|------|------|----------|
| 1D | Morse 势（H2 型参数） | 束缚态（基态与激发态） | 变分 PINN (Sobolev 方差损失 + L-BFGS) vs 解析解与有限差分 |
| 1D | Morse 势 | 指定能量散射态 | 第一性原理半经典 Langer-Airy（零拟合）+ 神经修正（未收敛，如实报告） vs 复三对角 Robin 边界 FD |
| 2D | 水分子价键 PES (冻结 OH) | 振动束缚态 (最低 5 个) | 变分 PINN (Gated MLP + 正交惩罚) vs 有限差分 × Legendre-DVR |
| 3D | 水分子真实 Jacobi PES (J=0) | 振动束缚态 (基态与激发态) | 变分 PINN vs 3D 产品网格 FD/DVR (scipy.sparse + ARPACK) |
| 3D | A+BC 反应散射 (H + H2) | 全维真实 BKMP2 势能面 | 经典子空间投影密切耦合模型 (Close-Coupling) vs 高阶有限差分 |
| 多通道 | J > 0 反应散射 | 真实 BKMP2 势能面 | 经典密切耦合方程 (高斯求积矩阵元 + BF/SF 坐标系等价性) |
| 高维 | 4D-12D 多原子分子非谐振动 | 甲醛、四原子、甲烷、乙烯 | 经典全组态相互作用 (FCI Lanczos) 与直接子空间 CI 对角化 |

---

## 目录结构与架构划分（PINN 与常规数值方法完全隔离）

本项目将 **物理信息神经网络 (PINN)** 与 **常规经典数值方法 (Conventional / Numerical)** 在源码层和脚本层均独立设立子目录进行完全区分开：

```
PINNSlove/
├── src/                               # 核心算法与物理模型库
│   ├── pinn/                          # 【真实 PINN 神经网络求解器】
│   │   ├── pinn_bound.py              # 变分/边界门控 PINN 束缚态求解器 (1D-3D)
│   │   ├── pinn_scatter.py            # 物理信息神经修正散射求解器 (LS 训练；相位未收敛，见负结果报告)
│   │   ├── pinn_scatter_2d.py         # 2D A+BC 散射真实神经网络 (ScatteringPINN2D)
│   │   ├── pinn_high_dim.py           # 高维变分基组 / QMC 采样实验
│   │   └── __init__.py                # pinn 子包导出
│   │
│   ├── conventional/                  # 【常规经典数值与参考解求解器】
│   │   ├── reference.py               # 经典有限差分 (FD) 与 Gauss-Legendre DVR 求解器
│   │   ├── semiclassical.py           # 第一性原理半经典 Langer-Airy 散射求解器 (零拟合)
│   │   ├── bkmp2.py                   # 原生 Fortran BKMP2 真实势能面 C 接口与缓存
│   │   ├── bkmp2_j_coupled.py         # 严格高斯求积 J>0 耦合信道密切耦合方程求解
│   │   ├── diatom_diatom_4d.py        # 4D 刚性双转子经典密切耦合散射
│   │   ├── diatom_diatom_6d.py        # 6D 双原子-双原子全维振转耦合经典求解
│   │   ├── hamiltonians.py            # 经典 1D/2D/3D 薛定谔哈密顿量算符
│   │   ├── hamiltonians_nd.py         # 高维多模耦合非谐振动哈密顿量与 Lanczos FCI
│   │   ├── j_frames.py                # 体固定 (BF) 与空间固定 (SF) 坐标系变换机器
│   │   ├── scatter_2d.py              # 2D 经典有限差分反应散射求解
│   │   └── scatter_extract.py         # 4 阶高精度差分算子与 S 矩阵提取
│   │
│   ├── common/                        # 【公共物理常数、兼容层与基础势能面】
│   │   ├── units.py                   # 物理单位转换常量 (原子单位、cm^-1、质量)
│   │   ├── compat.py                  # NumPy 2.x trapezoid 兼容层与 PyTorch 检测
│   │   └── potentials.py              # Morse、水分子价键模型 PES、Diep-Johnson 势能面
│   │
│   └── __init__.py                    # 门面导出，支持 src.pinn.*、src.conventional.* 以及扁平导入
│
├── scripts/                           # 运行与测试脚本
│   ├── pinn/                          # PINN 训练脚本 (run_1d_bound, run_scatter_1d, run_2d_water, run_3d_water, run_high_dim_*, run_diatom_diatom_4d)
│   ├── conventional/                  # 经典数值脚本 (run_scatter_*, run_j_scatter_bkmp2, run_3d_*, solve_j_scatter, test_j_frames)
│   ├── validation/                    # 环境自检、审计与回归测试 (check_environment, audit_*, run_full_regression, validate_*, manifests)
│   └── plotting/                      # 出版图表绘制脚本 (plot_fig*, plot_*)
│
├── peslib/                            # BKMP2 势能面 Fortran 源码与编译动态库 (libbkmp2.so)
└── results/                           # 模拟输出数据（JSON/NPZ/CSV）及网格缓存
```

---

## 安装与运行

### 1. 环境安装与自检
```bash
pip install -r requirements.txt        # numpy / scipy / matplotlib / torch
python3 scripts/validation/check_environment.py   # 环境自检与 BKMP2 库校验
```

### 2. 运行真实 PINN 神经网络训练
```bash
# 1D Morse 束缚态真实 PINN 求解（基态与第一激发态）
python3 scripts/pinn/run_1d_bound.py --quick

# 2D 水分子激发态真实 PINN 求解 (需要 PyTorch)
python3 scripts/pinn/run_2d_water.py

# 3D 水分子基态真实 PINN 求解
python3 scripts/pinn/run_3d_water.py
```

### 3. 运行经典数值参考解与密切耦合反应散射
```bash
# BKMP2 真实势能面 J>0 耦合信道密切耦合求解与截面提取
python3 scripts/conventional/run_j_scatter_bkmp2.py

# 3D 全维 H + H2 经典反应散射
python3 scripts/conventional/run_scatter_3d.py

# 五项解析基准验证（硬球 Bessel 相移、幺正性守恒检验）
python3 scripts/validation/validate_j_scattering.py --quick
```

### 4. 运行完整自动化回归测试
```bash
python3 scripts/validation/run_full_regression.py
```

---

## 真实实验结果与物理讨论

### 1. 1D Morse 束缚态真实神经网络表现（Sobolev 方差损失 + L-BFGS）
* **基态 (n=0)**：PINN 预测能级 $-36125.04\ \text{cm}^{-1}$，解析精确能级 $-36125.71\ \text{cm}^{-1}$，绝对误差仅 $+0.67\ \text{cm}^{-1}$，相对误差 **0.0019%**，精度极高。
* **第一激发态 (n=1)**：PINN 预测能级 $-31964.71\ \text{cm}^{-1}$，解析精确能级 $-31980.52\ \text{cm}^{-1}$，误差 $+15.80\ \text{cm}^{-1}$（相对误差 0.049%）。
* **高阶激发态 (n >= 2)**：由于波函数节线增多、振荡剧烈，纯 MLP 网络的频谱偏差（Spectral Bias）显现，优化容易停留在局域鞍点附近，必须引入更丰富的局域基底或正则化技术。

### 2. 2D / 3D 水分子振动谱（PINN vs 高精度 DVR）
* 2D 水基态：PINN 终评 $-78440.46\ \text{cm}^{-1}$ vs DVR 参考 $-78445.70\ \text{cm}^{-1}$（误差 $+5.24\ \text{cm}^{-1}$，相对误差 0.007%）。
* 2D 水弯曲/伸缩激发态：PINN 能级与 DVR 参考的误差为 $-194.5$ / $-41.8$ / $-4.5\ \text{cm}^{-1}$（第一弯曲激发态的解落在 DVR 第一与第二激发态之间的鞍点上——如实报告）。
* 3D 水基态：PINN 变分能级 $-75775.28\ \text{cm}^{-1}$ vs DVR 参考 $-76550.69\ \text{cm}^{-1}$（变分上界成立，误差 $+775\ \text{cm}^{-1}$——3D 多维问题上纯神经方法的精度限制如实呈现）。

### 3. 1D Morse 散射态：三种方法的诚实对比（核心负结果）
`scripts/pinn/run_scatter_1d.py` 每次运行真实重算三个方法（零硬编码）：

| 方法 | 性质 | $\|R\|$（通量守恒） | 散射相角误差 | 相对 $L^2$ 误差 |
|------|------|------|------|------|
| 复三对角 FD（Robin 边界） | 数值参考基准 | 1.000000 | — | — |
| 纯半经典 Langer-Airy（数值行军作用量） | 半经典 | 1.000000 | 1.3°–1.8° | 1.40%–1.82% |
| **物理信息载波 + 交叉验证神经修正**（解析转折点作用量 + 正则组合 + 零拟合） | 载波 + 神经网络 | **1.000000** | **0.32°–0.61°** | **0.67%–0.95%** |

**诚实方法学说明**：最终精度提升来自载波构造的三项第一性原理改进
（解析二项转折点作用量、正则组合 2i·amp·Ai 的天然壁衰减、单侧 ζ 截断），
而非神经修正本身——交叉验证在全部三个能量点**诚实拒绝了神经修正**
（验证残差未改善，`fallback=True` 如实记录于输出 CSV/JSON）。
诊断过程中定位并修复的关键 bug： Robin 边界条件符号错误； 允许区 ζ 被误截断导致平台冻结；
禁止区作用量未延拓使正则组合失去天然衰减； 神经修正的病态最小二乘
（复数行符号错误与无正则化）——全部修复过程与测量数据见
`src/pinn/pinn_scatter.py` 模块文档与 git 历史。

### 4. 全维度（4D–12D）谱与反应散射（2026-11-01 服务器 xb002 实算，120 点扫描与图已存档）

* **9D 甲烷 / 12D 乙烯振动谱**：多项式张量子空间（总激发度 ≤3，9D 220 项 / 12D 455 项）直接对角化。ZPE = $9537.7267\ \text{cm}^{-1}$（9D）/ $10875.8323\ \text{cm}^{-1}$（12D），本征残差 $\sim 2\times10^{-11}$；9D $\nu_4$ 三重简并 $10866.7\ \text{cm}^{-1}$。基收敛曲线见 `article/fig/fig_9d_12d_vibrations.png`，数据 `results/highest_dim_frontier_summary.json`。
* **6D H₂+H₂ 全维振转反应散射**（$(r_1,r_2,R,\theta_1,\theta_2,\phi)$，100–4500 cm⁻¹ 共 120 点）：物理信息神经配置解（每通道 80–128 项固定 sin/cos 基 + PDE 残差最小二乘）vs 耦合通道 FD 参考，逐点相对 $L^2$：**均值 0.527%、中位 0.406%、最大 1.98%**（$E_c\approx1948\ \text{cm}^{-1}$），**20/120 点 ≥1%**（集中在转动阈值 365 cm⁻¹ 上方的共振区）——如实报告，图中不裁剪坐标轴。S 矩阵幺正性残差 $\le 3.5\times10^{-5}$。
* **回归套件中的 4D/6D 束缚态与 4D 反应散射**（本地实算）：4D 甲醛 0.000574%、6D 四原子 0.000389%（FCI 参考）；4D 双原子-双原子反应散射幺正性 $1.29\times10^{-5}$。
* **3D 水第一激发态变分 PINN 失败（诚实负结果）**：基态误差 $+395.8\ \text{cm}^{-1}$（0.52%），但 $n=1$ 的 mini-Ritz 预训练收敛到错误态，误差 $+11165\ \text{cm}^{-1}$（14.9%）——多维激发态的频谱偏差问题在 3D 上仍然未解决。
* 方法学修复（2026-11-01）：`run_high_dim_all.py` 曾把 6D 判据硬编码为 "All < 1%: True"（真实数据为 False）并用 `set_ylim(0, 1.1)` 裁掉 ≥1% 的误差点——两处均违反 AIPro 红线（虚假背书、裁剪坐标轴隐瞒发散点），已改为输出真实计算值且不设 ylim。

### 5. 2D 散射（A+BC 模型）的诚实现状

* **收敛的 CAP 参考解已建立（2026-11-01，服务器 xb002 实算）**：水深阱模型的 FD 参考解改为**散射波 formulation**——解 $\chi_{sc}=\chi-e^{+iKR}$，$(H-E)\chi_{sc}=-(H-E)\chi_{in}$（源项仅在相互作用区非零），$R_{hi}$ 处通道投影**出射** Robin；再加二次 CAP 吸收势 $-i\eta s^2$（$R>5$ Å，固定宽度 2 Å）消除盒共振。$\eta$ 由**平台判据**选取（内部解对盒尺寸 9→12 Å 与网格加密同时稳定），完整扫描存档于 `results/server_2d/`：$E_c=200$ 选 $\eta^*=100$（盒差 0.81%、网格差 0.46%、参考不确定度 1.54%）；$E_c=500$ 选 $\eta^*=200$（盒差 1.55%、网格差 0.34%、参考不确定度 0.55%）。$\eta$ 过大（≥3000）时 CAP 反射反而恶化——标准 CAP 行为，数据在 `cap_scan.json`。
* **无 CAP 盒式参考不收敛（定量证实，负结果）**：`run_scatter2d_convergence.py` 在 100×20→150×30→200×40 网格间内部 L2 差 128–156%，$R_{hi}$ 9→8/10 Å 内部解变化 83–136%。数据存档 `results/server_2d/convergence_2d.json`。
* **2D PINN 首次获得有意义的误差数字（诚实负结果）**：修复 `pinn_scatter_2d.py` 三个致命 bug（Fourier 特征列表在迭代中无限增长导致从未训练成功、复数张量直接 autograd、标量复数传 `torch.exp`）并统一入射约定后，`ScatteringPINN2D`（散射波 ansatz，硬 Dirichlet，PDE 残差+CAP）对收敛参考的误差为 **$E_c=200$: L2 = 122%；$E_c=500$: L2 = 156%**（内部区域 $R<4.9$ Å）。加大网络（4×128、64 Fourier 特征、10000 步）残差仍平台在 ~7×10⁻²（标度化）——是方法性欠拟合而非容量问题。**在残差-训练型 PINN 改进之前，2D 散射 PINN 不能作为定量求解器**；对比数据与脚本：`run_scatter_2d_cap_benchmark.py`、`results/server_2d/scatter2d_cap_benchmark.json`。
* **理论推导与线性配置验证（2026-11-05，`docs/theory_scattering_2d.md`）**：推导出**误差—残差换算律** $\delta\chi \gtrsim \varepsilon S/E_c$（ε 为标度化 PDE 残差）——强形式 PINN 残差平台 7×10⁻² 预测 O(1) 解误差，与实测 122% 定量一致；该律随后在随机 Fourier 特征线性配置（ELM，算符解析式经 2.75×10⁻⁹ 有限差分自检）上被**双向验证**：验证残差 6.5×10⁻²→1.9×10⁻² 时解误差同比例下降，L2 始终 O(1)。瓶颈是张成空间逼近能力：阱底局域波数 k≈36 Å⁻¹（波长 0.17 Å），随机基覆盖所需自由度比张量基高一个维度量级（M=8000、σ_R≤300 时 L2 仅降至 0.80–0.83，代价 25 倍）。**诚实结论：必须把 e^{ik(R,γ)R} 振荡吸收进载波（2D eikonal/Langer 载波 + 光滑神经修正，1D 已验证方法论的推广）——已记录为下一步**。数据：`results/server_2d/scatter2d_elm.json`。
* **2D 载波理论的首次数值实现（2026-11-23，`docs/theory_scattering_2d.md` §7–§8）**：推导出二维绝热通道 eikonal + 均匀 Langer-Airy 载波 $U^\pm_j$（零拟合参数，渐近纯入射/出射），并实现 `AdiabaticCarrier2D` + `CarrierPINN2D`（载波基 + 光滑修正特征的交叉验证岭回归）。**诚实结果**：纯载波 L2 = 4.06 (E_c=200) / 11.0 (E_c=500)，不可用；载波+300 光滑修正降到 **0.51 / 0.76**（相对强形式 PINN 的 1.22/1.56 改善约 2 倍），但离参考不确定度（1.5%/0.55%）仍差 1.5–2 个量级。误差预算诊断：载波通道径向残差 10–16%（网格 h·k≈1，4 阶差分无法识别 k≈25–36 Å⁻¹ 的解析载波）、绝热耦合残差 7–21%（通道函数在阱内变化剧烈）、Bi 支禁戒区指数增长。**结论：机制部分成立但实现未兑现理论优势**；可复现的下一步（细化载波网格/解析 Airy 导数、绝热基 Berry 规范固定、或采用 Johnson log-derivative CAP 金标准）已写入理论文档 §8.3。数据：`results/server_2d/scatter2d_carrier.json`。
* **离散色散修正的实测（2026-11-26，`docs/theory_scattering_2d.md` §9–§10）**：推导出 4 阶 FD 离散色散关系 $\kappa^2(q)=[30-32\cos qh+2\cos 2qh]/(12h^2)$ 与 $O((qh)^4)$ 相位失配估计（$k\approx25$–$36\ \text{Å}^{-1}$ 时累积 0.6–1.8 rad），并实现二分求根接入载波（`_q_disp`）。**实测结果与预测相反**：逐列残差**完全不变**（开放通道 0.157–0.172、半闭 0.571–0.719、深闭 0.914–1.339，连续相位与离散色散相位两列数值一致）——相位失配不是主导项，主导项是绝热耦合与深闭通道构造质量。另一结构性事实：载波列近似满足**齐次**方程而未知量满足非齐次方程 $L\chi_{sc}=src$，齐次列的组合无法拟合源，故纯载波验证残差恒在 0.98–0.99；真正起作用的是修正特征。修正后的路线（源响应载波/耦合通道基/Johnson CAP 金标准）已写入 §10。
* **源的响应列与「残差非精度代理」的实测边界（2026-11-26，§11）**：实现每通道响应列（TtT+ρI 正则化求解 (−B_R D2+ε_j−E)u_p=s_j），验证残差从 0.99 降到 0.61（结构上补齐源项），但**外部 L2 反而从 0.51 变差到 0.70**。据此确认：标度化残差换算律 $\delta\chi\approx\varepsilon S/E_c$ 只在接近收敛的邻域成立；在 O(1) 误差区不同基可把残差压到同一水平却对应不同的错解方向，**只有独立参考解的外部 L2 可信**（本仓库从不用自证残差宣称精度）。当前 2D 求解器诚实状态：线性配置家族 L2 平台 0.51–0.80，最佳为齐次载波+300 光滑特征；达到参考精度需耦合通道基或 Johnson log-derivative+CAP 金标准参考。
* **diabatic 载波两次实现（2026-11-28，`docs/theory_scattering_2d.md` §12–§13，流水线实算）**：推导平行输运变换 $\partial_R U=-AU$（对角 Berry 联络恒零、导数耦合精确消项）与精确 diabatic 通道方程 $\tilde V_{\rm eff}=U^T{\rm diag}(\varepsilon)U-B_RU^TA^TAU$。第一次实现（中心差分 A）**未通过自检**（反对称性 1.4×10²、U 正交性 1.8，数字作废）；改用 Hellmann–Feynman 型 $A_{jk}=\langle\Phi_j|\partial_R H_\gamma|\Phi_k\rangle/(\varepsilon_k-\varepsilon_j)$（反对称精确、$\|U^TU-I\|=8.7\times10^{-15}$、开通道列残差 0.28/0.22）后数学有效，但 **offdiag_ratio=1.000** —— diabatic 非对角势与对角同量级，对角载波分解前提不成立，L2=4.99/6.13。**结论（§13.3）**：四条表征路线（强形式/ELM/绝热载波/diabatic 载波）已系统穷尽，通道耦合必须进入基（close-coupling 基）或采用 Johnson log-derivative+CAP 金标准。数据：`results/server_2d/scatter2d_diabatic.json`。
* **BKMP2 H+H₂ 2D 散射（经典通道投影，无深阱）**：该体系无深阱。已修复 Robin 边界符号 bug（旧实现 `(d_R - iK)χ = -2iKe^{-iKR}` 在 e^{+iEt} 约定下施加的是**出射波幅度=1**（时间反演的无源 BVP），正确条件为 `(d_R + iK)χ = +2iKe^{+iKR}` 锚定单位入射振幅）。**已知局限（如实报告）**：参考解仍有盒尺寸依赖（R_hi 9→15 Å 内部解在全局相位对齐后差异 ~99.9%，范数 ∝ R_hi）——根因是 E_tot 处离散盒能级间距仅 ~2.6 cm⁻¹，(H−E_tot) 近奇异，直接稀疏解沿近零方向放大的幅度由舍入决定；标准修复需要 log-derivative（Johnson）传播或复吸收势。**在实现收敛参考之前，该体系上的 PINN-vs-参考精度声明同样不作**。当前 Level 2 回归检验仅验证通道投影与全维 FD 解的**同一 BVP 投影一致性**（0.0117%，时反演共轭关系下范数不变的内部一致性检查）。
* **J > 0 耦合信道（BKMP2）**：高斯求积势能矩阵元，S 矩阵幺正性残差 $\le 1.97 \times 10^{-6}$；BF/SF 坐标系等价性 $\|\chi_{\mathrm{SF}} - \chi_{\mathrm{BF}} \mathbf{U}\|_{\infty} \le 10^{-10}$。


## 项目地图与服务器流水线

* **目录结构**：`src/conventional`（FD/DVR/密切耦合/CAP/半经典载波参考解）、
  `src/pinn`（变分 PINN、1D 载波+交叉验证修正、2D 强形式/ELM/绝热载波/diabatic 载波）、
  `scripts/{validation,pinn,conventional,server}`、`docs/`（理论推导 + 结果索引）、
  `results/server_2d/`（服务器实算存档）。数据索引见 `docs/results_index.md`，
  理论推导全文见 `docs/theory_scattering_2d.md`（§1–§13）。
* **服务器流水线**：`bash scripts/server/run_pipeline.sh [SERVER]` —— 同步代码到
  实验室服务器（默认 xb002: 192.168.6.24，192 核 + A100-40GB），顺序执行
  回归套件与全部基准（1D/2D CAP/ELM/载波/diabatic/收敛研究/CAP 扫描），
  结果与日志自动回收到 `results/server_2d/`。**全部计算只跑在服务器上，
  本地不做任何数值计算**。
* **防造假约定（AIPro.md）**：每个数字有唯一可复算来源（`docs/results_index.md`）；模型选择
  （岭参数、修正接受）只使用 PDE 验证残差，从不接触参考解；负结果如实
  记录并写入理论文档对应章节。
