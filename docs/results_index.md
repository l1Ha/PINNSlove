# 结果索引（单一事实来源地图）

每个数字的可复算来源。JSON 为唯一权威数据；`results/server_2d/*.log` 是
对应的完整运行日志；`*.npz` 是逐点解数据（供绘图与外部复核）。

## 2D 散射（水深阱模型，H+OH 冻结键）

| 文件 | 内容 | 生成脚本 |
|---|---|---|
| `server_2d/convergence_2d.json` | 无 CAP 盒式参考不收敛的定量证据（网格差 128–156%、盒差 83–136%） | `scripts/conventional/run_scatter2d_convergence.py` |
| `server_2d/cap_scan.json` | CAP 强度平台扫描（η* 判据的完整数据） | `scripts/conventional/run_cap_scan.py` |
| `server_2d/scatter2d_cap_benchmark.json` | 收敛 CAP 参考 + 强形式 PINN 误差（122%/156%） | `scripts/pinn/run_scatter_2d_cap_benchmark.py` |
| `server_2d/scatter2d_elm.json` | 随机特征线性配置（ELM）7 配置扫描 | `scripts/pinn/run_scatter_2d_elm.py` |
| `server_2d/scatter2d_carrier.json` | 绝热载波 + 光滑修正基准 | `scripts/pinn/run_scatter_2d_carrier.py` |
| `server_2d/scatter2d_diabatic.json` | diabatic 载波基准（§12 变换的实算检验） | `scripts/pinn/run_scatter_2d_diabatic.py` |

## 全维度（束缚态与散射）

| 文件 | 内容 | 生成脚本 |
|---|---|---|
| `highest_dim_frontier_summary.json` | 9D 甲烷 / 12D 乙烯谱 + 6D H₂+H₂ 反应散射 120 点逐点误差 | `scripts/pinn/run_high_dim_all.py` |
| `3d_water/summary.json` | 3D 水变分 PINN vs DVR（含 n=1 失败的诚实记录） | `scripts/pinn/run_3d_water.py` |
| `full_regression_report.json` | 8 级回归套件报告 | `scripts/validation/run_full_regression.py` |
| `scatter2d_summary.json` / `_scatter2d_*.npz` | 历史无 CAP 对比（已被 CAP 基准取代，仅存档） | `scripts/pinn/run_scatter_2d.py`（兼容转发） |

## 1D

| 文件 | 内容 | 生成脚本 |
|---|---|---|
| `scatter_1d_*.csv/json`（脚本输出） | FD / 半经典 / 载波+交叉验证修正 三方法对比 | `scripts/pinn/run_scatter_1d.py` |

运行方式见 README「服务器流水线」一节：`bash scripts/server/run_pipeline.sh`。
