# pyCFRAM 技术说明文档

## 1. 项目概述

pyCFRAM 是 CFRAM（Climate Feedback–Response Analysis Method，气候反馈-响应分析方法）的 Python 封装实现。该框架以 Fortran RRTMG 辐射传输模型为核心引擎，通过 Python 完成数据预处理、矩阵求解、并行调度和可视化，实现了对温度变化的定量归因分解。

当前版本（v0.01）已成功复现 Wu et al. (2025, *J. Climate*) 论文中 Fig.3（温度分解空间图）、Fig.4（PAP 柱状图）和 Fig.5（气溶胶分种类空间图）的结果。

### 1.1 CFRAM 方法原理

CFRAM 的核心思想是将观测到的温度变化 ΔT 分解为各物理过程的贡献：

```
ΔT = -(∂R/∂T)⁻¹ × [Σ ΔF_radiative + Σ ΔF_non-radiative]
```

其中：
- `(∂R/∂T)` 是 Planck 反馈矩阵（38×38，37 个大气层 + 地表），描述温度对辐射通量的偏导数
- `ΔF_radiative` 包括水汽、云、气溶胶、CO₂、O₃、太阳辐射、地表反照率、地表温度等 9 项辐射强迫
- `ΔF_non-radiative` 包括潜热通量、感热通量、地表过程等非辐射强迫
- 大气动力项通过残差法获得：`ΔT_atmdyn = ΔT_observed - Σ(其他所有项)`

### 1.2 参考文献

- Lu, J., and M. Cai, 2009: A new framework for isolating individual feedback processes in coupled general circulation climate models. *Climate Dynamics*.
- Cai, M., and J. Lu, 2009: A new framework for isolating individual feedback processes in coupled general circulation climate models. Part II. *Climate Dynamics*.
- Wu, Q., et al., 2025: Quantitative attribution of 2013 and 2022 extremely hot events in China. *Journal of Climate*, 38(17), 4331–4349.

---

## 2. 系统架构

### 2.1 Fortran-Python 分工

```
┌─────────────────────────────────────────────────────────┐
│ Fortran RRTMG（辐射引擎，cfram_rrtmg_1col）            │
│                                                         │
│  输入：单柱大气廓线 + 地表状态（base 和 perturbed）     │
│  计算：11 次 RRTMG 辐射传输调用                         │
│        + ~37 次 Planck 矩阵构建（逐层 +1K 扰动）       │
│  输出：9 项辐射 forcing (W/m²) + Planck 逆矩阵          │
└──────────────────────┬──────────────────────────────────┘
                       │ (临时二进制文件)
┌──────────────────────▼──────────────────────────────────┐
│ Python（分解求解 + 并行调度）                           │
│                                                         │
│  辐射项：  dT_i = -(∂R/∂T)⁻¹ × frc_i    (numpy 矩阵乘)│
│  非辐射项：dT_lhflx/shflx/sfcdyn         (同一 Planck)  │
│  气溶胶分种类：dT_bc/oc/sulf/seas/dust   (同一 Planck)  │
│  大气动力：dT_atmdyn = dT_obs - Σ(所有其他项)  (残差法) │
│  并行：multiprocessing Pool, 每格点独立进程              │
└─────────────────────────────────────────────────────────┘
```

### 2.2 为什么不用 OpenMP？

RRTMG 内部的 Fortran module 全局变量在线程间共享，导致数据竞争。尝试 OpenMP 并行后出现静默错误（无输出但不崩溃）。

采用 Python multiprocessing 方案：每个 worker 进程独立运行一个 Fortran 可执行文件（cfram_rrtmg_1col），通过进程隔离天然避免了 RRTMG 的线程安全问题。在 256 核的 hqlx204 上，80 进程并行，每个 case（9801 格点）约 19 分钟完成。

### 2.3 关键设计决策

| 决策 | 选择 | 理由 |
|------|------|------|
| dT 求解位置 | Python（非 Fortran） | 新增分解项不需要改 Fortran，只需提供 forcing |
| 并行方式 | multiprocessing（非 OpenMP） | RRTMG 线程不安全 |
| 非辐射项 | 用户提供 forcing，Planck 矩阵求解 | 保持通用性 |
| 气溶胶分种类 | per-species forcing × 自计算 Planck 矩阵 | 半自计算 |
| 大气动力项 | 残差法 | 无需额外数据 |
| 输入格式 | 标准 NetCDF（base + perturbed 两个状态） | 不预设 event/climatology |

---

## 3. 目录结构

```
pyCFRAM/
├── run_case.py              # 一键入口：python3 run_case.py eh13
├── requirements.txt         # Python 依赖
├── configs/
│   └── defaults.yaml        # 全局配置（辐射参数、气压层、气溶胶映射、出图设置）
├── cases/                   # 每个 case 一个子目录
│   ├── eh13/
│   │   ├── case.yaml        # case 配置（输入路径、出图区域）
│   │   ├── input/           # 输入数据（NetCDF，符合 input_spec）
│   │   ├── output/          # CFRAM 计算结果
│   │   └── figures/         # 出图结果
│   └── eh22/
├── core/                    # Python 核心模块
│   ├── config.py            # 统一配置加载器
│   ├── constants.py         # 物理常量 + RRTMG 固有参数
│   ├── radiation.py         # 辐射预处理（坐标变换、柱密度）
│   ├── aerosol_optics.py    # 气溶胶光学属性（混合比 → AOD/SSA/g）
│   ├── planck_matrix.py     # Planck 矩阵构建与求解
│   ├── decomposition.py     # CFRAM 分解流程编排
│   ├── cfram_runner.py      # Fortran 子进程调用封装
│   └── fortran_io.py        # Fortran 直接访问文件读写
├── fortran/                 # Fortran RRTMG 辐射引擎
│   ├── cfram_rrtmg.f90      # 主程序（只输出 forcing + Planck 矩阵）
│   ├── rad_driver.f90       # RRTMG 调用驱动
│   ├── rad_driver_lw.f90    # LW Planck 矩阵构建
│   ├── drdt.f90             # ∂R/∂T 计算
│   ├── math.f90             # 矩阵求逆
│   ├── writeout.f90         # 输出工具
│   ├── makefile             # 编译（同时生成 cfram_rrtmg 和 cfram_rrtmg_1col）
│   ├── data_prep/aerosol/   # GOCART 气溶胶查找表（opticsBands_*.nc）
│   ├── RRTMG_LW-master/     # 上游 RRTMG 长波库（含 rrtmg_lw.nc）
│   └── RRTMG_SW-master/     # 上游 RRTMG 短波库（含 rrtmg_sw.nc）
├── scripts/                 # 工作流脚本
│   ├── build_case_input.py         # ERA5 + MERRA-2 → 标准 NetCDF 输入（case.yaml 驱动）
│   ├── download_era5_flux.py       # ERA5 PL+SL 下载（CDS API）
│   ├── download_merra2_aerosol.py  # MERRA-2 M2I3NVAER 下载（NASA GES DISC）
│   ├── extract_full_field.py       # 输入 NetCDF → Fortran 二进制（含 GOCART 气溶胶光学）
│   ├── run_parallel_python.py      # 并行 CFRAM 计算（multiprocessing）
│   ├── plot_fig3_independent.py    # Fig.3 双列对比图（paper | independent）
│   ├── plot_fig3_self.py           # Fig.3 温度分解空间图（基于 paper_data 复跑）
│   ├── plot_fig3.py                # Fig.3 paper_data 直出（验证用）
│   ├── plot_fig4.py                # Fig.4 PAP 柱状图
│   ├── plot_fig5.py                # Fig.5 气溶胶分种类
│   └── validate_vs_paper.py        # 对比 surface dT 与 Wu et al. 结果
├── data/                    # 数据源模块
│   ├── source_base.py              # DataSource 抽象基类 + 工厂注册器
│   ├── era5_source.py              # ERA5DailySource：Wu et al. warm-period 方法全场版
│   └── merra2_aerosol.py           # MERRA-2 气溶胶加载器（垂直 log-p + 水平双线性插值）
├── plotting/                # 可视化模块
├── tests/                   # 单元测试
├── docs/                    # 文档
│   ├── algorithm_spec.md    # 算法规格说明
│   ├── input_spec.md        # 输入数据格式规范
│   └── technical_notes_zh.md # 本文档
└── validation/              # 验证输出
```

---

## 4. 安装与编译

### 4.1 环境要求

- Python 3.8+
- gfortran（用于编译 RRTMG）
- LAPACK/BLAS 库（通常系统自带）

### 4.2 安装步骤

```bash
git clone git@github.com:lzhenn/pyCFRAM.git
cd pyCFRAM
pip install -r requirements.txt

# 编译 Fortran（同时生成全场版和单柱版）
cd fortran
make
cd ..
```

`make` 会自动生成两个可执行文件：
- `cfram_rrtmg`：全场版（nlat=81, nlon=121），用于调试
- `cfram_rrtmg_1col`：单柱版（nlat=1, nlon=1），用于 Python 并行调用

---

## 5. 运行流程

### 5.1 使用 ERA5 + MERRA-2 独立驱动（标准流程）

#### Step 0：获取源数据

- **ERA5 PL+SL**：`scripts/download_era5_flux.py`（需 CDS 账号）；放到 `era5_data/daily/`
- **MERRA-2 气溶胶**：`scripts/download_merra2_aerosol.py`（需 NASA Earthdata）；放到 `era5_data/merra2/`

目录结构详见 README "Setup - Step 3"。

#### Step 1：生成标准 NetCDF 输入

```bash
python3 scripts/build_case_input.py --case eh13
python3 scripts/build_case_input.py --case eh22
```

`build_case_input.py` 从 `cases/<case>/case.yaml` 的 `source:` 段读取配置，生成：
- `cases/<case>/input/{base,perturbed}_{pres,surf}.nc` — state 变量
- `cases/<case>/input/nonrad_forcing.nc` — 非辐射 forcing（lhflx, shflx）

#### Step 2：一键运行

```bash
python3 run_case.py eh13
python3 run_case.py eh22
```

`run_case.py` 依次执行三步：
1. **extract**：读取输入 NetCDF，计算气溶胶光学属性，写入 Fortran 二进制
2. **run**：multiprocessing 并行运行 CFRAM（默认用满所有 CPU）
3. **plot**：生成温度分解空间图

也可分步执行：

```bash
python3 run_case.py eh13 --step extract
python3 run_case.py eh13 --step run --nproc 40
python3 run_case.py eh13 --step plot
```

#### Step 3：查看结果

- 计算结果：`cases/eh13/output/cfram_result.nc`
- 图形结果：`cases/eh13/figures/fig3_decomposition.png`

额外出图：

```bash
python3 scripts/plot_fig5.py --case eh13 eh22   # Fig.5 气溶胶分种类
```

### 5.2 添加新 case

1. 创建目录和配置：

```bash
mkdir -p cases/my_case/input
```

2. 编写 `cases/my_case/case.yaml`：

```yaml
case_name: MY_CASE
description: "My analysis: scenario A vs scenario B"

input:
  base_pres: input/base_pres.nc
  base_surf: input/base_surf.nc
  perturbed_pres: input/perturbed_pres.nc
  perturbed_surf: input/perturbed_surf.nc
  # nonrad_forcing: input/nonrad_forcing.nc  # 可选

plot:
  key_region:
    lon: [100, 120]
    lat: [25, 40]
```

3. 准备输入 NetCDF 文件（格式见 `docs/input_spec.md`）

4. 运行：

```bash
python3 run_case.py my_case
```

---

## 6. 输入数据格式

详见 [input_spec.md](input_spec.md)。核心要求：

- **base_pres.nc / perturbed_pres.nc**：大气廓线（37层），包含温度、比湿、臭氧、云量/水含量、CO₂、6种气溶胶
- **base_surf.nc / perturbed_surf.nc**：地表状态（温度、气压、太阳辐射、反照率）
- **nonrad_forcing.nc**（可选）：非辐射 forcing（潜热、感热、地表过程），以及气溶胶分种类 forcing（bc, oc, sulf, seas, dust）

pyCFRAM 不限定两个状态的定义方式——可以是事件期 vs 气候态、未来 vs 历史、敏感性实验 vs 控制实验等。

---

## 7. 配置系统

### 7.1 全局配置 `configs/defaults.yaml`

包含所有与 case 无关的参数：

| 配置项 | 说明 |
|--------|------|
| `radiation.scon` | 太阳常数 (1360.98 W/m²) |
| `radiation.co2_ppmv` | CO₂ 默认浓度 (395 ppmv，输入数据有则覆盖) |
| `radiation.ch4_ppmv` | CH₄ 浓度 (1.6 ppmv) |
| `radiation.n2o_ppmv` | N₂O 浓度 (0.28 ppmv) |
| `radiation.cloud_re_ice/liq` | 云粒子有效半径 (5/20 microns) |
| `grid.pressure_levels` | 37 层标准气压层 (hPa) |
| `aerosol.species` | 气溶胶种类 → GOCART 查找表映射 |
| `run.nproc` | 并行进程数 (auto = 所有 CPU) |
| `plotting.levels_dt` | 温度分解图色阶 |

### 7.2 Case 配置 `cases/<name>/case.yaml`

每个 case 独立的配置：输入文件路径、出图参数等。

### 7.3 配置加载

所有脚本通过 `core/config.py` 统一加载配置：

```python
from core.config import load_case, defaults, get_plev, get_aerosol_map
```

---

## 8. 输出格式

### 8.1 cfram_result.nc

每个 case 的 CFRAM 计算结果保存为 NetCDF：

| 变量 | 维度 | 说明 |
|------|------|------|
| `dT_q` | (lev, lat, lon) | 水汽贡献的温度变化 |
| `dT_cloud` | (lev, lat, lon) | 云贡献 |
| `dT_aerosol` | (lev, lat, lon) | 气溶胶总贡献（RRTMG one-at-a-time） |
| `dT_co2` | (lev, lat, lon) | CO₂ 贡献 |
| `dT_o3` | (lev, lat, lon) | 臭氧贡献 |
| `dT_ts` | (lev, lat, lon) | 地表温度辐射效应 |
| `dT_solar` | (lev, lat, lon) | 太阳辐射贡献 |
| `dT_albedo` | (lev, lat, lon) | 地表反照率贡献 |
| `dT_warm` | (lev, lat, lon) | 总辐射变化 |
| `dT_lhflx` | (lev, lat, lon) | 潜热通量贡献（需 nonrad_forcing） |
| `dT_shflx` | (lev, lat, lon) | 感热通量贡献 |
| `dT_sfcdyn` | (lev, lat, lon) | 地表过程贡献 |
| `dT_bc/oc/sulf/seas/dust` | (lev, lat, lon) | 气溶胶分种类贡献（需 nonrad_forcing） |
| `dT_observed` | (lev, lat, lon) | 观测温度变化 (T_perturbed - T_base) |
| `dT_atmdyn` | (lev, lat, lon) | 大气动力贡献（残差） |
| `frc_*` | (lev, lat, lon) | 对应的辐射 forcing (W/m²) |

### 8.2 level 维度

38 层 = 37 个大气层（TOA→地表）+ 1 个地表层。地表层为最后一个 index。

---

## 9. 验证结果

### 9.1 单柱验证 (115°E, 32°N)

| 项目 | Fortran 输出 | paper_data | 差值 |
|------|-------------|------------|------|
| Total dT (surface) | 4.650 K | 4.659 K | 0.009 K |
| WV dT (surface) | -3.181 K | -3.143 K | 0.038 K |

### 9.2 全场验证（地表层空间相关系数）

| 分解项 | EH13 | EH22 | 说明 |
|--------|------|------|------|
| Water Vapor | 0.997 | 0.998 | RRTMG 自计算 |
| Cloud | 0.945 | 0.956 | RRTMG 自计算 |
| Aerosol (total) | 0.885 | 0.842 | RRTMG 自计算 |
| Surface Process | 0.999 | 0.998 | Planck 矩阵 × paper forcing |
| Atm. Dynamics | 0.621 | 0.514 | 残差法 |
| **Total (observed)** | **0.981** | **0.981** | ERA5 观测 |

### 9.3 已知差异原因

1. **Cloud/Aerosol 互补偏差**：气溶胶从 6 种聚合物种近似映射到单一 GOCART bin，导致气溶胶偏大、云偏小，但两者之和一致
2. **大气动力项精度较低**：因 cloud/aerosol 偏差传播到残差
3. **论文的 `dyn` = `atmdyn + sfcdyn + lhflx + shflx`**（全部非辐射项之和）

---

## 10. 关键代码说明

### 10.1 Fortran 辐射引擎 (`fortran/cfram_rrtmg.f90`)

主程序流程：

```fortran
! 对每个格点 (ilat, ilon)：
! 1. 确定有效大气层数 nlayer (基于地表气压)
! 2. 调用 rad_driver × 11 次：
!    - base 态 + warm 态 → 总变化
!    - 逐项替换 (CO₂, q, ts, O₃, solar, albedo, cloud, aerosol) → 各项 forcing
! 3. 调用 calc_drdt → 构建 Planck 矩阵 (∂R/∂T)
! 4. 矩阵求逆 → drdt_inv
! 5. 输出 9 项 forcing + drdt_inv
```

### 10.2 Python 并行运行器 (`scripts/run_parallel_python.py`)

```python
def process_column(args):
    ilat, ilon = args
    # 1. 提取单柱数据
    # 2. 计算气溶胶光学属性
    # 3. 写临时二进制文件到 tmpdir
    # 4. 运行 cfram_rrtmg_1col
    # 5. 读取 forcing + drdt_inv
    # 6. dT = -drdt_inv @ frc  (所有项统一)
    # 7. 清理 tmpdir，返回结果

with Pool(nproc) as pool:
    for ilat, ilon, result in pool.imap_unordered(process_column, tasks):
        # 收集结果到全场数组
```

### 10.3 气溶胶光学属性计算 (`core/aerosol_optics.py`)

从 MERRA-2 气溶胶混合比 → RRTMG 波段分辨的 AOD/SSA/g：

```
mixing_ratio × air_density × layer_thickness × mass_extinction_coeff → AOD per band
```

查找表来自 GOCART（opticsBands_*.nc），包含 30 个波段（14 SW + 16 LW）的消光系数，依赖相对湿度。

当前限制：paper_data 提供 6 种聚合物种（bc, ocphi, ocpho, sulf, ss, dust），但 GOCART 有 15 个 bins（DU001-005, SS001-005 等）。每种聚合物种映射到单一代表性 bin，这是气溶胶精度的主要瓶颈。

---

## 11. 性能指标

| 指标 | 值 |
|------|------|
| 格点数 | 81 × 121 = 9,801 |
| 并行度 | 80 进程 (hqlx204, 256 核) |
| 单 case 耗时 | ~19 分钟 |
| 单格点耗时 | ~1 秒 |
| 吞吐率 | 8.6 pts/s |
| 内存 | ~2 GB/进程（主进程加载全场数据 ~500 MB） |

---

## 12. 后续开发方向

1. **气溶胶 15-bin 精确映射**：替代当前 6→1 bin 近似，改善 cloud/aerosol 分配精度
2. **ERA5 数据自动预处理**：从原始 ERA5 日数据计算 base/perturbed 态，完全脱离 paper_data
3. **f2py 封装 RRTMG**：避免临时文件 I/O，提升单格点性能
4. **支持不同气压层数**：通过 Fortran namelist 运行时指定（当前固定 37 层）
5. **DAMIP/CMIP6 接口**：支持模式输出作为输入
6. **Energy Gain Kernel 框架集成**
