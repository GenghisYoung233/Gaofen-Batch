<img src="https://raw.githubusercontent.com/GenghisYoung233/Gaofen-Batch/main/assets/app_icon.png" alt="logo" align="left" height="100"/>

# Gaofen Batch

基于rasterio和Orfeo Toolbox开发的国产卫星影像预处理工具，前端由Electron驱动，支持批量处理。

### 特性:
- 😉 ~~GDAL~~ rasterio 驱动，省时高效，结果可靠。
- 😋 解压即可运行，可视化界面，小白友好。

### 支持的卫星传感器与功能

| 卫星传感器 | RPC正射校正 | 大气校正 | 融合 |
|------------|--------------|----------|------|
| CB04A_WPM | ✅ | ❌ | ✅ |
| GF1_PMS | ✅ | ✅ | ✅ |
| GF1_WFV | ✅ | ✅ | / |
| GF1B/C/D-PMS | ✅ | ✅ | ✅ |
| GF2_PMS | ✅ | ✅ | ✅ |
| GF4_PMI | ✅ | ❌ | ✅ |
| GF5_AHSI | ✅ | ❌ | / |
| GF5B_AHSI | ✅ | ❌ | / |
| GF5B_VIMI | ✅ | ❌ | / |
| GF6_PMS | ✅ | ✅ | ✅ |
| GF6_WFV | ✅ | ✅ | / |
| GF7_BWD | ✅ | ✅ | ✅ |
| GF7_DLC | ✅ | ✅ | ✅ |
| HJ2A_CCD | ✅ | ❌ | / |
| HJ2B_CCD | ✅ | ❌ | / |
| ZY1E_VNIC | ✅ | ❌ | / |
| ZY1F_AHSI | ✅ | ❌ | / |
| ZY303_TMS | ✅ | ✅ | ✅ |

> **注**: "✅" 表示支持，"❌" 表示不支持（缺少太阳辐照度数据）, "/" 表示不需要。大气校正 = 辐射定标 + 大气表观反射率计算 + 地表反射率计算，使用暗像元法。

## 如何使用

1. 从项目的[Release](https://github.com/GenghisYoung233/Gaofen-Batch/releases)下载程序，解压到本地，路径避免中文和空格。
2. 双击GaofenBatch.exe，添加待处理数据（原始.tar.gz压缩包）。
3. 点击“运行”按钮，选择输出文件夹，程序开始逐个解压并处理数据。

## 演示Gif（×10倍速）

<img src="/assets/GaofenBatch.gif" alt="demo" width="500"/>

## 从源码编译

### 环境要求

- **Python**: 3.14+ (推荐使用 [uv](https://github.com/astral-sh/uv) 管理 Python 版本)
- **Node.js**: 用于 Electron 前端打包
- **包管理**: 使用 [uv](https://github.com/astral-sh/uv) 管理 Python 依赖

### 开始编译

1. 克隆本仓库：
    ```bash
    git clone https://github.com/GenghisYoung233/Gaofen-Batch.git
    ```

2. 安装并打包Electron：
    ```bash
    npm install
    npm run pack
    ```

3. 安装Python依赖（使用 uv）：

    ```bash
    uv sync
    ```

4. 编译Python：

    ```bash
    pip install pyinstaller
    pyinstaller main.py
    ```

    编译完成后会在`dist`目录中生成`main`目录, 请`main`重命名为`bin`。

5. 配置Orfeo Toolbox

    下载[Orfeo Toolbox](https://www.orfeo-toolbox.org/download/)，并解压。

    并将解压好的`OTB-*-Win64`文件夹与本项目附带的`data`文件夹一起放入上一步的`bin`文件夹中。

6. 将`bin`文件夹放入打包好的Electron文件夹中

## 常见问题

1. 可用的OTB版本有哪些？

    目前仅使用 `OTB-8.1.2-Win64` 以及 `OTB-9.0.0-Win64` 进行了测试，如果在其他版本中出现了问题，请提交issue。

2. 我看不懂编译说明

    ~~请移步 `Github Action`查看具体执行流程 ，本项目的自动编译由Github Action驱动。~~ 还没调试好。
    
    **请不要在生产环境使用action内编译的版本，他们不一定是稳定的！**

3. 有关跨平台支持

    由于Electron的特性，本程序理论上支持跨平台运行，但尚未进行测试。

    可以自行安装对应平台的OTB并修改源码进行尝试。

    欢迎提交PR！

## TODO

- [ ] 解决action编译版无法使用的问题

- [ ] Linux支持

- [ ] 支持高分系列15米级数据正射校正