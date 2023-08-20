<img src="https://raw.githubusercontent.com/GenghisYoung233/Gaofen-Batch/main/assets/app_icon.png" alt="logo" align="left" height="100"/>

# Gaofen Batch

基于GDAL和Electron的国产卫星影像预处理工具

### 特性:
- 😉 GDAL驱动，省时高效，结果可靠。
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

如果你没有数据，可从[这里](https://gaofen-batch-r2.remotesensing.top)获取测试数据。

## 演示Gif（×10倍速）

<!-- <img src="/assets/GaofenBatch.gif" alt="demo" width="500"/> -->
<img src="/assets/GaofenBatch.gif" alt="demo" width="500" onerror="this.onerror=null;this.src='410445dd94fe2c91201cdc6ce03d7006.r2.cloudflarestorage.com/gaofen-batch/assets/GaofenBatch.gif';" />


## 从源码编译

1. 克隆本仓库：
    ```bash
    git clone https://github.com/GenghisYoung233/Gaofen-Batch.git
    ```
2. 安装并打包Electron：
    ```bash
    npm install
    npm run pack
    ```
3. 编译Python：
    ```bash
    pyinstaller main.py
    ```
    编译后的文件夹改名为'bin'。
4. 下载[Orfeo Toolbox](https://www.orfeo-toolbox.org/download/)，并将它与data文件夹一起放入'bin'中。
