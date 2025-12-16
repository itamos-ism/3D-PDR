import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
from matplotlib import gridspec
from scipy.ndimage import gaussian_filter1d
import sys

def read_pdr_columns(chem_path, file_path, columns=['Tgas', 'rho', 'CO', 'CI', 'CII', 'HI', 'HII']):
    """
    从PDR输出文件中读取指定列的数据
    
    参数:
        chem_path: network文件路径
        file_path: fin文件路径
        columns: 要读取的列名列表，可选['Tgas', 'rho', 'abundance']
    
    返回:
        dict: 包含指定列数据的字典
    """

    # 读取化学物种文件建立列名映射
    species_names = {}
    with open(chem_path, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 2:
                species = parts[1].strip()
                species_id = int(parts[0]) + 8
                species_names[species] = species_id

    other_names = {
        'p':0, 'x':1, 'y':2, 'z':3,
        'Tgas':4, 'Tdust':5, 'etype':6,
        'rho':7, 'UV':8
    }

    # 创建列名到列索引的映射
    col_map = {**other_names, **species_names}
    
    # 读取所有数据
    data = np.loadtxt(file_path)
    
    # 提取指定列
    result = {}
    for col in columns:
        if col in col_map:
            result[col] = data[:, col_map[col]]
        else:
            raise ValueError(f"无效列名: {col}")
    
    return result

import numpy as np

def plot_bins(ax, rho, values, species, scatter_color, line_color, ls, alpha, label=None, num_bins=100, sigma=1.0, scatter=True, marker='o'):
    """
    对单个 level 的数据，按密度分箱，绘制物种丰度的均值线与范围阴影。

    Parameters:
        ax: matplotlib axis
        rho: array-like, 氢密度 (nH)
        values: array-like, 物种相对丰度
        species: str, 物种名称（用于调试）
        color: str, 线条/阴影颜色
        ls: str, 线型
        alpha: float, 阴影透明度
        label: str, 图例标签
        num_bins: int, 分箱数量
        sigma: float, 高斯平滑参数
    """
    # 过滤有效数据
    valid = (rho > 0) & (values > 0)
    if not np.any(valid):
        return
    
    # === 绘制原始散点图（可选）===
    if scatter:
        ax.scatter(rho[valid], values[valid], edgecolors=scatter_color, facecolors='none', s=20, alpha=alpha, marker=marker, rasterized=True, linewidths=1.2)    
    
    log_rho = np.log10(rho[valid])
    values_valid = values[valid]

    # 创建分箱
    bins = np.linspace(log_rho.min(), log_rho.max(), num_bins)
    bin_indices = np.digitize(log_rho, bins)

    mean_vals, centers = [], []

    for i in range(1, len(bins)):
        mask = (bin_indices == i)
        if np.sum(mask) == 0:
            continue
        subset = values_valid[mask]
        mean_vals.append(np.mean(subset))
        # 计算对数空间 bin 中心，转回线性
        center_log = 0.5 * (bins[i-1] + bins[i])
        centers.append(10**center_log)

    if len(centers) == 0:
        return

    # 平滑
    mean_smooth = gaussian_filter1d(mean_vals, sigma)

    # 绘图
    plot_kwargs = {
        'color': line_color,
        'linestyle': ls,
        'linewidth': 6
    }
    if label is not None:
        plot_kwargs['label'] = label

    ax.plot(centers, mean_smooth, **plot_kwargs)

# 输入参数
chem_file = "species_reduced.d"
level0_file = "level0/diffuse_32_base.pdr.fin"
level1_file = "level1/test.pdr.fin"
col = 2 # column number of the figure
row = 3 # row number of the figure
num_bins = 100
sigma = 1.0 # 平滑强度参数（可根据需要调整）
cnlev_list = [41,5,5,5]
ilevel = 1

# 设置全局样式
# plt.style.use('default')
plt.rcParams.update({
    'font.size': 48,
    'axes.labelsize': 32,
    'xtick.labelsize': 32,
    'ytick.labelsize': 32,
    'legend.fontsize': 32,
    'axes.linewidth': 2.0
})

# 读取数据
l0_data = read_pdr_columns(chem_file,level0_file,
                           columns = ['rho', 'CO', 'CI', 'CII', 'HI', 'H2', 'Tgas'])
l1_data = read_pdr_columns(chem_file,level1_file,
                           columns = ['rho', 'CO', 'CI', 'CII', 'HI', 'H2', 'Tgas'])

# 创建画布和子图布局
fig = plt.figure(figsize=(col*16, row*8))
gs = gridspec.GridSpec(row, col, height_ratios=[1]*row, hspace=0.05)
ax1_list = []
ax2_list = []

# 设置画的物种数
plot_config = [
    {'species_1': 'CII', 'species_2': 'HI'},
    {'species_1': 'CI', 'species_2': 'H2'},
    {'species_1': 'CO', 'species_2': 'H2'}
]

# 设置颜色和样式配置
level_styles = {
    'level0': {'scatter_color': '#2E86AB', 'line_color': '#1A4D69', 'ls': '--', 'alpha': 0.2, 'marker': 'o', 'label': r'$Base$'},
    'level1': {'scatter_color': '#A23B72', 'line_color': '#5F253B', 'ls': '--', 'alpha': 0.2, 'marker': '^', 'label': r'$DDA$'}
}

for i, config in enumerate(plot_config):
    species_1 = config['species_1']
    species_2 = config['species_2']

    # --------------------------
    # 左列：CI, CII, CO
    # --------------------------
    # 创建子图（共享x轴）
    ax1 = fig.add_subplot(gs[i,0]) if i == 0 else fig.add_subplot(gs[i,0], sharex=ax1_list[0])
    ax1_list.append(ax1)

    # 处理所有level数据
    for level, data in zip(['level0', 'level1'], [l0_data, l1_data]):
        style = level_styles[level]
        plot_bins(
            ax=ax1,
            rho=data['rho'],
            values=data[species_1],
            species=species_1,
            scatter_color=style['scatter_color'],
            line_color=style['line_color'],
            ls=style['ls'],
            alpha=style['alpha'],
            label=style['label'],
            num_bins=num_bins,
            sigma=sigma,
            marker=style['marker']
        )

    # 设置子图属性
    ax1.set_ylabel(f'Relative abundance of {species_1}', labelpad=15)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_xlim(10,100000)
    # 设置图例位置（根据物种动态调整）
    legend_loc = 'lower left' if config['species_1'] == 'CII' else 'lower right'
    ax1.legend(loc=legend_loc, frameon=False, ncol=1)
    
    # 只在最下方子图显示x轴标签
    if i < row - 1:
        plt.setp(ax1.get_xticklabels(), visible=False)
    else:
        ax1.set_xlabel(r'$n_\mathrm{H}$ [cm$^{-3}$]')

    # --------------------------
    # 右列：HI, HII
    # --------------------------
    ax2 = fig.add_subplot(gs[i, 1]) if i == 0 else fig.add_subplot(gs[i, 1], sharex=ax1_list[0])
    ax2_list.append(ax2)

    for level, data in zip(['level0', 'level1'], [l0_data, l1_data]):
        style = level_styles[level]
        plot_bins(
            ax=ax2,
            rho=data['rho'],
            values=data[species_2],
            species=species_2,
            scatter_color=style['scatter_color'],
            line_color=style['line_color'],
            ls=style['ls'],
            alpha=style['alpha'],
            label=style['label'],
            num_bins=num_bins,
            sigma=sigma,
            marker=style['marker']
        )

    if config['species_2'] == 'H2':
        ax2.set_ylabel(r'Relative abundance of $\mathrm{H_2}$', labelpad=15)
        ax2.legend(loc='lower right', frameon=False, ncol=1)
    elif config['species_2'] == 'HI':
        ax2.set_ylabel(f'Relative abundance of {species_2}', labelpad=15)
        ax2.legend(loc='lower left', frameon=False, ncol=1)  # HI的图例放左下角
    else:
        ax2.set_ylabel(f'Relative abundance of {species_2}', labelpad=15)
        ax2.legend(loc='upper right', frameon=False, ncol=1)
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_xlim(10,100000)
    if i < row - 1:
        plt.setp(ax2.get_xticklabels(), visible=False)
    else:
        ax2.set_xlabel(r'$n_\mathrm{H}$ [cm$^{-3}$]') 

# ========================
# 单独绘制 Tgas 到最后一格
# =======================
# 移除最后一个子图（右下角），用于绘制 Tgas
ax2_list[-1].remove()

ax_tgas = fig.add_subplot(gs[row-1, col-1])
for level, data in zip(['level0', 'level1'], [l0_data, l1_data]):
    style = level_styles[level]
    rho = data['rho']
    Tgas = data['Tgas']

    plot_bins(
        ax=ax_tgas,
        rho=rho,
        values=Tgas,
        species='Tgas',
        scatter_color=style['scatter_color'],
        line_color=style['line_color'],
        ls=style['ls'],
        alpha=style['alpha'],
        label=style['label'],
        num_bins=num_bins,
        sigma=sigma,
        marker=style['marker']
    )

ax_tgas.set_xlabel(r'$n_\mathrm{H}$ [cm$^{-3}$]')
ax_tgas.set_ylabel('Gas Temperature [K]')
ax_tgas.set_xscale('log')
ax_tgas.set_yscale('log')
ax_tgas.set_xlim(10,100000)
ax_tgas.legend(loc='upper right', frameon=False, ncol=1)    

plt.tight_layout()
plt.savefig('abundances.png', dpi=300, bbox_inches='tight')

