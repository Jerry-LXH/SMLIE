import matplotlib.pyplot as plt
import platform

def prob(data_results_clean, comp, fontsize=12):
    # 设置中文字体
    if platform.system() == 'Windows':
        plt.rcParams['font.family'] = 'SimHei'
    elif platform.system() == 'Darwin':
        plt.rcParams['font.family'] = 'Arial Unicode MS'
    else:
        plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['axes.unicode_minus'] = False

    # 构造频数表
    triples_counter_1 = (
        data_results_clean
        .groupby(['故障情况'])['故障方式']
        .value_counts()
        .unstack(fill_value=0)
    )
    triples_counter_2 = (
        data_results_clean
        .groupby(['故障情况'])['解决方案']
        .value_counts()
        .unstack(fill_value=0)
    )

    # 获取 comp 对应的频率分布
    try:
        result1 = triples_counter_1.loc[comp]
        result1 = result1[result1 > 0].sort_values(ascending=False)
        prob_result_1 = (result1 / result1.sum()).sort_values(ascending=False)
    except KeyError:
        print(f"'{comp}' 不在 triples_counter_1 中")
        prob_result_1 = None

    try:
        result2 = triples_counter_2.loc[comp]
        result2 = result2[result2 > 0].sort_values(ascending=False)
        prob_result_2 = (result2 / result2.sum()).sort_values(ascending=False)
    except KeyError:
        print(f"'{comp}' 不在 triples_counter_2 中")
        prob_result_2 = None

    # 画饼图（如果数据存在）
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    if prob_result_1 is not None:
        wedges, texts, autotexts = axes[0].pie(
            prob_result_1,
            labels=prob_result_1.index,
            autopct='%1.1f%%',
            startangle=140,
            textprops={'fontsize': fontsize}
        )
        axes[0].set_title(f"{comp} - 故障方式分布", fontsize=fontsize + 2)
        axes[0].axis('equal')
    else:
        axes[0].text(0.5, 0.5, "无数据", ha='center', va='center', fontsize=fontsize)
        axes[0].axis('off')

    if prob_result_2 is not None:
        wedges, texts, autotexts = axes[1].pie(
            prob_result_2,
            labels=prob_result_2.index,
            autopct='%1.1f%%',
            startangle=140,
            textprops={'fontsize': fontsize}
        )
        axes[1].set_title(f"{comp} - 解决方案分布", fontsize=fontsize + 2)
        axes[1].axis('equal')
    else:
        axes[1].text(0.5, 0.5, "无数据", ha='center', va='center', fontsize=fontsize)
        axes[1].axis('off')

    plt.suptitle(f"{comp} 的问题分布分析", fontsize=fontsize + 4)
    plt.tight_layout()
    plt.show()

    return fig, prob_result_1, prob_result_2