import streamlit as st
import pandas as pd
from prob import prob

data_results_clean = pd.read_pickle("data_results_clean.pkl")

st.set_page_config(page_title="故障问题分析", layout="wide")
st.title("🔧 故障部件 - 问题分布分析")
st.markdown("请输入一个故障部件名称，查看其故障方式和解决方案的概率分布饼图。")

# 用户输入部件名称
comp = st.text_input("🔍 请输入部件名称：", value="门锁")

# 字体大小滑块
fontsize = st.slider("图表字体大小：", min_value=8, max_value=20, value=12)

# 按钮触发
if st.button("生成图表"):
    # 创建双饼图
    with st.spinner("正在绘图，请稍候..."):
        fig, prob_result_1, prob_result_2 = prob(data_results_clean, comp, fontsize=fontsize)
        st.pyplot(fig)
    
    # 还可选输出数据表格
    if prob_result_1 is not None:
        st.subheader("📊 故障方式概率分布")
        st.dataframe(prob_result_1.reset_index().rename(columns={0: "概率"}))
    if prob_result_2 is not None:
        st.subheader("📊 解决方案概率分布")
        st.dataframe(prob_result_2.reset_index().rename(columns={0: "概率"}))