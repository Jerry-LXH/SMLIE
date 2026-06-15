import nltk
from nltk.corpus import brown

# 下载语料库
nltk.download('brown')

# 获取 Brown 语料库中的单词
words = brown.words()

# 使用 FreqDist 统计词频
freq_dist = nltk.FreqDist(words)

# 获取前 50 个高频词汇
most_common_words = freq_dist.most_common(500) # 返回一个元组列表

# 保存
with open('Toolkits/Recognize_pdf/Data/freq_500.txt', 'w') as file:
    for word in most_common_words:
        file.write(f"{word[0]}\n")