import PyPDF2
import nltk
import string
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize
from collections import Counter


# First use:
# nltk.download('punkt')
# nltk.download('stopwords')

def extract_text_from_pdf(pdf_file):
    with open(pdf_file, 'rb') as file: # 读取二进制文件使用“rb”
        reader = PyPDF2.PdfReader(file) # 建立了一个针对file的reader器（实例）
        text = ""
        # print(reader.pages[0])
        for page_num in range(len(reader.pages)): # reader的pages属性是以一页为一个元素的数组
            page = reader.pages[page_num] # 将所有数据读入一个字符串！
            text += page.extract_text()
    return text

def tokenize_and_clean(text):
    # 分词
    words = word_tokenize(text.lower()) # 小写+分成列表

    # 清理数据：去除标点符号和停用词
    stop_words = set(stopwords.words('english')) # 英语的停用列表的集合
    cleaned_words = [word for word in words if word.isalpha() and word not in stop_words]
    
    return cleaned_words #返回一个列表

def get_word_frequencies(words):
    return Counter(words) # 统计词频为一个词典

def filter_new_words(word_freq, known_words_file): # word_freq 是一个词典

    # 从已知词汇文件known_words_file中加载已知词汇
    with open(known_words_file, 'r') as file:
        known_words = set(word.strip().lower() for word in file.readlines()) # 单词为元素的字符串列表

    # 找到不在【已知词汇】表中的词汇
    new_words = {word: freq for word, freq in word_freq.items() if ((word not in known_words) and (len(word)>5))} # 返回一个词典

    new_words_sorted=sorted(new_words.items(), key=lambda item: item[1],reverse=True)

    # 过滤掉重复的单词
    for new_word in new_words_sorted: 
        new_word_s_list=new_word[0]+"s"
    for new_word in new_words_sorted:
        if new_word[0] in new_word_s_list:
            new_words_sorted.remove(new_word)

    return new_words_sorted # 返回排序的元组列表

def save_word_frequencies(word_freq, output_file): # 保存词典
    with open(output_file, 'w') as file:
        for item in word_freq:
            if item[1]>4:
                file.write(f"{item[0]}\n")      
    
def main(pdf_file, known_words_file, output_file):
    # 从PDF提取文本
    text = extract_text_from_pdf(pdf_file)
    
    # 分词和清理
    words = tokenize_and_clean(text)
    
    # 统计词频
    word_freq = get_word_frequencies(words)
    
    # 过滤生词
    new_words = filter_new_words(word_freq, known_words_file)
    
    # 保存生词频率
    save_word_frequencies(new_words, output_file)

main("Toolkits/Recognize_pdf/pdf/Tacikowski et al. - 2024 - Human hippocampal and entorhinal neurons encode the temporal structure of experience.pdf","Toolkits/Recognize_pdf/Data/freq_5000.txt","Toolkits/Recognize_pdf/Results/test_output.txt")