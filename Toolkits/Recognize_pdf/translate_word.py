from googletrans import Translator

def translate_word_into_chinese(word):
    # 创建翻译对象
    translator = Translator()
    # 翻译为中文
    translated = translator.translate(word, src='en', dest='zh-cn')
    # 输出翻译后的中文字符串
    return translated.text


def save_translated_words(input_file, output_file): # 保存词典
    with open(input_file, 'r') as file:
        words = file.read()
        word_list=words.splitlines()
    with open(output_file, 'w') as file:
        chns=translate_word_into_chinese(words)
        chns_list= chns.splitlines()
        i=0
        for (chn) in chns_list:
            file.write(f"{word_list[i]}:{chn}\n") 
            i+=1

save_translated_words('Toolkits/Recognize_pdf/Results/test_output.txt','Toolkits/Recognize_pdf/Results/translated.txt')