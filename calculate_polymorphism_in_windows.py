#!/usr/bin/env python3
# coding=utf-8

"""
For details, see https://github.com/shelkmike/Polymorphism_suite

Notes:
1) All comments, except this one, are in Russian. Sorry, but it's somewhat easier for me to write in Russian than in English. To understand some comment, you can use Google Translate. Names of variables are usually self-explanatory, so it is often possible to understand the meaning of a piece of code without comments. In case of trouble understanding code, ask a question at https://github.com/shelkmike/Polymorphism_suite/issues .
2) Throughout the code I use a Hungarian notation, which means I denote meaning of a word by using special prefixes. In particular:
s_ - variables that contain strings
n_ - variables that contain numbers
l_ - lists
d_ - dictionaries
f_ - file handlers
o_ - more complex objects
Nested data structures are denoted by several letters. For example, dl_ are dictionaries of lists and ll_ are lists of lists.
"""

import sys
import os
import re
import natsort
import statistics

f_infile = open(sys.argv[1], "r")
n_window_size = int(sys.argv[2])
n_window_step = int(sys.argv[3])
s_alignment_circular_or_not_circular = sys.argv[4]
s_position_to_start = sys.argv[5]
n_number_of_cpu_threads_to_use = int(sys.argv[6])
s_path_to_the_output_folder = sys.argv[7]

#Определяю путь к текущей папке. Нужен, чтобы воспользоваться скриптом calculate_polymorphism_using_entire_alignment.py , который необходим для работы данного скрипта
s_path_to_the_folder_where_this_script_lies = os.path.dirname(os.path.realpath( __file__ ))#Использую os.path.realpath, чтобы если скрипт запускается по мягкой ссылке на исполняемый файл, всё равно удалось обнаружить папку со вторым скриптом.

#Проверяю, что выходной папки не существует, либо она существует и пустая. В противном случае, говорю пользователю, что это ошибка.
if os.path.exists(s_path_to_the_output_folder):
	if len(os.listdir(s_path_to_the_output_folder)) != 0:
		print("Error: the output folder already exists and is not empty.")
		sys.exit()
else:
	os.mkdir(s_path_to_the_output_folder)

#f_log = open("log.txt", "w")

#Делаю три словаря, в которых ключи это позиции центра окна, а значения это, соответственно, MPPI, MPMLD или pi. Если на вход скрипту было дано не множественное выравнивание, а парное, то словарь про MPMLD не заполняется (см. заголовочный комментарий).
d_position_of_the_window_center__to__MPPI = {}
d_position_of_the_window_center__to__MPMLD = {}
d_position_of_the_window_center__to__pi = {}

"""
Величины n_window_left_spread и n_window_right_spread равны тому, сколько нуклеотидов налево и направо нужно брать от данного (не включая его), чтобы получить окно. 
Например, если размер окна 100, то n_window_left_spread будет 50, n_window_right_spread будет 49. 
А если размер окна 99, то n_window_left_spread будет 49 и n_window_right_spread будет 49. 
"""
n_window_left_spread = int(n_window_size / 2)
n_window_right_spread = n_window_size - n_window_left_spread - 1

#загружаю все последовательности в словарь. Ключи в словаре это заголовки последовательностей без ">" (но только часть до первого пробельного символа), а значение это сами последовательности. Беру только часть заголовка до первого пробельного символа, чтобы у IQ-TREE не было проблем.
d_sequence_title_to_sequence = {}
s_sequence_title = "" #заголовок последовательности, на которую я сейчас смотрю.

for s_line in f_infile:
	if re.search(r"^>(\S+)", s_line):
		o_regular_expression_results = re.search(r"^>(\S+)", s_line)
		s_sequence_title = o_regular_expression_results.group(1)
		#удаляю символ переноса строки
		s_sequence_title.rstrip('\n')
	elif re.search(r"^(.+)", s_line): #если это не строка с заголовком, то считаю, что строка с последовательностью
		o_regular_expression_results = re.search(r"^(.+)", s_line)
		#если для этого контига последовательность ещё не инициализирована, то инициализирую её
		if s_sequence_title not in d_sequence_title_to_sequence:
			d_sequence_title_to_sequence[s_sequence_title] = ""
		
		s_sequence_from_this_string = o_regular_expression_results.group(1)
		#удаляю всякие пробельные символы, в том числе символ переноса строки.
		s_sequence_from_this_string = re.sub(r"\s","",s_sequence_from_this_string)

		d_sequence_title_to_sequence[s_sequence_title] += s_sequence_from_this_string

#Длину выравнивания считаю просто по последней последовательности, которую видел. Всё равно длина всех последовательностей во множественном выравнивании одинаковая.
n_alignment_length = len(d_sequence_title_to_sequence[s_sequence_title])

#print("For " + sys.argv[1] + " - " + str(n_alignment_length) + "\n")

#Если длина выравнивания меньше, чем размер окна, то выхожу с ошибкой.
if n_alignment_length < n_window_size:
	sys.exit("Error: alignment length is smaller than the window size\n")

#Смотрю в соответствие с использованной опцией, должен ли быть центр первого окна в позиции 1 или в такой позиции, чтобы левый край окна был в позиции 1. Переменная n_position_of_the_window_center характеризует центр окна; эта переменная считается one-based.
if s_position_to_start == "start_from_1":
	n_position_of_the_window_center = 1
elif s_position_to_start == "start_from_window_center":
	n_position_of_the_window_center = (1 + int(n_window_size / 2))

f_table_with_polymorphism_in_windows = open(s_path_to_the_output_folder + "/polymorphism_in_windows.csv", "w", buffering = 1)
#Если на вход было дано парное выравнивание, то столбец про MPMLD не нужен (см. заголовочный комментарий)
if len(d_sequence_title_to_sequence) > 2:
	f_table_with_polymorphism_in_windows.write("Position;MPPI;MPMLD;pi; ;\n")
else:
	f_table_with_polymorphism_in_windows.write("Position;MPPI;pi; ;\n")

while n_position_of_the_window_center <= n_alignment_length:
	
	#Теперь мне нужно сделать словарь, в котором ключ это заголовок последовательности, а значение это последовательность конкретно этого окна.
	d_sequence_title_to_sequence_of_the_window = {}
	
	"""
	В следующих трёх словарях ключ это название последовательности, а значения части последовательности окна, которые, соответственно:
	1) Выпирают за левый край и потому должны быть взяты с правого края. Значения этого словаря будут непустые только для кольцевого выравнивания.
	2) Не выпирают ни за левый, ни за правый край.
	3) Выпирают за правый край и потому должны быть взяты с левого края. Значения этого словаря будут непустые только для кольцевого выравнивания.
	
	Значения словаря d_sequence_title_to_sequence_of_the_window будут получаться конкатенацией значений этих трёх словарей.
	"""
	
	d_sequence_title_to_sequence_of_the_left_overhang_of_the_window = {}
	d_sequence_title_to_sequence_of_the_non_overhanging_part_of_the_window = {}
	d_sequence_title_to_sequence_of_the_right_overhang_of_the_window = {}
	
	
	for s_sequence_title in d_sequence_title_to_sequence:
		
		d_sequence_title_to_sequence_of_the_left_overhang_of_the_window[s_sequence_title] = ""
		d_sequence_title_to_sequence_of_the_non_overhanging_part_of_the_window[s_sequence_title] = ""
		d_sequence_title_to_sequence_of_the_right_overhang_of_the_window[s_sequence_title] = ""
		
		#Если выравнивание кольцевое, а окно расположено так, что оно высовывается за левый край.
		if (s_alignment_circular_or_not_circular == "circular") and (n_position_of_the_window_center - n_window_left_spread < 1):
			d_sequence_title_to_sequence_of_the_left_overhang_of_the_window[s_sequence_title] = d_sequence_title_to_sequence[s_sequence_title][((n_alignment_length - (n_window_left_spread - n_position_of_the_window_center)) - 1):]
		
		#Если выравнивание кольцевое, а окно расположено так, что оно высовывается за правый край.
		if (s_alignment_circular_or_not_circular == "circular") and (n_position_of_the_window_center + n_window_right_spread > n_alignment_length):
			d_sequence_title_to_sequence_of_the_right_overhang_of_the_window[s_sequence_title] = d_sequence_title_to_sequence[s_sequence_title][:(n_position_of_the_window_center + n_window_right_spread - n_alignment_length)]
		
		#Та часть выравнивания, которая не выпирает ни за левый, и на за правый край. max и min тут используются, чтобы учесть возможное выпирание за край.
		d_sequence_title_to_sequence_of_the_non_overhanging_part_of_the_window[s_sequence_title] = d_sequence_title_to_sequence[s_sequence_title][(max(1, n_position_of_the_window_center - n_window_left_spread) - 1):(min(n_position_of_the_window_center + n_window_right_spread, n_alignment_length))]
		
		d_sequence_title_to_sequence_of_the_window[s_sequence_title] = d_sequence_title_to_sequence_of_the_left_overhang_of_the_window[s_sequence_title] + d_sequence_title_to_sequence_of_the_non_overhanging_part_of_the_window[s_sequence_title] + d_sequence_title_to_sequence_of_the_right_overhang_of_the_window[s_sequence_title]
	
	#Выписываю последовательность окна в отдельный файл и запускаю для него скрипт calculate_several_distance_and_polymorphism_metrics.py
	os.mkdir(s_path_to_the_output_folder + "/Temp")
	f_outfile = open(s_path_to_the_output_folder + "/Temp/alignment_for_the_window.fasta", "w")
	for s_sequence_title in d_sequence_title_to_sequence_of_the_window:
		f_outfile.write(">" + s_sequence_title + "\n" + d_sequence_title_to_sequence_of_the_window[s_sequence_title] + "\n")
	f_outfile.close()
	
	os.system("python3 " + s_path_to_the_folder_where_this_script_lies + "/calculate_polymorphism_using_entire_alignment.py " + s_path_to_the_output_folder + "/Temp/alignment_for_the_window.fasta " + str(n_number_of_cpu_threads_to_use) + " " + s_path_to_the_output_folder + "/Temp/Analysis_results_of_window_polymorphism")
	
	#Беру значения трёх метрик полиморфизма из созданного файла.
	f_infile = open(s_path_to_the_output_folder + "/Temp/Analysis_results_of_window_polymorphism/metrics_of_polymorphism.txt", "r")
	for s_line in f_infile:
		"""
		MPPI is 68.70790903720112
		MPMLD is 0.5236608295047619
		pi is 0.21222624462081333
		"""
		
		o_regular_expression_results = re.search(r"MPPI is (.+)$", s_line)
		if o_regular_expression_results:
			n_MPPI = float(o_regular_expression_results.group(1))
			d_position_of_the_window_center__to__MPPI[n_position_of_the_window_center] = n_MPPI
		
		o_regular_expression_results = re.search(r"MPMLD is (.+)$", s_line)
		if o_regular_expression_results:
			n_MPMLD = float(o_regular_expression_results.group(1))
			d_position_of_the_window_center__to__MPMLD[n_position_of_the_window_center] = n_MPMLD
		
		o_regular_expression_results = re.search(r"pi is (.+)$", s_line)
		if o_regular_expression_results:
			n_pi = float(o_regular_expression_results.group(1))
			d_position_of_the_window_center__to__pi[n_position_of_the_window_center] = n_pi
	
	
	f_infile.close()
	
	#Если на вход было дано парное выравнивание, то значения MPMLD нет (см. заголовочный комментарий)
	if len(d_sequence_title_to_sequence) > 2:
		f_table_with_polymorphism_in_windows.write(str(n_position_of_the_window_center) + ";" + str(n_MPPI) + ";" + str(n_MPMLD) + ";" + str(n_pi) + "; ;\n")
	else:
		f_table_with_polymorphism_in_windows.write(str(n_position_of_the_window_center) + ";" + str(n_MPPI) + ";" + str(n_pi) + "; ;\n")
	
	"""
	f_log.write("For window with the center " + str(n_position_of_the_window_center) + "\n")
	for s_sequence_title in d_sequence_title_to_sequence_of_the_window:
		f_log.write("The sequence for " + s_sequence_title + " is " + d_sequence_title_to_sequence_of_the_window[s_sequence_title] + "\n")
	f_log.write("The list of divergences is " + ", ".join(map(str, l_pairwise_divergences_in_the_window)) + "\n")
	f_log.write("The diversity is " + str(n_diversity_in_the_window) + "\n")
	"""
	
	os.system("rm -rf " + s_path_to_the_output_folder + "/Temp")
	
	n_position_of_the_window_center += n_window_step	
		





