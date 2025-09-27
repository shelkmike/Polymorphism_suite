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
import statistics
import natsort
import dendropy

s_path_to_alignment = sys.argv[1]
n_number_of_cpu_threads_to_use = int(sys.argv[2])
s_path_to_the_output_folder = sys.argv[3]

#Определяю путь к текущей папке. Нужен, чтобы воспользоваться IQ-TREE, который прилагается к этому скрипту.
s_path_to_the_folder_where_this_script_lies = os.path.dirname(os.path.realpath( __file__ )) #Путь к папке, где лежит Mabs. Использую os.path.realpath, чтобы если скрипт запускается по мягкой ссылке на исполняемый файл, всё равно удалось обнаружить папку с IQ-TREE и папку.


#Проверяю, что выходной папки не существует, либо она существует и пустая. В противном случае, говорю пользователю, что это ошибка.
if os.path.exists(s_path_to_the_output_folder):
	if len(os.listdir(s_path_to_the_output_folder)) != 0:
		print("Error: the output folder already exists and is not empty.")
		sys.exit()
else:
	os.mkdir(s_path_to_the_output_folder)

f_log = open(s_path_to_the_output_folder + "/log.txt", "w", buffering = 1)

#загружаю все последовательности в словарь. Ключи в словаре это заголовки последовательностей без ">", а значение это сами последовательности. Для простоты, называю последовательности "контигами". Чтобы у IQ-TREE не было проблем, удаляю всё, что в заголовке идёт после первого пробельного символа.
#При загрузке последовательностей делаю все буквы большими, чтобы маленькая и такая же большая буква не создавали мисматчей.
f_infile = open(s_path_to_alignment, "r")

d_contig_title_to_contig_sequence = {}
s_contig_title = "" #заголовок контига, на который я сейчас смотрю.

for s_line in f_infile:
	if re.search(r"^>(\S+)", s_line):
		o_regular_expression_results = re.search(r"^>(\S+)", s_line)
		s_contig_title = o_regular_expression_results.group(1)
	elif re.search(r"^(.+)", s_line): #если это не строка с заголовком, то считаю, что строка с последовательностью
		o_regular_expression_results = re.search(r"^(.+)", s_line)
		#если для этого контига последовательность ещё не инициализирована, то инициализирую её
		if s_contig_title not in d_contig_title_to_contig_sequence:
			d_contig_title_to_contig_sequence[s_contig_title] = ""
		
		s_sequence_from_this_string = o_regular_expression_results.group(1)
		#удаляю всякие пробельные символы, в том числе символ переноса строки.
		s_sequence_from_this_string = re.sub(r"\s","",s_sequence_from_this_string)
		#Делаю маленькие буквы большими.
		s_sequence_from_this_string = s_sequence_from_this_string.upper()

		d_contig_title_to_contig_sequence[s_contig_title] += s_sequence_from_this_string

#Список заголовков контигов.
l_contig_titles = list(d_contig_title_to_contig_sequence.keys())

#Делаю естественную сортировку заголовков контигов.
l_contig_titles = natsort.natsorted(l_contig_titles)

if len(l_contig_titles) == 1:
	f_log.write("There is " + str(len(l_contig_titles)) + " sequence in the input alignment.\n")
else:
	f_log.write("There are " + str(len(l_contig_titles)) + " sequences in the input alignment.\n")

#Если в выравнивании меньше двух последовательностей, то завершаю работу скрипта.
if len(l_contig_titles) < 2:
	f_log.write("Stopping calculations because the script requires at least 2 sequences in the alignment.\n")
	sys.exit()

f_metrics_of_polymorphism = open(s_path_to_the_output_folder + "/metrics_of_polymorphism.txt", "w") #Файл, в котором несколько метрик полиморфизма.

#Теперь считаю парные сходства и заполняю соответствующую таблицу. Кроме этого, запишу все значения из таблицы (кроме диагональных) в список l_pairwise_comparison_values__except_diagonal_values . По ней потом посчитаю среднее значение.

l_pairwise_comparison_values__except_diagonal_values = []

f_table = open(s_path_to_the_output_folder + "/percent_identity_table.csv", "w")
f_table.write("Sequence title")
for s_contig_title in l_contig_titles:
	f_table.write(";" + s_contig_title)
f_table.write("\n")

for s_first_contig_title in l_contig_titles:
	f_table.write(s_first_contig_title)
	for s_second_contig_title in l_contig_titles:
		
		#Если я сравниваю контиг сам с собой, то просто ставлю в ячейке пробел, без чисел.
		if s_first_contig_title == s_second_contig_title:
			f_table.write("; ")
		else:
			s_first_contig_sequence = d_contig_title_to_contig_sequence[s_first_contig_title]
			s_second_contig_sequence = d_contig_title_to_contig_sequence[s_second_contig_title]
			
			#Длина у последовательностей одинаковая, потому что это последовательности из множественного выравнивания.
			
			n_alignment_length = 0 #Длина выравнивания конкретно этих двух последовательностей. Она определяется после удаления тех столбцов, где гэп у обоих последовательностей.
			n_number_of_matches = 0 #Количество матчей между двумя последовательностями.
			for s_position in range(1, (len(s_first_contig_sequence) + 1)):
				s_character_in_the_first_contig = s_first_contig_sequence[s_position - 1]
				s_character_in_the_second_contig = s_second_contig_sequence[s_position - 1]
				
				#Если в обеих последовательностях гэп
				if (s_character_in_the_first_contig == "-") and (s_character_in_the_second_contig == "-"):
					pass
				#Если между последовательностями матч
				elif s_character_in_the_first_contig == s_character_in_the_second_contig:
					n_number_of_matches += 1
					n_alignment_length += 1
				#Если между последовательностями мисматч
				else:
					n_alignment_length += 1
			
			#Если во множественном выравнивании две последовательности вообще не пересекаются
			if n_alignment_length == 0:
				n_percent_identity = 0
			else:
				n_percent_identity = 100 * n_number_of_matches / n_alignment_length
			
			l_pairwise_comparison_values__except_diagonal_values.append(n_percent_identity)
			
			f_table.write(";" + str(n_percent_identity))
	
	f_table.write("\n")

#Это потом удалить
#if len(l_pairwise_comparison_values__except_diagonal_values) == 0:
#	print("Error: for the alignment " + s_path_to_alignment + " len(l_pairwise_comparison_values__except_diagonal_values) == 0")

if len(l_pairwise_comparison_values__except_diagonal_values) != 1:
	#Пополам делю потому, что сравнение A с B и сравнение B с A это, фактически, одно сравнение. В том смысле, что наличие двух сравнений, а не одного, никак не влияет на результат.
	f_log.write("Calculation of MPPI done based on " + str(int(len(l_pairwise_comparison_values__except_diagonal_values) / 2)) + " pairwise comparisons.\n")
else:
	f_log.write("Calculation of MPPI done based on 1 pairwise comparison.\n")

n_MPPI = statistics.mean(l_pairwise_comparison_values__except_diagonal_values)
f_metrics_of_polymorphism.write("MPPI is " + str(n_MPPI) + "\n")

#Теперь, если последовательностей как минимум 3, то делаю анализы методом максимального правдоподобия непосредственно для них. Если же последовательности всего 2, то добавляю 3-ю, которая идентична одной из этих двух, но при этом имеет заголовок ">dummy". После построения дерева удалю её из вычислений.

if len(l_contig_titles) >= 3:
	#Все значения из таблицы (кроме диагональных) в список l_pairwise_comparison_values__except_diagonal_values . По ней потом посчитаю среднее значение.
	l_pairwise_comparison_values__except_diagonal_values = []

	os.mkdir(s_path_to_the_output_folder + "/IQ-TREE_results")
	os.system(s_path_to_the_folder_where_this_script_lies + "/Additional/IQ-TREE1/iqtree -s " + s_path_to_alignment + " -nt " + str(n_number_of_cpu_threads_to_use) + " -seed 12345 -pre " + s_path_to_the_output_folder + "/IQ-TREE_results/iqtree_results")
	
	#Открываю дерево и смотрю на матрицу расстояний.
	s_path_to_the_tree = s_path_to_the_output_folder + "/IQ-TREE_results/iqtree_results.treefile"
	o_tree = dendropy.Tree.get(path = s_path_to_the_tree, schema = "newick", preserve_underscores = True) #"preserve_underscores = True" чтобы он не заменял подчёркивания на пробелы. По умолчанию заменяет.
	o_distance_matrix = o_tree.phylogenetic_distance_matrix()
	
	#Сделаю двойной словарь, в котором первый ключ это название первого контига, второй ключ это название второго контига, а значение это MLD. В случае, если ключи одинаковые, значения не будет.
	dd_first_contig_title__and__second_contig_title__to__MLD = {}
	for o_first_contig_object in o_tree.taxon_namespace:
		for o_second_contig_object in o_tree.taxon_namespace:
			s_first_contig_title = o_first_contig_object.label
			s_second_contig_title = o_second_contig_object.label
						
			if s_first_contig_title not in dd_first_contig_title__and__second_contig_title__to__MLD:
				dd_first_contig_title__and__second_contig_title__to__MLD[s_first_contig_title] = {}
			
			dd_first_contig_title__and__second_contig_title__to__MLD[s_first_contig_title][s_second_contig_title] = o_distance_matrix.patristic_distance(o_first_contig_object, o_second_contig_object)
	
	f_table = open(s_path_to_the_output_folder + "/MLD_table.csv", "w")
	f_table.write("Sequence title")
	for s_contig_title in l_contig_titles:
		f_table.write(";" + s_contig_title)
	f_table.write("\n")

	for s_first_contig_title in l_contig_titles:
		f_table.write(s_first_contig_title)
		for s_second_contig_title in l_contig_titles:
			
			#Если я сравниваю контиг сам с собой, то просто ставлю в ячейке пробел, без чисел.
			if s_first_contig_title == s_second_contig_title:
				f_table.write("; ")
			else:
				n_MLD = dd_first_contig_title__and__second_contig_title__to__MLD[s_first_contig_title][s_second_contig_title]
				
				l_pairwise_comparison_values__except_diagonal_values.append(n_MLD)
				
				f_table.write(";" + str(n_MLD))
		
		f_table.write("\n")
	
	#Пополам делю потому, что сравнение A с B и сравнение B с A это, фактически, одно сравнение. В том смысле, что наличие двух сравнений, а не одного, никак не влияет на результат.
	f_log.write("Calculation of MPMLD done based on " + str(int(len(l_pairwise_comparison_values__except_diagonal_values) / 2)) + " pairwise comparisons.\n")
	
	n_MPMLD = statistics.mean(l_pairwise_comparison_values__except_diagonal_values)
	f_metrics_of_polymorphism.write("MPMLD is " + str(n_MPMLD) + "\n")
	

#Код для случая со всего двумя контигами частично основан на коде для множества контигов. В частности, поэтому я использую "таблицу сходства". Но на результаты это не влияет.
elif len(l_contig_titles) >= 2:
	
	#Все значения из таблицы (кроме диагональных) в список l_pairwise_comparison_values__except_diagonal_values . По ней потом посчитаю среднее значение.
	l_pairwise_comparison_values__except_diagonal_values = []
	
	#Как я и пишу выше, я создаю выравнивание, в которое в качестве третьей последовательности добавлена вторая, в которой заголовок поменян на "dummy".
	f_outfile = open(s_path_to_the_output_folder + "/alignment_with_dummy_sequence.fasta", "w")
	for s_contig_title in d_contig_title_to_contig_sequence:
		f_outfile.write(">" + s_contig_title + "\n")
		f_outfile.write(d_contig_title_to_contig_sequence[s_contig_title] + "\n")
	
	#Поскольку в переменной s_contig_title сейчас записан заголовок второй последовательности, то чтобы сделать "dummy", достаточно записать эту последовательность ещё раз.
	f_outfile.write(">dummy\n")
	f_outfile.write(d_contig_title_to_contig_sequence[s_contig_title] + "\n")
	
	f_outfile.close()
	
	os.mkdir(s_path_to_the_output_folder + "/IQ-TREE_results")
	os.system(s_path_to_the_folder_where_this_script_lies + "/Additional/IQ-TREE1/iqtree -s " + s_path_to_the_output_folder + "/alignment_with_dummy_sequence.fasta -nt " + str(n_number_of_cpu_threads_to_use) + " -seed 12345 -pre " + s_path_to_the_output_folder + "/IQ-TREE_results/iqtree_results")
	
	#Открываю дерево и смотрю на матрицу расстояний.
	s_path_to_the_tree = s_path_to_the_output_folder + "/IQ-TREE_results/iqtree_results.treefile"
	o_tree = dendropy.Tree.get(path = s_path_to_the_tree, schema = "newick", preserve_underscores = True) #"preserve_underscores = True" чтобы он не заменял подчёркивания на пробелы. По умолчанию заменяет.
	#Удаляю лист "dummy"
	o_tree.prune_taxa_with_labels(["dummy"], update_bipartitions = True)
	#Удаляю название таксона "dummy" из списка имеющихся таксонов. Иначе он выдаётся при использовании ".taxon_namespace" ниже.
	o_tree.purge_taxon_namespace()
	o_distance_matrix = o_tree.phylogenetic_distance_matrix()
	
	#Сделаю двойной словарь, в котором первый ключ это название первого контига, второй ключ это название второго контига, а значение это MLD. В случае, если ключи одинаковые, значения не будет.
	dd_first_contig_title__and__second_contig_title__to__MLD = {}
	for o_first_contig_object in o_tree.taxon_namespace:
		for o_second_contig_object in o_tree.taxon_namespace:
			s_first_contig_title = o_first_contig_object.label
			s_second_contig_title = o_second_contig_object.label
						
			if s_first_contig_title not in dd_first_contig_title__and__second_contig_title__to__MLD:
				dd_first_contig_title__and__second_contig_title__to__MLD[s_first_contig_title] = {}
			
			dd_first_contig_title__and__second_contig_title__to__MLD[s_first_contig_title][s_second_contig_title] = o_distance_matrix.patristic_distance(o_first_contig_object, o_second_contig_object)
	
	f_table = open(s_path_to_the_output_folder + "/MLD_table.csv", "w")
	f_table.write("Sequence title")
	for s_contig_title in l_contig_titles:
		f_table.write(";" + s_contig_title)
	f_table.write("\n")

	for s_first_contig_title in l_contig_titles:
		f_table.write(s_first_contig_title)
		for s_second_contig_title in l_contig_titles:
			
			#Если я сравниваю контиг сам с собой, то просто ставлю в ячейке пробел, без чисел.
			if s_first_contig_title == s_second_contig_title:
				f_table.write("; ")
			else:
				n_MLD = dd_first_contig_title__and__second_contig_title__to__MLD[s_first_contig_title][s_second_contig_title]
				
				l_pairwise_comparison_values__except_diagonal_values.append(n_MLD)
				
				f_table.write(";" + str(n_MLD))
		
		f_table.write("\n")
	
	f_log.write("Calculation of MPMLD done based on 1 pairwise comparison.\n")
	
	n_MPMLD = statistics.mean(l_pairwise_comparison_values__except_diagonal_values)
	f_metrics_of_polymorphism.write("MPMLD is " + str(n_MPMLD) + "\n")
	
#Теперь считаю p-distance и заполняю соответствующую таблицу. Кроме этого, запишу все значения из таблицы (кроме диагональных) в список l_pairwise_comparison_values__except_diagonal_values . По ней потом посчитаю среднее значение.

l_pairwise_comparison_values__except_diagonal_values = []

f_table = open(s_path_to_the_output_folder + "/p-distance_table.csv", "w")
f_table.write("Sequence title")
for s_contig_title in l_contig_titles:
	f_table.write(";" + s_contig_title)
f_table.write("\n")

for s_first_contig_title in l_contig_titles:
	f_table.write(s_first_contig_title)
	for s_second_contig_title in l_contig_titles:
		
		#Если я сравниваю контиг сам с собой, то просто ставлю в ячейке пробел, без чисел.
		if s_first_contig_title == s_second_contig_title:
			f_table.write("; ")
		else:
			s_first_contig_sequence = d_contig_title_to_contig_sequence[s_first_contig_title]
			s_second_contig_sequence = d_contig_title_to_contig_sequence[s_second_contig_title]
			
			#Длина у последовательностей одинаковая, потому что это последовательности из множественного выравнивания.
			
			n_alignment_length = 0 #Длина выравнивания конкретно этих двух последовательностей. Она определяется после удаления тех столбцов, где гэп у обоих последовательностей.
			n_number_of_matches = 0 #Количество матчей между двумя последовательностями.
			for s_position in range(1, (len(s_first_contig_sequence) + 1)):
				s_character_in_the_first_contig = s_first_contig_sequence[s_position - 1]
				s_character_in_the_second_contig = s_second_contig_sequence[s_position - 1]
				
				#Если хотя бы в одной последовательности гэп
				if (s_character_in_the_first_contig == "-") or (s_character_in_the_second_contig == "-"):
					pass
				#Если между последовательностями матч
				elif s_character_in_the_first_contig == s_character_in_the_second_contig:
					n_number_of_matches += 1
					n_alignment_length += 1
				#Если между последовательностями мисматч
				else:
					n_alignment_length += 1
			
			#Если во множественном выравнивании две последовательности вообще не пересекаются
			if n_alignment_length == 0:
				n_p_distance = 1
			else:
				n_p_distance = (n_alignment_length - n_number_of_matches) / n_alignment_length
			
			l_pairwise_comparison_values__except_diagonal_values.append(n_p_distance)
			
			f_table.write(";" + str(n_p_distance))
	
	f_table.write("\n")

if len(l_pairwise_comparison_values__except_diagonal_values) != 1:
	#Пополам делю потому, что сравнение A с B и сравнение B с A это, фактически, одно сравнение. В том смысле, что наличие двух сравнений, а не одного, никак не влияет на результат.
	f_log.write("Calculation of pi done based on " + str(int(len(l_pairwise_comparison_values__except_diagonal_values) / 2)) + " pairwise comparisons.\n")
else:
	f_log.write("Calculation of pi done based on 1 pairwise comparison.\n")
	
n_pi = statistics.mean(l_pairwise_comparison_values__except_diagonal_values)
f_metrics_of_polymorphism.write("pi is " + str(n_pi) + "\n")
