#python




def take_tsv_to_reactions_file(filename, out_filepath):

    f = open(filename, "r")
    file_str = f.read()
    f.close()
    file_lines = file_str.split('\n')
    file_table = []
    for line in file_lines:
        file_table.append(line.split('\t'))

    reaction_file_str = ''
    file_table = file_table[1:]
    for line in file_table:
        if len(line)> 9:
            reaction_file_str += line[0] + ": " + line[8] + '\n'
        else:
            print(line)

    g = open(out_filepath, "w")
    g.write(reaction_file_str)
    g.close()

    return 0

take_tsv_to_reactions_file('/Users/omreeg/Programs/Arkin_Lab_Research_Home/Current_Projects/Mika/Ecoli_added_vio_complete4.TSV/BiGG_Ecoli_MG1665-reactions.tsv', "Ec_reactions.txt") 

