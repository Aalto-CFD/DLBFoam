import numpy as np







def write_sp_list(sp_names, fileName):

    with open(fileName, 'a') as output:

        output.write('species\n(\n')

        for sp_i in sp_names:
            output.write('\t')
            output.write(sp_i)
            output.write('\n')

        output.write(');\n\n')





