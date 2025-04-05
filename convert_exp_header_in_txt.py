import os
import sys
import csv

dir = sys.argv[1]
inclusive = True

def convert_file(pattern, labels):
    for root, dirs, files in os.walk(dir, topdown=False):
        for name in files:
            if name == pattern:
                path = os.path.join(root, name)
                print(path)
                f = open(path)
                lines_comments = []
                lines_data = []
                for line in f:
                    line1 = line.strip(" \t")            
                    line = line.strip(" \t\n")            
                    if len(line1) == 0:
                        pass
                    elif line1[0] == "#":
                        lines_comments.append(line)
                    elif line1[0:5] =="depth" or line1[0:9] == "air_depth" or line1 == '\n':
                        pass
                    else:                    
                        lines_data.append(line)
                        # if path == 'paleochrono/AICC2023-Hulu1/VK-TALDICE/airair_synchro_horizons.txt':
                        #     print(repr(line))
                f.close()
                if len(lines_data) > 0:
                    # print(repr(lines_data[0]))
                    dialect = csv.Sniffer().sniff(lines_data[0])
                    sep = dialect.delimiter
                else:
                    sep = "\t"
                if name == 'isotopes.txt':
                    nb_cols = len(lines_data[0].split(sep))
                    if nb_cols == 4:
                        header = "depth"+sep+"d18O"+sep+"deut"+sep+"d18Osw\n"
                        print('4 columns file.')
                    elif nb_cols == 2:
                        header = "depth"+sep+"deut\n"
                        print("2 columns file.")
                    else:
                        print('incorrect number of columns in isotopes.txt file')
                else:
                    header = ''                
                    for label in labels:
                        header += label+sep
                    header += 'comment\n'
                f = open(path, 'w')
                for line in lines_comments:
                    f.write(line+'\n')
                f.write(header)
                for line in lines_data:
                    # if path == 'paleochrono/AICC2023-Hulu1/EDC-EDML/iceice_synchro_horizons.txt':
                    #     print(repr(line), len(line.split(sep)), len(header))
                    if len(line.split(sep)) < len(labels):
                        if sep == ' ':
                            line = line.replace('\t', ' ')
                        elif sep == '\t':
                            line = line.replace(' ', '\t')
                    line1 = sep
                    line = line1.join(line.split())
                        # print(repr(line))
                    f.write(line+'\n')
                    # line_split = line.replace('\t', ' ').split()
                    # for elem in line_split[-1]:
                    #     f.write(elem+sep)
                    # f.write(line_split[-1]+'|n')
                f.close()
            
def convert_file_iso():
    pattern = 'isotopes.txt'
    for root, dirs, files in os.walk(dir, topdown=False):
        for name in files:
            if name == pattern:
                path = os.path.join(root, name)
                print(path)
                f = open(path)
                lines_comments = []
                lines_data = []
                for line in f:
                    line = line.strip(" \t\n")            
                    if len(line) == 0:
                        pass
                    elif line[0] == "#":
                        lines_comments.append(line)
                    elif line[0:5] =="depth" or line == '\n':
                        pass
                    else:                    
                        lines_data.append(line)
                f.close()
                if len(lines_data) > 0:
                    dialect = csv.Sniffer().sniff(lines_data[0])
                    sep = dialect.delimiter
                else:
                    sep = "\t"
                nb_cols = len(lines_data[0].split(sep))
                if nb_cols == 4:
                    header = "depth"+sep+"d18O"+sep+"deut"+sep+"d18Osw\n"
                    print('4 columns file.')
                elif nb_cols == 2:
                    header = "depth"+sep+"deut\n"
                    print("2 columns file.", header)
                else:
                    print('incorrect number of columns in isotopes.txt file')
                f = open(path, 'w')
                for line in lines_comments:
                    f.write(line+'\n')
                f.write(header)
                for line in lines_data:
                    if len(line.split(sep)) < len(labels):
                        if sep == ' ':
                            line = line.replace('\t', ' ')
                        elif sep == '\t':
                            line = line.replace(' ', '\t')
                    elif len(line.split(sep)) > len(labels)+1:
                        line1 = sep
                        line = line1.join(line.split())
                        # print(repr(line))
                    f.write(line+'\n')
                f.close()
                
labels = ["depth","age","age_unc"]
for name in ['age_horizons.txt', 'ice_age_horizons.txt', 'air_age_horizons.txt']:
    convert_file(name, labels)
convert_file('density.txt', ["depth", "rel_dens"])
convert_file('deposition.txt', ["depth", "deporate", "rel_unc"])
convert_file('thinning.txt', ["depth", "thinning", "rel_unc"])
convert_file('lock_in_depth.txt', ["depth", "LID", "rel_unc"])
labels = ["depth_top", "depth_bot", "duration", "dur_unc"]
for name in ['ice_age_intervals.txt', 'air_age_intervals.txt', 'age_intervals.txt']:
    convert_file(name, labels)
labels = ["depth1", "depth2", "age_unc"]
for name in ['synchro_horizons.txt', 'iceice_synchro_horizons.txt', 'airair_synchro_horizons.txt',
             "iceair_synchro_horizons.txt", "airice_synchro_horizons.txt",
             "ice_synchro_horizons.txt", "air_synchro_horizons.txt"]:
    convert_file(name, labels)
convert_file('delta_depths.txt', ["air_depth", "Ddepth", "Ddepth_unc"])
convert_file('isotopes.txt', ['depth', 'd18O', 'deut', 'd18Osw'])

