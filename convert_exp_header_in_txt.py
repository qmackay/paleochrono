import os
import re
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
                    line = line.strip(" \t")            
                    if line[0] == "#":
                        lines_comments.append(line)
                    elif line[0:5] =="depth" or line == '\n':
                        pass
                    else:                    
                        lines_data.append(line)
                if len(lines_data) > 0:
                    dialect = csv.Sniffer().sniff(lines_data[0])
                    sep = dialect.delimiter
                else:
                    sep = "\t" 
                header = ''
                for label in labels[:-1]:
                    header += label+sep
                header += labels[-1]+"\n"
                f.close()
                f = open(path, 'w')
                for line in lines_comments:
                    f.write(line)
                f.write(header)
                for line in lines_data:
                    f.write(line)
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
                    line = line.strip(" \t")            
                    if line[0] == "#":
                        lines_comments.append(line)
                    elif line[0:5] =="depth" or line == '\n':
                        pass
                    else:                    
                        lines_data.append(line)
                if len(lines_data) > 0:
                    dialect = csv.Sniffer().sniff(lines_data[0])
                    sep = dialect.delimiter
                else:
                    sep = "\t"
                nb_cols = len(lines_data[0].split(sep))
                if nb_cols == 4:
                    header = "depth"+sep+"d18O"+sep+"deut"+sep+"d18Osw\n"
                elif nb_cols == 2:
                    header = "depth"+sep+"deut\n"
                header += labels[-1]+"\n"
                f.close()
                f = open(path, 'w')
                for line in lines_comments:
                    f.write(line)
                f.write(header)
                for line in lines_data:
                    f.write(line)
                f.close()
                
labels = ["depth","age","age_unc"]
for name in ['age_horizons.txt', 'ice_age_horizons.txt', 'air_age_horizons.txt']:
    convert_file(name, labels)
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
convert_file_iso()

