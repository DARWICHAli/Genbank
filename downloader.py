from Bio import Entrez
import pandas as pd
import os
import shutil
from Bio import SeqIO
import pickle
import random
import string
from ftplib import FTP
import time
import os.path
import datetime

save_pickle = False
DEBUG = False
VERBOSE = False

################################################################################
#def reset_tree(root, progressbar, progress = None, window = None):
def reset_tree(root):

    # Reset the tree stored locally and download new one from ftp server

    if os.path.exists('../GENOME_REPORTS/overview.txt') and os.path.isfile("../pickle/organism_df"):
        time = os.path.getmtime("../GENOME_REPORTS/overview.txt") #get the time of last modification
        pred_time = datetime.datetime.fromtimestamp(time)
        now = datetime.datetime.now()

       # if pred_time.month == now.month and pred_time.day == now.day and pred_time.year == now.year : # No update in this case
           # root.destroy() # Exit
           # return

    try: shutil.rmtree('../GENOME_REPORTS') # remove all files/subdirectories
    except: pass
    os.mkdir("../GENOME_REPORTS")
    os.chdir('../GENOME_REPORTS')
    os.mkdir('IDS')
    Kingdom_ID_list = ["Archaea.ids", "Bacteria.ids", "Eukaryota.ids", "Viruses.ids"]

    with FTP('ftp.ncbi.nlm.nih.gov') as ftp:
        ftp.login()  # Connexion to the FTP server

        # Retrieve Genome overview and IDS
        ftp.cwd('genomes/GENOME_REPORTS')
        ftp.retrbinary('RETR overview.txt', open("overview.txt",'wb').write)
        ftp.cwd('IDS')
        for Kingdom_ID in Kingdom_ID_list:
            ftp.retrbinary('RETR '+ Kingdom_ID, open("IDS/" + Kingdom_ID, 'wb').write)
    if VERBOSE:
        print("Genome overview and IDS downloaded !")


    # Parsing of "overview.txt"
    organism_names = []
    organism_paths = []

    os.chdir('../script')
    with open('../GENOME_REPORTS/overview.txt') as f:
        first_row = True
        count_rows = 1
        for row in f:
            if DEBUG:
                print(count_rows, " / 59674")
            count_rows += 1
            if first_row:
                first_row=False
                continue
            parsed_row = row.split('\t')

            # Extraction of the tree components
            try :
                organism = parsed_row[0].replace(' ','_').replace('/','_')
                kingdom = parsed_row[1].replace(' ','_').replace('/','_')
                group = parsed_row[2].replace(' ','_').replace('/','_')
                subgroup = parsed_row[3].replace(' ','_').replace('/','_')
                path = '../Results/' + kingdom +'/' + group +'/' + subgroup +'/' + organism
                organism_names.append(parsed_row[0])
                organism_paths.append('../Results/' + kingdom +'/' + group +'/' + subgroup +'/')
            except IndexError : pass

    # Parsing of the IDS
    ids_files = os.listdir('../GENOME_REPORTS/IDS/')
    if VERBOSE:
        print('overview done !')
    organism_names_ids = []
    organism_paths_ids = []
    organism_NC_ids = []
    i = 0

    for ids in ids_files:
        i += 1
        #progressbar['value'] = 0
        #root.update_idletasks()

        with open('../GENOME_REPORTS/IDS/' + ids) as f:
            n_line = sum(1 for _ in f)
        with open('../GENOME_REPORTS/IDS/' + ids) as f:
            if VERBOSE:
                print("ids")
            for row in f:
                parsed_row = row.replace('\n', '').split('\t')

                if (parsed_row[1][0:2] != 'NC'): # We need only the NC
                    continue
                try:
                    index = organism_names.index(parsed_row[5]) # Retrieve Kingdom
                except ValueError:
                    parsed_name = parsed_row[5].split(' ')[::-1]
                    try_name = parsed_row[5]
                    for word in parsed_name :
                        try_name = try_name.replace(' '+word, '')
                        try:
                            index = organism_names.index(try_name)
                            break
                        except : pass

                #progressbar['value'] += 1/n_line*100
                #root.update_idletasks()

                try:
                    organism_NC_ids[organism_names_ids.index(organism_names[index])].append(parsed_row[1])
                except ValueError:
                    organism_names_ids.append(organism_names[index])
                    organism_paths_ids.append(organism_paths[index])
                    organism_NC_ids.append([parsed_row[1]])
                    name = organism_names[index].replace(" ", "_")
                    name = name.replace("[", "_")
                    name = name.replace("]", "_")
                    name = name.replace(":", "_")
                    path = organism_paths[index] + name + "/"
                    if not os.path.exists(path):
                        os.makedirs(path)

    # Store the organisms in pandas DataFrame
    organism_df = pd.DataFrame({
                "name":organism_names_ids,
                "path":organism_paths_ids,
                "NC":organism_NC_ids})

    # Create a pickle file to save the dataframe in local
    if not os.path.exists("../pickle"):
        os.makedirs("../pickle")
    with open("../pickle/organism_df", 'wb') as f:
        pickle.dump(organism_df, f)

    #root.destroy() # EXIT

################################################################################

def load_df_from_pickle():

    # Loading pickle dataframe

    try:
        with open("../pickle/organism_df", 'rb') as f:
            organism_df = pickle.load(f)
    except IOError:
        print("Pickle file not accessible")
        return reset_tree()

    for i in range(len(organism_df)):
        name = organism_df["name"][i].replace(" ", "_")
        name = name.replace("[", "_")
        name = name.replace("]", "_")
        name = name.replace(":", "_")
        path = organism_df["path"][i] + name + "/"
        if not os.path.exists(path):
            os.makedirs(path)
    return organism_df

################################################################################
################################################################################

def main() :
    start_time = time.time()
    directory = os.getcwd()
    reset_tree(directory)
    load_df_from_pickle()
    print("Arborescence DONE !")
    print("--- %s seconds ---" % (time.time() - start_time))
if __name__ == "__main__":
    main()
