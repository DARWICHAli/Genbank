from Bio import *
import pandas as pd
import os
import pickle
import random
import string
import ftplib
import time
import os.path
import datetime
from multiprocessing import Pool

from genome_reports import download_ftp

save_pickle = False
DEBUG = False
VERBOSE = False

################################################################################
def reset_tree():

    # Reset the tree stored locally and download new one from ftp server
    download_ftp()
    if VERBOSE:
        print("Genome overview and IDS downloaded !")

    # Parsing of "overview.txt"
    organism_names = []
    organism_paths = []

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
                        except: pass

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

################################################################################

def load_df_from_pickle():

    # Loading the Tree from the pickle without Downloading IDS and overview from FTP
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

def load_tree():
    files = ["Bacteria.ids", "Eukaryota.ids", "Archaea.ids", "Viruses.ids"]

    if os.path.isdir("../pickle") and os.path.isfile("../pickle/organism_df"):
        # Initialization
        ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
        ftp.login()
        ftp.cwd('genomes/GENOME_REPORTS/IDS')

        last_ftp_change = 0
        last_local_change = 1e50

        # ftp files last modification timestamp
        for f in files:
            remote_datetime = ftp.voidcmd("MDTM " + f)[4:].strip()
            remote_timestamp = time.mktime(time.strptime(remote_datetime, '%Y%m%d%H%M%S'))
            if int(remote_timestamp) > int(last_ftp_change):
                last_ftp_change = remote_timestamp
        #print(last_ftp_change)

        # local files timestamp
        last_local_change = os.path.getmtime("../pickle/organism_df")
        #print(last_local_change)

        # Download the newest files and create the tree
        if int(last_ftp_change) > int(last_local_change):
            reset_tree()
            #print("loaded from ftp")
        else:
            load_df_from_pickle()
            #print("loaded from pickle")
    else:
        reset_tree()
        #print("loaded from ftp (file or directory doesn't exist)")


################################################################################
################################################################################

def main() :
    start_time = time.time()

    load_tree()

    #for (index, name, path, NC_list) in organism_df.itertuples():
    #    print(NC_list)

    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    main()
