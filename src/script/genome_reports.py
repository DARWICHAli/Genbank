import os
import sys
import ftplib
import logging
import shutil
from multiprocessing import Pool
import time
import datetime


########################################################################################################################
def download_ftp_file(arg):
    dst, dir, file = arg
    GENOME_PATH = "genomes/GENOME_REPORTS"

    ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
    ftp.login()
    if len(dir):
        ftp.cwd(GENOME_PATH + "/" + dir)
    else:
        ftp.cwd(GENOME_PATH)

    with open(os.path.join(dst, dir, file), "wb") as f:
        ftp.retrbinary(f"RETR {file}", f.write)

########################################################################################################################

def download_ftp():
    try: shutil.rmtree('../GENOME_REPORTS') # remove all files/subdirectories
    except: pass
    os.mkdir("../GENOME_REPORTS")
    os.chdir('../GENOME_REPORTS')

    directory = os.getcwd()

    if os.path.exists(os.path.join(directory, "IDS")):
        shutil.remove(os.path.join(directory, "IDS"))
    os.mkdir(os.path.join(directory, "IDS"))

    files = [("", "overview.txt"), ("IDS", "Bacteria.ids"),
             ("IDS", "Eukaryota.ids"), ("IDS", "Archaea.ids"),
             ("IDS", "Viruses.ids")]

    # Adding destination
    files = [(directory,) + f for f in files]

    with Pool(5) as p:
        p.map(download_ftp_file, files)


#def main() :
#    start_time = time.time()
#    directory = os.getcwd()
#    download_ftp()
#    print("--- %s seconds ---" % (time.time() - start_time))
#if __name__ == "__main__":
#    main()
