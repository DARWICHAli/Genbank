import os
import sys
import ftplib
import logging
import shutil
from multiprocessing import Pool
import time
import datetime


def download_ftp_file(arg):

    # msg = "Start Fetch"
    # print(msg)
    # self.any_signal.emit(msg) 

    dst, dir, file = arg
    GENOME_PATH = "genomes/GENOME_REPORTS"

    # msg = "Logging in to FTP server"
    # print(msg)
    # self.any_signal.emit(msg)

    ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
    ftp.login()


    if len(dir):
        ftp.cwd(GENOME_PATH + "/" + dir)
    else:
        ftp.cwd(GENOME_PATH)

    with open(os.path.join(dst, dir, file), "wb") as f:
        ftp.retrbinary(f"RETR {file}", f.write)