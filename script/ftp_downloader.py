import os
import ftplib
import socket
from ftplib import error_temp
from time import sleep 

def download_ftp_file(arg):

    dst, dir, file = arg
    GENOME_PATH = "genomes/GENOME_REPORTS"
    
    try_count = 0
    ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
    ftp.set_pasv(True)

    while(True):
        try:

            f= open(os.path.join(dst, dir, file), "wb")
            ftp.login()
            # IDs
            if len(dir):
                ftp.cwd(GENOME_PATH + "/" + dir)
            # Overview
            else:
                ftp.cwd(GENOME_PATH)
            ftp.retrbinary(f"RETR {file}", f.write)

            break
        except (error_temp, BrokenPipeError, socket.timeout) as e:    
            try_count += 1
            if try_count > 3:
                raise e
            else:
                print(e)
                sleep(try_count * 2)
            try: 
                f.close()
            except: 
                print("error file")
                pass
            try:
                reply = ftp.quit()
                #print(reply)
            except: pass
    
    f.close()
    ftp.quit()