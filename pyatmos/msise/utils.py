import requests
from tqdm import tqdm
from colorama import Fore
# ------------------------------------------------------------------- #
# ----------------------------- utilities --------------------------- #
# ------------------------------------------------------------------- #

'''
pyatmos utils

This submodule defines the following functions:

wraplon - Wrap a longitude in range of [0,360] to [-180,180]

hms2s   - Convert hour/minute/second to seconds

hms2h   - Convert hour/minute/second to hours        
'''

# ========================= convert position ======================== #

def wraplon(lon):
    if lon > 180:
        lonwrap = lon - 360
    else:
        lonwrap = lon
    return lonwrap 

# =========================== convert time ========================== #

def hms2s(h,m,s):
    return h*3.6E3 + m*60 + s

def hms2h(h,m,s):
    return h + m/60 + s/3.6E3

def tqdm_request(url,dir_to,file,desc):
    block_size = 1024
    bar_format = "{l_bar}%s{bar}%s{r_bar}" % (Fore.BLUE, Fore.RESET)
    for idownload in range(5):
        try:
            local_file = open(dir_to + file, 'ab')
            pos = local_file.tell()
            resume_header = {'Range': f'bytes={pos}-'}
            res = requests.get(url,stream=True,timeout=100,headers=resume_header)
            total_size = int(res.headers.get('content-length', 0))
            pbar = tqdm(desc = desc,total=total_size,unit='B',unit_scale=True,bar_format = bar_format,position=0,initial=pos)
            for chunk in res.iter_content(block_size):
                pbar.update(len(chunk))
                local_file.write(chunk)  
            break
        except: 
            sleep(2)
            if idownload == 4:
                remove(dir_to + file)
                print('No response, skip this file.') 
        finally:
            pbar.close()    
            local_file.close()    

