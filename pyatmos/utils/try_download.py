import requests
from tqdm import tqdm
from colorama import Fore

def tqdm_request(url,dir_to,file,desc):
    '''
    Try to download files from a remote server by request with a colored progress bar.
    '''
    block_size = 1024*10
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
            pbar.close()  
            res.close()  
            break
        except: 
            sleep(2)
            if idownload == 4:
                remove(dir_to + file)
                print('No response, skip this file.') 
        finally:    
            local_file.close() 

