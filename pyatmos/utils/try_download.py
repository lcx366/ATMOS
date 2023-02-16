import wget

def wget_download(url,dir_file,desc=None):
    """
    Download files by wget command

    Inputs:
        url -> [str]
        dir_file -> [str] output filename or directory
    Parameters:
        desc -> [str] description of the downloading   
    Outpits:
        wget_out -> [str] path and filename where URL is downloaded to   

    """
    if desc: print(desc)
    wget_out = wget.download(url,dir_file)
    print()

    return wget_out

