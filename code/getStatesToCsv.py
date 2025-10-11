# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 00:01:33 2024

@author: micha
"""
import requests
import json
from datetime import datetime
import sys

# ------------------------------------
fromT = datetime(2024,1,17,16,44).timestamp()
toT = datetime(2024,1,17,18,50).timestamp()

fromT = datetime(2024,4,17,5,44).timestamp()
toT = datetime(2024,4,18,21,15).timestamp()

# me home-b7-home-net-home
fromT = datetime(2024,5,12,6,00).timestamp()
toT = datetime(2024,5,13,6,00).timestamp()

# dora home-b7-home
fromT = datetime(2024,5,13,6,00).timestamp()
toT = datetime(2024,5,13,13,00).timestamp()

# me home-b7-home-ash-whatever-home (aome fast!!!)
fromT = datetime(2024,5,13,13,00).timestamp()
toT = datetime(2024,5,14,15,00).timestamp()

# me econ
fromT = datetime(2024,5,14,7,00).timestamp()
toT = datetime(2024,5,14,12,5).timestamp()
# ------------------------------------

def getStatesToCsv(fromT, toT, fileprefix=""):
    if not isinstance(fromT, str):
        fromT = f"{fromT:.0f}"
    if not isinstance(toT, str):
        toT = f"{toT:.0f}"   
    url = "https://api.tessie.com/LRW3E7FS8PC835950/states?from=" + fromT + "&to=" + toT + "&interval=1&condense=false&timezone=UTC&distance_format=km&temperature_format=c&format=csv"

    headers = {
        "accept": "application/json",
        "authorization": "Bearer QG7GgITOHCKxIZxv3laKRBDTuPMDEBOz"
    }
    response = requests.get(url, headers=headers)

    str1 = response.text
    str1 = str1.replace("(", "");
    str1 = str1.replace(")", "");
    str1 = str1.replace(" ", "_");
    str1 = str1.replace("%", "prc");
    str1 = str1.replace("/", "_per_");
    str1 = str1.replace("°", "");
    
    #tstamp0 = datetime.fromtimestamp(int(fromT)).strftime("%d%m%Y_%H%M")
    #tstamp1 = datetime.fromtimestamp(int(toT)).strftime("%d%m%Y_%H%M")
    
    tstamp0 = datetime.fromtimestamp(int(fromT)).strftime("%Y%m%d_%H%M")
    tstamp1 = datetime.fromtimestamp(int(toT)).strftime("%Y%m%d_%H%M")
    if fileprefix:
        filename = fileprefix + "_f" + tstamp0 + "_t" +  tstamp1 + ".csv"
    else:
        filename = "f" + tstamp0 + "_t" +  tstamp1 + ".csv"
    
    f = open(filename, "w")
    f.write(str1)
    f.close()
    
    
if __name__ == '__main__':
    args = sys.argv[1:]
    if args.__len__()>0:
        print(args)
        fromT = datetime.fromisoformat(args[0]).timestamp()
        toT = datetime.fromisoformat(args[1]).timestamp()
        getStatesToCsv(fromT, toT)
    else:
        print("enter start and end times in ISO format (YYYY-MM-DD hh:mm)")
    
    
