from Bio import SeqIO
from Bio.Seq import Seq
import re
import pandas as pd
import datetime
from multiprocessing import Pool

def search_pcr(args):
    index, pcr = args 
    found = False
    p = re.compile(pcr)
    results = [] 
    for i in SeqIO.parse("caminho/do/arquivo", "fasta"):
        seqstr = str(i.seq)
        matches = p.finditer(seqstr)
        for m in matches:
            if m.group():
                found = True
                start, end = m.span()
                id = i.id  # Get the sequence ID
                results.append((index, start, end, id))  # If need more information to put future in table add here
                print(i)
                print("{start} {match} {end} da linha {index}".format(start=start, match=m.group(), end=end, index=index))

    print(f"{index}, {found}")  
    return results  

def main():
    table = pd.read_csv('arquivo.csv')
    table['find'] = 0
    table['id'] = ""
    table['start'] = 0  
    table['end'] = 0
    start_clock = datetime.datetime.now()
    
    # Define the number of processes to be used (adjust as needed)
    num_processes = 4
    
    with Pool(num_processes) as pool:
        indices_geral = range(len(table['PCR Product']))
        pcr_list = table['PCR Product'].tolist()
        results_list = pool.map(search_pcr, zip(indices_geral, pcr_list))
        
        for aux, results in enumerate(results_list):
            if results:
                found = 1
                index, start, end, id = results[0] 
                table.at[aux, 'find'] = found
                table.at[aux, 'id'] = id
                table.at[aux, 'start'] = start
                table.at[aux, 'end'] = end

    table.to_csv('output.csv', index=False)
    end_clock = datetime.datetime.now()
    print(end_clock - start_clock)

if __name__ == "__main__":
    main()
