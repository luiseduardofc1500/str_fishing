from Bio import SeqIO
from Bio.Seq import Seq
import re
import pandas as pd
import datetime
from multiprocessing import Pool

def search_pcr(args):
    index, pcr = args 
    achou = False
    p = re.compile(pcr)
    results = []  # Store the results for this PCR product
    for i in SeqIO.parse("arraia.fasta", "fasta"):
        seqstr = str(i.seq)
        matches = p.finditer(seqstr)
        for m in matches:
            if m.group():
                achou = True
                start, end = m.span()
                id = i.id  # Get the sequence ID
                results.append((index, start, end, id))  # Append (index, start, end, id) to results
                print(i)
                print("{start} {match} {end} da linha {index}".format(start=start, match=m.group(), end=end, index=index))

    print(f"{index}, {achou}")  
    return results  # Return a list of results

def main():
    tabela = pd.read_csv('cortadinho1.csv', sep=";")
    tabela['find'] = 0
    tabela['id'] = ""
    tabela['start'] = 0  
    tabela['end'] = 0
    start_clock = datetime.datetime.now()
    
    # Define the number of processes to be used (adjust as needed)
    num_processes = 4
    
    with Pool(num_processes) as pool:
        indices_geral = range(len(tabela['PCR Product']))
        pcr_list = tabela['PCR Product'].tolist()
        results_list = pool.map(search_pcr, zip(indices_geral, pcr_list))
        
        for aux, results in enumerate(results_list):
            if results:
                achou = 1
                index, start, end, id = results[0]  # Take the first result
                tabela.at[aux, 'find'] = achou
                tabela.at[aux, 'id'] = id
                tabela.at[aux, 'start'] = start
                tabela.at[aux, 'end'] = end

    tabela.to_csv('output.csv', index=False)
    end_clock = datetime.datetime.now()
    print(end_clock - start_clock)

if __name__ == "__main__":
    main()
