# Metodo seqIO
from Bio import SeqIO
from Bio.Seq import Seq
import re
import pandas
import datetime

aux = 0
tabela = pandas.read_csv('STRs.csv')
tabela = pandas.DataFrame(tabela)
tabela['find'] = 0


start_clock = datetime.datetime.now()
for pcr in tabela['PCR Product']:

    print(aux)
    achou = False
    for i in SeqIO.parse("arraia.fasta", "fasta"):
        p = re.compile(pcr)
        seqstr = str(i.seq)
        matches = p.finditer(seqstr)
        for m in matches:
            if m.group():
                achou = True
                start, end = m.span()
                print(i)
                print(
                    "{start} {match} {end}".format(
                        start=start, match=m.group(), end=end)
                )
    print(achou)
    tabela['find'][aux] = 1 if achou else 0
    aux += 1
    if aux == 4:
        tabela.to_csv('output.csv', index=False)
        break
end_clock = datetime.datetime.now()
print(end_clock - start_clock)
