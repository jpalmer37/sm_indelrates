#%%
import os 
import sys
import pandas as pd 
from subprocess import Popen, PIPE, check_output
import re
from Bio import Entrez
#import snakemake
import calendar

def padded(num):
    if len(str(num)) == 1:
        return "0"+str(num)
    else:
        return str(num)

#dictionary for ABBREVIATED MONTH NAMES
cal_abbrev = {}
for k,v in enumerate(calendar.month_abbr):
    cal_abbrev[v] = padded(k)

print(cal_abbrev)
#Dictionary for FULL MONTH NAMES
cal_full = {}
for i, x in enumerate(calendar.month_name):
    cal_full[x] = padded(i)
print(cal_full)


# class snake:
#     config = {'delimiter':"_", "accno_pos":1, "date_pos":2, }
#     input=["aligned/A1_filter.fasta"]
#     output=["date/A1_dates.tsv"]
# snakemake = snake()
#%%
def extract_headers(filepath):
    headers, _ = Popen(f'grep ">" {filepath}', shell=True, stdout=PIPE).communicate()
    #headers = check_output('grep ">" ../data/conserved/A1_CR.fasta',shell=True)

    headers = headers.decode('utf-8').strip().split("\n")
    headers = [x.split(snakemake.config['delimiter']) for x in headers]
    header_dict = {x[snakemake.config['accno_pos']]: x[snakemake.config['date_pos']] for x in headers}

    return header_dict
#%%
def extract_date(accno_list):

    handle = Entrez.efetch(db='nucleotide', id=','.join(accno_list), retmode='xml')
    
    expr = re.compile(">(.+)<\/")

    fetched = []
    dates = []
    date_next = False
    found = True

    for line in handle:
        #print(line.decode('utf-8').strip())
        if date_next:
            dates.append(line.decode('utf-8').strip())
            date_next = False
            found = True

        if b'<GBSeq_locus>' in line:
            print(line)
            fetched.append(line.decode('utf-8').strip())
            if not found:
                dates.append("NA")
            found = False

        if b'>collection_date<' in line:
            date_next = True
            
    if not found:
        dates.append("NA")
    fetched = [expr.search(x).group(1) if expr.search(x) else x for x in fetched]
    dates = [expr.search(x).group(1) if expr.search(x) else x for x in dates]

    # missing = set(accno_list) - set(fetched)
    # missing = list(missing)
    # #print(missing)
    # if len(missing) > 0: 
    #     print("WARNING: There are one or more missing dates.")
        
    #     for m in missing:
    #         fetched.append(m)
    #         dates.append("NA")
    print(len(accno_list))
    print(len(dates))
    assert len(accno_list) == len(dates), print("ERROR: There are a differing number of dates found.")
    return dates
#%%

expr0 = re.compile("\d{4}-\d{2}-\d{2}")        # 2008-09-28   ISO FORMAT
expr1 = re.compile("\d{2}-[a-zA-Z]+-\d{4}")    # 28-Sept-2008
expr2 = re.compile("[a-zA-Z]+-\d{4}")          # Sept-2008
expr3 = re.compile("\d{4}")                    # 2008
expr4 = re.compile("[a-zA-Z]+\s\d+,\s\d{4}")   # September 28, 2008

def modify_dates(dates, accnos, header_dict):   
    final = []
    for date, acc in zip(dates, accnos):
        
        if date == " " or date == "NA":
            #print("Retrieving from header dictionary...")
            #print(acc)
            year = header_dict[acc]

            final.append(f"b({year}-01-01,{year}-12-31)")

        elif expr0.match(date):
            final.append(date)

        #Handles the rare September 28, 2008 format
        elif expr1.match(date):
            fields = date.split("-")

            # convert into ISO format 
            fields = [fields[2], cal_abbrev[fields[1]], fields[0]]
            final.append("-".join(fields))

        #Handles all normal dates with
        elif expr2.match(date):
            fields = date.split("-")
            _, days_in_month = calendar.monthrange(int(fields[1]),int(cal_abbrev[fields[0]]))
            fields = "-".join([fields[1], cal_abbrev[fields[0]]])
            final.append(f"b({fields}-01,{fields}-{days_in_month})")

        elif expr3.match(date):
            final.append(f"b({date}-01-01,{date}-12-31)")


        elif expr4.match(date):
            fields = date.split("-")
            fields = [x.strip(", ") for x in fields]
            fields = [fields[2], cal_full[fields[0]], padded(fields[1])]
            final.append("-".join(fields)) 

        else:
            print("CASE NOT FOUND")
            print(date)
            print(acc)
    return final

def write_dates(filename, accnos, dates):
    
    with open(filename,'w') as outfile:
        for acc, date in zip(accnos, dates):
            outfile.write(acc + "\t" + date + "\n")

#%%
if __name__ == '__main__':

    # headers = extract_accnos('../data/conserved/A1_CR.fasta')
    header_dict = extract_headers(snakemake.input[0])
    accnos = list(header_dict.keys())
    dates = extract_date(accnos)
    final_dates = modify_dates(dates, accnos, header_dict)
    write_dates(snakemake.output[0], accnos, final_dates)
    