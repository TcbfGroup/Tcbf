import os.path

from pandas import read_table
def parse_block(workdir,reference,query):
    query_boundries = os.path.join(workdir,"Step1",f"{query}")