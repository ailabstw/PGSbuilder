#!/usr/bin/python3
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys, itertools
from sqlalchemy.pool import NullPool
from sqlalchemy import create_engine, MetaData, Table, Column, Index
from sqlalchemy.types import Integer, BigInteger, String, Text, DateTime, Float
from sqlalchemy.sql import and_, or_


class RSID:
    """rsid
    """
    def __init__(self, _id, name, chrom, start, end, ref, observed):  # pylint: disable=too-many-arguments
        self._id = _id
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.ref = ref
        self.observed = observed


MYSQL_HOST = 'mysql.group-taigenomics.svc.cluster.local'
MYSQL_PORT = 3306
MYSQL_ROOT_PASSWORD = 'mysql-root'
host = MYSQL_HOST
port = MYSQL_PORT
passwd = MYSQL_ROOT_PASSWORD
engine = create_engine(f'mysql+pymysql://root:{passwd}@{host}:{port}/gene',
                       poolclass=NullPool)
metadata = MetaData()


rsid = Table(
    'rsid', metadata,
    Column('_id', Integer(), primary_key=True, autoincrement=True),
    Column('name', BigInteger(), nullable=False),
    Column('chrom', String(50), nullable=False),
    Column('start', BigInteger(), nullable=False),
    Column('end', BigInteger(), nullable=False),
    Column('ref', Text(), nullable=False),
    Column('observed', Text(), nullable=False),
    mysql_engine='MyISAM',
)



bim_file = sys.argv[1]
bim = pd.read_csv(bim_file, sep='\t', names=['chr', 'rsid', '_', 'pos', 'a1', 'a2'])


rsid2result = list()
n = 5000
for i in tqdm(range(int(np.ceil(bim.shape[0] / n)))):
    cond_list = list()
    for j in range(i*n, (i+1)*n):
        if j >= bim.shape[0]:
            break
        req = bim.iloc[j].rsid.split('rs')[1]
        cond_list.append(and_(rsid.c.name == str(req)))
    cond = or_(*cond_list)
    with engine.connect() as conn:
        query = rsid.select().where(cond)
        rows = conn.execute(query)
        results = list(itertools.starmap(RSID, rows))
    for result in results:
        rsid2result.append({'rsid':result.name, 'chr':result.__dict__['chrom'].split('chr')[1], 'pos':result.__dict__['end']})

map_df = pd.DataFrame(rsid2result)
map_df['rsid'] = map_df['rsid'].apply(lambda x: 'rs' + str(x))
map_df = map_df[~map_df['chr'].str.contains(r'[A-Za-z]')]
merge_df = bim.merge(map_df, on='rsid', how='left')
merge_df['chr_y'].fillna(merge_df['chr_x'], inplace=True)
merge_df['pos_y'].fillna(merge_df['pos_x'], inplace=True)
merge_df['pos_y'] = merge_df['pos_y'].astype('int')
merge_df[['chr_y', 'rsid', '_', 'pos_y', 'a1', 'a2']].to_csv(sys.argv[2], sep='\t', header=False, index=False)