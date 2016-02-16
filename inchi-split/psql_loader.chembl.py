import os,gzip,psycopg2

import splitter

def run():
    cn = psycopg2.connect(dbname='inchi_split')
    curs = cn.cursor()
    fields = 'chemblid text,'+','.join('%s text'%x for x in splitter.Layers._fields)
    qs = ','.join(['%s']*(len(splitter.Layers._fields)+1))
    curs.execute('drop table if exists chembl_export')
    curs.execute('drop table if exists chembl_export_standard')
    curs.execute('drop table if exists chembl_export_nonstandard')
    curs.execute('create table chembl_export (zinc_id text, standard_inchi text, nonstandard_inchi text, smiles text)')
    curs.execute('create table chembl_export_standard (%s)'%fields)
    curs.execute('create table chembl_export_nonstandard (%s)'%fields)
    cn.commit()

    filename = os.path.join('.','test_data','chembl_export.full.txt.gz')  # Note: this file isn't checked in because it's huge
    reader=gzip.open(filename,'r')
    inch_rows=[]
    sinch_rows=[]
    struct_rows=[]
    for i,line in enumerate(reader):
        if not i:
            continue
        if not (i+1)%100000:
            print("Done:",i+1)
            cn.commit()
        # if i>1000000:
        #     break
        line = line.decode('utf-8').strip().split()
        if len(line)<6:
            continue
        sinchi = line[1]
        inchi = line[3]
        smiles = line[5]
        struct_rows.append((line[0],sinchi,inchi,smiles))
        match = splitter.coreExpr.match(sinchi)
        sinch_rows.append((line[0],)+match.groups())
        match = splitter.coreExpr.match(inchi)
        inch_rows.append((line[0],)+match.groups())
        if len(sinch_rows)==10000:
            curs.executemany('insert into chembl_export values (%s, %s, %s, %s)',struct_rows)
            curs.executemany('insert into chembl_export_standard values (%s)'%qs,sinch_rows)
            curs.executemany('insert into chembl_export_nonstandard values (%s)'%qs,inch_rows)
            struct_rows=[]
            sinch_rows=[]
            inch_rows=[]
    if len(sinch_rows):
        curs.executemany('insert into chembl_export values (%s, %s, %s, %s)',struct_rows)
        curs.executemany('insert into chembl_export_standard values (%s)'%qs,sinch_rows)
        curs.executemany('insert into chembl_export_nonstandard values (%s)'%qs,inch_rows)
    cn.commit()

if __name__ == '__main__':
    run()
