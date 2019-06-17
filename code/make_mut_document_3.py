# -*- coding: utf-8 -*-                                                                                                                     
                                
#!/usr/bin/env python                                                                                                                       
                                

#スクリプト                                                                                                                                 
                                

#python make_mut_document_2.py test_CosmicNCV.csv hg38.fa    

import pandas as pd
from Bio import SeqIO
import re
import sys

mut_file=sys.argv[1]
ref_fasta=sys.argv[2]

df_cosmic=pd.read_table(mut_file)
df_cut=df_cosmic[['ID_SAMPLE', 'ID_tumour','genome position','WT_SEQ', 'MUT_SEQ']]

del df_cosmic

#レファレンスファイルについて                                                                                                               
                                
list_id=[]
list_desc=[]
list_seq=[]

#biopythonを用いて、idやdesc、seqを抽出                                                                                                     
                                
for record in SeqIO.parse(ref_fasta, 'fasta'):
    id_part = record.id
    desc_part = record.description
    seq = record.seq

    list_id.append(id_part)
    list_desc.append(desc_part)
    list_seq.append(seq)

#mutについて、genome positiionから、chrとsingle positionに分ける
list_cut_chr=[]
list_cut_pos=[]
for i in range(len(df_cut)):
    list_cut_chr.append((re.split('[:-]',df_cut['genome position'][i]))[0])
    list_cut_pos.append((re.split('[:-]',df_cut['genome position'][i]))[1])
df_cut['chr']=pd.Series(list_cut_chr)
df_cut['single position']=pd.Series(list_cut_pos)

#mutのdfからchr23、24、MUT_SEQのnanを除く

df_cut_chr=df_cut[(df_cut['chr'] != '23') & (df_cut['chr'] != '24') & (df_cut['chr'] != '25')]
df_cut_chr.dropna(subset=['MUT_SEQ'],inplace=True)
df_cut_chr.dropna(subset=['WT_SEQ'],inplace=True)

#indexを振り直す
df_cut_chr_i=df_cut_chr.reset_index()

del df_cut_chr

#refから変異位置の前後の配列を抽出
seq_before=[]
seq_after=[]
for i in range(len(df_cut_chr_i)):
    seq_before.append(list_seq[list_id.index('chr'+df_cut_chr_i['chr'][i])][int(df_cut_chr_i['single position'][i])-2])
    seq_after.append(list_seq[list_id.index('chr'+df_cut_chr_i['chr'][i])][int(df_cut_chr_i['single position'][i])])
    
#df_cut_chr_iにつなげる
df_cut_chr_i['seq_before']=pd.Series(seq_before)
df_cut_chr_i['seq_after']=pd.Series(seq_after)

#大文字に変更
df_cut_chr_i['seq_before']=df_cut_chr_i['seq_before'].str.upper()
df_cut_chr_i['seq_after']=df_cut_chr_i['seq_after'].str.upper()

#3文字の配列にする（seq_forwardが鋳型）
df_cut_chr_i['seq_forward']=df_cut_chr_i['seq_before']+df_cut_chr_i['MUT_SEQ']+df_cut_chr_i['seq_after']
df_cut_chr_i['seq_forward_r']=df_cut_chr_i['seq_after']+df_cut_chr_i['MUT_SEQ']+df_cut_chr_i['seq_before']

#辞書型で相補鎖に変換
dict_base={"A":"T","G":"C","C":"G","T":"A"}

list_seq_reverse=[]
list_wt_seq_reverse=[]
for i in range(len(df_cut_chr_i)):
    list_seq_reverse.append((df_cut_chr_i['seq_forward_r'][i]).translate(str.maketrans(dict_base)))
    list_wt_seq_reverse.append((str(df_cut_chr_i['WT_SEQ'][i])).translate(str.maketrans(dict_base)))
df_cut_chr_i['seq_reverse']=pd.Series(list_seq_reverse)
df_cut_chr_i['WT_seq_reverse']=pd.Series(list_wt_seq_reverse)

#ID_tumorのカウント
dict_tumor_id=df_cut_chr_i['ID_tumour'].value_counts().to_dict()
list_tumor_id = list(dict_tumor_id.keys())

#以下、ID_tumorごとに取り出し変異をテキストに書き込む
df_index_for_id=pd.DataFrame()
for i in list_tumor_id:
    df_index_for_id['index_{}'.format(i)]=pd.Series((df_cut_chr_i[df_cut_chr_i['ID_tumour']==i]).index)
    for j in df_index_for_id['index_{}'.format(i)]:
        if pd.isnull(j)==True:
            continue
        else:
            with open("text_2/{}.txt".format(i), "a") as f:
                if (df_cut_chr_i.iloc[int(j)])['WT_SEQ']=='C'or (df_cut_chr_i.iloc[int(j)])['WT_SEQ']=='T':
                    f.write((df_cut_chr_i.iloc[int(j)])['WT_SEQ']+'>'+(df_cut_chr_i.iloc[int(j)])['MUT_SEQ']+\
                    '_'+(df_cut_chr_i.iloc[int(j)])['seq_forward'][0:1]+\
                    (df_cut_chr_i.iloc[int(j)])['seq_forward'][2:3]+' ')
                elif(df_cut_chr_i.iloc[int(j)])['WT_SEQ']=='A'or (df_cut_chr_i.iloc[int(j)])['WT_SEQ']=='G':
                    f.write((df_cut_chr_i.iloc[int(j)])['WT_seq_reverse']+'>'+(df_cut_chr_i.iloc[int(j)])['seq_reverse'][1:2]+\
                    '_'+(df_cut_chr_i.iloc[int(j)])['seq_reverse'][0:1]+\
                    (df_cut_chr_i.iloc[int(j)])['seq_reverse'][2:3]+' ')
                else:
                    continue


