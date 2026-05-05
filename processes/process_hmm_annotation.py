# process_hmm_annotation.py
import pandas as pd
import re
import sys

def process_hmm_annotations(file_path, sample):
    df = pd.read_csv(file_path, sep=' ')
    df = df.iloc[1:, :]

    transcript_dict = {}
    for scaffold, group in df.groupby('#'):
        group['subdomain'] = group['Generated'].str.split('_clustalo', expand=True)[0]
        group['start'] = group['one'].str.split('-', expand=True)[0].astype(int)

        size = scaffold.split('size')[1]
        if int(size) > 400:
            count = len(re.findall(r'CIDRa', ','.join(group['subdomain'])))
            if count > 1:
                subset = group[group['Generated'].str.contains('CIDRa')]
                subset['of'] = subset['of'].astype(float)
                selected = subset[subset['of'] == subset['of'].min()].head(1)
                group = group[~group['Generated'].str.contains('CIDRa')]
                group = pd.concat([group, selected]).sort_values('start')

            count = len(re.findall(r'DBLa', ','.join(group['subdomain'])))
            if count > 1:
                subset = group[group['Generated'].str.contains('DBLa')]
                subset['of'] = subset['of'].astype(float)
                selected = subset[subset['of'] == subset['of'].min()].head(1)
                group = group[~group['Generated'].str.contains('DBLa')]
                group = pd.concat([group, selected]).sort_values('start')

            transcript = '-'.join(group['subdomain'])
            if len(transcript) >= 3:
                transcript_dict[scaffold] = transcript

    sig_ids = list(transcript_dict.keys())
    sig_ann = list(transcript_dict.values())

    sig_id_ann = pd.DataFrame({'Assembled_id': sig_ids, 'Annotation': sig_ann})
    sig_id_ann.to_csv(f"{sample}_assembled_id_and_var_annotation.csv", index=False)

    sig_id_only = pd.DataFrame({'Assembled_id': sig_ids})
    sig_id_only.to_csv("assembled_id_sig_annotation.txt", index=False)

if __name__ == "__main__":
    file_path = sys.argv[1]
    sample_name = sys.argv[2]
    process_hmm_annotations(file_path, sample_name)
