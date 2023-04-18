import json
from aligner import smith_waterman, needleman_wunsch

FILENAME_DRUGS = 'configs/drug_params.json'
FILENAME_DRUG_GROUPS = 'configs/drug_group_params.json'


def read_json_from(filename: str) -> dict:
    with open(filename) as f:
        return json.load(f)


def get_drug_params(filename=FILENAME_DRUGS) -> dict:
    return read_json_from(filename)


def get_drug_groups(filename=FILENAME_DRUG_GROUPS) -> dict:
    return read_json_from(filename)


def get_drug_group_names(filename=FILENAME_DRUG_GROUPS) -> list:
    drug_group_params = get_drug_groups(filename)
    return [{'name': name, 'fullname': drug_group_params[name]['fullname']} for name in drug_group_params]

def seq_len_check(seq_len, group, filename=FILENAME_DRUG_GROUPS) -> bool:
    drug_group_params = get_drug_groups(filename)
    return seq_len == drug_group_params[group]['seq_len']


def aligning(seq2, group, filename=FILENAME_DRUG_GROUPS):
    drug_group_params = get_drug_groups(filename)
    seq1 = drug_group_params[group]['seq']
    aligned_seq1, aligned_seq2 = smith_waterman(seq1, seq2)
    align2 = needleman_wunsch(seq1, aligned_seq2)
    new_str_len = len(align2)
    return new_str_len, align2

def get_bin_seq_from(seq: str, drug: dict) -> list:
    bin_seq = []

    pep_len = drug['len']
    # Разделение строки на пептиди заданной длины с перекрытием 2
    for start in range(0, len(seq) - (pep_len-1), 2):
        pep = seq[start:start + pep_len]
        bin_seq.append(1 if pep in drug['peps'] else 0)
    return bin_seq
