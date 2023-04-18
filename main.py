import drug_functions
import model_functions

# NNRTI and NRTI
# SEQ = 'PISPIETVPVKLKPGMDGPRVKQWPLTEEKIKALMEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSDKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKQKKSVTVLDVGDAYFSVPLDKDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQCSMTKILEPFRKQNPDLVIYQYMDDLYVGSDLEIGQHRTKIEELREHLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWT'
# PI
# SEQ = 'PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF'


def get_predict_from(seq: str, drug_group: str) -> dict:
    """
    Amino acid sequence based HIV drug resistance prediction
    :param seq: sequence
    :param drug_group: drug group (PI, NRTI, NNRTI)
    :return: response [dict]
    """

    response = {}
    seq_len = len(seq)
    drugs = drug_functions.get_drug_params()
    for drug_name in drugs:
        drug = drugs[drug_name]
        group = drugs[drug_name]['group']

        if group != drug_group:
            continue

        if not drug_functions.seq_len_check(seq_len, group):
            seq_len, seq = drug_functions.aligning(seq,group)

        model = model_functions.get_model(drug_name)
        bin_seq = drug_functions.get_bin_seq_from(seq, drug)
        response[drug_name] = model_functions.prediction_from(bin_seq, model)

    return response


# print(get_predict_from(SEQ, 'NNRTI'))
# print(drug_functions.get_drug_group_names())