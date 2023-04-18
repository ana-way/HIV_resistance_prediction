from joblib import load as joblib_load


def get_model(drug_name: str):
    return joblib_load(f'models/{drug_name}.joblib')


def prediction_from(bin_seq: list, model) -> str:
    predict = model.predict([bin_seq])
    return 'sensitive' if predict[0] == 'N' else 'resistance'
