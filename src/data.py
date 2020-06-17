import os
import pandas as pd

data_root = os.path.abspath(
        os.path.join(
            os.path.dirname(
                os.path.dirname(__file__)
            ),
            'data'
        )
)

google_country_codes = {
        'france': 'FR',
        'germany': 'DE',
        'denmark': 'DK',
}

apple_country_codes = {
        'france': 'France',
        'germany': 'Germany',
        'denmark': 'Denmark',
}

datahub_country_codes = {
        'france': 'France',
        'germany': 'Germany',
        'denmark': 'Denmark',
}

ic_country_codes = {
        'france': 'France',
        'germany': 'Germany',
        'denmark': 'Denmark',
}

epiforecast_country_codes = {
        'france': 'France',
        'germany': 'Germany',
        'denmark': 'Denmark',
}

def load_mobility(country='france', provider='google', agregate=''):
    fname = f'{country}_{provider}_mobility_report{"_" + agregate if agregate else ""}.csv'
    return pd.read_csv(
            os.path.join( data_root, 'processed',fname ),
            parse_dates=['date']
    ).set_index('date')

def load_covid19(country='france'):
    return pd.read_csv(
            os.path.join(
                data_root, 'processed', f'{country}_covid19.csv'),
            parse_dates=['Date']
    ).set_index('Date')

def load_imperial_college_results(country='france'):
    return pd.read_csv(
            os.path.join(
                data_root, 'processed', f'{country}_imperial_college_sir.csv'),
            parse_dates=['date']
    ).set_index('date')

def load_cosmo_phi():
    return pd.read_csv(
            os.path.join(
                data_root, 'processed', f'germany_cosmo_phi.csv'),
            parse_dates=['date']
    ).set_index('date')

def load_rki_nowcasting():
    return pd.read_csv(
            os.path.join(
                data_root, 'processed', 'germany_rki_nowcasting.csv'),
            parse_dates=['date_of_disease_onset']
    ).set_index('date_of_disease_onset')

def date_index_to_days(df, day0):
    """
    Recalculates index of a pandas DataFrame to days from given day0.
    
    Parameters
    ----------
    df : instance of pandas.DataFrame
        Dataframe with DatetimeIndex index.
    day0: instance of pands.Timestamp 
        Date of the day_0.
    """
    return df.set_index((df.index - day0).days)

def load_imperial_college_interventions(country):
    df = pd.read_csv(
        f'{data_root}/external/imperial_college_interventions/interventions.csv',
        parse_dates=['Date effective'],
        dayfirst=True
    )
    return df.loc[df.Country==ic_country_codes[country] ]
    
def load_epiforecast_rt(country='france'):
    return pd.read_csv(
            os.path.join(
                data_root, 'processed', f'{country}_epiforecast_rt.csv'),
            parse_dates=['date']
    ).set_index('date')

def load_epiforecast_cases(country='france'):
    return pd.read_csv(
            os.path.join(
                data_root, 'processed', f'{country}_epiforecast_cases.csv'),
            parse_dates=['date']
    ).set_index('date')

