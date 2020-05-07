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
        'germany': 'DE'
}

apple_country_codes = {
        'france': 'France',
        'germany': 'Germany'
}

datahub_country_codes = {
        'france': 'France',
        'germany': 'Germany'
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
