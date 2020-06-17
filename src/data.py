import os
import pandas as pd
import scipy.io as sio 

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


cosmo_dates = {
        'germany':[
            '2020-03-03',
            '2020-03-10',
            '2020-03-17',
            '2020-03-24',
            '2020-03-31',
            '2020-04-07',
            '2020-04-14',
            '2020-04-21',
            '2020-04-28',
            '2020-05-05',
            '2020-05-12',
            '2020-05-19',
            '2020-05-26',
            #'2020-06-09', # not yet
        ],
        'denmark': [ # let's go with mid-point wednesdays
            '2020-03-25', #Week 13: 23/03 - 29/03/2020
            '2020-04-01', #Week 14: 30/03 - 05/04/2020
            '2020-04-08', #Week 15: 06/04 - 12/04/2020
            '2020-04-15', #Week 16 - 13/04 - 19/04/2020
            '2020-04-22', #Week 17 - 20/04 - 26/04/2020
            '2020-04-29', #Week 18 - 27/04 - 03/05/2020
            '2020-05-06', #Week 19 - 04/05 - 10/05/2020
            '2020-05-13', #Week 20 - 11/05 - 17/05/2020
        ]
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

def load_cosmo_phi(country='germany'):
    return pd.read_csv(
            os.path.join(
                data_root, 'processed', f'{country}_cosmo_phi.csv'),
            parse_dates=['date']
    ).set_index('date')

def load_raw_cosmo_de(file_path):
        return pd.read_spss(file_path, convert_categoricals=False)

def load_raw_cosmo_dk(file_path):
        dkmat = sio.loadmat( file_path )
        colnames = [ name[0] for name in dkmat["COSMO_DK_varnames"].squeeze() ]
        df = pd.DataFrame( 
            data=dkmat["COSMO_DK"], 
            columns=colnames
        )
        return df.rename(
                axis='columns',
                mapper={ 'wave':'TIME' } 
        )

def phi_t(df):
    df_aff = df.filter(
            regex='^TIME|^AFF_FEAR|^AFF_THINK|^AFF_WORRY',
            axis=1
    ).groupby(
            'TIME'
    ).mean()
    df_aff = (df_aff -1 ) / 6 # data range 1-7 -> 0-1    

    return  df_aff.mean(axis=1).to_frame(name='phi')

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

