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

def load_google_mobility_france():
    return pd.read_csv(
            os.path.join(data_root, 'processed','france_mobility_report.csv'),
            parse_dates=['date']
    ).set_index('date')

def load_covid19_france():
    return pd.read_csv(
            os.path.join(
                data_root, 'processed', 'france_covid19.csv'),
            parse_dates=['Date']
    ).set_index('Date')
