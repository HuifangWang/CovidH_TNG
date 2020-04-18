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
    df  = pd.read_csv(
            os.path.join(data_root, 'external','Global_Mobility_Report.csv'),
            parse_dates=['date']
    )

    df_france = df.loc[
        (df['country_region_code']=='FR') &
        pd.isna(df['sub_region_1'])
        , [
            'date',
            'retail_and_recreation_percent_change_from_baseline',
            'grocery_and_pharmacy_percent_change_from_baseline',
            'parks_percent_change_from_baseline',
            'transit_stations_percent_change_from_baseline',
            'workplaces_percent_change_from_baseline',
            'residential_percent_change_from_baseline'
        ]
    ].set_index('date')

    return df_france

def load_covid19_france():
    df  = pd.read_csv(
            os.path.join(
                data.data_root, 'external', 'datahub-covid-19', 'time-series-19-covid-combined_csv.csv'),
            parse_dates=['Date']
    )

    df_france = df.loc[
            (df['Country/Region']=='France') &
            pd.isna(df['Province/State']),
            [
                'Date',
                'Confirmed',
                'Recovered',
                'Deaths'
            ]
    ].set_index('Date')


