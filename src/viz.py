import matplotlib.dates as mdates
import matplotlib.pylab as plt
import matplotlib
from src import data

def format_xaxis_dates(ax, rotation=45):
    ax.xaxis.set_major_locator(
            mdates.WeekdayLocator( byweekday=mdates.MO )
    )
    ax.xaxis.set_major_formatter(
            mdates.DateFormatter('%b %d')
    )
    ax.xaxis.set_tick_params(rotation=rotation)

def annotate_interventions(ax,country, legend=False, bbox_offset=-0.3):
    types = [
            'Schools + Universities',
            'Public events',
            'Lockdown',
            'Social distancing encouraged',
            'Self-isolating if ill',
            'Closure of cultural institutions',
            'Closure of restaurants',
            'Advice to work from home',
            'Keep distance from others',
            'Work from home',
            'Stay at home',
            'Public gatherings',
            ]
    markers = ['.', '+', 'x', '*', 'v', '^', '<', '>',  'o', '1', '2', '3',
            '4', '8', 's', 'p',  'h', 'H', 'D', 'd', '|', '_', 'P', 'X', 0,
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, ','] 
    type_markers = dict(zip(types, markers))
    df = data.load_imperial_college_interventions(country)
    for i,(_, (date, t)) in enumerate(df[['Date effective','Type']].iterrows()):
        ax.axvline(date, ymax=0.9+0.01*i, marker=type_markers[t], markevery=(1,1), lw=1, c='k', ls=':')

    if legend:
        legend_elements = [
            plt.Line2D([0], [0], lw=0,marker=type_markers[t], color='k',
                label=t, markerfacecolor='k', markersize=10)
            for _,(date,t,event) in df[['Date effective','Type','Event']].iterrows()
        ]

        legends = [c for c in ax.get_children() if isinstance(c, matplotlib.legend.Legend)]

        ax.legend(handles=legend_elements, loc='center', bbox_to_anchor=(0.5, bbox_offset), ncol=2)

        for legend in legends:
            ax.add_artist(legend)


