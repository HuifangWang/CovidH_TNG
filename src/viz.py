import matplotlib.dates as mdates

def format_xaxis_dates(ax):
    ax.xaxis.set_major_locator(
            mdates.WeekdayLocator( byweekday=mdates.MO )
    )
    ax.xaxis.set_major_formatter(
            mdates.DateFormatter('%b %d')
    )
