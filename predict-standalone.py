from skyfield.api import load, Topos
from skyfield.almanac import dark_twilight_day, risings_and_settings
from datetime import datetime, timedelta
import pandas as pd

# Configuration
LAT, LON = 38.7083, -79.5308  # Spruce Knob, WV
DAYS_AHEAD = 14
MIN_MW_ALT_DEGREES = 15  # minimum altitude for Milky Way core

# Load ephemeris
ts = load.timescale()
eph = load('de421.bsp')
obs = eph['earth'] + Topos(latitude_degrees=LAT, longitude_degrees=LON)

# Date range
start_date = datetime.utcnow().date()
end_date = start_date + timedelta(days=DAYS_AHEAD)
t0 = ts.utc(start_date.year, start_date.month, start_date.day)
t1 = ts.utc(end_date.year, end_date.month, end_date.day)

# Get darkness periods
eph_f = dark_twilight_day(eph, Topos(latitude_degrees=LAT, longitude_degrees=LON))
times, states = risings_and_settings(eph, eph['Moon'], Topos(latitude_degrees=LAT, longitude_degrees=LON))
dark_times, dark_states = load('de421.bsp'), eph_f

# Find darkness periods
from skyfield.almanac import find_discrete

dark_times, dark_states = find_discrete(t0, t1, eph_f)
moon_times, moon_events = find_discrete(t0, t1, risings_and_settings(eph, eph['Moon'], Topos(latitude_degrees=LAT, longitude_degrees=LON)))

# Filter for periods when it's dark and the moon is down
results = []
for i in range(0, len(dark_states)-1):
    if dark_states[i] == 2:  # Astronomical night
        dark_start = dark_times[i]
        dark_end = dark_times[i+1]
        
        for j in range(0, len(moon_events)-1, 2):
            if moon_events[j] == 1:  # Moon set
                moon_down_start = moon_times[j]
                moon_down_end = moon_times[j+1]
                
                # Overlap of moon-down and astronomical night
                start = max(dark_start.utc_datetime(), moon_down_start.utc_datetime())
                end = min(dark_end.utc_datetime(), moon_down_end.utc_datetime())
                
                if start < end:
                    # Check if MW core is visible
                    check_time = ts.utc(start.year, start.month, start.day, start.hour, start.minute)
                    mw = obs.at(check_time).observe(eph['galactic center']).apparent()
                    alt, az, _ = mw.altaz()
                    if alt.degrees > MIN_MW_ALT_DEGREES:
                        results.append({
                            'Date': start.date(),
                            'Start (UTC)': start.strftime('%H:%M'),
                            'End (UTC)': end.strftime('%H:%M'),
                            'MW Alt (°)': round(alt.degrees, 1)
                        })

# Display results
df = pd.DataFrame(results)
print("\nBest Milky Way Viewing Windows (Offline, No Weather):\n")
print(df.to_string(index=False))
