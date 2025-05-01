#!/usr/bin/env python
#
from skyfield.api import load, Topos
from skyfield.almanac import dark_twilight_day, risings_and_settings, find_discrete
from datetime import datetime, timedelta
import pandas as pd

# Configuration
LAT, LON = 38.7083, -79.5308  # Example: Spruce Knob, WV
DAYS_AHEAD = 30
MIN_MW_ALT_DEGREES = 15  # Minimum altitude for Milky Way core to be visible

# Load ephemeris
ts = load.timescale()
eph = load('de421.bsp')
observer = Topos(latitude_degrees=LAT, longitude_degrees=LON)
obs = eph['earth'] + observer

# Date range
start_date = datetime.utcnow().date()
end_date = start_date + timedelta(days=DAYS_AHEAD)
t0 = ts.utc(start_date.year, start_date.month, start_date.day)
t1 = ts.utc(end_date.year, end_date.month, end_date.day)

# Almanac functions
dark_func = dark_twilight_day(eph, observer)
moon_func = risings_and_settings(eph, eph['Moon'], observer)

# Evaluate events
dark_times, dark_states = find_discrete(t0, t1, dark_func)
moon_times, moon_events = find_discrete(t0, t1, moon_func)

# Find overlaps between moon-down and astronomical darkness
results = []
for i in range(0, len(dark_states) - 1):
    print(f"dark_states[{i}]: {dark_states[i]}\n")
    if dark_states[i] == 2:  # Astronomical night
        print("Astronomical night\n")
        dark_start = dark_times[i].utc_datetime()
        dark_end = dark_times[i + 1].utc_datetime()

        print(f"dark_start: {dark_start}\n")
        print(f"dark_end: {dark_end}\n")

        for j in range(0, len(moon_events) - 1, 2):
            if moon_events[j] == 1:  # Moon set
                print("moon set\n")
                moon_down_start = moon_times[j].utc_datetime()
                moon_down_end = moon_times[j + 1].utc_datetime()

                # Overlap
                overlap_start = max(dark_start, moon_down_start)
                overlap_end = min(dark_end, moon_down_end)

                if overlap_start < overlap_end:
                    check_time = ts.utc(overlap_start.year, overlap_start.month, overlap_start.day,
                                        overlap_start.hour, overlap_start.minute)

                    # Check Milky Way core elevation (~Sagittarius A*, Galactic center)
                    galactic_ra_hours = 17 + (45 / 60)  # 17h45m
                    galactic_dec_deg = -29.0
                    from skyfield.positionlib import ICRF
                    from skyfield.units import Angle
                    galactic_center = eph['earth'].at(check_time).observe(ICRF(
                        ra=Angle(hours=galactic_ra_hours),
                        dec=Angle(degrees=galactic_dec_deg),
                        distance=1.0  # Arbitrary nonzero
                    )).apparent()

                    alt, az, _ = galactic_center.altaz()

                    if alt.degrees >= MIN_MW_ALT_DEGREES:
                        results.append({
                            'Date': overlap_start.date(),
                            'Start (UTC)': overlap_start.strftime('%H:%M'),
                            'End (UTC)': overlap_end.strftime('%H:%M'),
                            'MW Altitude (°)': round(alt.degrees, 1)
                        })

# Output results
df = pd.DataFrame(results)
print("\nBest Milky Way Viewing Windows (Offline, No Weather):\n")
if df.empty:
    print("No suitable windows found in the next", DAYS_AHEAD, "days.")
else:
    print(df.to_string(index=False))
